"""Rank annotated variants by pathogenicity evidence."""
from __future__ import annotations
from dataclasses import dataclass
from track1_pipeline.annotate_vep import AnnotationResult
from track1_pipeline.clingen_validity import get_clingen_score, get_clingen_validity


@dataclass
class VariantScore:
    annotation: AnnotationResult
    score: float
    rank: int
    is_compound_het_partner: bool = False
    compound_het_partner_pos: int | None = None
    rationale: str = ""
    acmg_classification: str = ""


def rank_variants(
    annotations: list[AnnotationResult],
    max_gnomad_af: float = 0.01,
    min_quality_filter: bool = True,
) -> list[VariantScore]:
    # Filter: keep variants that are (a) rare OR (b) ClinVar pathogenic/LP.
    # Never filter out ClinVar P/LP even if gnomAD AF is missing.
    filtered = []
    for ann in annotations:
        if min_quality_filter:
            if ann.impact in ("LOW", "MODIFIER", None):
                # Allow MODIFIER/LOW only if ClinVar pathogenic (e.g., splice region)
                if not ann.is_pathogenic:
                    continue
            if not ann.is_rare and not ann.is_pathogenic:
                continue
        filtered.append(ann)

    # Score
    scored = []
    for ann in filtered:
        score, rationale, acmg = _score_variant(ann)
        scored.append((score, ann, rationale, acmg))
    scored.sort(key=lambda x: -x[0])

    # Compound-het detection
    gene_variants: dict[str, list[AnnotationResult]] = {}
    for score, ann, _, _ in scored:
        if ann.gene:
            gene_variants.setdefault(ann.gene, []).append(ann)

    compound_het_positions: set[tuple[int, int]] = set()
    for gene, vars_in_gene in gene_variants.items():
        if len(vars_in_gene) >= 2:
            p1 = vars_in_gene[0].pos
            p2 = vars_in_gene[1].pos
            compound_het_positions.add((p1, p2))

    # Build results
    results = []
    for i, (score, ann, rationale, acmg) in enumerate(scored):
        ch_partner = None
        is_ch = False
        for p1, p2 in compound_het_positions:
            if ann.pos == p1:
                is_ch = True
                ch_partner = p2
                break
            elif ann.pos == p2:
                is_ch = True
                ch_partner = p1
                break
        results.append(VariantScore(
            annotation=ann, score=score, rank=i + 1,
            is_compound_het_partner=is_ch,
            compound_het_partner_pos=ch_partner,
            rationale=rationale, acmg_classification=acmg,
        ))
    return results


def _score_variant(ann: AnnotationResult) -> tuple[float, str, str]:
    clingen_score = get_clingen_score(ann.gene) if ann.gene else 0.0
    base_score, rationale, acmg = _base_score(ann)
    if clingen_score == 0.0:
        final_score = base_score * 0.5
    else:
        final_score = base_score * clingen_score
    if ann.gene:
        cv = get_clingen_validity(ann.gene)
        if cv["classification"]:
            rationale += f" ClinGen gene-disease validity: {cv['classification']}."
    return final_score, rationale, acmg


def _base_score(ann: AnnotationResult) -> tuple[float, str, str]:
    if ann.is_lof:
        nmd_note = ""
        nmd_bonus = 0
        if ann.nmd_likely is True:
            nmd_note = " NMD likely (stop before last exon)."
            nmd_bonus = 5
        elif ann.nmd_likely is False:
            nmd_note = " NMD unlikely (stop in/near last exon)."

        # AlphaGenome splicing bonus
        splicing_note = ""
        splicing_bonus = 0
        if ann.splicing_merged_score is not None and ann.splicing_merged_score > 1.0:
            splicing_note = " AlphaGenome predicts large splicing disruption."
            splicing_bonus = 5

        if ann.is_pathogenic:
            score = 100 + nmd_bonus + splicing_bonus
            clinvar_str = ",".join(ann.clinvar_sig)
            return score, f"LoF ({ann.consequence}) + ClinVar {clinvar_str}.{nmd_note}{splicing_note}", "Likely Pathogenic"
        if ann.is_novel:
            score = 85 + nmd_bonus + splicing_bonus
            return score, f"LoF ({ann.consequence}), novel (no gnomAD/ClinVar).{nmd_note}{splicing_note}", "VUS"
        score = 80 + nmd_bonus + splicing_bonus
        af = ann.gnomad_exome_af if ann.gnomad_exome_af else "unknown"
        return score, f"LoF ({ann.consequence}), rare (gnomAD AF={af}).{nmd_note}{splicing_note}", "VUS"

    if ann.is_damaging_missense:
        domain_note = " In functional domain." if ann.domains else ""
        splicing_note = ""
        splicing_bonus = 0
        if ann.splicing_merged_score is not None and ann.splicing_merged_score > 1.0:
            splicing_note = " AlphaGenome predicts large splicing disruption."
            splicing_bonus = 5

        if ann.is_pathogenic:
            score = 70 + splicing_bonus
            return score, f"Damaging missense + ClinVar {','.join(ann.clinvar_sig)}.{domain_note}{splicing_note}", "Likely Pathogenic"
        if ann.is_novel or (ann.gnomad_exome_af is not None and ann.gnomad_exome_af < 1e-5):
            score = (60 if ann.domains else 55) + splicing_bonus
            return score, f"Damaging missense, novel/ultra-rare (SIFT={ann.sift}, PolyPhen={ann.polyphen}).{domain_note}{splicing_note}", "VUS"
        score = 50 + splicing_bonus
        return score, f"Damaging missense, rare (gnomAD AF={ann.gnomad_exome_af}).{domain_note}{splicing_note}", "VUS"

    if ann.is_missense:
        return 40, f"Missense, rare (SIFT={ann.sift}, PolyPhen={ann.polyphen}).", "VUS"

    return 10, f"Coding variant ({ann.consequence}), rare.", "VUS"
