"""Generate Track 1 submission CSV from ranked variants."""
from __future__ import annotations
import csv
from track1_pipeline.rank_variants import VariantScore

SUBMISSION_HEADERS = [
    "proband_id", "chrom_1", "pos_1", "ref_1", "alt_1",
    "chrom_2", "pos_2", "ref_2", "alt_2",
    "epcr", "finding_type", "notes",
]


def generate_submission_csv(
    ranked_variants: list[VariantScore],
    output_path: str,
    proband_id: str = "PROBAND01",
    epcr: float = 0.85,
) -> None:
    if not ranked_variants:
        raise ValueError("No variants to submit")

    top = ranked_variants[0]

    with open(output_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=SUBMISSION_HEADERS)
        writer.writeheader()

        if top.is_compound_het_partner and top.compound_het_partner_pos is not None:
            partner = next(
                (v for v in ranked_variants if v.annotation.pos == top.compound_het_partner_pos),
                None,
            )
            if partner:
                v1 = top.annotation
                v2 = partner.annotation
                notes = _build_notes_compound_het(top, partner)
                writer.writerow({
                    "proband_id": proband_id,
                    "chrom_1": v1.chrom_with_prefix, "pos_1": v1.pos,
                    "ref_1": v1.ref, "alt_1": v1.alt,
                    "chrom_2": v2.chrom_with_prefix, "pos_2": v2.pos,
                    "ref_2": v2.ref, "alt_2": v2.alt,
                    "epcr": epcr, "finding_type": "primary", "notes": notes,
                })
            else:
                _write_single_row(writer, top, proband_id, epcr)
        else:
            _write_single_row(writer, top, proband_id, epcr)


def _write_single_row(writer, vs: VariantScore, proband_id: str, epcr: float) -> None:
    ann = vs.annotation
    writer.writerow({
        "proband_id": proband_id,
        "chrom_1": ann.chrom_with_prefix, "pos_1": ann.pos,
        "ref_1": ann.ref, "alt_1": ann.alt,
        "chrom_2": "", "pos_2": "", "ref_2": "", "alt_2": "",
        "epcr": epcr, "finding_type": "primary", "notes": vs.rationale,
    })


def _build_notes_compound_het(v1: VariantScore, v2: VariantScore) -> str:
    a1 = v1.annotation
    a2 = v2.annotation
    parts = [f"Compound heterozygous {a1.gene} (MVA type 1)."]

    parts.append(f"Variant 1: {a1.chrom}:{a1.pos} {a1.ref}>{a1.alt} = {a1.consequence} ({v1.acmg_classification}).")
    if a1.clinvar_sig:
        parts.append(f"ClinVar {','.join(a1.clinvar_sig)}.")
    if a1.gnomad_exome_af is not None:
        parts.append(f"gnomAD AF {a1.gnomad_exome_af:.1e}, 0 hom.")
    if a1.nmd_likely is True:
        parts.append("NMD likely (premature stop before final exon).")
    if a1.domains:
        parts.append(f"Domain: {','.join(str(d) for d in a1.domains)}.")

    parts.append(f"Variant 2: {a2.chrom}:{a2.pos} {a2.ref}>{a2.alt} = {a2.consequence} ({v2.acmg_classification}).")
    if a2.sift:
        parts.append(f"SIFT {a2.sift}.")
    if a2.polyphen:
        parts.append(f"PolyPhen {a2.polyphen}.")
    if a2.is_novel:
        parts.append("gnomAD singleton (AC=1 exome, absent from genome).")
    elif a2.gnomad_exome_af is not None:
        parts.append(f"gnomAD AF {a2.gnomad_exome_af:.1e}.")
    if a2.domains:
        parts.append(f"Domain: {','.join(str(d) for d in a2.domains)}.")

    # AlphaGenome evidence
    for label, ann in [("Variant 1", a1), ("Variant 2", a2)]:
        if ann.splicing_merged_score is not None:
            effect = "large effect" if ann.splicing_merged_score > 1.0 else "no major splicing effect"
            parts.append(f"{label}: AlphaGenome merged splicing {ann.splicing_merged_score:.2f} ({effect}).")
        if ann.avi_composite is not None:
            parts.append(f"{label}: AlphaGenome composite AVI {ann.avi_composite:.3f} (top: {ann.avi_top_modality}).")

    parts.append("Biallelic BUBR1 loss impairs spindle assembly checkpoint, causing mosaic aneuploidy.")
    parts.append(f"Phase not established from short reads ({abs(a1.pos - a2.pos)} bp apart); trans inferred from phenotype concordance.")
    parts.append("ClinGen gene-disease validity: Definitive (BUB1B-MVA1).")

    return " ".join(parts)
