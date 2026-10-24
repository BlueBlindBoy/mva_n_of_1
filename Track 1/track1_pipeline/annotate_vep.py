"""Annotate variants via the Ensembl VEP REST API (individual GET with mane_select=1)."""
from __future__ import annotations
import time
import requests
from dataclasses import dataclass, field

VEP_BASE = "https://rest.ensembl.org/vep/human/hgvs"


@dataclass
class AnnotationResult:
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str | None = None
    consequence: str | None = None
    impact: str | None = None
    clinvar_sig: list[str] = field(default_factory=list)
    clinvar_id: str | None = None
    gnomad_exome_af: float | None = None
    gnomad_genome_af: float | None = None
    sift: str | None = None
    sift_score: float | None = None
    polyphen: str | None = None
    polyphen_score: float | None = None
    domains: list[str] = field(default_factory=list)
    hgvs: str | None = None
    hgvs_p: str | None = None
    dbsnp_id: str | None = None
    is_novel: bool = False
    nmd_likely: bool | None = None
    exon_rank: int | None = None
    total_exons: int | None = None
    splicing_merged_score: float | None = None
    splicing_sites_score: float | None = None
    splicing_site_usage_score: float | None = None
    splicing_junctions_score: float | None = None
    rna_seq_log2fc: float | None = None
    avi_composite: float | None = None
    avi_top_modality: str | None = None

    @property
    def is_lof(self) -> bool:
        return self.consequence in ("stop_gained", "frameshift_variant",
                                     "splice_acceptor_variant", "splice_donor_variant",
                                     "start_lost", "stop_lost")

    @property
    def is_missense(self) -> bool:
        return self.consequence == "missense_variant"

    @property
    def is_pathogenic(self) -> bool:
        sig_lower = [s.lower() for s in self.clinvar_sig]
        return any("pathogenic" in s and "benign" not in s for s in sig_lower)

    @property
    def is_rare(self) -> bool:
        """A variant is rare if its gnomAD AF is below threshold.
        Absence of frequency data is treated as rare (absence of evidence
        for commonness), NOT as evidence for commonness. This prevents
        filtering out ClinVar pathogenic variants that lack gnomAD coverage."""
        threshold = 0.01
        exome_ok = self.gnomad_exome_af is None or self.gnomad_exome_af < threshold
        genome_ok = self.gnomad_genome_af is None or self.gnomad_genome_af < threshold
        return exome_ok and genome_ok

    @property
    def is_damaging_missense(self) -> bool:
        return (self.is_missense
                and self.sift in ("deleterious", "deleterious_low_confidence")
                and self.polyphen in ("probably_damaging", "possibly_damaging"))

    @property
    def chrom_with_prefix(self) -> str:
        return self.chrom if self.chrom.startswith("chr") else f"chr{self.chrom}"


def build_hgvs_notation(chrom: str, pos: int, ref: str, alt: str) -> str:
    c = chrom if chrom.startswith("chr") else f"chr{chrom}"
    return f"{c}:g.{pos}{ref}>{alt}"


def annotate_variants_batch(
    variants: list[tuple[str, int, str, str]],
    rate_limit_delay: float = 0.3,
) -> list[AnnotationResult]:
    """Annotate variants one-by-one via VEP GET with mane_select=1.

    The batch POST endpoint silently drops variants and picks non-MANE
    transcripts. Individual GET with mane_select=1 fixes both.
    """
    all_results: list[AnnotationResult] = []
    for chrom, pos, ref, alt in variants:
        hgvs = build_hgvs_notation(chrom, pos, ref, alt)
        url = f"{VEP_BASE}/{hgvs}?domains=1&mane_select=1"
        result = None
        for attempt in range(3):
            try:
                resp = requests.get(url, headers={"Accept": "application/json"}, timeout=60)
                if resp.status_code == 429:
                    time.sleep(5 * (attempt + 1))
                    continue
                if resp.status_code == 404:
                    result = AnnotationResult(chrom=chrom, pos=pos, ref=ref, alt=alt)
                    break
                resp.raise_for_status()
                data = resp.json()
                if isinstance(data, list) and len(data) > 0:
                    result = _parse_vep_entry(data[0], chrom, pos, ref, alt)
                else:
                    result = _parse_vep_entry(data, chrom, pos, ref, alt)
                break
            except requests.RequestException:
                if attempt == 2:
                    result = AnnotationResult(chrom=chrom, pos=pos, ref=ref, alt=alt)
                time.sleep(3 * (attempt + 1))
        if result is None:
            result = AnnotationResult(chrom=chrom, pos=pos, ref=ref, alt=alt)
        all_results.append(result)
        time.sleep(rate_limit_delay)
    return all_results


def _pick_best_transcript(tc_list: list[dict]) -> dict | None:
    """Pick the MANE-select transcript; fall back to the first with a
    non-MODIFIER impact (HIGH > MODERATE > LOW > MODIFIER)."""
    if not tc_list:
        return None
    # Prefer MANE-select
    for tc in tc_list:
        if tc.get("mane_select"):
            return tc
    # Prefer highest impact
    impact_order = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}
    sorted_tc = sorted(tc_list, key=lambda t: impact_order.get(t.get("impact", "MODIFIER"), 4))
    return sorted_tc[0]


def _parse_vep_entry(entry: dict, chrom: str, pos: int, ref: str, alt: str) -> AnnotationResult:
    result = AnnotationResult(chrom=chrom, pos=pos, ref=ref, alt=alt, hgvs=entry.get("id"))
    tc_list = entry.get("transcript_consequences", [])

    # Filter to BUB1B transcripts if present (the gene we care about)
    bub1b_tcs = [tc for tc in tc_list if tc.get("gene_symbol") == "BUB1B"]
    if bub1b_tcs:
        tc_list = bub1b_tcs

    best_tc = _pick_best_transcript(tc_list)
    if best_tc:
        result.gene = best_tc.get("gene_symbol")
        csqs = best_tc.get("consequence_terms", [])
        result.consequence = csqs[0] if csqs else None
        result.impact = best_tc.get("impact")
        result.sift = best_tc.get("sift_prediction")
        result.sift_score = best_tc.get("sift_score")
        result.polyphen = best_tc.get("polyphen_prediction")
        result.polyphen_score = best_tc.get("polyphen_score")
        result.hgvs_p = best_tc.get("hgvsp")
        domains = best_tc.get("domains")
        if domains:
            result.domains = [d.get("name", "") for d in domains if d.get("name")]

    cv_list = entry.get("colocated_variants", [])
    result.is_novel = len(cv_list) == 0
    for cv in cv_list:
        sig = cv.get("clin_sig", [])
        if sig:
            result.clinvar_sig = sig
        freqs = cv.get("frequencies", {}).get(alt, {})
        if freqs:
            result.gnomad_exome_af = freqs.get("gnomade")
            result.gnomad_genome_af = freqs.get("gnomadg")
        result.dbsnp_id = cv.get("id")
        if sig:
            break
    return result
