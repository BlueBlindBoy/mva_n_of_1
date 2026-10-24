"""Full pipeline orchestration: VCF -> HPO genes -> ClinGen -> VEP annotation ->
NMD prediction -> AlphaGenome multi-modal scoring -> ranking -> CSV."""
from __future__ import annotations
import os
import json
import sys
from pathlib import Path

# Ensure the Track 1 folder is on the path (the package is track1_pipeline)
sys.path.insert(0, str(Path(__file__).parent.parent))

from track1_pipeline.parse_vcf import parse_vcf
from track1_pipeline.hpo_genes import get_candidate_genes, HPO_TERMS
from track1_pipeline.clingen_validity import get_clingen_validity
from track1_pipeline.annotate_vep import annotate_variants_batch
from track1_pipeline.nmd_predict import predict_nmd
from track1_pipeline.rank_variants import rank_variants
from track1_pipeline.generate_submission import generate_submission_csv

MVA_KNOWN_GENES = ["BUB1B", "CEP57", "TRIP13"]
BUB1B_CANONICAL_TRANSCRIPT = "ENST00000287598"


def run_pipeline(
    vcf_path: str,
    output_path: str = "output/bit-by-bit_hpo-clingen-nmd.csv",
    proband_id: str = "PROBAND01",
    hpo_terms: list[str] | None = None,
    epcr: float = 0.85,
) -> dict:
    if hpo_terms is None:
        hpo_terms = HPO_TERMS

    # Stage 1: Get candidate genes
    print("=" * 60)
    print("Stage 1: Resolving candidate genes from HPO terms")
    print("=" * 60)
    candidate_genes = []
    try:
        hpo_candidates = get_candidate_genes(hpo_terms)
        candidate_genes = [g for g, _ in hpo_candidates[:50]]
        print(f"  HPO API returned {len(candidate_genes)} candidate genes")
        print(f"  Top 10: {[g for g, c in hpo_candidates[:10]]}")
    except Exception as e:
        print(f"  Warning: HPO API query failed: {e}")

    for g in MVA_KNOWN_GENES:
        if g not in candidate_genes:
            candidate_genes.append(g)
            print(f"  Added known MVA gene: {g}")

    print(f"\n  ClinGen gene-disease validity:")
    for g in candidate_genes[:10]:
        cv = get_clingen_validity(g)
        if cv["classification"]:
            print(f"    {g}: {cv['classification']} ({cv['disease']})")

    # Stage 2: Resolve gene coordinates
    print("\n" + "=" * 60)
    print("Stage 2: Resolving gene coordinates")
    print("=" * 60)
    gene_regions = _resolve_gene_regions(candidate_genes)
    print(f"  Resolved {len(gene_regions)} gene regions")

    # Stage 3: Parse VCF variants in candidate gene regions
    print("\n" + "=" * 60)
    print("Stage 3: Parsing VCF variants in candidate gene regions")
    print("=" * 60)
    all_variants = []
    for gene, (chrom, start, end) in gene_regions.items():
        variants = parse_vcf(vcf_path, chrom=chrom, start=start, end=end, het_only=True)
        print(f"  {gene} ({chrom}:{start}-{end}): {len(variants)} het PASS variants")
        all_variants.extend(variants)
    print(f"\n  Total variants to annotate: {len(all_variants)}")

    # Stage 4: Annotate via Ensembl VEP (with domains=1)
    print("\n" + "=" * 60)
    print("Stage 4: Annotating variants via Ensembl VEP (domains=1)")
    print("=" * 60)
    variant_tuples = [(v.chrom, v.pos, v.ref, v.alt) for v in all_variants]
    annotations = annotate_variants_batch(variant_tuples)
    print(f"  Annotated {len(annotations)} variants")

    # Stage 5: Predict NMD for LoF variants
    print("\n" + "=" * 60)
    print("Stage 5: Predicting NMD for LoF variants")
    print("=" * 60)
    for ann in annotations:
        if ann.is_lof and ann.gene == "BUB1B":
            nmd_result = predict_nmd(
                variant_pos=ann.pos,
                transcript_id=BUB1B_CANONICAL_TRANSCRIPT,
            )
            ann.nmd_likely = nmd_result.nmd_likely
            ann.exon_rank = nmd_result.exon_rank
            ann.total_exons = nmd_result.total_exons
            print(f"  {ann.chrom}:{ann.pos} {ann.consequence}: "
                  f"NMD={nmd_result.nmd_likely}, "
                  f"exon {nmd_result.exon_rank}/{nmd_result.total_exons}")

    # Stage 5b: AlphaGenome multi-modal scoring
    print("\n" + "=" * 60)
    print("Stage 5b: AlphaGenome multi-modal scoring")
    print("=" * 60)
    api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if api_key:
        from track1_pipeline.alphagenome_annotate import score_variant_full
        for ann in annotations:
            if ann.gene == "BUB1B" and ann.impact in ("HIGH", "MODERATE"):
                try:
                    ag = score_variant_full(
                        ann.chrom, ann.pos, ann.ref, ann.alt, api_key=api_key,
                    )
                    ann.splicing_merged_score = ag.splicing_merged
                    ann.splicing_sites_score = ag.splicing_sites
                    ann.splicing_site_usage_score = ag.splicing_site_usage
                    ann.splicing_junctions_score = ag.splicing_junctions
                    ann.rna_seq_log2fc = ag.rna_seq_log2fc
                    ann.avi_composite = ag.composite_avi
                    ann.avi_top_modality = ag.top_modality
                    print(f"  {ann.chrom}:{ann.pos} {ann.consequence}: "
                          f"merged splicing={ag.splicing_merged:.3f} "
                          f"(large: {ag.is_large_splicing_effect}), "
                          f"composite AVI={ag.composite_avi:.3f} "
                          f"(top: {ag.top_modality})")
                except Exception as e:
                    print(f"  Warning: AlphaGenome query failed for {ann.pos}: {e}")
    else:
        print("  ALPHAGENOME_API_KEY not set; skipping AlphaGenome annotation.")

    # Stage 6: Filter and rank
    print("\n" + "=" * 60)
    print("Stage 6: Filtering and ranking variants")
    print("=" * 60)
    ranked = rank_variants(annotations)
    print(f"  {len(ranked)} variants survived filtering")
    for vs in ranked[:10]:
        ch_info = ""
        if vs.is_compound_het_partner:
            ch_info = f" [compound-het partner at {vs.compound_het_partner_pos}]"
        print(f"  Rank {vs.rank}: {vs.annotation.chrom}:{vs.annotation.pos} "
              f"{vs.annotation.ref}>{vs.annotation.alt} "
              f"gene={vs.annotation.gene} "
              f"score={vs.score:.1f} "
              f"ACMG={vs.acmg_classification} "
              f"{vs.rationale[:100]}{ch_info}")

    # Stage 7: Generate submission CSV
    print("\n" + "=" * 60)
    print("Stage 7: Generating submission CSV")
    print("=" * 60)
    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    generate_submission_csv(ranked, output_path, proband_id=proband_id, epcr=epcr)
    print(f"  Written to: {output_path}")

    print("\n--- Submission CSV ---")
    print(Path(output_path).read_text())

    return {
        "n_candidate_genes": len(candidate_genes),
        "n_variants_parsed": len(all_variants),
        "n_variants_annotated": len(annotations),
        "n_variants_after_filter": len(ranked),
        "top_variant": f"{ranked[0].annotation.chrom}:{ranked[0].annotation.pos}" if ranked else None,
        "is_compound_het": ranked[0].is_compound_het_partner if ranked else False,
        "top_acmg": ranked[0].acmg_classification if ranked else None,
        "output_path": output_path,
    }


def _resolve_gene_regions(genes: list[str]) -> dict[str, tuple[str, int, int]]:
    import requests
    regions = {}
    for gene in genes:
        url = f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}?content-type=application/json"
        try:
            resp = requests.get(url, timeout=30)
            if resp.status_code == 200:
                data = resp.json()
                chrom = data.get("seq_region_name")
                start = data.get("start")
                end = data.get("end")
                if chrom and start and end:
                    regions[gene] = (chrom, start, end)
        except Exception as e:
            print(f"  Error resolving {gene}: {e}")
    return regions


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Run Track 1 variant prediction pipeline")
    parser.add_argument("--vcf", required=True, help="Path to gzipped VCF")
    parser.add_argument("--output", default="output/bit-by-bit_hpo-clingen-nmd.csv")
    parser.add_argument("--proband-id", default="PROBAND01")
    parser.add_argument("--epcr", type=float, default=0.85)
    args = parser.parse_args()

    summary = run_pipeline(vcf_path=args.vcf, output_path=args.output,
                           proband_id=args.proband_id, epcr=args.epcr)
    print(f"\nPipeline complete: {json.dumps(summary, indent=2)}")
