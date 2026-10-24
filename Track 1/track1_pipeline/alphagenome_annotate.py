"""Query AlphaGenome for variant-level splicing, expression, and multi-modal impact.

Uses the AlphaGenome Python SDK to score variants across all available modalities.
Produces a composite AVI-like score by taking the max absolute score across all
modalities, plus a merged splicing score (splice_sites + splice_site_usage +
splice_junctions/5).
"""
from __future__ import annotations
import os
from dataclasses import dataclass

MERGED_SPLICING_THRESHOLD = 1.0


@dataclass
class AlphaGenomeResult:
    variant_chrom: str
    variant_pos: int
    ref: str
    alt: str
    splicing_merged: float | None
    splicing_sites: float | None
    splicing_site_usage: float | None
    splicing_junctions: float | None
    rna_seq_log2fc: float | None
    cage_max: float | None
    atac_max: float | None
    dnase_max: float | None
    chip_tf_max: float | None
    chip_histone_max: float | None
    procap_max: float | None
    composite_avi: float | None
    top_modality: str | None

    @property
    def is_large_splicing_effect(self) -> bool:
        if self.splicing_merged is None:
            return False
        return self.splicing_merged > MERGED_SPLICING_THRESHOLD


def score_variant_full(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    api_key: str | None = None,
) -> AlphaGenomeResult:
    """Score a variant across all AlphaGenome modalities in a single API call."""
    from alphagenome.data import genome
    from alphagenome.models import dna_client, variant_scorers

    if api_key is None:
        api_key = os.environ.get("ALPHAGENOME_API_KEY")
    if api_key is None:
        raise RuntimeError("ALPHAGENOME_API_KEY not set")

    model = dna_client.create(api_key)
    c = chrom if chrom.startswith("chr") else f"chr{chrom}"
    variant = genome.Variant(chromosome=c, position=pos,
                             reference_bases=ref, alternate_bases=alt)
    interval = variant.reference_interval.resize(dna_client.SEQUENCE_LENGTH_1MB)

    all_scorers = [
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["RNA_SEQ"],
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITES"],
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_SITE_USAGE"],
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["SPLICE_JUNCTIONS"],
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["CAGE"],
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["ATAC"],
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["DNASE"],
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["CHIP_TF"],
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["CHIP_HISTONE"],
        variant_scorers.RECOMMENDED_VARIANT_SCORERS["PROCAP"],
    ]

    scores = model.score_variant(
        interval=interval, variant=variant,
        variant_scorers=all_scorers,
        organism=dna_client.Organism.HOMO_SAPIENS,
    )
    df = variant_scorers.tidy_scores([scores])

    def _max_for(output_type: str) -> float:
        subset = df[df["output_type"] == output_type]
        return float(subset["raw_score"].max()) if len(subset) > 0 else 0.0

    ss = _max_for("SPLICE_SITES")
    su = _max_for("SPLICE_SITE_USAGE")
    sj = _max_for("SPLICE_JUNCTIONS")
    merged = ss + su + sj / 5.0

    rna = _max_for("RNA_SEQ")
    cage = _max_for("CAGE")
    atac = _max_for("ATAC")
    dnase = _max_for("DNASE")
    chip_tf = _max_for("CHIP_TF")
    chip_h = _max_for("CHIP_HISTONE")
    procap = _max_for("PROCAP")

    modality_maxes = {
        "RNA_SEQ": rna, "SPLICE_SITES": ss, "SPLICE_SITE_USAGE": su,
        "SPLICE_JUNCTIONS": sj, "CAGE": cage, "ATAC": atac,
        "DNASE": dnase, "CHIP_TF": chip_tf, "CHIP_HISTONE": chip_h,
        "PROCAP": procap,
    }
    composite = max(abs(v) for v in modality_maxes.values()) if modality_maxes else 0.0
    top_mod = max(modality_maxes, key=lambda k: abs(modality_maxes[k])) if modality_maxes else None

    return AlphaGenomeResult(
        variant_chrom=c, variant_pos=pos, ref=ref, alt=alt,
        splicing_merged=merged, splicing_sites=ss,
        splicing_site_usage=su, splicing_junctions=sj,
        rna_seq_log2fc=rna, cage_max=cage, atac_max=atac,
        dnase_max=dnase, chip_tf_max=chip_tf,
        chip_histone_max=chip_h, procap_max=procap,
        composite_avi=composite, top_modality=top_mod,
    )
