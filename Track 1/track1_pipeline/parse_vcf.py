"""Parse a gzipped VCF file using only stdlib."""
from __future__ import annotations
import gzip
from dataclasses import dataclass


@dataclass
class Variant:
    chrom: str
    pos: int
    id: str
    ref: str
    alt: str
    qual: float
    filter: str
    info: str
    gt: str
    ad: tuple[int, ...]
    dp: int
    gq: int

    @property
    def is_het(self) -> bool:
        return self.gt in ("0/1", "0|1", "1/0", "1|0")

    @property
    def is_hom_alt(self) -> bool:
        return self.gt in ("1/1", "1|1")

    @property
    def is_pass(self) -> bool:
        return self.filter == "PASS"

    @property
    def chrom_with_prefix(self) -> str:
        return self.chrom if self.chrom.startswith("chr") else f"chr{self.chrom}"


def parse_vcf(
    path: str,
    chrom: str | None = None,
    start: int | None = None,
    end: int | None = None,
    include_nonpass: bool = False,
    het_only: bool = True,
) -> list[Variant]:
    """Parse a gzipped VCF and return filtered variants."""
    variants: list[Variant] = []
    fmt_keys: list[str] = []

    with gzip.open(path, "rt") as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                continue

            cols = line.strip().split("\t")
            if len(cols) < 10:
                continue

            v_chrom = cols[0]
            v_pos = int(cols[1])

            # Region filter
            if chrom is not None and v_chrom != chrom:
                continue
            if start is not None and v_pos < start:
                continue
            if end is not None and v_pos > end:
                continue

            v_filter = cols[6]
            if not include_nonpass and v_filter != "PASS":
                continue

            # Parse FORMAT/sample columns - format keys are per-data-line in this VCF
            data_fmt_keys = cols[8].split(":")
            sample_fields = cols[9].split(":")
            sample = dict(zip(data_fmt_keys, sample_fields))

            gt = sample.get("GT", "./.")
            if het_only and gt not in ("0/1", "0|1", "1/0", "1|0"):
                continue

            ad_str = sample.get("AD", "0,0")
            ad = tuple(int(x) for x in ad_str.split(","))
            dp = int(sample.get("DP", "0"))
            gq = int(sample.get("GQ", "0"))

            variants.append(Variant(
                chrom=v_chrom, pos=v_pos, id=cols[2], ref=cols[3], alt=cols[4],
                qual=float(cols[5]) if cols[5] != "." else 0.0,
                filter=v_filter, info=cols[7], gt=gt, ad=ad, dp=dp, gq=gq,
            ))

    return variants
