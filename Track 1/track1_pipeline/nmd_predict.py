"""Predict nonsense-mediated decay (NMD) for LoF variants."""
from __future__ import annotations
import requests
from dataclasses import dataclass

ENSEMBL_OVERLAP_URL = "https://rest.ensembl.org/overlap/id/{transcript_id}"
NMD_THRESHOLD_BP = 50


@dataclass
class NMDResult:
    nmd_likely: bool | None
    exon_rank: int | None
    total_exons: int | None
    distance_to_last_junction: int | None


def predict_nmd(variant_pos: int, transcript_id: str, timeout: int = 30) -> NMDResult:
    url = ENSEMBL_OVERLAP_URL.format(transcript_id=transcript_id)
    try:
        resp = requests.get(url, params={"feature": "exon"},
                            headers={"Accept": "application/json"}, timeout=timeout)
        if resp.status_code != 200:
            return NMDResult(nmd_likely=None, exon_rank=None, total_exons=None,
                             distance_to_last_junction=None)
        exons_raw = resp.json()
    except Exception:
        return NMDResult(nmd_likely=None, exon_rank=None, total_exons=None,
                         distance_to_last_junction=None)

    seen = set()
    exons = []
    for e in exons_raw:
        key = (e["start"], e["end"])
        if key not in seen:
            seen.add(key)
            exons.append(e)
    exons.sort(key=lambda e: e["start"])
    total_exons = len(exons)
    if total_exons == 0:
        return NMDResult(nmd_likely=None, exon_rank=None, total_exons=0,
                         distance_to_last_junction=None)

    variant_exon_idx = None
    for i, e in enumerate(exons):
        if e["start"] <= variant_pos <= e["end"]:
            variant_exon_idx = i
            break

    if variant_exon_idx is None:
        return NMDResult(nmd_likely=None, exon_rank=None, total_exons=total_exons,
                         distance_to_last_junction=None)

    exon_rank = variant_exon_idx + 1
    if variant_exon_idx == total_exons - 1:
        return NMDResult(nmd_likely=False, exon_rank=exon_rank,
                         total_exons=total_exons, distance_to_last_junction=0)

    last_exon_start = exons[-1]["start"]
    distance = last_exon_start - variant_pos
    if distance < NMD_THRESHOLD_BP:
        return NMDResult(nmd_likely=False, exon_rank=exon_rank,
                         total_exons=total_exons, distance_to_last_junction=distance)
    return NMDResult(nmd_likely=True, exon_rank=exon_rank,
                     total_exons=total_exons, distance_to_last_junction=distance)
