"""Map HPO phenotype terms to candidate disease genes via the JAX HPO REST API."""
from __future__ import annotations
import requests
from collections import Counter

HPO_TERMS = [
    "HP:0002859",  # Rhabdomyosarcoma
    "HP:0000121",  # Nephrocalcinosis
    "HP:0004322",  # Short stature
    "HP:0001508",  # Failure to thrive
    "HP:0003202",  # Skeletal muscle atrophy
    "HP:0001622",  # Premature birth
    "HP:0001518",  # Small for gestational age
    "HP:0200067",  # Recurrent spontaneous abortion
]

HPO_API_BASE = "https://api.hpo.jax.org/api/hpo"


def get_genes_for_hpo_term(hpo_id: str, timeout: int = 30) -> set[str]:
    url = f"{HPO_API_BASE}/term/{hpo_id}/genes"
    resp = requests.get(url, timeout=timeout)
    resp.raise_for_status()
    data = resp.json()
    genes: set[str] = set()
    for g in data.get("genes", []):
        symbol = g.get("gene_symbol") or g.get("gene")
        if symbol:
            genes.add(symbol)
    return genes


def get_candidate_genes(hpo_terms: list[str], min_hits: int = 1) -> list[tuple[str, int]]:
    gene_counter: Counter[str] = Counter()
    for hpo_id in hpo_terms:
        try:
            genes = get_genes_for_hpo_term(hpo_id)
            for g in genes:
                gene_counter[g] += 1
        except Exception as e:
            print(f"Warning: failed to query HPO for {hpo_id}: {e}")
    return [(gene, count) for gene, count in gene_counter.most_common() if count >= min_hits]
