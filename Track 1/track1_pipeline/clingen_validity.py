"""ClinGen gene-disease validity classifications for MVA-associated genes."""
from __future__ import annotations

CLINGEN_VALIDITY_SCORES: dict[str, float] = {
    "Definitive": 1.0, "Strong": 0.8, "Moderate": 0.6,
    "Limited": 0.4, "Disputed": 0.1, "No classification": 0.0,
}

MVA_GENE_CLINGEN: dict[str, dict] = {
    "BUB1B": {"classification": "Definitive", "disease": "Mosaic variegated aneuploidy syndrome 1",
              "url": "https://search.clinicalgenome.org/kb/genes/HGNC:1149"},
    "CEP57": {"classification": "Definitive", "disease": "Mosaic variegated aneuploidy syndrome 2",
              "url": "https://search.clinicalgenome.org/kb/genes/HGNC:30794"},
    "TRIP13": {"classification": "Strong", "disease": "Mosaic variegated aneuploidy syndrome 3",
               "url": "https://search.clinicalgenome.org/kb/genes/HGNC:12307"},
}


def get_clingen_validity(gene: str) -> dict:
    if gene in MVA_GENE_CLINGEN:
        return MVA_GENE_CLINGEN[gene]
    return {"classification": None, "disease": None, "url": None}


def get_clingen_score(gene: str) -> float:
    result = get_clingen_validity(gene)
    cls = result["classification"]
    if cls is None:
        return CLINGEN_VALIDITY_SCORES["No classification"]
    return CLINGEN_VALIDITY_SCORES.get(cls, 0.0)
