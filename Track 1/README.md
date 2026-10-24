# Track 1: Variant Prediction (reproduction kit)

**Challenge:** Rare Disease, Real Kid: The MVA Hackathon 2026, Track 1.

This folder holds the full pipeline code and the steps to reproduce the run.
The submission file `bit_by_bit_hpo-clingen-nmd.csv` is deliberately kept
out of this repository, because it contains the proband's variant data. It
stays on the local machine and was submitted to the challenge site
separately. The long-form methods preprint is out of this repository for
the same reason. The final Track 1 report (the challenge deliverable, which
names the two BUB1B variants) and the methods description form (repo root)
are in the repository.

## Contents

| Path | What it is |
|---|---|
| `track1_pipeline/` | The full pipeline (10 Python modules): VCF parsing, HPO mapping, ClinGen gene-disease validity, VEP annotation, NMD prediction, optional AlphaGenome scoring, ranking, CSV generation. |
| `bit_by_bit_track1_report.pdf` | The final Track 1 report (challenge deliverable). |
| `requirements.txt` | Python dependencies. |

## Reproduce the run

```bash
pip install -r requirements.txt
python track1_pipeline/run_pipeline.py \
    --vcf <path-to-gated-WGS-VCF> \
    --output bit_by_bit_hpo-clingen-nmd.csv \
    --proband-id <PROBAND_ID> \
    --epcr 0.85
```

- The output filename follows the challenge convention
  `<username>_<model-or-approach>.csv` (username `bit_by_bit`; the approach
  slug covers the HPO gene mapping, ClinGen validity, and NMD-aware
  prioritization core).

- The VCF is the gated SageBio challenge dataset (request access via the
  Hugging Face dataset page). It is not redistributable and it is not in
  this repo.
- All other inputs come from public, documented APIs (Ensembl VEP, HPO JAX,
  ClinGen KB, optional AlphaGenome). No other private data is needed.
- The output is a ranked variant CSV. The pinned run's output (one primary
  row: a compound-heterozygous pair in the causal MVA gene, epcr 0.85) is
  kept locally and was submitted to the challenge site. It is not committed
  here, because it contains the proband's variant data.

## Privacy

The final report and the methods description form (repo root) name the
proband's two BUB1B variants; both are challenge deliverables. The
submission CSV is kept out of the repository for the same reason. The
pipeline code and every other file in this folder contain no proband
variant data. See the repo-root `README.md` and `.gitignore`.
