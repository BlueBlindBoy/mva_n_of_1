# Track 2: Drug Repurposing (reproduction kit)

**Challenge:** Rare Disease, Real Kid: The MVA Hackathon 2026, Track 2.

This folder holds the computational scripts and the steps to reproduce the
run. The final Track 2 report (challenge deliverable) is in this folder as
`bit_by_bit_track2_report.pdf`, and the submitted methods description form
is in the repo root (`Methods_Description_Form_Track1_Track2.pdf`). The
narrative framework and audit documents, along with any output containing
the proband's variant data, are kept locally under the hackathon data-use
agreement and are out of this repository on purpose.

## Contents

| Path | What it is |
|---|---|
| `run_pangolin_bub1b.py` | Local Pangolin splice-site scoring (CPU). |
| `requirements.txt` | Python dependencies for the Pangolin step. |
| `variants.example.json` | Format of the local-only `variants.json` config (coordinates, kept out of git). |
| `drug_search/` | Approved-drug repurposing screen and verdict table (section 3). |

## 1. Pangolin splice-site scoring

```bash
pip install -r requirements.txt   # numpy, torch, pangolin
python run_pangolin_bub1b.py
```

- Variants to score are read from `variants.json` (format:
  `variants.example.json`). That file is local-only and gitignored, because
  it holds the proband's variant coordinates, which must not be shared.
- The `<name>_ref.fa` files are 10001 bp windows of public GRCh38 sequence
  centered on each variant. The script reads the sequence from line 2 of
  each FASTA. Regenerate them from UCSC or Ensembl using the coordinates in
  your local `variants.json`. They are public reference data, but they
  encode the variant positions, so they stay local.
- For each variant the script scores P(splice) over a 10001 bp window
  centered on the variant (public GRCh38 reference sequence, supplied as
  `<name>_ref.fa` alongside) and reports the max gain/loss within ±50 bp.
  CPU-only.
- The output `pangolin_results.json` contains scores only (variant name +
  tissue scores, no coordinates) and is gitignored. The submitted copy was
  produced the same way.

## 2. Structure predictions

See `../track2-structures/README.md`. Wild-type reference models (kinase
domain with BUB3/CDC20 complexes; full-length protein) were generated with
JapanFold (Boltz-2), ESMFold-2 and Protenix from the wild-type FASTA inputs
kept there. Models of the proband's variant allele were generated the same
way and are kept locally only.

## 3. Approved-drug repurposing screen

See `drug_search/README.md`. The reproducible 22-gene approved-drug screen
(`drug_search/run_approved_screen.py`, ChEMBL-based) and the
model-to-approved-drug verdict table
(`drug_search/Model_to_Approved_Drugs.csv`) live in `Track 2/drug_search/`.

## 4. Henze SAC model and attachment-error module (chromosome mis-segregation)

The disease model behind the drug repurposing analysis. It is a 1:1
reconstruction of the published spindle-assembly-checkpoint network (Henze,
Dittrich and Ibrahim 2017, Scientific Reports 7, 3865; BioModels
MODEL1812210002, original SBML deposit included), extended with a BUBR1
deficiency perturbation layer and a kinetochore attachment-error module
that predicts chromosome mis-segregation risk. Every kinetic value is read
from the model deposit (`models/henze2017_original/henze_model.json`). No
parameters were invented. The perturbation values are literature-informed
scenarios, not patient measurements.

### Files

| Path | What it is |
|---|---|
| `models/henze2017_original/henze2017.py` | The Henze 2017 network (ODE integration, conservation checks). |
| `models/henze2017_original/henze_model.json` | Machine-readable parameter extraction from the SBML deposit. |
| `models/henze2017_original/Henze2017_original_SBML.xml` | The original BioModels deposit, kept for provenance. |
| `models/henze2017_original/test_structural.py` | Self-check of the network structure. |
| `models/henze2017_bubr1_extension/bubr1_perturbation.py` | BUBR1 pool scaling plus activity scenarios (A1 MCC assembly, A2 BCC assembly, A3 APC/C inhibition) in cis and trans configurations. |
| `simulations/attachment_module/attachment_module.py` | Kinetochore attachment module (unattached, merotelic, correct states; 46 chromosome pairs). |
| `simulations/shared_eval.py` | Shared simulation and scoring source used by the run scripts. |
| `simulations/attachment_module/run_*.py` | Run scripts: depletion panel, chromosome sensitivity, intervention screen, setting robustness, validation. |

### Reproduce

```bash
pip install -r requirements.txt
python models/henze2017_original/test_structural.py
python simulations/attachment_module/run_depletion_panel.py
```

`run_depletion_panel.py` must run first, because it writes the
`Attachment_Depletion_Panel.csv` that the validation and sensitivity
scripts read. The remaining scripts then run independently:

```bash
python simulations/attachment_module/run_chromosome_sensitivity.py
python simulations/attachment_module/run_intervention_screen.py
python simulations/attachment_module/run_setting_robustness.py
python simulations/attachment_module/run_validation.py
```

All outputs are written to `Track 2/results/` (gitignored, kept on the
local machine). `run_validation.py` checks the module against the BUBR1
depletion observations of Suijkerbuijk et al. (2010) in HeLa cells.

### Notes

- The BUBR1 perturbation layer is parameterized by a pool fraction and an
  activity scenario. It contains no patient variant data. The two patient
  variants appear only as de-identified cis/trans configurations.
- The pinned run outputs (depletion panel, validation summary, atlas
  tables) are kept locally. The drug-feasibility atlas built on these
  scripts lives in the main project repository, not here.

## Privacy

This repository contains code and steps. It has no proband variant data of
any kind. See the repo-root `README.md` and `.gitignore`.
