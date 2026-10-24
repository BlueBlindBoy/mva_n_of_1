# drug_search: Approved-drug repurposing screen (reproduction kit)

**Challenge:** Rare Disease, Real Kid: The MVA Hackathon 2026, Track 2
(drug repurposing for the MVA disease model).

This folder holds the reproducible approved-drug screen and the
model-to-approved-drug verdict table. The long-form curated report (FDA and
EMA label details with primary citations) is kept locally under the
hackathon data-use agreement and is out of this repository on purpose.

## Contents

| Path | What it is |
|---|---|
| `run_approved_screen.py` | Reproducible 22-gene approved-drug screen (ChEMBL REST API; stdlib only, CPU, no API key). |
| `approved_screen_results.csv` | Pinned full-run output: one row per (panel gene, target, drug-mechanism) record with phase flags. |
| `approved_screen_rows.json` | Same rows plus the resolved target map (for auditing). |
| `Model_to_Approved_Drugs.csv` | Verdict table: model intervention lever to candidate approved drug to verdict (scrubbed; see Privacy). |

## 1. The 22-gene approved-drug screen (automated)

```bash
python run_approved_screen.py --out approved_screen_results.csv --json approved_screen_rows.json
```

- No dependencies beyond the Python standard library.
- For each gene in the 22-gene mitotic-checkpoint and spindle-apparatus
  panel (BUB1B, BUB3, MAD2L1, TRIP13, ANAPC1/2, APC11, CDC20, CDC25C,
  UBE2C, PDS5A/B, PTTG1, CDK1, BUB1, AURKA/B, PLK1, KIF11, TUBB, CEP57,
  CDCA7), the script resolves the ChEMBL target(s). It then pulls every
  documented drug mechanism against them and flags records with
  `max_phase >= 4` (approved in at least one jurisdiction).
- Non-approved rows are kept in the output, so the per-target "no approved
  modulator" conclusions are auditable from the data, not just the summary.
- The pinned outputs in this folder are from a run on 2026-10-24 against
  ChEMBL 37. Re-running against a future ChEMBL release will change the rows
  (database drift). Review the diff before re-committing.
- The original 2026-10-23 screen used DGIdb and Open Targets (which
  aggregate ChEMBL) plus manual FDA/EMA label curation. Those API endpoints
  are not resolvable from the team network, so the reproducible core runs on
  ChEMBL directly. The conclusion is unchanged: all approved direct hits are
  TUBB (beta-tubulin) microtubule agents (taxanes, epothilone, eribulin,
  vinca alkaloids, colchicine, maytansinoid/auristatin ADC payloads), and no
  approved drug targets any of the other 21 panel checkpoint proteins.
- Off-target interpretation: a drug mapped to a panel gene may carry a
  different primary mechanism in the `mechanism_of_action` column (e.g. a
  multi-kinase inhibitor). The curated layer counts a drug as an approved
  *direct* modulator of a panel target only when that target is the drug's
  approved mechanism. ADCs are judged on the payload's direct target.

## 2. Model-to-approved-drug verdict table (curated)

`Model_to_Approved_Drugs.csv` maps each intervention lever of the disease
model (BubR1:Bub3 pool restoration, nonsense-allele readthrough, open-MAD2
pool enlargement, MCC assembly, APC/C:Cdc20 association, Cdc20 abundance)
to candidate approved medications with a verdict
(conditional / quantitatively insufficient / negative / excluded). It is a
curated table built from the model evidence and the screen results, with
supporting literature. It is not script output. Its `evidence_from_atlas`
column cites the drug-feasibility atlas, which lives in the main project
repository (outside this track-deliverables repo).

## Privacy

This folder contains code and curated results. It has no proband variant
data of any kind. The verdict table is scrubbed: proband variant
nomenclature is replaced with generic "proband's missense / nonsense
allele" labels. See the repo-root `README.md` and `.gitignore`.
