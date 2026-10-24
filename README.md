# mva_n_of_1

Track 1 and Track 2 deliverables for the MVA Hackathon 2026
("Rare Disease, Real Kid").

## Privacy rule

This repository holds code, reproduction steps, the methods description
form, and the final track reports. The methods form and the two reports are
challenge deliverables, and they name the proband's two BUB1B variants.
The Track 1 submission CSV is deliberately kept out of this repository,
because it contains the proband's variant rows; it stays on the local
machine and was submitted through the challenge site. Per the challenge
rules the repository is private during the hackathon and must be public
after it ends. Everything else proband-related stays on the local machine
under the hackathon data-use agreement.

## Layout

| Path | What it is |
|---|---|
| `Track 1/` | Variant-prediction pipeline (10 Python modules) + the final Track 1 report + requirements + README with reproduction steps. |
| `Track 2/` | Pangolin splice-scoring script, the drug repurposing screen in `drug_search/`, the Henze SAC model and attachment-error module, the final Track 2 report, requirements, and README with reproduction steps. |
| `track2-structures/` | Wild-type reference structure predictions (Boltz-2, ESMFold-2, Protenix) + README documenting which patient-variant models were kept out. |
| `Methods_Description_Form_Track1_Track2.pdf` | The submitted methods description form, covering both tracks. |

Each folder's README has the exact steps to reproduce its part of the run.
