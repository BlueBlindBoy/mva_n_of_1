# track2-structures

Wild-type reference structure predictions for the Track 2 drug-repurposing
analysis, plus the public reference structures and sequences used as inputs.

## What's here

- `BUBR1_AF_O60566.pdb`: public AlphaFold reference structure (UniProt O60566).
- `*_wt*.fa`, `bub3.fa`, `cdc20.fa`: wild-type and partner-protein sequences (public reference).
- `BUBR1_kinase_WT.cif`, `BUBR1_kinase_WT_results.json`: predicted kinase-domain model (Boltz-2 via JapanFold).
- `kinase_wt_bub3/`, `kinase_wt_cdc20/`, `kinase_wt_protenix/`: predicted kinase-domain plus partner complexes (JapanFold/Protenix).
- `bubr1_full_wt/`, `bubr1_full_wt_esmfold2_target_1.cif`, `bubr1_full_wt_esmfold2_results.json`: full-length protein models (JapanFold, ESMFold-2).
- `japanfold_jobs*.json`: job registries for the wild-type model runs.

## Patient-derived models (not in this repo)

Models predicted from the proband's variant allele are not in this
repository and must never be committed here. This includes its FASTA
sequences, the predicted structures, the per-variant analysis summaries, the
wild-type-vs-variant comparison files, the local-regulatory scan, and the
interactive structure viewer annotated with the proband's variant positions.
Those files are retained on the local machine under the hackathon data-use
agreement. The `.gitignore` in the repo root enforces this.
