#!/usr/bin/env python3
"""Chromosome-count sensitivity of the segregation-risk readout (review
response, 2026-10-23).

HeLa is aneuploid (reported 76-80 chromosomes; Theunissen et al. 2016),
while the module simulates 46 sister-chromosome pairs (the diploid
convention of the published 92-kinetochore network). The readout is a
binomial in the per-pair error probability p:
    P(N) = 1 - (1 - p)^N
    P_pred(N) = floor + P(N),  floor = 0.02
The empirical floor is the measured WT mis-segregation level in HeLa, so
it is not rescaled with N.

This script recomputes the readout for N in {46, 50, 60, 70, 80} from the
existing panel (p_pair per arm) - no new simulations. A full coupled
re-run at a different kinetochore count is not defined within the
published network (calibrated at 92 kinetochores; the switch threshold
91.98 is a model constant), so this binomial readout sensitivity at fixed
per-pair error is the closest reproducible implementation.
"""
from __future__ import annotations

import csv
import os

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
PANEL = os.path.join(ROOT, "results", "Attachment_Depletion_Panel.csv")
OUT = os.path.join(ROOT, "results", "Attachment_Chromosome_Sensitivity.csv")

OBS = {0.13: 0.03, 0.06: 0.50}
FLOOR = 0.02
N_LIST = [46, 50, 60, 70, 80]
BEST = "w0.5_g0.01_k1_e0.05"
PASS = (1.0 / 3.0, 3.0)

rows = list(csv.DictReader(open(PANEL)))
for r in rows:
    r["mu"] = float(r["mu"])
    r["p_pair"] = None if r["p_pair"] in ("", "None") else float(r["p_pair"])

# per setting: p_pair per held-out arm
pp = {}
for r in rows:
    if r["mu"] in (0.13, 0.06):
        pp.setdefault(r["setting"], {})[f"p_{r['mu']:g}"] = r["p_pair"]


def perr(p, N):
    return None if p is None else 1.0 - (1.0 - p) ** N


def pred_at(p, N):
    pe = perr(p, N)
    if pe is None:
        return None
    return min(1.0, FLOOR + pe)


# --- best setting, per N ----------------------------------------------------
print(f"best setting {BEST} (floor {FLOOR}, phi 1):")
best_rows = []
for N in N_LIST:
    rec: dict = {"N": N}
    for mu in (0.13, 0.06):
        p = pp[BEST].get(f"p_{mu:g}")
        pv = pred_at(p, N)
        ratio = None if pv is None else pv / OBS[mu]
        rec[mu] = (pv, ratio)
    best_rows.append(rec)
    p13, r13 = rec[0.13]
    p6, r6 = rec[0.06]
    print(f"  N={N}: 13% P_pred={p13:.4f} (ratio {r13:.2f}) | "
          f"6% P_pred={p6:.4f} (ratio {r6:.2f})")

# --- ensemble: settings passing both held-out arms within 3x, per N ---------
ens = {}
print("\nensemble (54 settings), both held-out arms within 3x:")
for N in N_LIST:
    n = 0
    for s, d in pp.items():
        ok = True
        for mu in (0.13, 0.06):
            pv = pred_at(d.get(f"p_{mu:g}"), N)  # None (arrest) -> None
            if pv is None:
                ok = False
                continue
            if not (PASS[0] <= pv / OBS[mu] <= PASS[1]):
                ok = False
                break
        if ok:
            n += 1
    ens[N] = n
    print(f"  N={N}: {n}/54 settings pass")

with open(OUT, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["N_pairs", f"P_pred_13pct_{BEST}", "ratio_13pct",
                f"P_pred_6pct_{BEST}", "ratio_6pct",
                "ensemble_both_heldout_arms_within_3x"])
    for rec in best_rows:
        (p13, r13), (p6, r6) = rec[0.13], rec[0.06]
        w.writerow([rec["N"], round(p13, 4), round(r13, 3),
                    round(p6, 4), round(r6, 3), ens[rec["N"]]])
print(f"\nwrote {OUT}")
