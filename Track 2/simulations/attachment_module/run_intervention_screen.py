#!/usr/bin/env python3
"""Dual-outcome intervention screen (item 10, compact version).

Question: with the attachment-error output added, does restoring BUBR1
abundance remain the leading intervention for BOTH checkpoint signaling
(canonical score S) and predicted segregation-error risk (P_pred)?

Interventions (same levers as the canonical atlas, restricted to the
small set needed for the item-10 question):
  - BUBR1 pool x1 (disease baseline), x2.0, x3.231 (canonical
    partial / near-reference rescue multipliers for D1)
  - O_Mad2 x2 (checkpoint signaling lever)
  - Cdc20 x0.25 (harmful-control lever: over-restrains APC/C)
  - BUBR1 pool x3.231 WITH attachment module off (checkpoint-only
    comparison)

Run for D1 (patient_trans a=0.5 gamma=0.5 A1) and D3 (abundance 0.2)
with the best-validated attachment setting
(w_a=0.5, gamma0=0.01, k=1, f_err=0.05).
"""
from __future__ import annotations

import csv
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "models", "henze2017_bubr1_extension"))
sys.path.insert(0, os.path.join(ROOT, "models", "henze2017_original"))
sys.path.insert(0, os.path.join(ROOT, "simulations"))
sys.path.insert(0, HERE)

import numpy as np
import attachment_module as AM
import shared_eval as SE

BEST = {"w_a": 0.5, "gamma0": 0.01, "k": 1, "f_err": 0.05}
FLOOR = 0.02


def disease_pools():
    """(label, params, y0_19, mu_equivalent) for D1, D3, and S6.

    mu_equivalent: fraction of the WT pool (0.22 uM) represented by
    the scenario's BUBR1 pool, used to drive the attachment module's
    BUBR1-dependent rates (the same pool-conserving logic as
    bubr1_perturbation).

    S6: the severe depletion scenario, BUBR1 pool at 6% of WT - the 6%
    arm of Suijkerbuijk et al. (2010) - where the segregation-risk
    readout is well above its floor (review response, 2026-10-23).
    """
    p1, y01, _ = SE.base_scenario()                 # D1
    p3, y03, _ = AM.B.abundance_only(0.2)           # D3
    p6, y06, _ = AM.B.abundance_only(0.06)          # S6
    mu1 = float(y01[AM.H.IDX["BubR1_Bub3"]]) / AM.B.POOL_WT
    mu3 = float(y03[AM.H.IDX["BubR1_Bub3"]]) / AM.B.POOL_WT
    mu6 = float(y06[AM.H.IDX["BubR1_Bub3"]]) / AM.B.POOL_WT
    return [("D1", p1, y01, mu1), ("D3", p3, y03, mu3), ("S6", p6, y06, mu6)]


def main():
    wt = SE.wt_ref()
    rows = []
    interventions = [
        ("baseline", lambda pi, y0i: None),
        ("BubR1_x2", None),
        ("BubR1_x3.231", None),
        ("BubR1_x3.231_module_off", None),
        ("OMad2_x2", lambda pi, y0i: y0i.__setitem__(AM.H.IDX["O_Mad2"],
                                                       y0i[AM.H.IDX["O_Mad2"]] * 2.0)),
        ("Cdc20_x0.25", lambda pi, y0i: y0i.__setitem__(AM.H.IDX["Cdc20"],
                                                         y0i[AM.H.IDX["Cdc20"]] * 0.25)),
    ]
    for label, p, y0, mu in disease_pools():
        for iv, fn in interventions:
            pi = dict(p)
            y0i = y0.copy()
            mu_i = mu
            if iv.startswith("BubR1"):
                mult = float(iv.split("_x")[1].split("_")[0])
                y0i[AM.H.IDX["BubR1_Bub3"]] *= mult
                mu_i = mu * mult
            if fn is not None:
                fn(pi, y0i)
            if iv == "BubR1_x3.231_module_off":
                res = AM.simulate_extended(
                    mu_i, {"w_a": 0.0, "gamma0": 0.005, "k": 1,
                           "f_err": 0.0},
                    p=pi, y0_19=y0i)
                eo = None
            else:
                res = AM.simulate_extended(mu_i, BEST, p=pi, y0_19=y0i)
                eo = AM.error_output(res, floor=FLOOR)
            cr = AM.checkpoint_readouts(res)
            # canonical checkpoint score needs leak/mcc/t_commit vs WT
            r = {"leak": cr["leak_uM"], "mcc_peak": cr["mcc_peak_uM"],
                 "t_commit": cr["t_commit_s"]}
            S, arrest, parts = SE.score(r, wt)
            rows.append({
                "scenario": label, "intervention": iv,
                "pool_uM": round(float(y0i[AM.H.IDX["BubR1_Bub3"]]), 4),
                "S_checkpoint": round(S, 4),
                "arrest": arrest,
                "t_commit_s": round(cr["t_commit_s"], 1)
                             if cr["t_commit_s"] else None,
                "mcc_peak_uM": round(cr["mcc_peak_uM"], 5),
                "leak_uM": round(cr["leak_uM"], 5),
                "P_error": None if eo is None or eo["P_error"] is None
                          else round(eo["P_error"], 5),
                "P_pred": None if eo is None or eo["P_pred"] is None
                        else round(eo["P_pred"], 5),
            })

    out_csv = os.path.join(ROOT, "results",
                           "Attachment_Intervention_Screen.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out_csv}")
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
