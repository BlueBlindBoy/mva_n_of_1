#!/usr/bin/env python3
"""Intervention-ranking robustness across the settings that pass both
held-out arms (review response, 2026-10-23).

Question: do the 3 of 54 settings that pass both held-out arms
(13%/6% within 3x) agree on the intervention ranking, and does the
severe depletion scenario S6 (BUBR1 pool at 6% of WT, where the
segregation-risk readout is well above its floor) change the ranking or
the required restoration level?

Scenarios: D1 (representative, SE.base_scenario), D3 (abundance 0.2),
S6 (abundance 0.06 - the 6% arm of Suijkerbuijk et al. 2010).
Settings: w0.5_g0.01_k1_e0.02, w0.5_g0.01_k1_e0.05 (best),
w0.5_g0.005_k1_e0.01 (the three with heldout_arms_within_3x == 2).
Interventions: baseline, BUBR1 x2.0, BUBR1 x3.231, O_Mad2 x2,
Cdc20 x0.25.
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

import attachment_module as AM
import shared_eval as SE

FLOOR = 0.02
SETTINGS = [
    {"w_a": 0.5, "gamma0": 0.01, "k": 1, "f_err": 0.02},
    {"w_a": 0.5, "gamma0": 0.01, "k": 1, "f_err": 0.05},
    {"w_a": 0.5, "gamma0": 0.005, "k": 1, "f_err": 0.01},
]
SET_LABELS = ["w0.5_g0.01_k1_e0.02", "w0.5_g0.01_k1_e0.05",
              "w0.5_g0.005_k1_e0.01"]

INTS = [
    ("baseline", None),
    ("BubR1_x2", 2.0),
    ("BubR1_x3.231", 3.231),
    ("OMad2_x2", None),
    ("Cdc20_x0.25", None),
]


def pools():
    p1, y01, _ = SE.base_scenario()
    p3, y03, _ = AM.B.abundance_only(0.2)
    p6, y06, _ = AM.B.abundance_only(0.06)
    out = []
    for lab, p, y0 in (("D1", p1, y01), ("D3", p3, y03), ("S6", p6, y06)):
        mu = float(y0[AM.H.IDX["BubR1_Bub3"]]) / AM.B.POOL_WT
        out.append((lab, p, y0, mu))
    return out


def main():
    wt = SE.wt_ref()
    rows = []
    for lab, p, y0, mu in pools():
        for mech, slabel in zip(SETTINGS, SET_LABELS):
            for iv, mult in INTS:
                y0i = y0.copy()
                mu_i = mu
                if mult is not None:
                    y0i[AM.H.IDX["BubR1_Bub3"]] *= mult
                    mu_i = mu * mult
                if iv == "OMad2_x2":
                    y0i[AM.H.IDX["O_Mad2"]] *= 2.0
                if iv == "Cdc20_x0.25":
                    y0i[AM.H.IDX["Cdc20"]] *= 0.25
                res = AM.simulate_extended(mu_i, mech, p=p, y0_19=y0i)
                eo = AM.error_output(res, floor=FLOOR)
                cr = AM.checkpoint_readouts(res)
                r = {"leak": cr["leak_uM"], "mcc_peak": cr["mcc_peak_uM"],
                     "t_commit": cr["t_commit_s"]}
                S, arrest, _ = SE.score(r, wt)
                rows.append({
                    "scenario": lab, "setting": slabel, "intervention": iv,
                    "S_checkpoint": round(S, 4), "arrest": arrest,
                    "P_error": None if eo["P_error"] is None else round(eo["P_error"], 5),
                    "P_pred": None if eo["P_pred"] is None else round(eo["P_pred"], 5),
                })
    out_csv = os.path.join(ROOT, "results",
                           "Attachment_Screen_Setting_Robustness.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out_csv} ({len(rows)} rows)")
    for r in rows:
        if r["scenario"] == "S6":
            print(r)


if __name__ == "__main__":
    main()
