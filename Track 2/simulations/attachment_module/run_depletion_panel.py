#!/usr/bin/env python3
"""Depletion panel for the attachment-error extension (items 5, 8).

Runs the integrated Henze+attachment model at the Suijkerbuijk et al.
(2010) residual-BUBR1 levels (100%, 13%, 6%; HeLa doxycycline-inducible
shRNA) across a small bounded mechanism grid, plus the item-5 scenario
comparison (checkpoint-only / attachment-only / combined deficiency).

Mechanism grid (all values literature-derived or bounded assumptions;
see docs/Attachment_Error_Extension.md for the parameter-evidence table):
  w_a   in {0.5, 1.0}      BUBR1 contribution to congression speed
  gamma0 in {0.002,0.005,0.01} s^-1  WT error-correction rate
  k     in {1, 2, 3}       order of BUBR1 dependence of correction
  f_err in {0.01, 0.02, 0.05}        merotelic fraction of captures

Outputs
-------
results/Attachment_Depletion_Panel.csv   one row per (setting, mu)
results/Attachment_Scenario_Comparison.csv  item-5 scenarios
"""
from __future__ import annotations

import csv
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "models", "henze2017_bubr1_extension"))
sys.path.insert(0, os.path.join(ROOT, "models", "henze2017_original"))
sys.path.insert(0, HERE)

import numpy as np
import attachment_module as AM

MU_GRID = [1.0, 0.13, 0.06]          # Suijkerbuijk residual BUBR1
OBS = {1.0: 0.02, 0.13: 0.03, 0.06: 0.50}   # observed segregation defects
FLOOR = 0.02                           # empirical WT background (100% arm)

W_A_GRID = [0.5, 1.0]
GAMMA0_GRID = [0.002, 0.005, 0.01]
K_GRID = [1, 2, 3]
F_ERR_GRID = [0.01, 0.02, 0.05]


def main():
    t0 = time.time()
    rows = []
    n = 0
    for w_a in W_A_GRID:
        for gamma0 in GAMMA0_GRID:
            for k in K_GRID:
                for f_err in F_ERR_GRID:
                    sid = (f"w{w_a:g}_g{gamma0:g}_k{k}_e{f_err:g}")
                    mech = {"w_a": w_a, "gamma0": gamma0, "k": k,
                            "f_err": f_err}
                    for mu in MU_GRID:
                        res = AM.simulate_extended(mu, mech)
                        eo = AM.error_output(res, floor=FLOOR)
                        cr = AM.checkpoint_readouts(res)
                        pred = eo["P_pred"]
                        obs = OBS[mu]
                        ratio = (pred / obs) if (pred is not None and obs > 0) else None
                        rows.append({
                            "setting": sid, "w_a": w_a, "gamma0": gamma0,
                            "k": k, "f_err": f_err, "mu": mu,
                            "pool_uM": round(AM.B.POOL_WT * mu, 4),
                            "t_commit_s": None if eo["t_commit_s"] is None
                                            else round(eo["t_commit_s"], 1),
                            "arrest": eo["arrest"],
                            "kinu_at_commit": None if eo["kinu_at_commit"] is None
                                            else round(eo["kinu_at_commit"], 6),
                            "e_at_commit_pairs": None if eo["e_at_commit_pairs"] is None
                                            else round(eo["e_at_commit_pairs"], 4),
                            "p_pair": None if eo["p_pair"] is None
                                        else round(eo["p_pair"], 6),
                            "P_error": None if eo["P_error"] is None
                                      else round(eo["P_error"], 5),
                            "P_pred": None if pred is None else round(pred, 5),
                            "obs_missegregation": obs,
                            "pred_over_obs": None if ratio is None else round(ratio, 3),
                            "mcc_peak_uM": round(cr["mcc_peak_uM"], 5),
                            "leak_uM": round(cr["leak_uM"], 5),
                            "securin_frac_120min": round(cr["securin_frac_120min"], 5),
                        })
                        n += 1
                    print(f"done {sid}  ({n}/{len(W_A_GRID)*len(GAMMA0_GRID)*len(K_GRID)*len(F_ERR_GRID)*len(MU_GRID)} runs, {time.time()-t0:.0f}s)", flush=True)

    os.makedirs(os.path.join(ROOT, "results"), exist_ok=True)
    panel_csv = os.path.join(ROOT, "results", "Attachment_Depletion_Panel.csv")
    with open(panel_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {panel_csv} ({len(rows)} rows, {time.time()-t0:.0f}s)")

    # ---- item-5 scenario comparison at the 6% and 13% arms ---------------
    base = {"w_a": 0.5, "gamma0": 0.005, "k": 1, "f_err": 0.02}
    srows = []
    for mu in (0.13, 0.06):
        # checkpoint-only: original network + BUBR1 pool perturbation,
        # attachment module off (f_err=0, w_a=0) -> cannot represent
        # attachment errors at all
        res_co = AM.simulate_extended(mu, {"w_a": 0.0, "gamma0": 0.005,
                                           "k": 1, "f_err": 0.0})
        cr = AM.checkpoint_readouts(res_co)
        srows.append({
            "scenario": "checkpoint_only", "mu": mu,
            "t_commit_s": round(cr["t_commit_s"], 1) if cr["t_commit_s"] else None,
            "mcc_peak_uM": round(cr["mcc_peak_uM"], 5),
            "leak_uM": round(cr["leak_uM"], 5),
            "securin_frac_120min": round(cr["securin_frac_120min"], 5),
            "P_error": None, "P_pred": None,
            "note": "error output undefined: no attachment states in original network",
        })
        # attachment-only: WT MCC (pool intact), attachment module at mu
        p_wt = dict(AM.H.PARAMS)
        res_ao = AM.simulate_extended(mu, base, p=p_wt,
                                      y0_19=AM.H.Y0_VEC.copy())
        eo = AM.error_output(res_ao, floor=FLOOR)
        cr = AM.checkpoint_readouts(res_ao)
        srows.append({
            "scenario": "attachment_only", "mu": mu,
            "t_commit_s": round(eo["t_commit_s"], 1) if eo["t_commit_s"] else None,
            "mcc_peak_uM": round(cr["mcc_peak_uM"], 5),
            "leak_uM": round(cr["leak_uM"], 5),
            "securin_frac_120min": round(cr["securin_frac_120min"], 5),
            "P_error": None if eo["P_error"] is None else round(eo["P_error"], 5),
            "P_pred": None if eo["P_pred"] is None else round(eo["P_pred"], 5),
            "note": "WT MCC pool; module at mu",
        })
        # combined: BUBR1 pool perturbation + module at mu (main panel)
        res_cb = AM.simulate_extended(mu, base)
        eo = AM.error_output(res_cb, floor=FLOOR)
        cr = AM.checkpoint_readouts(res_cb)
        srows.append({
            "scenario": "combined", "mu": mu,
            "t_commit_s": round(eo["t_commit_s"], 1) if eo["t_commit_s"] else None,
            "mcc_peak_uM": round(cr["mcc_peak_uM"], 5),
            "leak_uM": round(cr["leak_uM"], 5),
            "securin_frac_120min": round(cr["securin_frac_120min"], 5),
            "P_error": None if eo["P_error"] is None else round(eo["P_error"], 5),
            "P_pred": None if eo["P_pred"] is None else round(eo["P_pred"], 5),
            "note": "abundance perturbation + module at mu (base setting)",
        })
    comp_csv = os.path.join(ROOT, "results", "Attachment_Scenario_Comparison.csv")
    with open(comp_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(srows[0].keys()))
        w.writeheader()
        w.writerows(srows)
    print(f"wrote {comp_csv}")

    # ---- limiting-case test (item 4) -------------------------------------
    lt = AM.limit_test(mu=0.13)
    print("limit test:", lt)
    assert lt["ok"], "extension does not reduce to the original model"


if __name__ == "__main__":
    main()
