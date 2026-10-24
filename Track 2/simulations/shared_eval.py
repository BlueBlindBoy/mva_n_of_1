"""Canonical evaluation of the Henze et al. (2017) network: ONE shared
simulation + scoring source for the intervention atlas, the rescue-threshold
search, and the depletion/combination scripts (reconciliation pass,
2026-10-23, review response item 1).

Why this module exists
----------------------
The threshold search originally ran on a 6 s output grid (1201 points)
while the atlas and combination scripts used a 2 s grid (3601 points).
The coarse grid undersampled the early MCC peak and the checkpoint-leak
transient, which changed the disease-baseline rescue score (0.640 vs
0.572 for D1) and shifted the computed BUBR1-pool rescue thresholds
(e.g. D1 partial rescue 1.765x on the coarse grid vs ~2.0x on the
canonical grid). All downstream scores, rescue thresholds, and pool
requirements now come from this single module.

Canonical choices
-----------------
T_END   = 7200 s
N_PTS   = 3601 (2 s grid) - required to resolve the early MCC peak
          (~5-8 s) and the leak transient
score   = (s1 + s2 + s3) / 3 with
    s1 = restraint:   max(0, 1 - (leak - leak_WT)/leak_WT)
         (may exceed 1 when a perturbation over-restricts APC/C; the
          >1 value is retained so over-suppression is visible in s1)
    s2 = signaling:   min(1, mcc_peak/mcc_peak_WT)
    s3 = completion:  1 if t_commit <= 1.5 x t_commit_WT else 0
    arrest = t_commit is None or t_commit > 1.5 x t_commit_WT
Model-defined targets (NOT clinical thresholds):
    S >= 0.75  model-predicted partial checkpoint recovery
    S >= 0.90  model-predicted near-reference checkpoint behavior
"""
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(ROOT, "models", "henze2017_bubr1_extension"))
sys.path.insert(0, os.path.join(ROOT, "models", "henze2017_original"))
import henze2017 as H            # noqa: E402
import bubr1_perturbation as B   # noqa: E402

T_END = 7200.0
N_PTS = 3601          # 2 s grid (canonical; see module docstring)
S_PARTIAL = 0.75      # model-predicted partial checkpoint recovery
S_FULL = 0.90         # model-predicted near-reference checkpoint behavior


def readouts(res):
    """Canonical readouts from a simulate() result."""
    ro = H.readouts(res)
    leak = B.checkpoint_leak_readout(res)["leak_max_uM"]
    y = res["y"]
    t = res["t"]
    kinu = y[:, H.IDX["KinU"]]
    tc = ro["t_anaphase_commitment_s"]
    kinu_at_commit = None
    if tc is not None:
        j = int(np.searchsorted(t, tc))
        kinu_at_commit = float(kinu[min(j, len(kinu) - 1)]) / 92.0
    below = np.where(kinu < 0.1)[0]
    t999 = float(t[below[0]]) if len(below) else None
    return {
        "leak": leak,
        "mcc_peak": ro["mcc_peak_uM"],
        "t_commit": tc,
        "securin_frac_120min": ro["securin_frac_at_30min"],  # value is at end of the 7200 s window (=120 min); model key name refers to the deposition's original window
        "apcc_cdc20_max": ro["apcc_cdc20_max_uM"],
        "kinu_at_commit": kinu_at_commit,
        "t_attach_999": t999,
    }


def simulate_readouts(y0=None, p=None, t_end=T_END, n_points=N_PTS):
    """Run the network and return canonical readouts."""
    if y0 is None:
        y0 = H.Y0_VEC.copy()
    else:
        y0 = y0.copy()
    if p is None:
        p = dict(H.PARAMS)
    else:
        p = dict(p)
    return readouts(H.simulate(t_end=t_end, y0=y0, p=p, n_points=n_points))


def score(r, ref):
    """Canonical multi-objective rescue score vs a reference (WT) readout.

    Returns (S, arrest, (s1, s2, s3)).
    """
    s1 = max(0.0, 1.0 - (r["leak"] - ref["leak"]) / ref["leak"])
    s2 = min(1.0, r["mcc_peak"] / ref["mcc_peak"])
    s3 = (1.0 if (r["t_commit"] is not None
                  and r["t_commit"] <= 1.5 * ref["t_commit"]) else 0.0)
    arrest = (r["t_commit"] is None
              or r["t_commit"] > 1.5 * ref["t_commit"])
    return (s1 + s2 + s3) / 3.0, arrest, (s1, s2, s3)


def wt_ref(k5_dissoc=None):
    """WT baseline readout on the canonical grid.

    k5_dissoc: optional override of the MCC-dissociation rate p["k_5"]
    (the documented 10x model-version discrepancy: 0.2 in the deposited
    SBML vs 0.02 in the paper's Table S1).
    """
    p = None
    if k5_dissoc is not None:
        p = dict(H.PARAMS)
        p["k_5"] = k5_dissoc
    return simulate_readouts(p=p)


def disease_scenarios():
    """D1 (trans a=0.5 gamma=0.5 A1), D2 (missense 0.25/0.25 A1),
    D3 (abundance 0.2, activity intact)."""
    d1 = base_scenario()
    p2, y02, lab2 = B.patient_trans(0.25, 0.25, "A1_mcc_assembly")
    p3, y03, lab3 = B.abundance_only(0.2)
    return [("D1",) + d1, ("D2",) + (p2, y02, lab2),
            ("D3",) + (p3, y03, lab3)]


def base_scenario():
    """Reference disease scenario D1: patient_trans a=0.5 gamma=0.5 A1."""
    p, y0, lab = B.patient_trans(0.5, 0.5, "A1_mcc_assembly")
    return p, y0, lab
