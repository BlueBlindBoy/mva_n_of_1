"""Deterministic biochemical spindle-assembly-checkpoint model of
Henze, Dittrich & Ibrahim (2017), Scientific Reports 7, 3865.

Reconstructed 1:1 from the original SBML deposit:
  BioModels MODEL1812210002 (v2), file Henze2017.xml
  (local copy: Henze2017_original_SBML.xml in this directory;
   machine-readable extraction: henze_model.json).

Conventions
-----------
- One compartment of size 1.0; time in seconds.
- Soluble species in uM; kinetochores (KinU/KinA) are counts.
- Reversible reaction pairs are stored as net mass-action rates.
- Three assignment rules override the constant parameters at all times:
      k3(t) = 0.015 * KinU(t)    (Mad2 turnover, per unattached kinetochore)
      k8(t) = 0.014 * KinA(t)    (MCC silencing, per attached kinetochore)
      k9(t) = 0.014 * KinA(t)    (same for the APC/C:Cdc20:MCC form)
- Reactions 9-10 (Mad1 'template') are present in the SBML but dormant
  (Mad1_Mad2(0) = 0 and never produced); kept for fidelity.

No parameters were invented: every value is read from the deposit
(henze_model.json). See docs/Published_SAC_Model_Reproduction.md for
provenance and the one flagged discrepancy (k_5: 0.2 SBML vs 0.02
Table S1).
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp

# ---------------------------------------------------------------------------
# Model definition (order of the state vector)
# ---------------------------------------------------------------------------
SPECIES = [
    "KinU", "KinA",                       # kinetochore counts
    "APCC", "APCC_Cdc20",                 # APC/C forms
    "Securin", "CyclinB",                 # APC/C substrates
    "Cdc20",                              # free activator
    "O_Mad2", "C_Mad2",                   # open / closed MAD2
    "Cdc20_C_Mad2",                       # Cdc20:closed-MAD2 (pre-MCC)
    "MCC", "BCC",                         # Mitotic Checkpoint Complex / BCC
    "APCC_MCC", "APCC_Cdc20_MCC", "APCC_BCC", "APCC_Cdc20_C_Mad2",  # APC/C-bound
    "BubR1_Bub3",                         # free BubR1:Bub3 (BUBR1 pool)
    "Mad1_Mad2", "Mad1_Mad2_C_Mad2",      # dormant template species
]
IDX = {s: i for i, s in enumerate(SPECIES)}

PARAMS = {
    "kattach": 0.0035,
    "k1": 5.0, "k_1": 0.08,
    "kD": 0.05,
    "k2": 100.0, "k_2": 0.08,
    "k3": 1.38, "k_3": 0.2,
    "k4": 10.0, "k_4": 0.1,
    "k5": 10.0, "k_5": 0.2,
    "k6": 0.001, "k_6": 0.01,
    "k7": 0.01, "k_7": 0.2,
    "k8": 0.0, "k9": 0.0,
    "k10": 0.2,
    "kT1": 0.01, "k_T1": 0.02, "kT2": 10.0,
    # assignment-rule coefficients (per-kinetochore rates, s^-1 per kinetochore)
    "k3_per_kinU": 0.015,
    "k8_per_kinA": 0.014,
    "k9_per_kinA": 0.014,
}

Y0 = {
    "KinU": 92.0, "KinA": 0.0,
    "APCC": 0.09, "APCC_Cdc20": 0.0,
    "Securin": 0.1, "CyclinB": 0.0,
    "Cdc20": 0.13,
    "O_Mad2": 0.15, "C_Mad2": 0.01875,
    "Cdc20_C_Mad2": 0.0,
    "MCC": 0.0, "BCC": 0.0,
    "APCC_MCC": 0.0, "APCC_Cdc20_MCC": 0.0, "APCC_BCC": 0.0,
    "APCC_Cdc20_C_Mad2": 0.0,
    "BubR1_Bub3": 0.22,
    "Mad1_Mad2": 0.0, "Mad1_Mad2_C_Mad2": 0.0,
}
Y0_VEC = np.array([Y0[s] for s in SPECIES])


def rhs(t: float, y: np.ndarray, p: dict | None = None) -> np.ndarray:
    """ODE right-hand side. 19 reactions + 3 assignment rules."""
    p = PARAMS if p is None else p
    g = lambda s: max(y[IDX[s]], 0.0)  # guard tiny negatives from the solver

    # --- assignment rules (SBML) -------------------------------------------
    k3 = p["k3_per_kinU"] * g("KinU")
    k8 = p["k8_per_kinA"] * g("KinA")
    k9 = p["k9_per_kinA"] * g("KinA")

    # --- 19 net reaction rates (compartment = 1) ---------------------------
    r1 = p["kattach"] * g("KinU")
    r2 = p["k1"] * g("APCC") * g("Cdc20") - p["k_1"] * g("APCC_Cdc20")
    r3 = p["kD"] * g("APCC_Cdc20") * g("Securin")
    r4 = p["kD"] * g("APCC_Cdc20") * g("CyclinB")
    r5 = p["k2"] * g("APCC") * g("MCC") - p["k_2"] * g("APCC_MCC")
    r6 = k3 * g("O_Mad2") - p["k_3"] * g("C_Mad2")
    r7 = p["k4"] * g("Cdc20") * g("C_Mad2")
    r8 = p["k5"] * g("Cdc20_C_Mad2") * g("BubR1_Bub3") - p["k_5"] * g("MCC")
    r9 = p["kT1"] * g("Mad1_Mad2") * g("O_Mad2") - p["k_T1"] * g("Mad1_Mad2_C_Mad2")
    r10 = p["kT2"] * g("Mad1_Mad2_C_Mad2") * g("Cdc20")
    r11 = p["k6"] * g("O_Mad2") * g("Cdc20") - p["k_6"] * g("Cdc20_C_Mad2")
    r12 = p["k7"] * g("Cdc20") * g("BubR1_Bub3") - p["k_7"] * g("BCC")
    r13 = p["kattach"] * g("O_Mad2")
    r14 = p["k4"] * g("APCC_Cdc20") * g("C_Mad2") - p["k_4"] * g("APCC_Cdc20_C_Mad2")
    r15 = p["k5"] * g("APCC_Cdc20_C_Mad2") * g("BubR1_Bub3") - p["k_5"] * g("APCC_MCC")
    r16 = p["k2"] * g("APCC_Cdc20") * g("MCC") - p["k_2"] * g("APCC_Cdc20_MCC")
    r17 = k8 * g("APCC_MCC")
    r18 = k9 * g("APCC_Cdc20_MCC")
    r19 = p["k10"] * g("APCC_BCC")

    d = np.zeros_like(y)
    d[IDX["KinU"]] = -r1
    d[IDX["KinA"]] = r1
    d[IDX["APCC"]] = -r2 - r5
    d[IDX["APCC_Cdc20"]] = r2 - r14 - r16 + r19
    d[IDX["Securin"]] = -r3
    d[IDX["CyclinB"]] = -r4
    d[IDX["Cdc20"]] = -r2 - r7 - r11 - r12 - r10
    d[IDX["O_Mad2"]] = -r6 - r9 - r11 - r13 + r17 + r18
    d[IDX["C_Mad2"]] = r6 - r7 - r14
    d[IDX["Cdc20_C_Mad2"]] = r7 + r11 + r10 - r8
    d[IDX["MCC"]] = r8 - r5 - r16
    d[IDX["BCC"]] = r12
    d[IDX["APCC_MCC"]] = r5 + r15 - r17
    d[IDX["APCC_Cdc20_MCC"]] = r16 - r18
    d[IDX["APCC_BCC"]] = r17 + r18 - r19
    d[IDX["APCC_Cdc20_C_Mad2"]] = r14 - r15
    d[IDX["BubR1_Bub3"]] = -r8 - r12 - r15 + r19
    d[IDX["Mad1_Mad2"]] = -r9 + r10
    d[IDX["Mad1_Mad2_C_Mad2"]] = r9 - r10
    return d


def simulate(t_end: float = 3600.0, y0: np.ndarray | None = None,
             p: dict | None = None, rtol: float = 1e-9, atol: float = 1e-12,
             n_points: int = 7201) -> dict:
    """Integrate the wild-type (or perturbed) model. Returns dict with
    't' (s) and 'y' (n_points x n_species) plus species names."""
    y0 = Y0_VEC if y0 is None else np.asarray(y0, dtype=float)
    pp = PARAMS if p is None else p
    sol = solve_ivp(lambda t, y: rhs(t, y, pp), (0.0, t_end), y0, method="LSODA",
                    rtol=rtol, atol=atol, dense_output=False,
                    t_eval=np.linspace(0.0, t_end, n_points))
    if not sol.success:
        raise RuntimeError(f"integration failed: {sol.message}")
    return {"t": sol.t, "y": sol.y.T, "species": SPECIES, "params": p or PARAMS}


# ---------------------------------------------------------------------------
# Biological readouts (model-defined; NOT validated clinical observables)
# ---------------------------------------------------------------------------
def _col(res, name):
    return res["y"][:, IDX[name]]


def readouts(res: dict) -> dict:
    """Model-defined readouts. 'Anaphase commitment' is a model construct
    (securin < 10% of initial), not a measured cell event."""
    t = res["t"]
    sec0 = res["y"][:, IDX["Securin"]][0] or 1e-30
    mcc_peak_i = int(np.argmax(_col(res, "MCC")))
    out = {
        "t_mcc_peak_s": float(t[mcc_peak_i]),
        "mcc_peak_uM": float(_col(res, "MCC")[mcc_peak_i]),
        "mcc_at_30min_uM": float(_col(res, "MCC")[-1]),
        "apcc_cdc20_max_uM": float(np.max(_col(res, "APCC_Cdc20"))),
        "apcc_cdc20_at_30min_uM": float(_col(res, "APCC_Cdc20")[-1]),
        "securin_frac_at_30min": float(_col(res, "Securin")[-1] / sec0),
        "kinA_at_30min": float(_col(res, "KinA")[-1]),
        "bubr1_free_at_30min_uM": float(_col(res, "BubR1_Bub3")[-1]),
    }
    # model-defined anaphase commitment time: securin first drops below 10%
    below = np.where(_col(res, "Securin") < 0.10 * sec0)[0]
    out["t_anaphase_commitment_s"] = float(t[below[0]]) if below.size else None
    return out


# ---------------------------------------------------------------------------
# Conserved quantities
# ---------------------------------------------------------------------------
# NOTE: the deposited network is NOT Cdc20-mass-conserving (a lumped-
# representation artifact): reaction 16 consumes two Cdc20-bearing species
# (APCC_Cdc20 + MCC) while reaction 18 releases one (APCC_BCC). No content
# weighting conserves both, so free/total Cdc20 is not an invariant. The
# true invariants are the 4 left-null-space vectors of the stoichiometric
# matrix, computed here from henze_model.json (single source of truth).
#
# The kinetochore and BubR1 pool checks are kept as *approximate* sanity
# checks (they hold exactly in this network); the Mad2 pool is not conserved
# (reaction 13 is the published sequestration convention).
def _invariants():
    import os, json
    here = os.path.dirname(os.path.abspath(__file__))
    M = json.load(open(os.path.join(here, "henze_model.json")))
    sp = list(M["species_initial_conditions"].keys())
    S = np.zeros((len(sp), len(M["reactions"])))
    for j, r in enumerate(M["reactions"]):
        for s, st in r["reactants"].items():
            S[sp.index(s), j] -= st
        for s, st in r["products"].items():
            S[sp.index(s), j] += st
    _, svals, _ = np.linalg.svd(S)
    tol = 1e-9 * svals[0]
    U, _, _ = np.linalg.svd(S)
    C = U[:, svals < tol]
    # keep only components on our 19 species (json order == SPECIES order check)
    assert sp == SPECIES, "species order mismatch between json and module"
    return C  # shape (19, n_invariants)

INV = _invariants()


def conservation_report(res: dict, tol: float = 1e-6) -> dict:
    """Check the 4 exact invariants of the stoichiometric matrix, plus
    explicit sanity pools (kinetochores exact; BubR1 exact here)."""
    report = {}
    yT = res["y"]  # (n_t, n_species)
    vals = yT @ INV  # (n_t, 4)
    for k in range(INV.shape[1]):
        drift = float(np.max(np.abs(vals[:, k] - vals[0, k])))
        scale = max(abs(vals[0, k]), 1.0)
        report[f"invariant_{k+1}"] = {"value0": float(vals[0, k]),
                                      "max_abs_drift": drift,
                                      "ok": drift <= tol * scale}
    for name, (sp, _total_unused) in [("kinetochores", (["KinU", "KinA"], None)),
                                      ("bubr1_pool", (["BubR1_Bub3", "MCC", "BCC",
                                                       "APCC_MCC", "APCC_Cdc20_MCC",
                                                       "APCC_BCC"], None))]:
        v = sum(yT[:, IDX[s]] for s in sp)
        total0 = float(v[0])  # network conserves each pool exactly; t=0 is the reference
        drift = float(np.max(np.abs(v - total0)))
        report[name] = {"total0": total0, "max_abs_drift": drift,
                        "ok": drift <= max(tol * max(total0, 1e-9), 1e-9)}
    return report


if __name__ == "__main__":
    import json
    res = simulate(t_end=3600.0)
    ro = readouts(res)
    rep = conservation_report(res)
    print("== Wild-type baseline, t in [0, 3600] s ==")
    for k, v in ro.items():
        print(f"  {k:28s} {v}")
    print("conservation:")
    for k, v in rep.items():
        print(f"  {k:16s} drift={v['max_abs_drift']:.3e} ok={v['ok']}")
