"""Attachment-error module coupled to the Henze et al. (2017) SAC network.

Extends the published model with a minimal, literature-supported
kinetochore-attachment state module to predict chromosome-segregation
error risk (to-do list item 3-6, 2026-10-23).

Model (per sister-chromosome pair; N_PAIRS = 46, 92 kinetochores)
-----------------------------------------------------------------
States:
  U  unattached / incompletely attached (2 unattached kinetochores;
     CHECKPOINT-VISIBLE: each unattached kinetochore signals)
  E  erroneously (merotelically) attached (2 attached kinetochores;
     CHECKPOINT-SILENT: merotelic attachments evade the SAC)
  C  correctly (amphitelic) attached (checkpoint-silent)

Transitions (first order per pair):
  U -> C   rate alpha(mu)   correct capture / congression
  U -> E   rate beta        merotelic error capture (BUBR1-independent
                            in the base case; documented assumption)
  E -> C   rate gamma(mu)   error correction (Aurora B / PP2A-B56
                            surveillance; BUBR1-dependent)

  du/dt = -(alpha+beta) u
  de/dt = beta u - gamma e
  (c = 46 - u - e)

Coupling to Henze (item 4)
--------------------------
Henze's prescribed one-way attachment trajectory (reaction r1,
KinU -> KinA with rate kattach * KinU) is REPLACED by the module:
  KinU(t) = 2 u(t)          (only U pairs carry unattached kinetochores)
  KinA(t) = 92 - 2 u(t)     (E and C pairs both count as attached)
All 17 other molecular species, all 19 reactions, and the three
SBML assignment rules (k3 = 0.015*KinU, k8/k9 = 0.014*KinA) are
unchanged. r1 enters d[KinU] and d[KinA] only, so no other species
feels the replacement.

Limiting case (reduction to the original model):
  beta = 0 and alpha(mu) = alpha0 for all mu  =>  u(t) = 46 e^{-alpha0 t},
  KinU = 92 e^{-alpha0 t}, identical to Henze with kattach = alpha0.
Verified numerically in limit_test().

BUBR1 dependence (item 5)
-------------------------
mu = functional BUBR1 pool fraction vs WT (same mu as the existing
pool-conserving abundance perturbation in bubr1_perturbation.py).
  alpha(mu) = alpha0 * (1 - w_a * (1 - mu))
      w_a in [0,1]: BUBR1's contribution to congression speed.
      Evidence: BUBR1-depleted HeLa cells show chromosome-alignment
      defects that WT BUBR1 restores (Suijkerbuijk et al. 2010,
      Fig 4C); BUBR1 recruits PP2A-B56 which promotes congression
      (Xu et al. 2013). w_a is NOT assumed to be 1: congression has
      many drivers besides BUBR1 (dynein, CENP-E, etc.).
  gamma(mu) = gamma0 * mu**k
      k: effective order of BUBR1 dependence of error correction.
      Evidence: correction of erroneous attachments requires the
      BUBR1:B56:PP2A complex opposing Aurora B (Xu et al. 2013;
      Foley et al. 2011). k=1 (first order in complex availability)
      is the base case; k=2,3 tested as bounded scenarios.
  beta: fixed at f_err * alpha0 (fraction of captures that are
      merotelic); base f_err = 0.02 (1-5% scenario range, human
      metaphase merotelic frequencies).

Baseline values (provenance in docs/Attachment_Error_Extension.md):
  alpha0 = 0.0035 s^-1  (= Henze kattach; preserves the WT
      attachment trajectory exactly and reproduces human ~20-min
      congression completion)
  gamma0 = 0.002/0.005/0.01 s^-1 (minutes-scale Aurora B/PP2A
      correction; bounded ASSUMPTION, no direct human measurement
      of the merotelic-correction rate was found)
"""
from __future__ import annotations

import copy

import numpy as np

import henze2017 as H
import bubr1_perturbation as B

N_PAIRS = 46
KT_TOTAL = 92.0

ALPHA0 = 0.0035          # s^-1, = Henze kattach (WT congression rate)

MECH_DEFAULT = {
    "w_a": 0.5,          # BUBR1 contribution to congression speed
    "gamma0": 0.005,     # s^-1, WT error-correction rate
    "k": 1,              # order of BUBR1 dependence of correction
    "f_err": 0.02,       # merotelic fraction of captures (beta/alpha0)
}


def mech_rates(mu: float, mech: dict) -> dict:
    """Per-pair rates at BUBR1 pool fraction mu."""
    m = dict(MECH_DEFAULT, **mech)
    alpha = ALPHA0 * (1.0 - m["w_a"] * (1.0 - mu))
    beta = m["f_err"] * ALPHA0
    gamma = m["gamma0"] * mu ** m["k"]
    return {"alpha": alpha, "beta": beta, "gamma": gamma}


def rhs_extended(t: float, y: np.ndarray, p: dict, rates: dict) -> np.ndarray:
    """21-variable RHS: 19 Henze species + (u, e).

    KinU/KinA are slaved to the attachment module (KinU = 2u,
    KinA = 92 - 2u); the other 17 species follow the original
    network exactly.
    """
    u = max(y[19], 0.0)
    kinu = 2.0 * u
    # slave the checkpoint-visible/attached counts for the assignment rules
    yy = y.copy()
    yy[H.IDX["KinU"]] = kinu
    yy[H.IDX["KinA"]] = KT_TOTAL - kinu
    d = H.rhs(t, yy, p)[:19]  # H.rhs zeros_like's its input; keep 19 species
    a, b, g = rates["alpha"], rates["beta"], rates["gamma"]
    d[H.IDX["KinU"]] = -(a + b) * kinu
    d[H.IDX["KinA"]] = (a + b) * kinu
    d = list(d)
    d.append(-(a + b) * u)                      # du/dt
    d.append(b * u - g * max(y[20], 0.0))       # de/dt
    return np.array(d)


def y0_extended(p: dict | None = None, y0_19: np.ndarray | None = None) -> np.ndarray:
    """Initial condition: 19 species (BUBR1 pool already scaled by the
    caller via bubr1_perturbation) + u=46, e=0 (all pairs unattached
    at NEB)."""
    if y0_19 is None:
        y0_19 = H.Y0_VEC.copy()
    else:
        y0_19 = np.asarray(y0_19, dtype=float).copy()
    return np.concatenate([y0_19, [N_PAIRS, 0.0]])


def simulate_extended(mu: float, mech: dict, t_end: float = 7200.0,
                      n_points: int = 3601,
                      p: dict | None = None,
                      y0_19: np.ndarray | None = None) -> dict:
    """Integrated simulation.

    mu: BUBR1 pool fraction. If p/y0_19 are None, the BUBR1 pool is
    set to POOL_WT * mu (pool-conserving, as in B.abundance_only).
    """
    if p is None:
        p, y0_19, _ = B.abundance_only(mu)
    rates = mech_rates(mu, mech)
    y0 = y0_extended(p, y0_19)
    from scipy.integrate import solve_ivp
    sol = solve_ivp(lambda t, y: rhs_extended(t, y, p, rates),
                    (0.0, t_end), y0, method="LSODA",
                    rtol=1e-9, atol=1e-12,
                    t_eval=np.linspace(0.0, t_end, n_points))
    if not sol.success:
        raise RuntimeError(f"integration failed: {sol.message}")
    y = sol.y.T
    res = {
        "t": sol.t,
        "y": y[:, :19],                       # 19 Henze species (original order)
        "species": H.SPECIES,
        "u": y[:, 19],                        # unattached pairs
        "e": y[:, 20],                        # erroneous pairs
        "mu": mu,
        "mech": dict(MECH_DEFAULT, **mech),
        "rates": rates,
        "params": p,
    }
    return res


# ---------------------------------------------------------------------------
# Readouts (item 6)
# ---------------------------------------------------------------------------
def error_output(res: dict, floor: float = 0.02, phi: float = 1.0) -> dict:
    """Segregation-error output at anaphase commitment.

    t_commit: model-defined anaphase commitment (securin < 10% of
    initial), identical convention to the canonical readouts.
    e_tc: expected number of unresolved erroneous (merotelic) pairs
    at commitment. p_pair: per-pair probability (mean field).
    P_error: probability that >= 1 pair is unresolved at commitment
    (binomial over 46 independent pairs).
    P_pred: floor + phi * P_error, where floor = 0.02 is the
    empirical WT mis-segregation level (Suijkerbuijk 100% arm)
    attributed to model-unresolved background (late stochastic
    errors, tension noise), and phi = 1 assumes each unresolved
    merotelic pair produces a segregation defect in that cell
    (merotelics generate lagging chromosomes: CKAP5 2024).
    """
    t = res["t"]
    sec = res["y"][:, H.IDX["Securin"]]
    sec0 = sec[0] or 1e-30
    below = np.where(sec < 0.10 * sec0)[0]
    tc = float(t[below[0]]) if below.size else None
    kinu = res["y"][:, H.IDX["KinU"]]
    out = {
        "mu": res["mu"],
        "t_commit_s": tc,
        "arrest": tc is None,
        "kinu_at_commit": None,
        "e_at_commit_pairs": None,
        "p_pair": None,
        "P_error": None,
        "P_pred": None,
    }
    if tc is not None:
        j = int(np.searchsorted(t, tc))
        j = min(j, len(t) - 1)
        out["kinu_at_commit"] = float(kinu[j]) / 92.0
        e_tc = max(float(res["e"][j]), 0.0)
        p = e_tc / N_PAIRS
        out["e_at_commit_pairs"] = e_tc
        out["p_pair"] = p
        out["P_error"] = 1.0 - (1.0 - p) ** N_PAIRS
        out["P_pred"] = min(1.0, floor + phi * out["P_error"])
    return out


def checkpoint_readouts(res: dict) -> dict:
    """Canonical checkpoint-side readouts from an integrated run
    (same conventions as shared_eval; leak window uses KinU > 1)."""
    t = res["t"]
    ro = H.readouts({"t": t, "y": res["y"], "species": H.SPECIES})
    kinu = res["y"][:, H.IDX["KinU"]]
    apcc = res["y"][:, H.IDX["APCC_Cdc20"]]
    mask = kinu > 1.0
    leak = float(apcc[mask].max()) if mask.any() else float(apcc.max())
    return {
        "mcc_peak_uM": ro["mcc_peak_uM"],
        "t_commit_s": ro["t_anaphase_commitment_s"],
        "securin_frac_120min": ro["securin_frac_at_30min"],  # end of 7200 s window
        "apcc_cdc20_max_uM": ro["apcc_cdc20_max_uM"],
        "leak_uM": leak,
    }


# ---------------------------------------------------------------------------
# Limiting-case test (item 4): reduction to the original model
# ---------------------------------------------------------------------------
def limit_test(mu: float = 0.5, tol_traj: float = 1e-6,
               tol_rhs: float = 1e-12) -> dict:
    """beta = 0 (f_err = 0) and w_a = 0 (alpha = alpha0 for all mu)
    must reproduce the original Henze model (kattach = alpha0).

    Two checks:
    (1) ALGEBRAIC: at points along the extended trajectory, the
        extended RHS restricted to the 19 species must equal the
        original RHS evaluated at the slaved state to machine
        precision (structural identity of the networks).
    (2) NUMERIC: the full trajectories must agree to solver-noise
        level. Two systems of different dimension (19 vs 21) take
        different LSODA step paths, so exact agreement is not
        expected; tol_traj = 1e-6 is ~4 orders of magnitude below
        the smallest biologically relevant quantity (1e-10 uM
        species, 1e-4 kinetochore fraction).
    """
    mech = {"w_a": 0.0, "gamma0": 0.005, "k": 1, "f_err": 0.0}
    p, y0_19, _ = B.abundance_only(mu)
    ext = simulate_extended(mu, mech, t_end=3600.0, p=p, y0_19=y0_19)
    rates = ext["rates"]
    # (1) algebraic identity at sample points along the trajectory
    max_rhs = 0.0
    for i in range(0, len(ext["t"]), max(1, len(ext["t"]) // 20)):
        t = float(ext["t"][i])
        y21 = np.concatenate([ext["y"][i], [ext["u"][i], ext["e"][i]]])
        d_ext = rhs_extended(t, y21, p, rates)[:19]
        yy = y21.copy()
        yy[H.IDX["KinU"]] = 2.0 * max(ext["u"][i], 0.0)
        yy[H.IDX["KinA"]] = KT_TOTAL - 2.0 * max(ext["u"][i], 0.0)
        d_orig = H.rhs(t, yy, p)[:19]
        scale = np.maximum(np.abs(d_orig), 1e-30)
        max_rhs = max(max_rhs, float(np.max(np.abs(d_ext - d_orig) / scale)))
    # (2) numeric trajectory match
    orig = H.simulate(t_end=3600.0, p=p, y0=y0_19, n_points=len(ext["t"]))
    diff = float(np.max(np.abs(ext["y"] - orig["y"])))
    return {
        "mu": mu,
        "max_rel_rhs_diff": max_rhs,
        "max_abs_traj_diff": diff,
        "ok": max_rhs < tol_rhs and diff < tol_traj,
    }
