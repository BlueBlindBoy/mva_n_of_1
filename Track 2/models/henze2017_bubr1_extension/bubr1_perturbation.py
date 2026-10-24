#!/usr/bin/env python3
"""BUBR1-deficiency perturbation layer for the Henze et al. (2017) network.

Translates the patient's two BUB1B variants into explicitly labeled,
mechanistically distinct network perturbations (addendum Part II).

Design rules
------------
1. Pool-conserving abundance perturbation: every BUBR1-bearing species in the
   network is made from, and returns to, the free `BubR1_Bub3` pool, so
   scaling the initial `BubR1_Bub3` concentration by a factor f scales the
   ENTIRE functional BUBR1 pool (free + all bound forms) by f, exactly.
   (Verified: `bubr1_pool` conservation check in henze2017.conservation_report.)

2. No invented parameters: all perturbation values are SCENARIOS drawn from
   literature-informed ranges, recorded as such. No kinetic constant for this
   patient was measured.

3. Mechanistically distinct activity alternatives (none selected as 'the'
   mechanism, per addendum Sec. 4):
      A1: impaired MCC assembly        -> k5  (r8 forward)  x gamma
      A2: impaired BCC assembly        -> k7  (r12 forward) x gamma
      A3: impaired APC/C-inhibition    -> k2  (r5, r16 fwd) x gamma
          (MCC/BCC binding to APC/C)

4. Allele-to-pool mapping (documented limitation: the published network has a
   single BUBR1 species; allele-specific species would be new, poorly
   constrained biology, which we do NOT add). Reference: 0.22 uM pool =
   diploid WT = 2 x 0.11 uM per fully functional allele.

   patient_trans (V1 null allele in trans to V2 missense):
       pool = 0.11 * a           (a = missense abundance fraction)
       rate scale for mechanism A_i = gamma  (only BUBR1 is the mutant's)
   patient_cis (both variants on one chromosome; other allele assumed WT
       until proven otherwise - assumption recorded, NOT a measurement):
       pool = 0.11 * (1 + a)
       rate scale = (1 + a*gamma) / (1 + a)  (linear mole-fraction mixing of
       a WT population at full activity with a mutant population at activity
       gamma; documented simplification - the network cannot track two
       BUBR1 populations)
"""
from __future__ import annotations

import copy
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "models", "henze2017_original"))
import henze2017 as H  # noqa: E402

POOL_WT = 0.22          # uM, diploid WT reference (from the deposit)
PER_ALLELE = 0.11       # uM, one fully functional allele (documented split)

# scenario ranges (literature-informed; NOT patient measurements)
ABUNDANCE_FRACTIONS = [1.0, 0.5, 0.2, 0.1]        # deficiency pool fractions
MISSENSE_ABUNDANCE_A = [1.0, 0.5, 0.25]           # missense abundance fraction
ACTIVITY_GAMMA = [1.0, 0.5, 0.25, 0.1]            # specific-activity scale

MECHANISMS = {
    "A1_mcc_assembly": ["k5"],
    "A2_bcc_assembly": ["k7"],
    "A3_apcc_binding": ["k2"],
}


def _perturb_abundance(params: dict, y0: np.ndarray, f: float):
    """Scale the whole functional BUBR1 pool by f (pool-conserving)."""
    y0 = y0.copy()
    y0[H.IDX["BubR1_Bub3"]] = POOL_WT * f
    return params, y0


def _perturb_activity(params: dict, y0: np.ndarray, gamma: float, mech: str):
    p = copy.deepcopy(params)
    for key in MECHANISMS[mech]:
        p[key] = params[key] * gamma
    return p, y0


def patient_trans(a: float, gamma: float, mech: str):
    """Null allele in trans to missense allele.

    Returns (params, y0, label). The only BUBR1 present (beyond residual NMD
    escape, ignored) is the missense protein at abundance fraction a and
    specific activity gamma under mechanism `mech`.
    """
    p = dict(H.PARAMS)
    y0 = H.Y0_VEC.copy()
    y0[H.IDX["BubR1_Bub3"]] = PER_ALLELE * a
    p, y0 = _perturb_activity(p, y0, gamma, mech)
    label = (f"trans a={a:g} gamma={gamma:g} [{mech}] "
             f"pool={PER_ALLELE * a:.3f}uM")
    return p, y0, label


def patient_cis(a: float, gamma: float, mech: str):
    """Both variants on one chromosome; other allele assumed WT (assumption).

    Linear mole-fraction mixing of the WT allele's BUBR1 (full activity) with
    the mutant chromosome's BUBR1 (abundance a, activity gamma).
    """
    p = dict(H.PARAMS)
    y0 = H.Y0_VEC.copy()
    pool = PER_ALLELE * (1.0 + a)
    y0[H.IDX["BubR1_Bub3"]] = pool
    rate_scale = (1.0 + a * gamma) / (1.0 + a)
    p, _ = _perturb_activity(p, y0, rate_scale, mech)
    label = (f"cis a={a:g} gamma={gamma:g} [{mech}] "
             f"pool={pool:.3f}uM rate_scale={rate_scale:.3f}")
    return p, y0, label


def abundance_only(f: float):
    """BUBR1 abundance deficiency without activity change (scenario)."""
    p = dict(H.PARAMS)
    p, y0 = _perturb_abundance(p, H.Y0_VEC.copy(), f)
    label = f"abundance f={f:g} pool={POOL_WT * f:.3f}uM"
    return p, y0, label


def build_all_scenarios():
    """Full Part II scenario grid. Returns list of (name, params, y0, label)."""
    out = []
    for f in ABUNDANCE_FRACTIONS:
        p, y0, lab = abundance_only(f)
        out.append((f"abundance_f{f:g}", p, y0, lab))
    for a in MISSENSE_ABUNDANCE_A:
        for gamma in ACTIVITY_GAMMA:
            for mech in MECHANISMS:
                p, y0, lab = patient_trans(a, gamma, mech)
                out.append((f"trans_a{a:g}_g{gamma:g}_{mech.split('_')[0]}",
                            p, y0, lab))
                p, y0, lab = patient_cis(a, gamma, mech)
                out.append((f"cis_a{a:g}_g{gamma:g}_{mech.split('_')[0]}",
                            p, y0, lab))
    return out


def checkpoint_leak_readout(res: dict) -> dict:
    """Model-defined 'checkpoint competence under challenge':
    max APC/C:Cdc20 over the window in which >=1 kinetochore remains
    unattached. No external threshold; the window is defined by the state."""
    t = res["t"]
    kinu = res["y"][:, H.IDX["KinU"]]
    apcc = res["y"][:, H.IDX["APCC_Cdc20"]]
    mask = kinu > 1.0
    if not mask.any():
        return {"leak_max_uM": float(apcc.max()), "leak_window_s": 0.0}
    t_last = t[mask][-1]
    return {"leak_max_uM": float(apcc[mask].max()),
            "leak_window_s": float(t_last)}
