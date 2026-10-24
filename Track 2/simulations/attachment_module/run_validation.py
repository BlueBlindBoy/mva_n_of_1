#!/usr/bin/env python3
"""Validation of the attachment-error extension against the Suijkerbuijk
et al. (2010) HeLa depletion observations (item 8), with the item-9
keep-or-reject evidence summary.

Observed (held-out for the mechanism parameters; the WT arm supplies
the empirical floor only):
  residual BUBR1 100%  ->  mis-segregation ~0.02
  residual BUBR1  13%  ->  mis-segregation ~0.03
  residual BUBR1   6%  ->  mis-segregation ~0.50
"""
from __future__ import annotations

import csv
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "models", "henze2017_bubr1_extension"))
sys.path.insert(0, os.path.join(ROOT, "models", "henze2017_original"))
sys.path.insert(0, HERE)

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PANEL_CSV = os.path.join(ROOT, "results", "Attachment_Depletion_Panel.csv")
OBS = {1.0: 0.02, 0.13: 0.03, 0.06: 0.50}
MU_LIST = [1.0, 0.13, 0.06]
ARM_PASS = (1.0 / 3.0, 3.0)   # within 3x of observation


def load_panel():
    with open(PANEL_CSV) as f:
        rows = list(csv.DictReader(f))
    for r in rows:
        for k in ("w_a", "gamma0", "k", "f_err", "mu", "t_commit_s",
                  "e_at_commit_pairs", "p_pair", "P_error", "P_pred",
                  "obs_missegregation", "pred_over_obs", "mcc_peak_uM",
                  "leak_uM"):
            v = r[k]
            r[k] = None if v in ("", "None") else float(v)
        r["arrest"] = r["arrest"] == "True"
    return rows


def setting_summary(rows):
    """Per setting: pred/obs ratio for each arm, pass flags."""
    out = {}
    for r in rows:
        sid = r["setting"]
        d = out.setdefault(sid, {"setting": sid, "w_a": r["w_a"],
                                 "gamma0": r["gamma0"], "k": r["k"],
                                 "f_err": r["f_err"]})
        d[f"ratio_{r['mu']:g}"] = r["pred_over_obs"]
        d[f"P_pred_{r['mu']:g}"] = r["P_pred"]
    for d in out.values():
        ratios = [d[f"ratio_{m:g}"] for m in MU_LIST
                  if d[f"ratio_{m:g}"] is not None]
        d["arms_within_3x"] = sum(
            1 for m in MU_LIST
            if d[f"ratio_{m:g}"] is not None
            and ARM_PASS[0] <= d[f"ratio_{m:g}"] <= ARM_PASS[1])
        d["all_arms_within_3x"] = d["arms_within_3x"] == len(MU_LIST)
        # held-out arms only (WT arm is the floor calibration, not a
        # free prediction)
        d["heldout_arms_within_3x"] = sum(
            1 for m in (0.13, 0.06)
            if d[f"ratio_{m:g}"] is not None
            and ARM_PASS[0] <= d[f"ratio_{m:g}"] <= ARM_PASS[1])
        # missing (arrested) arms count as a large deviation so that
        # arrested settings rank last, not first
        dev = [max(np.log10(abs(x)), 0.0) for x in ratios]
        dev += [1.0] * (len(MU_LIST) - len(ratios))
        d["max_log10_dev"] = max(dev) if dev else None
    return list(out.values())


def main():
    rows = load_panel()
    ss = setting_summary(rows)
    ss.sort(key=lambda d: (d["max_log10_dev"] if d["max_log10_dev"] is not None
                           else 99))

    out_csv = os.path.join(ROOT, "results", "Attachment_Validation_Summary.csv")
    with open(out_csv, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(ss[0].keys()))
        w.writeheader()
        w.writerows(ss)
    print(f"wrote {out_csv}")

    n_pass = sum(1 for d in ss if d["all_arms_within_3x"])
    n_ho = sum(1 for d in ss if d["heldout_arms_within_3x"] == 2)
    n_2 = sum(1 for d in ss if d["arms_within_3x"] >= 2)
    print(f"settings: {len(ss)}; both held-out arms (13%,6%) within 3x: {n_ho}; "
          f">=2 arms: {n_2}; all 3 arms (incl. WT calibration): {n_pass}")
    print("\nbest settings (max log10 deviation):")
    for d in ss[:10]:
        print(f"  {d['setting']}: 100%={d['ratio_1']}, 13%={d['ratio_0.13']}, "
              f"6%={d['ratio_0.06']}  (max log10 dev "
              f"{d['max_log10_dev']:.3f})")

    # ---- figure (item 8 / 12): model vs experiment -----------------------
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.4))
    ax = axes[0]
    best = ss[0]
    # curve for the best setting
    ms = np.linspace(0.0, 1.0, 101)
    cur = [r for r in rows if r["setting"] == best["setting"]]
    cur = sorted(cur, key=lambda r: -r["mu"])
    ax.plot([r["mu"] for r in cur], [r["P_pred"] for r in cur], "k-o",
            lw=2, ms=6, label=f"model, {best['setting']}")
    ax.plot(MU_LIST, [OBS[m] for m in MU_LIST], "rs-", lw=2, ms=7,
            label="Suijkerbuijk 2010 (HeLa shRNA)")
    ax.axhline(0.02, color="gray", ls=":", lw=1)
    ax.text(0.99, 0.024, "WT empirical floor (0.02)", ha="right",
            fontsize=8, color="gray")
    ax.set_xlabel("residual BUBR1 (fraction of WT)")
    ax.set_ylabel("fraction of cells with segregation defects")
    ax.set_title("Model vs experiment (best setting)")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    ax = axes[1]
    x13 = [d["P_pred_0.13"] for d in ss if d["P_pred_0.13"] is not None]
    x6 = [d["P_pred_0.06"] for d in ss if d["P_pred_0.06"] is not None]
    for m, vals, col, lab in ((0.13, x13, "tab:blue", "13% BUBR1"),
                              (0.06, x6, "tab:red", "6% BUBR1")):
        ax.boxplot([vals], positions=[0 if m == 0.13 else 1], widths=0.45,
                   showfliers=False,
                   patch_artist=True, boxprops=dict(facecolor=col, alpha=0.35))
        ax.plot([0 if m == 0.13 else 1], [OBS[m]], "ks-", ms=8,
                label=f"observed {OBS[m]:.2f}")
    ax.axhline(OBS[0.13], color="tab:blue", ls=":", lw=1)
    ax.axhline(OBS[0.06], color="tab:red", ls=":", lw=1)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["13% BUBR1", "6% BUBR1"])
    ax.set_ylabel("predicted fraction with segregation defects")
    ax.set_title(f"All {len(ss)} mechanism settings (held-out arms)")
    ax.set_ylim(-0.02, 1.02)
    ax.legend(fontsize=8)
    ax.grid(alpha=0.3)

    fig.tight_layout()
    out_fig = os.path.join(ROOT, "results", "figures",
                           "Attachment_Model_vs_Experiment.png")
    os.makedirs(os.path.dirname(out_fig), exist_ok=True)
    fig.savefig(out_fig, dpi=150)
    print(f"wrote {out_fig}")


if __name__ == "__main__":
    main()
