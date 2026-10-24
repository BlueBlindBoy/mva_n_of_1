#!/usr/bin/env python3
"""Approved-drug screen for mitotic checkpoint / spindle-apparatus proteins.

Reproducible database screen (ChEMBL REST API, no API key required) that finds,
for a 22-gene panel of mitotic-checkpoint and spindle-apparatus proteins, every
drug with a documented molecular mechanism against a panel target, flagged by
development phase (max_phase >= 4 = approved in at least one jurisdiction).

This reproduces the database-search portion of the approved-drug checkpoint
screen (original report 2026-10-23, which used DGIdb + Open Targets + FDA/EMA
labels). DGIdb/Open Targets API endpoints are not resolvable from this network,
and ChEMBL is the upstream source for the ChEMBL-derived rows those tools
aggregate, so the screen runs on ChEMBL directly. FDA/EMA label details,
approval regions, and primary-citation curation remain a manual step on top of
the raw rows this script emits.

Usage:
    python run_approved_screen.py [--out screen_results.csv] [--json rows.json]

Output:
    One CSV row per (panel gene, target, drug mechanism) record:
    panel_gene, target_chembl_id, target_name, molecule_chembl_id,
    molecule_pref_name, max_phase, approved (1/0), mechanism_of_action,
    action_type, direct_interaction, chembl_version.
    The CSV includes non-approved rows (research-stage compounds) so the
    "no approved modulator" conclusion per target is auditable.

Interpretation notes:
- A drug mapped to a panel gene can be an off-target hit: the
  mechanism_of_action column shows the drug's primary documented mechanism.
  The curated screen counts a drug as an approved direct modulator of a panel
  target only when that panel target IS the drug's approved mechanism.
- ADCs: the payload is the direct target; ChEMBL annotations vary. The
  curated screen applies the inclusion rule (payload's direct target must be
  a panel protein) by hand.
"""
import argparse
import csv
import json
import sys
import time
import urllib.parse
import urllib.request

BASE = "https://www.ebi.ac.uk/chembl/api/data/"

# Panel gene -> ChEMBL target synonyms to try (first hit kept per target id).
# APC11/ANAPC11: the Cul1-family scaffold subunit of APC/C (NOT CUL1).
# PDS5A/PDS5B are cohesin regulators (NOT securin); PTTG1 is human securin.
PANEL = {
    "BUB1B": ["BUB1B", "BUBR1"],
    "BUB3": ["BUB3"],
    "MAD2L1": ["MAD2L1", "MAD2"],
    "TRIP13": ["TRIP13"],
    "ANAPC1": ["ANAPC1"],
    "ANAPC2": ["ANAPC2"],
    "APC11": ["APC11", "ANAPC11"],
    "CDC20": ["CDC20"],
    "CDC25C": ["CDC25C"],
    "UBE2C": ["UBE2C"],
    "PDS5A": ["PDS5A"],
    "PDS5B": ["PDS5B"],
    "PTTG1": ["PTTG1"],
    "CDK1": ["CDK1"],
    "BUB1": ["BUB1"],
    "AURKB": ["AURKB"],
    "AURKA": ["AURKA"],
    "PLK1": ["PLK1"],
    "KIF11": ["KIF11"],
    "TUBB": ["TUBB"],
    "CEP57": ["CEP57"],
    "CDCA7": ["CDCA7"],
}


def api(params, item_key, timeout=60, retries=4):
    """GET one ChEMBL API resource, retrying on errors and on empty-but-
    non-zero results (observed intermittent server behavior)."""
    url = BASE + params[0] + "?" + urllib.parse.urlencode(params[1])
    last = None
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "mva-screen/1.0 (research)"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                d = json.loads(r.read().decode())
            total = d.get("page_meta", {}).get("total_count", 0)
            items = d.get(item_key, [])
            if len(items) == 0 and total > 0:
                last = f"empty body vs total_count={total}"
                time.sleep(2 * (attempt + 1))
                continue
            return items, total
        except Exception as e:  # HTTP/URL/JSON errors
            last = f"{type(e).__name__}: {e}"
            time.sleep(2 * (attempt + 1))
    raise RuntimeError(f"ChEMBL query failed after {retries} attempts: {url} ({last})")


def chembl_version():
    try:
        req = urllib.request.Request(BASE + "status", headers={"User-Agent": "mva-screen/1.0 (research)"})
        with urllib.request.urlopen(req, timeout=30) as r:
            import xml.etree.ElementTree as ET
            root = ET.fromstring(r.read())
            return root.findtext("chembl_db_version", "?")
    except Exception:
        return "?"


def main():
    ap = argparse.ArgumentParser(
        description="Approved-drug screen for the 22-gene mitotic checkpoint panel (ChEMBL).")
    ap.add_argument("--out", default="approved_screen_results.csv")
    ap.add_argument("--json", dest="json_out", default="approved_screen_rows.json")
    args = ap.parse_args()

    version = chembl_version()
    print(f"ChEMBL version: {version}", flush=True)

    rows = []
    target_map = {}
    for gene, syns in PANEL.items():
        tids, seen = [], set()
        for s in syns:
            items, _ = api(("target.json", {"target_synonym__iexact": s, "limit": "20"}), "targets")
            for t in items:
                if t["target_chembl_id"] not in seen:
                    seen.add(t["target_chembl_id"])
                    tids.append((t["target_chembl_id"], t["pref_name"]))
            time.sleep(0.4)
        target_map[gene] = tids
        print(f"{gene:8s} -> {len(tids)} target(s): {[t for t, _ in tids]}", flush=True)
        for tid, tname in tids:
            items, _ = api(("mechanism.json", {"target_chembl_id": tid, "limit": "200"}), "mechanisms")
            for r in items:
                rows.append({
                    "panel_gene": gene,
                    "target_chembl_id": tid,
                    "target_name": tname,
                    "molecule_chembl_id": r.get("molecule_chembl_id"),
                    "max_phase": r.get("max_phase"),
                    "approved": 1 if (r.get("max_phase") or 0) >= 4 else 0,
                    "mechanism_of_action": r.get("mechanism_of_action"),
                    "action_type": r.get("action_type"),
                    "direct_interaction": r.get("direct_interaction"),
                })
            time.sleep(0.4)
        print(f"{gene:8s}    mechanism rows so far: {len(rows)}", flush=True)

    # resolve molecule names in batches
    mols = sorted({r["molecule_chembl_id"] for r in rows if r["molecule_chembl_id"]})
    names = {}
    for i in range(0, len(mols), 50):
        chunk = mols[i:i + 50]
        items, _ = api(("molecule.json", {"molecule_chembl_id__in": ",".join(chunk), "limit": "200"}), "molecules")
        for x in items:
            names[x["molecule_chembl_id"]] = x.get("pref_name") or x.get("molecule_id")
        time.sleep(0.4)

    for r in rows:
        r["molecule_pref_name"] = names.get(r["molecule_chembl_id"], "")
        r["chembl_version"] = version

    cols = ["panel_gene", "target_chembl_id", "target_name", "molecule_chembl_id",
            "molecule_pref_name", "max_phase", "approved", "mechanism_of_action",
            "action_type", "direct_interaction", "chembl_version"]
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    with open(args.json_out, "w") as f:
        json.dump({"chembl_version": version, "targets": target_map, "rows": rows}, f, indent=1)

    print(f"\nTotal mechanism rows: {len(rows)} -> {args.out}", flush=True)
    print("Per-gene summary (approved = max_phase >= 4):")
    for gene in PANEL:
        rs = [r for r in rows if r["panel_gene"] == gene]
        appr = sorted({r["molecule_pref_name"] or r["molecule_chembl_id"]
                       for r in rs if r["approved"]})
        print(f"  {gene:8s} approved={len(appr):2d}" + (f"  {appr}" if appr else ""))
    print("SCREEN_DONE")


if __name__ == "__main__":
    sys.exit(main())
