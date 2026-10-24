#!/usr/bin/env python3
"""Local Pangolin (Zeng & Li, Genome Biology 2022) splice-site scoring, run
on CPU over public GRCh38 sequence only.

For each variant listed in variants.json (10101 bp window, variant at base
5051, supplied as <name>_ref.fa next to this script):
  - score P(splice) at every position for REF and ALT sequence
  - report the largest increase/decrease within +/-50 bp of the variant
    (mirrors the Pangolin CLI default -d 50)

Models: heart P(splice) (num 0) and brain P(splice) (num 4), each averaged
over 5 independently trained checkpoints.

The output JSON contains scores only (variant name + tissue scores). It does
not record coordinates, which stay local-only in variants.json.
"""
import os, sys, json, numpy as np, torch
import pangolin as _pkg
from pangolin.model import Pangolin, L, W, AR
MODEL_DIR = os.path.join(os.path.dirname(_pkg.__file__), "models")

IN_MAP = np.asarray([[0, 0, 0, 0], [1, 0, 0, 0], [0, 1, 0, 0],
                     [0, 0, 1, 0], [0, 0, 0, 1]])
INDEX_MAP = {0: 1, 1: 2, 2: 4, 3: 5, 4: 7, 5: 8, 6: 10, 7: 11}

HERE = os.path.dirname(os.path.abspath(__file__))
MODEL_NUMS = [0, 4]  # heart P(splice), brain P(splice)
TISSUE = {0: "heart", 4: "brain"}


def one_hot_encode(seq, strand):
    seq = seq.upper().replace('A', '1').replace('C', '2')
    seq = seq.replace('G', '3').replace('T', '4').replace('N', '0')
    if strand == '+':
        arr = np.asarray(list(map(int, list(seq))))
    else:
        seq = seq[::-1]
        arr = np.asarray(list(map(int, list(seq))))
        arr = (5 - arr) % 5  # reverse complement
    return IN_MAP[arr.astype('int8')]


def load_models(nums):
    models = {}
    for i in nums:
        ms = []
        for j in range(1, 6):
            m = Pangolin(L, W, AR)
            w = torch.load(os.path.join(MODEL_DIR, f"final.{j}.{i}.3"),
                           map_location=torch.device('cpu'))
            m.load_state_dict(w)
            m.eval()
            ms.append(m)
        models[i] = ms
    return models


def score_profile(models, seq, strand, model_num):
    x = one_hot_encode(seq, strand).T
    x = torch.from_numpy(np.expand_dims(x, axis=0)).float()
    out = []
    for m in models[model_num]:
        with torch.no_grad():
            out.append(m(x)[0][INDEX_MAP[model_num], :].numpy())
    return np.mean(out, axis=0)


def main():
    cfg_path = os.path.join(HERE, "variants.json")
    if not os.path.exists(cfg_path):
        sys.exit("variants.json not found. See variants.example.json for the "
                 "format. It is intentionally not committed: it holds "
                 "patient-specific variant coordinates.")
    with open(cfg_path) as f:
        cfg = json.load(f)
    strand = cfg.get("strand", "+")
    models = load_models(MODEL_NUMS)
    results = {"strand": strand, "tissues": TISSUE, "variants": []}
    for v in cfg["variants"]:
        name, ref, alt = v["name"], v["ref"], v["alt"]
        seq = open(os.path.join(HERE, f"{name}_ref.fa")).read().split("\n")[1].strip()
        assert seq[5050] == ref, f"{name}: center base {seq[5050]} != ref {ref}"
        altseq = seq[:5050] + alt + seq[5051:]
        vres = {"name": name}
        for mn in MODEL_NUMS:
            prof_ref = score_profile(models, seq, strand, mn)
            prof_alt = score_profile(models, altseq, strand, mn)
            d = prof_alt - prof_ref
            # profile length 101; the variant is at profile index 50.
            prof_len = prof_ref.size
            assert prof_len == 101, f"profile length {prof_len}"
            win = d[0:101]
            gi, di = int(np.argmax(win)), int(np.argmin(win))
            vres.setdefault("scores", {})[TISSUE[mn]] = {
                "max_gain_in_50bp": float(win[gi]),
                "pos_gain_rel_variant": int(gi - 50),
                "max_loss_in_50bp": float(win[di]),
                "pos_loss_rel_variant": int(di - 50),
                "ref_score_at_variant": float(prof_ref[50]),
                "alt_score_at_variant": float(prof_alt[50]),
            }
        results["variants"].append(vres)
    out_path = os.path.join(HERE, "pangolin_results.json")
    with open(out_path, "w") as f:
        json.dump(results, f, indent=1)
    print(json.dumps(results, indent=1))
    print("PANGOLIN_DONE")


if __name__ == "__main__":
    main()
