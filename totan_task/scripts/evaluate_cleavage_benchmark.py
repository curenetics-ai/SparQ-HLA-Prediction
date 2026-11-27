#!/usr/bin/env python3
"""
Evaluate NetChop against the CALIS 2015 cleavage labels, stratified by
proteasome type.

Week 2, antigen-processing reference benchmark.

Three things this script is careful about, because each would otherwise produce
numbers that look fine and are wrong:

1. POSITION ALIGNMENT. NetChop scores every residue 1..L as "cleavage after this
   residue". CALIS annotates internal bonds 1..L-1 only. Position L is the
   peptide's own C-terminus - an artefact of synthesis, not an observed cut - so
   it is dropped. The join is then asserted to be exact: every label row must
   match exactly one prediction row, and the totals must agree.

2. I/C IS A PROPERTY OF THE LABELS, NOT THE PREDICTOR. NetChop has no
   immunoproteasome model; only Cterm-3.0 and 20S-3.0 exist, neither of which is
   immunoproteasome-specific. One prediction set is therefore evaluated twice,
   against the I ground truth and against the C ground truth. Any claim that
   NetChop "predicts immunoproteasome cleavage" would be fabricated.

3. UNCERTAINTY IS CLUSTERED. 368 positions come from only 16 peptides, and
   positions within a peptide are not independent. Confidence intervals use a
   cluster bootstrap resampling whole peptides, which is wider and honest, not
   a naive per-position bootstrap which would overstate precision.

Outputs:
    results/tables/netchop_cleavage_benchmark_by_proteasome_type.csv
    results/tables/netchop_cleavage_benchmark_by_peptide.csv
    results/tables/netchop_cleavage_joined_sites.csv

Run from the project root:
    python ./scripts/evaluate_cleavage_benchmark.py
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import numpy as np

TASK_ROOT = Path(__file__).resolve().parents[1]   # totan_task/
REPO_ROOT = Path(__file__).resolve().parents[2]   # repository root (Week 1 inputs)
ROOT = TASK_ROOT                                  # this task's outputs
LABELS = REPO_ROOT / "data" / "processed" / "cleavage" / "calis_2015_cleavage_site_labels.csv"
PREDS = ROOT / "data" / "processed" / "cleavage" / "netchop_predictions.csv"
ID_MAP = ROOT / "data" / "interim" / "calis_peptide_id_map.csv"
TABLES = ROOT / "results" / "tables"

THRESHOLD = 0.5
N_BOOT = 2000
SEED = 20261125


# ---------------------------------------------------------------- metrics ---

def auc(y: np.ndarray, s: np.ndarray) -> float:
    """ROC AUC via the rank (Mann-Whitney U) identity, ties averaged."""
    pos, neg = int(y.sum()), int((1 - y).sum())
    if pos == 0 or neg == 0:
        return float("nan")
    order = np.argsort(s, kind="mergesort")
    ranks = np.empty(len(s), dtype=float)
    ranks[order] = np.arange(1, len(s) + 1, dtype=float)
    # average ranks within ties
    srt = s[order]
    i = 0
    while i < len(srt):
        j = i
        while j + 1 < len(srt) and srt[j + 1] == srt[i]:
            j += 1
        if j > i:
            ranks[order[i:j + 1]] = ranks[order[i:j + 1]].mean()
        i = j + 1
    return (ranks[y == 1].sum() - pos * (pos + 1) / 2) / (pos * neg)


def confusion(y: np.ndarray, pred: np.ndarray) -> tuple[int, int, int, int]:
    tp = int(((y == 1) & (pred == 1)).sum())
    fp = int(((y == 0) & (pred == 1)).sum())
    tn = int(((y == 0) & (pred == 0)).sum())
    fn = int(((y == 1) & (pred == 0)).sum())
    return tp, fp, tn, fn


def mcc(tp: int, fp: int, tn: int, fn: int) -> float:
    denom = np.sqrt(float(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return float("nan") if denom == 0 else (tp * tn - fp * fn) / denom


def metrics(y: np.ndarray, s: np.ndarray) -> dict:
    pred = (s >= THRESHOLD).astype(int)
    tp, fp, tn, fn = confusion(y, pred)
    safe = lambda a, b: float("nan") if b == 0 else a / b
    return {
        "n_sites": len(y),
        "n_cleaved": int(y.sum()),
        "n_not_cleaved": int((1 - y).sum()),
        "prevalence": safe(int(y.sum()), len(y)),
        "auc": auc(y, s),
        "tp": tp, "fp": fp, "tn": tn, "fn": fn,
        "sensitivity": safe(tp, tp + fn),
        "specificity": safe(tn, tn + fp),
        "ppv": safe(tp, tp + fp),
        "npv": safe(tn, tn + fn),
        "accuracy": safe(tp + tn, len(y)),
        "f1": safe(2 * tp, 2 * tp + fp + fn),
        "mcc": mcc(tp, fp, tn, fn),
    }


def cluster_bootstrap_ci(y, s, groups, stat, n_boot=N_BOOT, seed=SEED):
    """Percentile CI, resampling whole peptides (clusters), not positions."""
    rng = np.random.default_rng(seed)
    uniq = np.unique(groups)
    vals = []
    for _ in range(n_boot):
        pick = rng.choice(uniq, size=len(uniq), replace=True)
        idx = np.concatenate([np.flatnonzero(groups == g) for g in pick])
        yb, sb = y[idx], s[idx]
        if yb.sum() == 0 or (1 - yb).sum() == 0:
            continue
        v = stat(yb, sb)
        if np.isfinite(v):
            vals.append(v)
    if len(vals) < 50:
        return float("nan"), float("nan"), len(vals)
    return (float(np.percentile(vals, 2.5)),
            float(np.percentile(vals, 97.5)), len(vals))


# ------------------------------------------------------------------- data ---

def read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def build_joined() -> list[dict]:
    for p in (LABELS, PREDS, ID_MAP):
        if not p.exists():
            sys.exit(f"missing input: {p}")

    name_to_id = {r["peptide_name"]: r["netchop_id"] for r in read_csv(ID_MAP)}

    # prediction lookup: (network, netchop_id, position) -> row, internal bonds only
    pred_idx: dict[tuple, dict] = {}
    dropped_terminal = 0
    for r in read_csv(PREDS):
        if r["is_internal_bond"] != "1":
            dropped_terminal += 1
            continue
        pred_idx[(r["network"], r["netchop_id"], int(r["position_1based"]))] = r
    print(f"dropped {dropped_terminal} terminal (position L) prediction rows")

    networks = sorted({k[0] for k in pred_idx})
    joined: list[dict] = []
    for lab in read_csv(LABELS):
        nid = name_to_id.get(lab["peptide_name"])
        if nid is None:
            sys.exit(f"label peptide not in id map: {lab['peptide_name']!r}")
        pos = int(lab["candidate_cleavage_after_residue_position_1based"])
        for net in networks:
            pr = pred_idx.get((net, nid, pos))
            if pr is None:
                sys.exit(f"no prediction for {net}/{lab['peptide_name']}/pos {pos}")
            if pr["residue"] != lab["left_residue"].upper():
                sys.exit(
                    f"residue mismatch {net}/{lab['peptide_name']} pos {pos}: "
                    f"NetChop {pr['residue']} vs CALIS {lab['left_residue']} "
                    "- position convention is misaligned"
                )
            joined.append({
                "peptide_name": lab["peptide_name"],
                "netchop_id": nid,
                "proteasome_type": lab["proteasome_type"],
                "proteasome_code": lab["proteasome_code"],
                "position_1based": pos,
                "left_residue": lab["left_residue"],
                "right_residue": lab["right_residue"],
                "network": net,
                "netchop_score": float(pr["netchop_score"]),
                "netchop_predicted_cleaved": int(pr["netchop_predicted_cleaved"]),
                "observed_cleaved": int(lab["is_cleaved"]),
            })

    expected = len(read_csv(LABELS)) * len(networks)
    if len(joined) != expected:
        sys.exit(f"join produced {len(joined)} rows, expected {expected}")
    print(f"joined {len(joined)} rows ({len(networks)} networks x {expected // len(networks)} labels)")
    return joined


# ------------------------------------------------------------------- main ---

def main() -> int:
    TABLES.mkdir(parents=True, exist_ok=True)
    joined = build_joined()

    xmod = {r["peptide_name"] for r in joined
            if any(c not in "ACDEFGHIKLMNPQRSTVWY" for c in r["left_residue"])}

    with (TABLES / "netchop_cleavage_joined_sites.csv").open(
            "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(joined[0].keys()))
        w.writeheader()
        w.writerows(joined)

    networks = sorted({r["network"] for r in joined})
    types = sorted({r["proteasome_type"] for r in joined})

    summary_rows = []
    for net in networks:
        for subset_name, keep in (("all_peptides", None),
                                  ("excluding_nonstandard_residue_peptides", xmod)):
            for ptype in types:
                sel = [r for r in joined
                       if r["network"] == net and r["proteasome_type"] == ptype
                       and (keep is None or r["peptide_name"] not in keep)]
                if not sel:
                    continue
                y = np.array([r["observed_cleaved"] for r in sel])
                s = np.array([r["netchop_score"] for r in sel])
                g = np.array([r["peptide_name"] for r in sel])
                m = metrics(y, s)
                lo, hi, nb = cluster_bootstrap_ci(y, s, g, auc)
                mlo, mhi, _ = cluster_bootstrap_ci(
                    y, s, g,
                    lambda yy, ss: mcc(*confusion(yy, (ss >= THRESHOLD).astype(int))))
                summary_rows.append({
                    "network": net,
                    "subset": subset_name,
                    "proteasome_type": ptype,
                    "n_peptides": len({r["peptide_name"] for r in sel}),
                    **m,
                    "auc_ci95_low": lo, "auc_ci95_high": hi,
                    "mcc_ci95_low": mlo, "mcc_ci95_high": mhi,
                    "bootstrap_replicates": nb,
                    "threshold": THRESHOLD,
                })

    fields = list(summary_rows[0].keys())
    out = TABLES / "netchop_cleavage_benchmark_by_proteasome_type.csv"
    with out.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        for r in summary_rows:
            w.writerow({k: (f"{v:.4f}" if isinstance(v, float) and np.isfinite(v)
                            else ("" if isinstance(v, float) else v))
                        for k, v in r.items()})

    # per-peptide, benchmark network only - shows how much the small sample moves
    per_pep = []
    for ptype in types:
        for pep in sorted({r["peptide_name"] for r in joined}):
            sel = [r for r in joined if r["network"] == "20S"
                   and r["proteasome_type"] == ptype and r["peptide_name"] == pep]
            if not sel:
                continue
            y = np.array([r["observed_cleaved"] for r in sel])
            s = np.array([r["netchop_score"] for r in sel])
            m = metrics(y, s)
            per_pep.append({"network": "20S", "proteasome_type": ptype,
                            "peptide_name": pep, **m})
    with (TABLES / "netchop_cleavage_benchmark_by_peptide.csv").open(
            "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(per_pep[0].keys()))
        w.writeheader()
        for r in per_pep:
            w.writerow({k: (f"{v:.4f}" if isinstance(v, float) and np.isfinite(v)
                            else ("" if isinstance(v, float) else v))
                        for k, v in r.items()})

    print(f"\nwrote {out.relative_to(ROOT)}")
    print("\n=== headline (all peptides, threshold 0.5) ===")
    hdr = f"{'network':7} {'proteasome':24} {'n':>4} {'pos':>4} {'AUC':>6} {'95% CI':>15} {'MCC':>7} {'sens':>6} {'spec':>6}"
    print(hdr)
    for r in summary_rows:
        if r["subset"] != "all_peptides":
            continue
        print(f"{r['network']:7} {r['proteasome_type']:24} {r['n_sites']:4d} "
              f"{r['n_cleaved']:4d} {r['auc']:6.3f} "
              f"[{r['auc_ci95_low']:.3f},{r['auc_ci95_high']:.3f}]".ljust(0) +
              f" {r['mcc']:7.3f} {r['sensitivity']:6.3f} {r['specificity']:6.3f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
