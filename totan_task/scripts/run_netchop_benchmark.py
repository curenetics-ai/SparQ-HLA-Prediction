#!/usr/bin/env python3
"""
Run NetChop 3.1 over the CALIS peptides and write parsed per-position scores.

Week 2, antigen-processing reference benchmark.

Network selection is the critical parameter and it is documented inconsistently
upstream, so it is pinned explicitly here:

    netChop -v 0  ->  Cterm-3.0   (THE DEFAULT)
    netChop -v 1  ->  20S-3.0

Verified from the binary's own man page and its runtime banner. IEDB's bundled
`predict.py` validation message reads "'0', '1', corresponding method name:
'20S', 'C term'", which implies the opposite mapping; it passes the value
straight through to `-v`, so trusting that message silently selects C-term while
believing 20S was chosen. See docs/week2_processing_reference_decision.md.

The benchmark uses 20S. C-term is also run, so the decision can be shown rather
than asserted, and so the two can be compared in the write-up.

Outputs (one row per peptide position per network):
    data/processed/cleavage/netchop_predictions.csv
    data/processed/cleavage/netchop_run_config.json

Run from the project root (Docker must be running):
    python ./scripts/run_netchop_benchmark.py
"""

from __future__ import annotations

import csv
import json
import platform
import re
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

TASK_ROOT = Path(__file__).resolve().parents[1]   # totan_task/
REPO_ROOT = Path(__file__).resolve().parents[2]   # repository root (Week 1 inputs)
ROOT = TASK_ROOT                                  # this task's outputs
FASTA = ROOT / "data" / "interim" / "calis_peptides.fasta"
ID_MAP = ROOT / "data" / "interim" / "calis_peptide_id_map.csv"
OUT_DIR = ROOT / "data" / "processed" / "cleavage"
PRED_CSV = OUT_DIR / "netchop_predictions.csv"
CONFIG_JSON = OUT_DIR / "netchop_run_config.json"

IMAGE = "sparq/netchop:3.1"
THRESHOLD = "0.500000"
NETWORKS = {"20S": "1", "Cterm": "0"}          # name -> -v value

# "   1   M  S   0.618935 CALIS01"
# NetChop truncates the identifier to a fixed width, so anything after the first
# token of the FASTA header can spill into this column. Inputs are written with
# id-only headers; the trailing `.*` keeps the parser robust regardless.
ROW_RE = re.compile(
    r"^\s*(?P<pos>\d+)\s+(?P<aa>[A-Za-z])\s+(?P<flag>[S.])\s+"
    r"(?P<score>[0-9.]+)\s+(?P<ident>\S+).*$"
)
BANNER_RE = re.compile(r"NetChop .* predictions using version (?P<net>\S+?)\.", re.I)


def run_netchop(v_value: str) -> str:
    """Run the container and return raw stdout."""
    cmd = [
        "docker", "run", "--rm",
        "-v", f"{ROOT}:/work",
        IMAGE,
        "netChop", "-t", THRESHOLD, "-v", v_value,
        "/work/data/interim/calis_peptides.fasta",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        sys.exit(f"netChop -v {v_value} failed (exit {proc.returncode}):\n{proc.stderr}")
    return proc.stdout


def parse(raw: str, network: str) -> tuple[list[dict], str]:
    """Parse NetChop long-format output into rows; also return the banner network."""
    banner = ""
    m = BANNER_RE.search(raw)
    if m:
        banner = m.group("net")

    rows = []
    for line in raw.splitlines():
        if line.startswith("#"):
            continue
        m = ROW_RE.match(line)
        if not m:
            continue
        rows.append({
            "netchop_id": m.group("ident"),
            "position_1based": int(m.group("pos")),
            "residue": m.group("aa").upper(),
            "netchop_cleavage_flag": m.group("flag"),
            "netchop_predicted_cleaved": int(m.group("flag") == "S"),
            "netchop_score": float(m.group("score")),
            "network": network,
        })
    if not rows:
        sys.exit(f"no prediction rows parsed for network {network}; raw output:\n{raw[:800]}")
    return rows, banner


def load_id_map() -> dict[str, dict]:
    with ID_MAP.open(encoding="utf-8-sig", newline="") as fh:
        return {r["netchop_id"]: r for r in csv.DictReader(fh)}


def main() -> int:
    if not FASTA.exists():
        sys.exit(f"missing {FASTA}\nRun scripts/prepare_calis_netchop_input.py first.")
    id_map = load_id_map()
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    image_id = subprocess.run(
        ["docker", "image", "inspect", IMAGE, "--format", "{{.Id}}"],
        capture_output=True, text=True,
    ).stdout.strip()

    all_rows: list[dict] = []
    banners: dict[str, str] = {}

    for network, v_value in NETWORKS.items():
        raw = run_netchop(v_value)
        rows, banner = parse(raw, network)
        banners[network] = banner
        print(f"{network:6} (-v {v_value}): {len(rows):4} rows, banner reports {banner!r}")

        unknown = {r["netchop_id"] for r in rows} - set(id_map)
        if unknown:
            sys.exit(f"unmapped NetChop identifiers (header truncation?): {sorted(unknown)[:5]}")

        for r in rows:
            meta = id_map[r["netchop_id"]]
            r["peptide_name"] = meta["peptide_name"]
            r["sequence_length"] = int(meta["sequence_length"])
            # NetChop scores every residue 1..L as "cleavage after this residue".
            # Position L is the peptide's own C-terminus: an artefact of how the
            # peptide was synthesised, not an observed internal cut, and CALIS
            # labels only bonds 1..L-1. Flag it here; the evaluator drops it.
            r["is_internal_bond"] = int(r["position_1based"] < r["sequence_length"])
        all_rows.extend(rows)

    # Sanity: the networks must genuinely differ, or -v was ignored.
    by_net = {n: {(r["netchop_id"], r["position_1based"]): r["netchop_score"]
                  for r in all_rows if r["network"] == n} for n in NETWORKS}
    shared = set(by_net["20S"]) & set(by_net["Cterm"])
    differing = sum(1 for k in shared if by_net["20S"][k] != by_net["Cterm"][k])
    if differing == 0:
        sys.exit("20S and Cterm produced identical scores: the -v flag was ignored.")
    print(f"networks differ at {differing}/{len(shared)} positions")

    fields = ["netchop_id", "peptide_name", "network", "position_1based", "residue",
              "is_internal_bond", "netchop_score", "netchop_cleavage_flag",
              "netchop_predicted_cleaved", "sequence_length"]
    all_rows.sort(key=lambda r: (r["network"], r["netchop_id"], r["position_1based"]))
    with PRED_CSV.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows({k: r[k] for k in fields} for r in all_rows)

    internal = sum(r["is_internal_bond"] for r in all_rows if r["network"] == "20S")
    config = {
        "generated_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "tool": "NetChop 3.1 (IEDB_NetChop-3.0 distribution)",
        "container_image": IMAGE,
        "container_image_id": image_id,
        "threshold": float(THRESHOLD),
        "networks_run": {n: f"-v {v}" for n, v in NETWORKS.items()},
        "benchmark_network": "20S",
        "benchmark_network_flag": "-v 1",
        "banner_reported_network": banners,
        "v_flag_mapping_source": "netchop.1 man page and runtime banner: -v 0 = Cterm-3.0, -v 1 = 20S-3.0",
        "v_flag_mapping_warning": (
            "IEDB predict.py validation text implies 0='20S', which contradicts the "
            "binary. The binary's mapping is authoritative and is used here."
        ),
        "input_fasta": str(FASTA.relative_to(ROOT)).replace("\\", "/"),
        "n_peptides": len({r["netchop_id"] for r in all_rows}),
        "n_rows_total": len(all_rows),
        "n_internal_bonds_per_network": internal,
        "host_platform": platform.platform(),
    }
    CONFIG_JSON.write_text(json.dumps(config, indent=2), encoding="utf-8")

    print(f"wrote {PRED_CSV.relative_to(ROOT)}  ({len(all_rows)} rows)")
    print(f"wrote {CONFIG_JSON.relative_to(ROOT)}")
    print(f"internal bonds per network: {internal}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
