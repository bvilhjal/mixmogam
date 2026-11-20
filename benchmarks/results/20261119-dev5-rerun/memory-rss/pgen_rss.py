"""Measure the same 500,000-sample step-2 scan in a fresh process per reader."""
import argparse
import json
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

ROOT = Path(__file__).resolve().parents[1]
PANEL = ROOT / "scalability"


def scan(mode):
    from mixmogam import HrattStep1, hratt_step2
    from mixmogam.io.pgen import iter_pgen

    fit = HrattStep1.load(PANEL / "highn-step1.npz")
    source = PANEL / "highn"
    chunk_bytes = 16 * 1024**2
    if mode != "native":
        source = iter_pgen(source, packed=mode == "packed", chunk_bytes=chunk_bytes)
    with TemporaryDirectory() as temp:
        result = hratt_step2(fit, source, chunk_bytes=chunk_bytes, out=Path(temp) / "scores.tsv")
        assert result["variants"] == 1024


def main():
    from _harness import digest, environment_manifest, measure, power_state, save_json, write_csv

    out = Path(__file__).resolve().parent
    archive = json.loads((PANEL / "manifest.json").read_text())
    hashes = {name: sha for name, sha in archive["input_sha256"].items() if name.startswith("highn")}
    assert len(hashes) == 4
    assert all(digest(PANEL / name) == sha for name, sha in hashes.items())
    save_json(out / "manifest.json", {
        "commit": (ROOT / "commit.txt").read_text().strip(),
        "environment": environment_manifest(),
        "inputs": hashes,
        "driver_sha256": digest(__file__),
        "power_state": power_state(),
        "scope": "One trait, saved HRATT step 1, step 2 writes TSV. RSS includes Python startup and input load.",
    })
    modes = ["native", "packed", "int8"]
    rows = []
    for rep in range(-1, 3):
        for j in range(3):
            mode = modes[(rep + j) % 3]
            label = f"{mode}-{rep + 1}"
            command = [sys.executable, str(Path(__file__).resolve()), "--scan", mode]
            row = measure(command, out / label, max_load=5)
            if rep >= 0:
                rows.append({"reader": mode, "rep": rep + 1, "peak_rss_bytes": row["peak_rss_bytes"],
                             "load_average": row["load_average_before"][0]})
                write_csv(out / "pgen-rss.csv", rows)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scan", choices=("native", "packed", "int8"))
    args = parser.parse_args()
    scan(args.scan) if args.scan else main()
