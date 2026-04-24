"""
download_camels.py — download CAMELS group catalog HDF5 files via HTTPS.

Prerequisites:
    python get_camels_token.py  (run once to get the token)

Usage:
    python download_camels.py [--sims CV_0 CV_1 ...] [--out-dir outputs/cache/camels]
    python download_camels.py --workers 8   # parallel connections (default: 6)

Downloads groups_066.hdf5 and groups_090.hdf5 for each requested simulation.
  snap 066 ≈ z=0.52 (early; descriptors measured here)
  snap 090 = z=0    (late;  targets measured here)

File layout:
    <out-dir>/<sim_id>/groups_066.hdf5
    <out-dir>/<sim_id>/groups_090.hdf5
"""
import argparse
import pathlib
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests

HTTPS_BASE  = "https://g-a77640.2c3d02.75bc.data.globus.org"
CAMELS_ROOT = "/FOF_Subfind/IllustrisTNG/CV"
TOKEN_FILE  = pathlib.Path.home() / ".globus" / "camels_https_token.txt"

SNAP_EARLY  = 66
SNAP_LATE   = 90

_print_lock = threading.Lock()


def get_token() -> str:
    if not TOKEN_FILE.exists():
        sys.exit(f"Token file not found: {TOKEN_FILE}\nRun:  python get_camels_token.py")
    return TOKEN_FILE.read_text().strip()


def download_file(url: str, dest: pathlib.Path, token: str, max_retries: int = 5) -> str:
    """Download one file. Returns a status string for printing."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    headers = {"Authorization": f"Bearer {token}"}

    head = requests.head(url, headers=headers, timeout=30)
    head.raise_for_status()
    expected = int(head.headers.get("content-length", 0))

    if dest.exists() and dest.stat().st_size == expected:
        return f"  skip  {dest.parent.name}/{dest.name}"

    if dest.exists():
        dest.unlink()

    for attempt in range(max_retries):
        tmp = dest.with_suffix(".part")
        try:
            with requests.get(url, headers=headers, stream=True, timeout=120) as r:
                r.raise_for_status()
                with open(tmp, "wb") as f:
                    for chunk in r.iter_content(chunk_size=1 << 20):
                        f.write(chunk)
            tmp.rename(dest)
            return f"  done  {dest.parent.name}/{dest.name}  ({expected/1e6:.1f} MB)"
        except Exception as e:
            if tmp.exists():
                tmp.unlink()
            if attempt == max_retries - 1:
                return f"  FAIL  {dest.parent.name}/{dest.name}  ({e})"

    return f"  FAIL  {dest.parent.name}/{dest.name}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sims", nargs="+", default=["CV_0"])
    ap.add_argument("--out-dir", default="outputs/cache/camels")
    ap.add_argument("--workers", type=int, default=6,
                    help="parallel download connections (default: 6)")
    args = ap.parse_args()

    token   = get_token()
    out_dir = pathlib.Path(args.out_dir)

    # Build full task list: (url, dest) for every (sim, snap) pair
    tasks = []
    for sim_id in args.sims:
        for snap in (SNAP_EARLY, SNAP_LATE):
            fname = f"groups_{snap:03d}.hdf5"
            url   = f"{HTTPS_BASE}{CAMELS_ROOT}/{sim_id}/{fname}"
            dest  = out_dir / sim_id / fname
            tasks.append((url, dest))

    n_total = len(tasks)
    n_done  = 0
    print(f"Downloading {n_total} files across {len(args.sims)} sims "
          f"({args.workers} parallel connections) ...")

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(download_file, url, dest, token): (url, dest)
                   for url, dest in tasks}
        for fut in as_completed(futures):
            n_done += 1
            msg = fut.result()
            with _print_lock:
                print(f"[{n_done:3d}/{n_total}] {msg}")

    print(f"\nDone.  Run:  python paper.py --sim-ids {' '.join(args.sims)} "
          f"--data-dir {args.out_dir}")


if __name__ == "__main__":
    main()
