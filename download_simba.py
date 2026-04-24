"""
download_simba.py — download CAMELS SIMBA group catalog HDF5 files via HTTPS.

CAMELS-SIMBA CV uses the same 91-snapshot scheme as TNG (000–090).
Confirmed snap numbers (same as config.py SIMBA_SNAP_EARLY/LATE):
    snap 066 = z=0.7747  (early epoch)
    snap 090 = z=0.0000  (z=0 final)

Prerequisites:
    python get_camels_token.py  (run once to get the token)

Usage:
    # Download snaps 066+090 for all 27 CV sims:
    python download_simba.py --sims $(python -c "print(' '.join(f'CV_{i}' for i in range(27)))")

    # Or just one sim for testing:
    python download_simba.py --sims CV_0

    # Probe snap→redshift table (slow, downloads all snaps for one sim):
    python download_simba.py --probe --sims CV_0

File layout:
    <out-dir>/<sim_id>/groups_<snap>.hdf5
"""
import argparse
import pathlib
import sys
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed

import requests

HTTPS_BASE   = "https://g-a77640.2c3d02.75bc.data.globus.org"
SIMBA_ROOT   = "/FOF_Subfind/SIMBA/CV"
TOKEN_FILE   = pathlib.Path.home() / ".globus" / "camels_https_token.txt"

# Confirmed snap numbers for CAMELS-SIMBA CV (same 91-snap scheme as TNG, 000–090).
# Verified from groups_066.hdf5 and groups_090.hdf5 Header/Redshift for CV_0.
#   snap 090 = z=0.0000  (confirmed)
#   snap 066 = z=0.7747  (confirmed)
SNAP_LATE_DEFAULT  = 90
SNAP_EARLY_DEFAULT = 66

_print_lock = threading.Lock()


def get_token() -> str:
    if not TOKEN_FILE.exists():
        sys.exit(f"Token file not found: {TOKEN_FILE}\nRun:  python get_camels_token.py")
    return TOKEN_FILE.read_text().strip()


def download_file(url: str, dest: pathlib.Path, token: str, max_retries: int = 5) -> str:
    """Download one file. Returns status string."""
    dest.parent.mkdir(parents=True, exist_ok=True)
    headers = {"Authorization": f"Bearer {token}"}

    try:
        head = requests.head(url, headers=headers, timeout=30)
        head.raise_for_status()
    except requests.HTTPError as e:
        return f"  FAIL  {dest.parent.name}/{dest.name}  ({e})"

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


def probe_snaps(token: str, sim_id: str = "CV_0", out_dir: pathlib.Path = pathlib.Path("outputs/cache/simba")) -> None:
    """
    Download snap 033 (z=0) and all snaps for one sim to map snap → redshift.
    Reads Header/Redshift from each downloaded HDF5 to build the table.

    Run once to find SNAP_EARLY for your desired target redshift.
    """
    try:
        import h5py
    except ImportError:
        sys.exit("h5py required: pip install h5py")

    print(f"\nProbing SIMBA snap schedule for {sim_id} ...")
    print(f"{'Snap':>6}  {'Redshift':>10}  {'Status'}")
    print("-" * 36)

    for snap in range(91):
        fname = f"groups_{snap:03d}.hdf5"
        url   = f"{HTTPS_BASE}{SIMBA_ROOT}/{sim_id}/{fname}"
        dest  = out_dir / sim_id / fname

        if not dest.exists():
            msg = download_file(url, dest, token)
            ok  = "done" in msg
        else:
            ok = True

        if ok and dest.exists():
            try:
                with h5py.File(dest, "r") as f:
                    z = float(f["Header"].attrs["Redshift"])
                print(f"  {snap:3d}    {z:10.4f}")
            except Exception as e:
                print(f"  {snap:3d}    <read error: {e}>")
        else:
            print(f"  {snap:3d}    <download failed>")


def main():
    ap = argparse.ArgumentParser(description="Download CAMELS-SIMBA CV group catalogs")
    ap.add_argument("--sims", nargs="+", default=["CV_0"])
    ap.add_argument("--out-dir",     default="outputs/cache/simba")
    ap.add_argument("--workers",     type=int, default=6)
    ap.add_argument("--snap-early",  type=int, default=SNAP_EARLY_DEFAULT,
                    help=f"Early snapshot number (default: {SNAP_EARLY_DEFAULT}; verify with --probe)")
    ap.add_argument("--snap-late",   type=int, default=SNAP_LATE_DEFAULT,
                    help=f"Late snapshot number (default: {SNAP_LATE_DEFAULT})")
    ap.add_argument("--probe",       action="store_true",
                    help="Download all snaps for first sim and print snap→z table")
    args = ap.parse_args()

    token   = get_token()
    out_dir = pathlib.Path(args.out_dir)

    if args.probe:
        probe_snaps(token, sim_id=args.sims[0], out_dir=out_dir)
        print(f"\nSIMBA uses same 91-snap scheme as TNG: snap 066=z≈0.77, snap 090=z=0")
        return

    tasks = []
    for sim_id in args.sims:
        for snap in (args.snap_early, args.snap_late):
            fname = f"groups_{snap:03d}.hdf5"
            url   = f"{HTTPS_BASE}{SIMBA_ROOT}/{sim_id}/{fname}"
            dest  = out_dir / sim_id / fname
            tasks.append((url, dest))

    n_total = len(tasks)
    n_done  = 0
    print(f"Downloading {n_total} SIMBA files ({args.workers} parallel connections) ...")
    print(f"  snap_early={args.snap_early}  snap_late={args.snap_late}")

    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(download_file, url, dest, token): (url, dest)
                   for url, dest in tasks}
        for fut in as_completed(futures):
            n_done += 1
            msg = fut.result()
            with _print_lock:
                print(f"[{n_done:3d}/{n_total}] {msg}")

    print(f"\nDone.  Run the SIMBA analysis with:")
    print(f"  python paper.py --suite simba --data-dir {args.out_dir} "
          f"--snap-early {args.snap_early} --snap-late {args.snap_late} "
          f"--run-label simba_CV --matching spatial")


if __name__ == "__main__":
    main()
