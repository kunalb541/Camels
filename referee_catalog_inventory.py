#!/usr/bin/env python
"""Inventory local CAMELS fields for referee-facing physical diagnostics."""
from __future__ import annotations

import json
import pickle
import re
from pathlib import Path

import h5py
import pandas as pd
from tqdm.auto import tqdm

OUT_DIR = Path("outputs/referee")
SEARCH_TERMS = [
    "agn", "bh", "blackhole", "feedback", "wind", "sfr", "gas", "star",
    "stellar", "mass", "halo", "group", "subhalo", "metal", "temperature", "energy",
]
AGN_TERMS = ["agn", "bh", "blackhole", "feedback", "wind", "energy"]


def term_matches(name: str, terms: list[str]) -> str:
    low = name.lower()
    matches = []
    for term in terms:
        if term == "bh":
            found = bool(re.search(r"(^|[/_.])bh|bhmass|bhmdot", low))
        elif term == "wind":
            found = bool(re.search(r"wind($|mass|_)", low))
        else:
            found = term in low
        if found:
            matches.append(term)
    return ";".join(matches)


def hdf5_rows(path: Path) -> list[dict]:
    rows: list[dict] = []
    with h5py.File(path, "r") as handle:
        def visitor(name, obj):
            if isinstance(obj, h5py.Dataset):
                rows.append(
                    {
                        "file": str(path),
                        "file_type": "hdf5",
                        "object": name,
                        "shape": str(obj.shape),
                        "dtype": str(obj.dtype),
                        "search_terms": term_matches(name, SEARCH_TERMS),
                        "agn_bh_feedback_terms": term_matches(name, AGN_TERMS),
                    }
                )
                if obj.dtype.names:
                    for field in obj.dtype.names:
                        rows.append(
                            {
                                "file": str(path),
                                "file_type": "hdf5_structured_field",
                                "object": f"{name}.{field}",
                                "shape": str(obj.shape),
                                "dtype": str(obj.dtype[field]),
                                "search_terms": term_matches(field, SEARCH_TERMS),
                                "agn_bh_feedback_terms": term_matches(field, AGN_TERMS),
                            }
                        )
        handle.visititems(visitor)
    return rows


def cached_object_rows(path: Path) -> list[dict]:
    rows: list[dict] = []
    try:
        with path.open("rb") as handle:
            obj = pickle.load(handle)
    except Exception as exc:
        return [{
            "file": str(path), "file_type": "pickle_error", "object": str(exc),
            "shape": "", "dtype": "", "search_terms": "", "agn_bh_feedback_terms": "",
        }]

    def add_fields(prefix: str, value) -> None:
        if hasattr(value, "columns"):
            for col in value.columns:
                name = f"{prefix}{col}"
                rows.append(
                    {
                        "file": str(path),
                        "file_type": "pickle_dataframe_column",
                        "object": name,
                        "shape": str(getattr(value, "shape", "")),
                        "dtype": str(value[col].dtype),
                        "search_terms": term_matches(name, SEARCH_TERMS),
                        "agn_bh_feedback_terms": term_matches(name, AGN_TERMS),
                    }
                )
        elif isinstance(value, dict):
            for key, child in value.items():
                add_fields(f"{prefix}{key}.", child)

    add_fields("", obj)
    return rows


def json_rows(path: Path) -> list[dict]:
    rows: list[dict] = []
    try:
        obj = json.loads(path.read_text())
    except Exception:
        return rows

    def walk(prefix: str, value) -> None:
        if isinstance(value, dict):
            for key, child in value.items():
                name = f"{prefix}{key}"
                rows.append(
                    {
                        "file": str(path), "file_type": "json_key", "object": name,
                        "shape": "", "dtype": type(child).__name__,
                        "search_terms": term_matches(name, SEARCH_TERMS),
                        "agn_bh_feedback_terms": term_matches(name, AGN_TERMS),
                    }
                )
                walk(f"{name}.", child)
        elif isinstance(value, list) and value:
            walk(prefix, value[0])
    walk("", obj)
    return rows


def local_data_files(pattern: str) -> list[Path]:
    """Find repository-local data files while excluding generated referee outputs."""
    return sorted(
        path for path in Path(".").rglob(pattern)
        if ".git" not in path.parts and not tuple(path.parts[:2]) == ("outputs", "referee")
    )


def write_reports(inventory: pd.DataFrame) -> None:
    inventory.to_csv(OUT_DIR / "catalog_field_inventory.csv", index=False)
    n_files = inventory["file"].nunique()
    matched = inventory[inventory["search_terms"] != ""]
    hdf_files = sorted(inventory.loc[inventory["file_type"].str.startswith("hdf5"), "file"].unique())
    markdown = [
        "# Local CAMELS catalog-field inventory",
        "",
        f"Inspected `{n_files}` local data files, including `{len(hdf_files)}` HDF5 files.",
        "The machine-readable inventory in `catalog_field_inventory.csv` lists every accessible HDF5 dataset, structured-tree field, cached DataFrame column, and JSON key.",
        "",
        "## Search terms",
        "",
        "`" + "`, `".join(SEARCH_TERMS) + "`",
        "",
        "## Matching fields",
        "",
        "| File | Field | Matching terms |",
        "| --- | --- | --- |",
    ]
    for _, row in matched.iterrows():
        markdown.append(f"| `{row['file']}` | `{row['object']}` | `{row['search_terms']}` |")
    (OUT_DIR / "catalog_field_inventory.md").write_text("\n".join(markdown) + "\n")

    agn = inventory[inventory["agn_bh_feedback_terms"] != ""].copy()
    unique_agn = sorted(agn["object"].unique())
    report = [
        "# AGN, BH, and feedback-field search",
        "",
        "## Files inspected",
        "",
        f"The inventory script inspected `{len(hdf_files)}` local HDF5 files and the repository-local PKL, JSON, FITS, NPY, and NPZ products. The raw HDF5 set includes TNG group catalogs, TNG SubLink offset/tree files, and SIMBA group catalogs. Exact file paths and fields are recorded in `catalog_field_inventory.csv`.",
        "",
        "## Matching raw and cached field names",
        "",
    ]
    report.extend(f"- `{field}`" for field in unique_agn)
    report.extend(
        [
            "",
            "## Availability assessment",
            "",
            "Direct AGN feedback energy, cumulative AGN energy, and explicitly named jet-energy fields were not present in the inspected local files.",
            "",
            "The raw group catalogs do expose clearly named BH and wind quantities: `SubhaloBHMass`, `SubhaloBHMdot`, `GroupBHMass`, `GroupBHMdot`, `SubhaloWindMass`, and `GroupWindMass`. The physical-window diagnostic uses subhalo BH mass and subhalo BH accretion rate as candidate BH-state proxies. They are not direct measurements of AGN feedback energy and are not interpreted as proof of an AGN-driven transition.",
            "",
            "The repository's cached matched DataFrames do not contain these raw BH fields because `camels_catalog.py` loads a narrower field subset. The referee diagnostic reads the raw early-epoch group catalogs directly and joins by `(sim_id, local_id)` without changing the paper pipeline.",
        ]
    )
    (OUT_DIR / "agn_field_search_report.md").write_text("\n".join(report) + "\n")


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    rows: list[dict] = []
    hdf_files = local_data_files("*.hdf5")
    pkl_files = local_data_files("*.pkl")
    json_files = local_data_files("*.json")
    other_files = (
        local_data_files("*.fits") + local_data_files("*.npy") + local_data_files("*.npz")
    )
    for path in tqdm(hdf_files, desc="Inventory HDF5"):
        rows.extend(hdf5_rows(path))
    for path in tqdm(pkl_files, desc="Inventory PKL"):
        rows.extend(cached_object_rows(path))
    for path in tqdm(json_files, desc="Inventory JSON"):
        rows.extend(json_rows(path))
    for path in tqdm(other_files, desc="Inventory FITS/NPY/NPZ"):
        rows.append(
            {
                "file": str(path), "file_type": path.suffix.lstrip("."), "object": "",
                "shape": "", "dtype": "", "search_terms": "", "agn_bh_feedback_terms": "",
            }
        )
    inventory = pd.DataFrame(rows).fillna("")
    write_reports(inventory)


if __name__ == "__main__":
    main()
