# AGN, BH, and feedback-field search

## Files inspected

The inventory script inspected `193` local HDF5 files and the repository-local PKL, JSON, FITS, NPY, and NPZ products. The raw HDF5 set includes TNG group catalogs, TNG SubLink offset/tree files, and SIMBA group catalogs. Exact file paths and fields are recorded in `catalog_field_inventory.csv`.

## Matching raw and cached field names

- `Group/GroupBHMass`
- `Group/GroupBHMdot`
- `Group/GroupWindMass`
- `Subhalo/SubhaloBHMass`
- `Subhalo/SubhaloBHMdot`
- `Subhalo/SubhaloWindMass`

## Availability assessment

Direct AGN feedback energy, cumulative AGN energy, and explicitly named jet-energy fields were not present in the inspected local files.

The raw group catalogs do expose clearly named BH and wind quantities: `SubhaloBHMass`, `SubhaloBHMdot`, `GroupBHMass`, `GroupBHMdot`, `SubhaloWindMass`, and `GroupWindMass`. The physical-window diagnostic uses subhalo BH mass and subhalo BH accretion rate as candidate BH-state proxies. They are not direct measurements of AGN feedback energy and are not interpreted as proof of an AGN-driven transition.

The repository's cached matched DataFrames do not contain these raw BH fields because `camels_catalog.py` loads a narrower field subset. The referee diagnostic reads the raw early-epoch group catalogs directly and joins by `(sim_id, local_id)` without changing the paper pipeline.
