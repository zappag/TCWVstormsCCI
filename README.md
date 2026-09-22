# ESACCI TCWV Masking and Track Analysis

This repository contains a workflow for preparing ESA CCI total column water vapour (TCWV) data, applying data-quality and surface masks, and then combining the masked fields with atmospheric track data (for example ERA5 cyclone tracks) to compute storm-centered statistics.

The project is built around two main processing stages:

1. `mask_tcwv_esa.py` — loads ESA CCI WV data, applies mask logic, and saves seasonal NetCDF files.
2. `add_tcwv_tracks.py` — loads the masked data and appends TCWV statistics to track files for downstream analysis.

---

## Overview

The workflow is designed for climate research applications where one wants to:

- work with ESA CCI WV products (`v3` or `v4`),
- apply quality filters and surface-type masks,
- keep only valid observations and complete seasonal periods,
- compare reanalysis fields (ERA5, MERRA2, JRA3Q) against ESA-derived products,
- attach TCWV information to cyclone/anticyclone tracks,
- calculate seasonal composites and summary statistics.

The processing is organized around seasonal batches for each meteorological season:

- DJF
- MAM
- JJA
- SON

---

## Main processing script: `mask_tcwv_esa.py`

This script handles the ESA CCI data preparation and masking workflow.

### What it does

- loads ESA CCI WV NetCDF files for a date range,
- normalizes longitudes to $0\! -\! 360^\circ$,
- fixes the dataset time coordinate to a daily noon-centered representation,
- applies version-specific surface masks and quality checks,
- saves seasonal files for downstream analysis,
- optionally creates cross-mask variants for sensitivity tests.

### Mask logic

The script derives masks from `surface_type_flag`, `atmospheric_conditions_flag`, and `tcwv_quality_flag`.

For ESA CCI v3:

- `land_coast_ice_mask` includes surface types `0, 2, 4, 5, 6, 7`
- `ocean_mask` includes `1, 3`
- `land_mask` includes `0, 2, 6`
- `land_clear_mask` includes `0`
- `valid_data_mask` keeps `tcwv_quality_flag == 0`

For ESA CCI v4:

- `land_coast_ice_mask` includes `0, 2, 3, 4, 5, 6`
- `ocean_mask` includes `1`
- `land_mask` includes `0`
- `land_clear_mask` uses `atmospheric_conditions_flag == 1`
- `valid_data_mask` keeps `tcwv_quality_flag == 0`

The script also sets values near the pole to a fixed ice/flag value to avoid invalid high-latitude observations.

### Seasonal output format

It splits the dataset into complete meteorological seasons and writes files like:

- `tcwv_ESA_2003_MAM.nc`
- `tcwv_ESA_LandCoastIce_2003_MAM.nc`
- `tcwv_ESA_Ocean_2003_MAM.nc`

Output directories are typically organized by dataset and resolution, for example:

- `.../cdr2/day_v2/`
- `.../cdr2/day_masked_v2/`
- `.../cdr2/masks/`

### Command-line usage

```bash
python mask_tcwv_esa.py --fdata ESA --year 2003 --version 4
python mask_tcwv_esa.py --fdata ERA5 --year 2003 --version 4 --fullrun
```

Key arguments:

- `--fdata`: dataset type (`ERA5`, `MERRA2`, `ESA`, `JRA3Q`)
- `--year`: single year to process; if omitted, the full period is used
- `--version`: dataset version (`3` or `4`)
- `--fullrun`: performs reanalysis remapping and seasonal extraction
- `--crossmask`: applies masks from the other ESA version for cross-check experiments

---

## Downstream analysis: `add_tcwv_tracks.py`

This second script takes track files and attaches TCWV fields from the processed data products.

### Purpose

It is used to answer questions such as:

- how much water vapour is associated with a cyclone at a given radius,
- whether the track passes over land, ocean, or clear-sky conditions,
- how the thermodynamic environment differs between dataset versions or masks.

### Supported datasets

The script supports:

- `ERA5`
- `MERRA2`
- `ESA`
- `JRA3Q`
- `ESAcross`

It also supports submasks such as:

- `LandCoastIce`
- `Ocean`
- `LandClear`

and different resolutions such as:

- `50km`
- `50km-day`
- `50km-day-masked`

### Processing steps

The script loops over seasons and years and then performs tasks such as:

- read track files for each season,
- optionally subsample time steps to daily 12UTC values,
- identify relevant TCWV files for that season,
- append mean/min/max TCWV fields around each track point,
- compute track-level summary statistics.

Main argument flags include:

- `--fdata`: dataset to add to the tracks
- `--fres`: resolution suffix
- `--fsubmask`: optional submask name
- `--radius`: radius for neighbourhood statistics
- `--dosubsample`: perform temporal subsampling
- `--doaddtcwv`: add TCWV field to tracks
- `--dostats`: compute statistics
- `--doradial`: run radial analysis
- `--testing`: restrict processing to a smaller test period
- `--version`: ESA version (`3` or `4`)
- `--feature`: feature type (`cyc` or `anticyc`)

### Example

```bash
python add_tcwv_tracks.py \
  --fdata ESA \
  --fres 50km-day-masked \
  --fsubmask LandCoastIce \
  --radius 5 \
  --doaddtcwv \
  --dostats \
  --version 4 \
  --feature cyc
```

---

## Data flow

The typical research workflow is:

1. Download ESA CCI WV files and reanalysis data.
2. Run `mask_tcwv_esa.py` to prepare and seasonally split masked TCWV data.
3. Run `add_tcwv_tracks.py` to attach TCWV to storm tracks.
4. Compute seasonal composites or summary statistics across track points.
5. Compare results across datasets, masks, and versions.

---

## Dependencies

This project relies on:

- Python 3
- `xarray`
- `numpy`
- `pandas`
- `matplotlib`
- `cdo` (Climate Data Operators)
- `netCDF4`
- `glob`-based file management and standard Python libraries

---

## Repository structure

```text
ESACCI/
├── mask_tcwv_esa.py
├── add_tcwv_tracks.py
├── split_to_seasons.py
├── readme
├── submit.job
├── submit_tcwv_mask.job
├── submit_advars.job
├── old/
└── ...
```

---

## Notes

- This project assumes a local data layout on large storage systems (for example `/mnt/naszappa/...`).
- Output filenames are heavily versioned and season-based, which makes the workflow suitable for batch climate analysis.
- The code is written for research workflows rather than a general-purpose public library; configuration paths and directories are largely hard-coded.

---

## License

This project does not currently include a formal license file. Please check with the repository owner before redistribution or broader reuse.
