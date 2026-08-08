# Local Scientific Data Inventory 001

## Purpose

Inventory the scientific data physically present in the runner's local repository before attempting another clean weak-lensing benchmark.

This audit is deliberately non-physical: it does not run PBUF, choose a source, construct a proxy, fit a parameter, normalize a field, or compare a prediction to lensing.

## What it scans

The script recursively scans the local repository for:

- FITS / FIT / FTS
- NPY / NPZ
- HDF5 / H5
- CSV / TSV

It excludes `.git/` and `runs/` so historical run artifacts are not mistaken for source data.

For FITS files it reports header metadata, dimensions and shape without loading the full cube into memory. Known kappa/gamma/shear products are identified from filenames. Other source-like files are flagged only as possible candidates for human review; nothing is automatically accepted.

## Run

```bash
PYTHONPATH=. python pbuf/labs/foundation/local_scientific_data_inventory001.py
```

## Interpretation

A candidate flag is not a source selection. The output is intended to tell us what independent source data actually exist locally so the frozen native weak-lensing lane can be benchmarked without inventing another source assumption.
