# Methanimine Hydrolysis Example

This folder provides inputs from two active-learning cycles for the hydrolysis of methanimine in solution.

## Contents
- `gen-2/` – second generation dataset and scripts
- `gen-2.1/` – follow-up generation with additional training data

Each generation contains:
- `interm.md/` – generation of intermediate MD configurations
- `relabel.md/` – selection of frames for DFT labeling
- `train/` (when present) – DeepMD training utilities

The scripts are written for a SLURM-based cluster and DeepMD-kit. Update resource requests and paths to suit your setup.

Refer to the [top-level README](../README.md) for a description of the overall RAL methodology.
