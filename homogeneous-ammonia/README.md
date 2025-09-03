# NH3 Synthesis Example

This directory contains scripts and data demonstrating one active-learning generation for gas-phase ammonia synthesis using the RAL workflow.

## Contents
- `gen-1/` – first active-learning generation
  - `interm.md/` – tools for generating intermediate MD trajectories
  - `relabel.md/` – frames selected for DFT relabeling
  - `train/` – DeepMD training scripts
  - `README.interm.md` – step-by-step instructions for the interm MD stage

Scripts assume a SLURM scheduler and a working DeepMD-kit installation. Adjust paths and submission commands for your environment.

See the [root README](../README.md) for an overview of RAL and citations.
