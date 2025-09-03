# CH$_4$ Activation on TiC

Example data and scripts for applying RAL to methane activation and C–C coupling on TiC(001).

## Notebooks
- `TiC-init-data.ipynb` – generate initial MD data
- `TiC-gen-struct.ipynb` – build reactant/product structures and initial NEB paths
- `TiC-adsorbate.ipynb` – create adsorbate geometries
- `TiC-filtering.ipynb` – filter NEB frames using committee disagreement

## Generations
- `gen-0/` – initial committee training and NEB exploration
- `gen-1/` – first active-learning cycle
  - `md/` – MD exploration scripts
  - `neb/` – NEB path generation utilities
  - `relabel.md/` and `relabel.neb/` – frames selected for DFT relabeling

## Utilities
- `neb_codes/` – helper scripts for processing NEB trajectories

All scripts assume a SLURM scheduler and DeepMD-kit. Adjust file paths and submission settings for your computing environment.
