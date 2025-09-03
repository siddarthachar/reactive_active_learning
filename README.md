# Reactive Active Learning (RAL)

[![DOI](https://img.shields.io/badge/DOI-10.1021%2Facs.jctc.5c00920-blue)](https://doi.org/10.1021/acs.jctc.5c00920)

**Reactive Active Learning (RAL)** is an automated framework for training machine-learned interatomic potentials for reactive systems. It combines automated reaction exploration, active learning based on model uncertainty, and retraining cycles to converge on accurate reactive MLIPs with minimal data. The method is validated across gas-phase, solution-phase, and heterogeneous catalysis examples.

---

## Repository Contents
```
reactive_active_learning/
├── homogeneous-ammonia/       # Example: NH3 synthesis (SE-GSM RAL)
├── homogeneous-methanimine/   # Example: Methanimine hydrolysis
├── heterogeneous-TiC/         # Example: CH4 activation/coupling on TiC
└── README.md
```

Each folder contains scripts and example data demonstrating the workflow described in the paper.

---

## Dependencies

Recommended installation using conda:

```bash
conda create -y -n ral python=3.10
conda activate ral
pip install numpy scipy pandas matplotlib tqdm pyyaml
pip install ase networkx rdkit-pypi
pip install jupyter ipykernel
```
Optional: PyTorch or another ML framework of choice for training committee models.

RAL requires access to a DFT package (e.g., VASP, Quantum ESPRESSO) for generating reference data.
