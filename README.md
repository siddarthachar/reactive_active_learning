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

## Quick Start

	1.	Train an initial MLIP committee on short DFT-MD of reactants.
	2.	Generate product candidates and reaction paths (SE-GSM for homogeneous systems; NEB + CatKit for surfaces).
	3.	Run path exploration with a committee member, collect frames.
	4.	Filter frames by model uncertainty (force deviation ε).
	5.	Relabel selected frames with DFT and retrain the committee.
	6.	Repeat until convergence criteria are satisfied (ε, UMAP convex hull, NN distance).

See the example folders for scripts reproducing workflows from the paper.

## Citation

If you use RAL, please cite:

> S. K. Achar, P. B. Shukla, C. V. Mhatre, L. Bernasconi, C. Y. Vinger, and J. K. Johnson.  
> *Reactive Active Learning: An Efficient Approach for Training Machine Learning Interatomic Potentials for Reacting Systems*.  
> *J. Chem. Theory Comput.* (2025). [https://doi.org/10.1021/acs.jctc.5c00920](https://doi.org/10.1021/acs.jctc.5c00920)

### BibTeX:
```
@article{achar2025ral,
  author  = {Achar, Siddarth K. and Shukla, Priyanka B. and Mhatre, Chinmay V. and Bernasconi, Leonardo and Vinger, Caitlyn Y. and Johnson, J. Karl},
  title   = {Reactive Active Learning: An Efficient Approach for Training Machine Learning Interatomic Potentials for Reacting Systems},
  journal = {Journal of Chemical Theory and Computation},
  year    = {2025},
  doi     = {10.1021/acs.jctc.5c00920}
}
```

## Disclaimer

This repository provides example workflows and scripts associated with the publication. Users are expected to adapt the code for their own systems and provide their own DFT software. The repository is getting tidied up.

For questions regarding the codes, contact ska31@pitt.edu.
