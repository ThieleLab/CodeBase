# Integrated PBPK-WBM Irinotecan Model

This repository contains the MATLAB code used for the simulations and analyses presented in:

**"An integrated whole-body metabolic and pharmacokinetic model predicts early irinotecan-induced metabolic perturbations driven by genetic and microbial variability"**

The framework integrates a physiologically based pharmacokinetic (PBPK) model of irinotecan with a whole-body metabolic (WBM) model containing irinotecan metabolism. It is used to simulate irinotecan and its metabolites, estimate pharmacokinetic parameters from clinical data, link drug disposition with organ-specific metabolism, and investigate genetic and gut microbial variability.

---

## Where to start

The workflow is organised into three MATLAB Live Scripts that should be run sequentially:

```text
MethodSectionIrin_part1.mlx
        ↓
MethodSectionIrin_part2.mlx
        ↓
MethodSectionIrin_part3.mlx
```

**Start with `MethodSectionIrin_part1.mlx`.** Parts 2 and 3 may use models or results generated in the preceding script.

---

## Repository structure

```text
Irinotecan/
│
├── PBPK_models/
│
├── src/
│   ├── Analysis/
│   ├── Cluster_Newton_Method/
│   ├── Integration_PBPK_WBM/
│   ├── PBPK_scripts/
│   ├── Personalisation/
│   └── Visualization/
│
├── MethodSectionIrin_part1.mlx
├── MethodSectionIrin_part2.mlx
├── MethodSectionIrin_part3.mlx
└── Readme.md
```

---

## Part 1 - Integrated PBPK-WBM model development

**File:** `MethodSectionIrin_part1.mlx`

This is the main starting point of the repository. It covers the development of the integrated model:

1. **PBPK model construction** – construction and parameterisation of the irinotecan PBPK model.
2. **Parameter estimation** – estimation of unknown PBPK parameters from clinical PK data using the Cluster Newton method.
3. **PBPK model validation** – comparison of simulated irinotecan and metabolite profiles with available PK data.
4. **WBM-irinotecan reconstruction** – incorporation of irinotecan metabolism into the sex-specific Harvey and Harvetta WBM models.
5. **PBPK-WBM integration** – coupling of time-dependent drug disposition with organ-specific whole-body metabolism.

Functions supporting these steps are located under `src/PBPK_scripts/`, `src/Cluster_Newton_Method/` and `src/Integration_PBPK_WBM/`.

---

## Part 2 - Genetic, microbial and dose variability

**File:** `MethodSectionIrin_part2.mlx`

Part 2 uses the integrated model from Part 1 to investigate variability in irinotecan pharmacokinetics and metabolism.

The simulations examine:

- different irinotecan doses;
- variability in UGT activity; and
- active or inactive intestinal microbial GUS.

The corresponding PBPK model variants are stored in `PBPK_models/`. The filenames identify the UGT condition, dose and GUS activity used in each simulation.

The results generated in this part are used for the analyses in Part 3.

---

## Part 3 - Flux and microbiome analysis

**File:** `MethodSectionIrin_part3.mlx`

Part 3 contains the downstream analysis of the integrated model simulations. This includes:

- organ-specific flux analysis;
- comparison of metabolic perturbations across UGT and dose conditions;
- analysis of drug and metabolite transport;
- microbiome-level analysis; and
- assessment of microbial deglucuronidation of SN-38G to SN-38.

Supporting analysis and plotting functions are available under `src/Analysis/` and `src/Visualization/`.

---

## Source code

The supporting MATLAB functions are organised under `src/`:

| Folder | Description |
|---|---|
| `Analysis/` | Processing and analysis of simulation and flux results |
| `Cluster_Newton_Method/` | Cluster Newton parameter-estimation functions |
| `Integration_PBPK_WBM/` | PBPK-WBM integration functions |
| `PBPK_scripts/` | PBPK model construction and simulation |
| `Personalisation/` | Generation of model-specific simulation conditions |
| `Visualization/` | Plotting and visualisation functions |

### `PBPK_models/`

Contains the PBPK model variants used for different irinotecan doses, UGT activity states and intestinal microbial GUS conditions.

---

## Requirements

The main requirements are:

- MATLAB
- COBRA Toolbox v3.0
- PSCM extension for the COBRA Toolbox
- Compatible LP/QP solver
- Symbolic Math Toolbox
- Parallel Computing Toolbox

Some parameter-estimation and large-scale flux simulations can be computationally intensive and may benefit from a high-performance computing system.

---

## Running the workflow

Clone or download the repository and set the `Irinotecan/` folder as the MATLAB working directory. Make sure the required dependencies and `src/` functions are available on the MATLAB path.

Run the Live Scripts in the following order:

```text
1. MethodSectionIrin_part1.mlx   → Model development and integration
2. MethodSectionIrin_part2.mlx   → Variability simulations
3. MethodSectionIrin_part3.mlx   → Flux and microbiome analyses
```

For a first-time user, **`MethodSectionIrin_part1.mlx` is the entry point.**

---

## Acknowledgements

This study was funded by the European Research Council (ERC) under the European Union’s
Horizon 2020 research and innovation programme (101125633) to IT and the Science Foundation
Ireland under Grant number 12/RC/2273-P2. Funding for this project was also provided through
NIA grant U19AG063744 Alzheimer’s Gut Microbiome Project (AGMP), PI Kaddurah-Daouk at
Duke University along with several academic institutions.