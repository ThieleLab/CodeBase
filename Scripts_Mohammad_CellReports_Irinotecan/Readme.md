# Irinotecan PBPK-WBM Modeling Project

This project develops an advanced computational framework that integrates a Whole-Body Drug Metabolic (WBM) model with a Physiologically Based Pharmacokinetic (PBPK) model for **Irinotecan**. The goal is to simulate the temporal distribution of Irinotecan, estimate individual-specific parameters using clinical pharmacokinetic data, and link drug transport with organ-specific metabolism. This framework is designed to model the effects of genetic, microbial, and dosage variability on the metabolism and pharmacokinetics of Irinotecan.

## Project Structure

The project is organized into the following folders and files:
wbm_modelingcode\WBM_PBK\Irinotecan
├── Irinotecan_WBMPBPK_modelBuilding.mlx
├── Irinotecan_study.mlx
├── PBPK_Models_Folder
├── Model_Source_Code
└── Irinotecan_Model_Analysis.mlx


### Files

- **Irinotecan_WBMPBPK_modelBuilding.mlx**: 
  - This MATLAB Live Script integrates the Whole-Body Metabolic (WBM) model with the Physiologically Based Pharmacokinetic (PBPK) model for Irinotecan. 
  - It replicates the simulations presented in the study titled *"An integrated whole-body metabolic and pharmacokinetic model predicts early irinotecan-induced metabolic perturbations driven by genetic and microbial variability"*. 
  - The script details the development of the generic PBPK model, parameter estimation via the Cluster Newton method, and validation with pharmacokinetic profiles.
  - A comprehensive WBM model incorporating irinotecan metabolism, including sex-specific versions (Harvey and Harvetta), is built. 
  - PBPK and WBM models are integrated using an indirect coupling approach.

- **Irinotecan_study.mlx**: 
  - This script delves into the framework’s capability to support various analyses focusing on the variability in irinotecan metabolism.
  - The main areas of study include:
    1. **High-dose Irinotecan Exposure and UGT/GUS Enzyme Variability**: Investigates how high-dose irinotecan impacts UGT enzyme activity and the potential for GUS enzyme inactivity at different dosage levels, influencing drug metabolism and response.
    2. **Microbiome-Specific Irinotecan Metabolism**: Examines the role of the microbiome in irinotecan metabolism using personalized host-microbial models, exploring the deglucuronidation potential of GUS in the human intestine.

- **Irinotecan_Model_Analysis.mlx**:
  - A detailed analysis of flux-based and microbiome-level dynamics across various organ systems is conducted.
  - Key aspects of the analysis include:
    - **Flux Distribution Comparisons**: Assess how UGT and GUS enzyme polymorphisms and dosage differences impact flux across organs (e.g., liver, urine, blood plasma, and feces). Visualization techniques like violin plots and pathway scores are used to evaluate metabolic perturbations due to enzyme activity changes.
    - **Microbiome-Level Analysis**: Integrates microbiome data from personalized host-microbial models, focusing on microbial deglucuronidation of SN38G to SN38, and how microbial factors contribute to irinotecan metabolism.

### Folders

- **PBPK_Models_Folder**: 
  - This folder contains the different PBPK model files related to Irinotecan, including various formulations, doses, and assumptions for pharmacokinetic simulations.
  
- **Model_Source_Code**: 
  - Contains source code that defines and simulates the integrated PBPK-WBM models. It includes functions, libraries, and configuration files necessary for running simulations.

## Getting Started

To run the project and generate results, follow these steps:

1. **Clone or Download** the repository to your local machine.

2. **Install Required Dependencies**:
    - **MATLAB**: Ensure you have the latest version of MATLAB.
    - **COBRA Toolbox (v3.0)**: You can download it [here](https://opencobra.github.io/cobratoolbox/stable/installation.html).
    - **PSCM Toolbox extension** for the COBRA Toolbox (required for advanced modeling).
    - **Violin Plot Functions** by Bastian Bechtold: Available [here](https://github.com/bastibe/Violinplot-Matlab).
    - **Linear and Quadratic Programming Solver**: The code has been tested using IBM ILog Cplex 12.1 (academic version, solver ILOGcomplex in the COBRA Toolbox).
    - **Symbolic Math Toolbox**: Required for symbolic computations in MATLAB.
    - **Deep Learning Toolbox**: For advanced model training and optimization.
    - **Parallel MATLAB Toolbox**: Required for parallel computation to speed up simulations.
  
3. **Optional**: If available, access to a high-performance computing facility can improve computation time for large-scale simulations.

4. Navigate to the folder containing the project files and open the main `.mlx` file:
    - Open `Irinotecan_WBMPBPK_modelBuilding.mlx` to start building the integrated model.
    - Open `Irinotecan_study.mlx` to explore the different analyses on drug metabolism variability.
    - Open `Irinotecan_Model_Analysis.mlx` for detailed pharmacokinetic analyses.

5. Run the scripts in MATLAB to execute model-building processes, simulations, and analyses.

## Prerequisites

- **MATLAB** (with SimBiology, Optimization Toolbox, and other necessary toolboxes).
- **COBRA Toolbox v3.0**.
- **PSCM Toolbox Extension** for COBRA Toolbox.
- **Violin Plot Functions** by Bastian Bechtold.
- **IBM ILog Cplex 12.1** (or another LP/QP solver compatible with the COBRA Toolbox).
- **Symbolic Math Toolbox** for symbolic computations.
- **Deep Learning Toolbox** for advanced model training and optimization.
- **Parallel MATLAB Toolbox** for parallel computations.
- **Optional**: Access to a high-performance computing facility to improve computation time.

## Acknowledgments

- The authors express their gratitude for the financial assistance provided by the European Research Council through the Horizon 2020 research and innovation program of the European Union (Grant \#757922). The authors thank all members of Thiele's lab for their support.
