# Koopman Operator Framework for Off-Road Vehicle Modeling and Control

[![Paper](https://img.shields.io/badge/arXiv-2603.28965-b31b1b.svg)](https://arxiv.org/abs/2603.28965)
[![MATLAB](https://img.shields.io/badge/MATLAB-R2023b-orange.svg)](https://www.mathworks.com/products/matlab.html)

MATLAB research code accompanying the paper **“Koopman Operator Framework for Modeling and Control of Off-Road Vehicle on Deformable Terrain”** by Kartik Loya and Phanindra Tallapragada.

The paper is available on [arXiv](https://arxiv.org/abs/2603.28965) and has been submitted to the *ASME Journal of Autonomous Vehicles and Systems* (JAVS-26-1012).

## About this repository

This repository contains the MATLAB simulation, data-processing, Koopman-model training, prediction-analysis, and model predictive control code used to study off-road vehicle dynamics on sandy-loam and clay terrain.

**Paper conclusion:** Terrain-specific Koopman operators provide accurate short-horizon predictions at a computational cost suitable for constrained MPC. The results also show that using an operator learned for the wrong terrain can cause large prediction errors and unsafe closed-loop behavior, motivating terrain-aware model selection or online adaptation.

## Overview

This repository implements a hybrid physics-informed and data-driven framework for predicting and controlling off-road vehicle motion on deformable terrain. The workflow combines:

- a five-degree-of-freedom vehicle simulation with Bekker–Wong terramechanics;
- simulated trajectories on sandy loam and clay;
- recursive subspace identification with Grassmannian-distance-based data selection;
- a finite-dimensional linear Koopman predictor with Gaussian-process lifting; and
- constrained Koopman model predictive control (KMPC) for trajectory tracking.

The learned models provide stable short-horizon predictions, can account for mild terrain-height variation, and support closed-loop tracking subject to steering and wheel-torque constraints.

## Prediction results

The Koopman predictors are evaluated separately on sandy-loam and clay test trajectories. The plots report the predicted body-frame longitudinal velocity $u$, lateral velocity $v$, and yaw rate $\dot{\psi}$.

### Sandy loam

The representative trajectory compares the simulated response (blue) with the K-SSID prediction (red) over 20 seconds.

<p align="center">
  <img src="result_analysis/Sandyloam_prediction/New_Images/Final/SL_rand_traj_6_150.png" alt="Sandy-loam trajectory and velocity prediction" width="650">
</p>

Prediction error generally increases as the interval between measurement refreshes grows and as the open-loop prediction advances in time.

| Error over the refreshing horizon | Error over prediction time |
|:---:|:---:|
| <img src="result_analysis/Sandyloam_prediction/New_Images/Final/SL_refreshing_time_NRMSE.png" alt="Sandy-loam prediction error versus refreshing horizon" width="80%"> | <img src="result_analysis/Sandyloam_prediction/New_Images/Final/SL_error_with_time_time_RMSE.png" alt="Sandy-loam prediction error over time" width="80%"> |

### Clay

The clay example likewise compares the simulated and Koopman-predicted trajectory and velocity states over 20 seconds.

<p align="center">
  <img src="result_analysis/Clay_prediction/clay_123_rand_trajectory_23_v2.png" alt="Clay trajectory and velocity prediction" width="650">
</p>

The clay error plots show the effect of extending the measurement-refresh interval and the growth of prediction error with time.

| Error over the refreshing horizon | Error over prediction time |
|:---:|:---:|
| <img src="result_analysis/Clay_prediction/clay_refreshing_time_RMSE_123_v2.png" alt="Clay prediction error versus refreshing horizon" width="80%"> | <img src="result_analysis/Clay_prediction/clay_error_with_time_123_v2.png" alt="Clay prediction error over time" width="80%"> |

<!-- ### Need for terrain-specific Koopman operators

A Koopman operator trained on sandy loam fails when it is used to control the vehicle on clay because the terrain mismatch causes prediction drift and loss of safe trajectory tracking. In contrast, the clay-specific operator remains stable on clay and successfully completes the maneuver, motivating terrain-aware operator selection or online adaptation.

| **Failure:** sandy-loam operator on clay | **Success:** clay operator on clay |
|:---:|:---:|
| <img src="docs/media/kmpc_sandy_operator_on_clay_failure.gif" alt="Failure of the sandy-loam Koopman operator during closed-loop control on clay" width="80%"> | <img src="docs/media/kmpc_clay_operator_on_clay_safe.gif" alt="Successful closed-loop control using the clay Koopman operator on clay" width="80%"> | -->

## Repository structure

```text
functions/
  simulation_model/       Vehicle and deformable-terrain dynamics
  utility/                Recursive SSID, Koopman fitting, and prediction
scripts/
  simulation/             Flat-terrain data generation
  simulation_elev/        Terrain-height-variation data generation
  data_handler/           Dataset assembly and train/validation/test splits
  koopman_training/       Model training and hyperparameter sweeps
result_analysis/
  Clay_prediction/        Clay prediction analysis and figures
  Sandyloam_prediction/   Sandy-loam prediction analysis and figures
  koopman_MPC/            Koopman MPC scripts and results
```

Large simulation datasets, trained model workspaces, and generated MATLAB result files are intentionally excluded from Git.

## Requirements

- MATLAB R2023b or a compatible release
- System Identification Toolbox
- Statistics and Machine Learning Toolbox
- Optimization Toolbox
- Parallel Computing Toolbox (for `parfor` training and simulation loops)
- [CasADi for MATLAB](https://web.casadi.org/get/) with IPOPT for the KMPC scripts

The MATLAB release shown above matches the supplied SLURM scripts. Other recent releases may work but have not been documented here.

## Workflow

Run MATLAB from the repository root and add the project functions to the path:

```matlab
addpath(genpath("functions"));
```

The main stages are:

1. Generate input signals with `scripts/offroad_Inputsignals.mlx`.
2. Generate terrain-specific trajectories with the scripts in `scripts/simulation/` or `scripts/simulation_elev/`.
3. Assemble the trajectories into training, validation, and test sets with `scripts/data_handler/gather_plot_data.m`.
4. Train Koopman models with `scripts/koopman_training/Offroad_Koopman_Training_main.m`. The example SLURM array scripts in `scripts/koopman_training/results/` pass the dataset path and the hyperparameters listed in `scripts/koopman_training/params.txt`.
5. Reproduce prediction plots with the scripts under `result_analysis/Clay_prediction/` and `result_analysis/Sandyloam_prediction/`.
6. Run constrained trajectory-tracking examples from `result_analysis/koopman_MPC/` after changing the CasADi and model paths for your machine.

Several scripts contain experiment-specific relative paths and are intended as research artifacts rather than a packaged MATLAB application. Update dataset, output, and CasADi paths before running them.

## Data and generated models

The full datasets and trained MATLAB workspaces are too large for standard Git hosting. The code expects terrain datasets under `datasets/` and model files under:

```text
scripts/koopman_training/results/clay_noelev_models/models_with_error/
scripts/koopman_training/results/sandyloam_noelev_models/models_with_error/
```

Contact the authors if you need access to the exact data or trained models used in the paper.

## Citation

If this repository is useful in your research, please cite:

```bibtex
@article{loya2026koopman,
  title   = {Koopman Operator Framework for Modeling and Control of Off-Road Vehicle on Deformable Terrain},
  author  = {Loya, Kartik and Tallapragada, Phanindra},
  journal = {arXiv preprint arXiv:2603.28965},
  year    = {2026},
  doi     = {10.48550/arXiv.2603.28965}
}
```

## Paper

K. Loya and P. Tallapragada, “Koopman Operator Framework for Modeling and Control of Off-Road Vehicle on Deformable Terrain,” arXiv:2603.28965, 2026. [Paper](https://arxiv.org/abs/2603.28965) · [PDF](https://arxiv.org/pdf/2603.28965)
