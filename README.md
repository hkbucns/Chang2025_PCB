
# Table of Content
- [Reconstruction Framework](#reconstruction-framework)
- [Getting Started](#getting-started)
  - [Quick start: run the reconstruction procedure](#quick-start-run-the-reconstruction-procedure)
  - [Exponential Scaling](#exponential-scaling)
  - [Hidden Inhibition](#hidden-inhibition)

# Reconstruction Framework

The reconstruction procedure is a inversion-problem designed to reconstruct the effective heterogeneity and asymmetric structural connectivity from observed neural activity. In brief, the reconstruction procedure consists of two part: Temporal Reconstruction: we use existed reconstruction framework: Dynamical Differential Covariance to infer the Jacobian matrix from neural activity under the assumption of stable stochastic process. Spatial Reconstruction: by incorporating the symmetric structural connectivity, we further seperate the directionality contributed from effective heterogeneity and asymmetric structural connectivity. 


# Getting Started

## Quick start: run the reconstruction procedure

The following provides a step-by-step instruction to simulate the neural mass model with heterogeneity and asymmetric connectivity and run the reconstruction procedure.
Refer to the Manuscript, the steps below can get the results from Fig. 1 and 2.
1. Open the main performing file, **HetergeneousMainTestScript1.m**.
2. Run the first 4 blocks in the main performing file to create the model parameters, run the simulation and calculate the ground-truth Jacobian matrix.
3. Run the 5th block for Temporal Reconstruction. This call function, **LinearReconst.m** to estimate the Jacobian Matrix.
5. Run the 6th block for Spatial Reconstruction, **RevealHHetero1.m** is called to further separate the Jacobian to effective heterogeneity and asymmetric structural connectivity. (Actually you can run **RevealHHetero2.m** which is a compact version of spatial reconstruction, only reconstruct asymmetric SC and effective heterogeneity.)
6. Run the following blocks until **Line 192** for validation and evaluation of detailed reconstruction of model A in **HetergeneousMainTestScript1.m**.
7. Run the following codes start from **Line 193: Replacement or tau & b from w & I** to Run the sample results at Fig. 4 and Fig. S3, illustrating the replacement of heterogeneity pairs. This will call function **dMFM_eq.m** to simulate the heterogeneous model B with tau_i and b_i.

## Exponential Scaling
8. Open the script testing the sampling effect, **ErrorSampling.m**.
9. Run this file, will automatically plot the main results in Fig. 3.

## Hidden Inhibition
10. Open the script to run the final section of hidden inhibition effect, **EI_HetergeneousMainTestScript1.m**.
11. This script focus on model C with local inhibitory populations, and only go through the reconstruction with excitatory population activity.
12. Running this file will generate results in Fig. 5.
