# Targeting cancer-associated fibroblasts for treatment of ER+ breast cancer: A mathematical modeling perspective and optimization of treatment strategies

Authors: Tuğba Akman, Kristian Pietras, A Köhn-Luque and Ahmet Acar.

Codes and data to reproduce the results and figures in the preprint published in bioRxiv ([https://doi.org/10.1007/s11538-023-01253-1](https://doi.org/10.64898/2026.03.27.714662) are provided in this repository. 

# Practical guide

## identifiability/ (Structural identifiability)
  runCAF2_with_denom_v2.m and runCAF2_wo_CAF_with_denom_v2.m must be run. These codes work with GenSSI 2.0 software toolbox (https://github.com/genssi-developer/GenSSI).

## Monolix_after_Dec2025/ (Parameter estimation)
  Monolix 2024R1 has been used for model calibration (https://monolixsuite.slp-software.com/getting-started/2024R1/download). Model_MCF7_wo_CAF2_high_med_low.txt and Data_MCF7_wo_CAF2_high_med_low_v2.txt, Model_MCF7_with_CAF2_high_med_low.txt and Data_MCF7_with_CAF2_high_med_low.txt must be uploaded for model 1 and 5, respectively. Adjust the parameter intervals and decide the fixed and random effects. Then, run the code.

## Boxplots/ (Boxplots)
  Boxplots can be reproduced by running the R file.

## OptimalTreatment/ (Optimal control problem)
  For both constant and optimal treatment for three different E2 levels, main_OCP_xxx.m must be run.
  


  
