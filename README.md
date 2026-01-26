# Doubly-Robust Machine Learning Causal Estimation of Pathway Activity Effects in Breast Cancer Stage Progression and Survival

- Authors: Amos Okutse, Ian Liu, Man-Fang Liang

## Introduction

This repository contains code and data for the paper "Doubly-Robust Machine Learning Causal Estimation of Pathway Activity Effects in Breast Cancer Stage Progression and Survival". The study focuses on estimating the causal effects of pathway activities on breast cancer progression and survival using advanced machine learning techniques.

## Repository Structure

- `data/`: Contains the datasets used in the analysis. These data are from the Cancer Genome Atlas (TCGA) and refer to breast cancer. Within this folder, we have data files for gene expression (`data_mrna_illumina_microarray.txt`), clinical information (`data_clinical_patient.txt`), patient sample information (`data_clinical_sample.txt`), gene mutation (`data_mutations_extended.txt`), and gene mutation data (`data_methylation_promoters_rrbs.txt`).
- `scripts/`: Contains the R and Python scripts for data preprocessing, model training, and evaluation. Within this file, the methods for data cleaning, feature engineering, and causal effect estimation are implemented separated for each task. `methods_survival.R` contains the code for survival analysis of the effect of the identified pathway on breast cancer survival, while `methods_tumor_stage.R` contains the code for stage progression analysis. All tasks employ a doubly-robust estimation procedure using machine learning models. The `survival_helpers.R` file contains helper functions for survival analysis.
- `results/`: Contains the output results, including figures and tables. 
- `notebooks/`: Jupyter notebooks for exploratory data analysis and visualization.
- `README.md`: This file. 

## Reproduction Steps

To reproduce the findings in this study, follow these steps:
1. **Data Preparation**: Download the necessary datasets from TCGA and place them in the `data/` directory. These data are available from the following link: https://www.cbioportal.org/
2. **Install Dependencies**: Ensure you have R and Python installed along with the required libraries listed in `requirements.txt` and `environment.yml`.
3. **Run Preprocessing Scripts**: Execute the scripts in the `scripts/` directory to preprocess the data.
4. **Model Training and Evaluation**: Run the main analysis scripts to train the models and evaluate the causal effects.
5. **Results Visualization**: Use the scripts provided in the `scripts/` directory to visualize the results.
6. **Review Results**: Check the `results/` directory for output figures and tables. 

## Contact
For any questions or issues, please contact:

(1) [Amos Okutse](amos_okutse@brown.edu)

(2) [Ian Liu](ian_y_liu@brown.edu)

(3) [Man-Fang Liang](man-fang_liang@brown.edu)
