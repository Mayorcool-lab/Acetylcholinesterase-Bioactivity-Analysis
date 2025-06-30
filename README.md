# 🧠 Acetylcholinesterase Bioactivity Prediction Pipeline

A fully reproducible cheminformatics pipeline designed to **predict the bioactivity (pIC₅₀)** of acetylcholinesterase inhibitors using **machine learning** and **molecular descriptors**.

> 🚀 Built for drug discovery teams seeking robust and interpretable compound screening tools.

---

## 🎯 Project Objectives

- Retrieve and curate **bioactivity data** for acetylcholinesterase inhibitors from **ChEMBL (v26)**
- Generate molecular descriptors using **RDKit** and **PaDEL-Descriptor**
- Apply **Lipinski’s Rule of Five** to assess drug-likeness
- Transform **IC50** values into **pIC₅₀** for modeling
- Train predictive **machine learning models** (e.g., Random Forest)
- Interpret structure–activity relationships
- Automate the workflow using **Snakemake**
- Ensure reproducibility with **Conda** and **Git**

---

## 🧪 Tools & Technologies

| Category           | Tools / Libraries                                                                 |
|--------------------|------------------------------------------------------------------------------------|
| Data & Chemistry   | [ChEMBL](https://www.ebi.ac.uk/chembl/), RDKit, PaDEL-Descriptor (Java)           |
| ML & Analysis      | scikit-learn, Pandas, NumPy, Matplotlib, Seaborn                                  |
| Workflow & Repro   | Snakemake, Conda, Git, GitHub                                                      |

---

## 🔁 Pipeline Overview

### 1. **Data Retrieval**
- Downloaded IC50 values and SMILES for inhibitors from ChEMBL using the Webresource Client.

### 2. **Data Preprocessing**
- Removed invalid/duplicate records
- Standardized IC50 values and converted to pIC₅₀
- Created activity classes: Active, Inactive, Intermediate

### 3. **Descriptor Generation**
- **RDKit**: MW, LogP, H-bond donors/acceptors
- **PaDEL-Descriptor**: 881 binary fingerprints
- Applied Lipinski's rules

### 4. **Exploratory Analysis**
- Distribution plots, MW vs LogP, pIC₅₀ boxplots
- Mann–Whitney U test for feature significance

### 5. **Machine Learning**
- Feature Set: RDKit + PaDEL fingerprints
- Target: pIC₅₀ (regression)
- Model: Random Forest (optimized)
- Best R² Score: **0.765**, RMSE: **0.738**
- Feature importance analysis and prediction scatter plots

### 6. **Snakemake Workflow**
- Modular pipeline for data merging, training, and evaluation
- Easily reproducible across environments

---

## 📁 Project Structure

```bash
Acetylcholinesterase_Project/
├── data/                 # Input datasets
├── padel/                # PaDEL-Descriptor setup
├── results/              # Model outputs and evaluation plots
├── scripts/              # Python scripts (train, merge, evaluate)
├── Snakefile             # Snakemake workflow
├── config.yaml           # File paths for pipeline
├── environment.yml       # Conda environment (optional)
├── .gitignore            # Git exclusion rules
└── Acetylcholinesterase_Bioactivity_Pipeline.ipynb # Main notebook
