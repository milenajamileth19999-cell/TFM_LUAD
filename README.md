# TFM_LUAD

## Comparative evaluation of transcriptomics, DNA methylation, proteomics and phosphoproteomics for prognosis prediction in lung adenocarcinoma

### Project description

This repository contains the R scripts developed for the Master's Thesis entitled:

**"Comparative evaluation of transcriptomics, DNA methylation, proteomics and phosphoproteomics for prognosis prediction in lung adenocarcinoma."**

The study compares the prognostic performance of four molecular modalities from the LinkedOmicsKB LUAD cohort using survival analysis and machine learning approaches.

---

## Data source

Data were obtained from **LinkedOmicsKB** for the LUAD (Lung Adenocarcinoma) cohort.

| Molecular modality | Patients | Variables |
|-------------------|---------:|----------:|
| Survival | 110 | 5 |
| Transcriptomics (RNA-seq) | 110 | 60,669 |
| DNA methylation | 105 | 13,178 |
| Proteomics | 110 | 12,433 |
| Phosphoproteomics | 110 | 61,705 |

After harmonization of patient identifiers across all molecular layers, a final integrated cohort of **105 patients** was obtained.

---

## Analytical workflow

✔ Cohort construction

✔ Multi-omics integration

✔ Data preprocessing

✔ Feature selection using Cox and LASSO

✔ Prognostic model development

✔ Repeated nested cross-validation

✔ Comparative evaluation of molecular modalities

✔ Biomarker prioritization

---

## Repository structure

```text
TFM_LUAD/
│
├── README.md
│
├── scripts/
├── 01_construccion_cohorte.R
├── 02_cohorte_comun.R
├── 03_procesamiento_matrices_omicas.R
├── 04_modelado.R
├── 05_evaluacion_aparente.R
├── 06_validacion_cruzada.R
├── 07_comparacion_modalidades.R
├── 08_biomarcadores.R
└── 09_anotacion_biomarcadores.R


```

---

## Software

- R 4.5.1
- RStudio
- Survival
- glmnet
- survminer
- randomForestSRC
- dplyr
- tidyr
- ggplot2

---

## Author

**Milena Vera García**

Master's Thesis  
Master's Degree in Bioinformatics
Master's Degree in Bioinformatics
