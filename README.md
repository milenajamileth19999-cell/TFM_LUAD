# TFM_LUAD

## Comparative evaluation of RNA-seq, methylation, proteomics and phosphoproteomics for prognosis prediction in lung adenocarcinoma

### Project description

This project evaluates the prognostic value of four molecular modalities available in LinkedOmicsKB:

- RNA sequencing (RNA-seq)
- DNA methylation
- Proteomics
- Phosphoproteomics

The objective is to compare the predictive performance of each modality independently using survival analysis and machine learning approaches.

---

## Data source

Data were obtained from LinkedOmicsKB for the LUAD (Lung Adenocarcinoma) cohort.

Available modalities:

| Modality | Patients | Variables |
|-----------|-----------|-----------|
| Survival | 110 | 5 |
| RNA-seq | 110 | 60,669 |
| Methylation | 105 | 13,178 |
| Proteomics | 110 | 12,433 |
| Phosphoproteomics | 110 | 61,705 |

After harmonization of patient identifiers and integration across modalities, a final cohort of 105 patients was obtained.

---

## Current workflow

- [x] Data acquisition
- [x] Patient identifier harmonization
- [x] Cohort construction
- [ ] Feature filtering
- [ ] Survival modeling
- [ ] Model validation
- [ ] Comparative evaluation

---

## Repository structure

```text
TFM_LUAD/
│
├── scripts/
│   └── 01_construccion_cohorte.R
│
├── README.md
└── .gitignore
```

---

## Author

Milena Vera García

Master's Degree in Bioinformatics
