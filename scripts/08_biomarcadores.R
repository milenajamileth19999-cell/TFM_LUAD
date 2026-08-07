###############################################################
# TFM LUAD
# SCRIPT 8: EXTRACCIÓN Y PRIORIZACIÓN DE BIOMARCADORES CANDIDATOS
#
# Autora: Milena Vera García
#
# Objetivo:
# Organizar las variables seleccionadas mediante LASSO-Cox en
# la cohorte completa y evaluar su estabilidad durante la
# validación cruzada anidada.
#
# Para cada variable se integran:
# - coeficiente LASSO;
# - dirección de la asociación;
# - hazard ratio del Cox univariado;
# - valor p del Cox univariado;
# - frecuencia de selección en los 15 pliegues externos.
#
# Este script no ajusta nuevos modelos.
###############################################################


# 1. Directorio del proyecto

# Se define la carpeta principal para cargar los resultados
# obtenidos en los Scripts 4 y 6.

project_directory <- "C:/Users/milen/OneDrive/Desktop/tesis/TFM_LUAD"

setwd(project_directory)


# 2. Cargar los resultados previos

# Se recuperan los modelos globales y los resultados de la
# validación cruzada anidada.

omic_models <- readRDS(
  "results/LUAD_omic_models.rds"
)

nested_validation_results <- readRDS(
  "results/LUAD_nested_validation.rds"
)


# 3. Función para organizar los biomarcadores candidatos

# Para cada modalidad se toman las variables seleccionadas por
# lambda.min en la cohorte completa y se integran sus resultados
# univariados y su frecuencia de selección en validación.

extract_candidate_biomarkers <- function(
    model_results,
    feature_frequency,
    modality
) {
  
  selected_features <-
    model_results$Selected_features_lambda_min
  
  selected_features <- selected_features[
    ,
    c(
      "Feature_ID",
      "Coefficient",
      "Absolute_coefficient",
      "Association"
    ),
    drop = FALSE
  ]
  
  colnames(selected_features) <- c(
    "Feature_ID",
    "LASSO_coefficient",
    "Absolute_LASSO_coefficient",
    "Association"
  )
  
  
  # Recuperar el HR y el valor p del análisis Cox univariado.
  
  univariate_results <- model_results$Univariate_results[
    ,
    c(
      "Feature_ID",
      "Hazard_ratio",
      "P_value"
    ),
    drop = FALSE
  ]
  
  colnames(univariate_results) <- c(
    "Feature_ID",
    "Univariate_HR",
    "Univariate_p_value"
  )
  
  
  # Integrar la información del modelo global.
  
  candidate_table <- merge(
    selected_features,
    univariate_results,
    by = "Feature_ID",
    all.x = TRUE,
    sort = FALSE
  )
  
  
  # Integrar la frecuencia observada durante la validación.
  
  candidate_table <- merge(
    candidate_table,
    feature_frequency[
      ,
      c(
        "Feature_ID",
        "Times_selected",
        "Selection_percentage"
      ),
      drop = FALSE
    ],
    by = "Feature_ID",
    all.x = TRUE,
    sort = FALSE
  )
  
  
  # Las variables que pertenecen a la firma global pero que
  # no reaparecieron en validación reciben frecuencia cero.
  
  candidate_table$Times_selected[
    is.na(candidate_table$Times_selected)
  ] <- 0
  
  candidate_table$Selection_percentage[
    is.na(candidate_table$Selection_percentage)
  ] <- 0
  
  
  candidate_table$Modality <- modality
  
  
  # Las variables se ordenan primero por estabilidad durante la
  # validación y después por magnitud del coeficiente LASSO.
  
  candidate_table <- candidate_table[
    order(
      -candidate_table$Selection_percentage,
      -candidate_table$Absolute_LASSO_coefficient
    ),
    ,
    drop = FALSE
  ]
  
  
  candidate_table <- candidate_table[
    ,
    c(
      "Modality",
      "Feature_ID",
      "LASSO_coefficient",
      "Association",
      "Univariate_HR",
      "Univariate_p_value",
      "Times_selected",
      "Selection_percentage"
    ),
    drop = FALSE
  ]
  
  rownames(candidate_table) <- NULL
  
  return(candidate_table)
}


# 4. Extraer candidatos de RNA-seq

rna_candidates <- extract_candidate_biomarkers(
  model_results =
    omic_models$rna,
  
  feature_frequency =
    nested_validation_results$rna$Feature_frequency,
  
  modality =
    "RNA-seq"
)


# 5. Extraer candidatos de metilación

methylation_candidates <- extract_candidate_biomarkers(
  model_results =
    omic_models$methylation,
  
  feature_frequency =
    nested_validation_results$methylation$Feature_frequency,
  
  modality =
    "Methylation"
)


# 6. Extraer candidatos de proteómica

proteomics_candidates <- extract_candidate_biomarkers(
  model_results =
    omic_models$proteomics,
  
  feature_frequency =
    nested_validation_results$proteomics$Feature_frequency,
  
  modality =
    "Proteomics"
)


# 7. Extraer candidatos de fosfoproteómica

phosphoproteomics_candidates <- extract_candidate_biomarkers(
  model_results =
    omic_models$phosphoproteomics,
  
  feature_frequency =
    nested_validation_results$phosphoproteomics$Feature_frequency,
  
  modality =
    "Phosphoproteomics"
)


# 8. Unir los biomarcadores candidatos

# Se construye una única tabla con las variables seleccionadas
# en el modelo global de las cuatro modalidades.

candidate_biomarkers <- rbind(
  rna_candidates,
  methylation_candidates,
  proteomics_candidates,
  phosphoproteomics_candidates
)

rownames(candidate_biomarkers) <- NULL

print(candidate_biomarkers)


# 9. Resumen por modalidad

# Se resume el número de variables de la firma global y cuántas
# de ellas reaparecieron al menos una vez durante la validación.

biomarker_summary <- data.frame(
  
  Modality = c(
    "RNA-seq",
    "Methylation",
    "Proteomics",
    "Phosphoproteomics"
  ),
  
  LASSO_candidates = c(
    nrow(rna_candidates),
    nrow(methylation_candidates),
    nrow(proteomics_candidates),
    nrow(phosphoproteomics_candidates)
  ),
  
  Candidates_seen_in_validation = c(
    sum(rna_candidates$Times_selected > 0),
    sum(methylation_candidates$Times_selected > 0),
    sum(proteomics_candidates$Times_selected > 0),
    sum(phosphoproteomics_candidates$Times_selected > 0)
  ),
  
  Maximum_selection_percentage = c(
    max(rna_candidates$Selection_percentage),
    max(methylation_candidates$Selection_percentage),
    max(proteomics_candidates$Selection_percentage),
    max(phosphoproteomics_candidates$Selection_percentage)
  ),
  
  stringsAsFactors = FALSE
)

print(biomarker_summary)


# 10. Seleccionar los candidatos prioritarios

# Para la interpretación biológica se priorizan las variables
# pertenecientes a la firma global que además reaparecieron
# durante la validación cruzada.

rna_priority <- rna_candidates[
  rna_candidates$Times_selected > 0,
  ,
  drop = FALSE
]

methylation_priority <- methylation_candidates[
  methylation_candidates$Times_selected > 0,
  ,
  drop = FALSE
]

proteomics_priority <- proteomics_candidates[
  proteomics_candidates$Times_selected > 0,
  ,
  drop = FALSE
]

phosphoproteomics_priority <- phosphoproteomics_candidates[
  phosphoproteomics_candidates$Times_selected > 0,
  ,
  drop = FALSE
]


# Se conservan como máximo las 10 variables más estables de
# cada modalidad para la posterior interpretación biológica.

rna_priority <- head(
  rna_priority,
  10
)

methylation_priority <- head(
  methylation_priority,
  10
)

proteomics_priority <- head(
  proteomics_priority,
  10
)

phosphoproteomics_priority <- head(
  phosphoproteomics_priority,
  10
)


priority_biomarkers <- rbind(
  rna_priority,
  methylation_priority,
  proteomics_priority,
  phosphoproteomics_priority
)

rownames(priority_biomarkers) <- NULL

print(priority_biomarkers)


# 11. Crear carpeta de resultados

# Las tablas utilizadas para la interpretación biológica se
# guardan juntas para facilitar su localización.

biomarker_directory <- file.path(
  "results",
  "tables",
  "Biomarkers"
)

dir.create(
  biomarker_directory,
  recursive = TRUE,
  showWarnings = FALSE
)


# 12. Guardar las tablas necesarias

# Se conserva la tabla completa de candidatos y la tabla reducida
# que se utilizará para la interpretación biológica.

write.csv(
  candidate_biomarkers,
  file.path(
    biomarker_directory,
    "candidate_biomarkers.csv"
  ),
  row.names = FALSE
)

write.csv(
  biomarker_summary,
  file.path(
    biomarker_directory,
    "biomarker_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  priority_biomarkers,
  file.path(
    biomarker_directory,
    "priority_biomarkers.csv"
  ),
  row.names = FALSE
)


# 13. Guardar el objeto necesario para el siguiente script

# El archivo RDS conserva las tablas completas por modalidad y
# los candidatos priorizados para su posterior anotación biológica.

biomarker_results <- list(
  rna = rna_candidates,
  methylation = methylation_candidates,
  proteomics = proteomics_candidates,
  phosphoproteomics = phosphoproteomics_candidates,
  summary = biomarker_summary,
  priority_biomarkers = priority_biomarkers
)

saveRDS(
  biomarker_results,
  "results/LUAD_biomarker_results.rds"
)
