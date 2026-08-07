###############################################################
# SCRIPT 3: PREPROCESAMIENTO DE LAS MATRICES ÓMICAS
#
# Objetivo:
# Preparar las matrices de RNA-seq, metilación, proteómica y
# fosfoproteómica para el modelado pronóstico.
###############################################################


# 1. Directorio del proyecto

# Se define la carpeta principal para localizar los datos cargados
# previamente y guardar las matrices preprocesadas.

project_directory <- "C:/Users/milen/OneDrive/Desktop/tesis/TFM_LUAD"

setwd(project_directory)


# 2. Preparar la información de supervivencia

# Se convierten el tiempo y el evento a formato numérico. Después
# se conservan únicamente los pacientes con tiempo de seguimiento
# positivo y estado de supervivencia definido como 0 o 1.

survival_common$OS_days <- as.numeric(
  survival_common$OS_days
)

survival_common$OS_event <- as.numeric(
  survival_common$OS_event
)

survival_complete <- survival_common[
  !is.na(survival_common$OS_days) &
    survival_common$OS_days > 0 &
    !is.na(survival_common$OS_event) &
    survival_common$OS_event %in% c(0, 1),
  ,
  drop = FALSE
]

analysis_ids <- rownames(survival_complete)


# 3. Resumen de supervivencia

# Esta tabla documenta el número de pacientes que serán incluidos
# en el modelado, así como la distribución de eventos y censuras.

survival_summary <- data.frame(
  Pacientes = nrow(survival_complete),
  Eventos = sum(survival_complete$OS_event == 1),
  Censurados = sum(survival_complete$OS_event == 0),
  stringsAsFactors = FALSE
)

print(survival_summary)


# 4. Umbrales de valores faltantes

# RNA-seq no presenta valores faltantes. Para metilación y
# proteómica se eliminan las variables con más del 20 % de datos
# ausentes. En fosfoproteómica se admite hasta un 50 % debido a
# su mayor proporción inicial de valores faltantes.

missing_thresholds <- c(
  rnaseq = 0,
  methylation = 0.20,
  proteomics = 0.20,
  phosphoproteomics = 0.50
)


# 5. Función de preprocesamiento

# La misma función se aplica a las cuatro modalidades para asegurar
# que todas sean procesadas mediante una estrategia comparable.

prepare_omic_matrix <- function(
    omic_data,
    patient_ids,
    modality,
    missing_threshold
) {
  
  feature_ids <- omic_data$idx
  
  omic_matrix <- as.matrix(
    omic_data[
      ,
      patient_ids,
      drop = FALSE
    ]
  )
  
  storage.mode(omic_matrix) <- "numeric"
  
  rownames(omic_matrix) <- feature_ids
  
  initial_missing <- sum(
    is.na(omic_matrix)
  )
  
  missing_fraction <- rowMeans(
    is.na(omic_matrix)
  )
  
  keep_features <- missing_fraction <= missing_threshold
  
  omic_filtered <- omic_matrix[
    keep_features,
    ,
    drop = FALSE
  ]
  
  missing_before_imputation <- sum(
    is.na(omic_filtered)
  )
  
  for (i in seq_len(nrow(omic_filtered))) {
    
    missing_positions <- is.na(
      omic_filtered[i, ]
    )
    
    if (any(missing_positions)) {
      
      omic_filtered[i, missing_positions] <- median(
        omic_filtered[i, ],
        na.rm = TRUE
      )
    }
  }
  
  omic_transposed <- t(
    omic_filtered
  )
  
  feature_variance <- apply(
    omic_transposed,
    2,
    var
  )
  
  omic_final <- omic_transposed[
    ,
    feature_variance > 0,
    drop = FALSE
  ]
  
  summary_table <- data.frame(
    Modalidad = modality,
    Pacientes = nrow(omic_final),
    Variables_iniciales = nrow(omic_matrix),
    Valores_faltantes_iniciales = initial_missing,
    Umbral_faltantes_porcentaje = missing_threshold * 100,
    Variables_eliminadas_por_faltantes = sum(!keep_features),
    Valores_imputados = missing_before_imputation,
    Variables_eliminadas_por_varianza_cero =
      sum(feature_variance == 0),
    Variables_finales = ncol(omic_final),
    stringsAsFactors = FALSE
  )
  
  list(
    matrix = omic_final,
    summary = summary_table
  )
}


# 6. Preprocesar las cuatro modalidades

# Cada matriz queda organizada con pacientes en las filas y
# variables moleculares en las columnas, formato requerido para
# los modelos de supervivencia posteriores.

rna_result <- prepare_omic_matrix(
  omic_data = rnaseq,
  patient_ids = analysis_ids,
  modality = "RNA-seq",
  missing_threshold = missing_thresholds["rnaseq"]
)

methylation_result <- prepare_omic_matrix(
  omic_data = methylation,
  patient_ids = analysis_ids,
  modality = "Metilación",
  missing_threshold = missing_thresholds["methylation"]
)

proteomics_result <- prepare_omic_matrix(
  omic_data = proteomics,
  patient_ids = analysis_ids,
  modality = "Proteómica",
  missing_threshold = missing_thresholds["proteomics"]
)

phosphoproteomics_result <- prepare_omic_matrix(
  omic_data = phosphoproteomics,
  patient_ids = analysis_ids,
  modality = "Fosfoproteómica",
  missing_threshold = missing_thresholds["phosphoproteomics"]
)


# 7. Matrices finales

rna_complete <- rna_result$matrix

methylation_complete <- methylation_result$matrix

proteomics_complete <- proteomics_result$matrix

phosphoproteomics_complete <-
  phosphoproteomics_result$matrix


# 8. Resumen del preprocesamiento

# La tabla reúne las variables iniciales, las eliminadas y las
# conservadas en cada modalidad. Esta será la principal evidencia
# del proceso de control de calidad aplicado a los datos.

preprocessing_summary <- rbind(
  rna_result$summary,
  methylation_result$summary,
  proteomics_result$summary,
  phosphoproteomics_result$summary
)

rownames(preprocessing_summary) <- NULL

print(preprocessing_summary)


# 9. Guardar los datos preprocesados

# Las matrices y la información de supervivencia se agrupan en un
# único archivo porque serán utilizadas conjuntamente en los
# scripts de selección de variables y modelado pronóstico.

processed_data <- list(
  survival = survival_complete,
  rna = rna_complete,
  methylation = methylation_complete,
  proteomics = proteomics_complete,
  phosphoproteomics = phosphoproteomics_complete
)

saveRDS(
  processed_data,
  "results/LUAD_processed_data.rds"
)

write.csv(
  survival_summary,
  "results/survival_summary.csv",
  row.names = FALSE
)

write.csv(
  preprocessing_summary,
  "results/preprocessing_summary.csv",
  row.names = FALSE
)
