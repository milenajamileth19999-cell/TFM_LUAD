###############################################################
# TFM LUAD
# SCRIPT 7: COMPARACIÓN FINAL DE LAS MODALIDADES ÓMICAS
#
# Autora: Milena Vera García
#
# Objetivo:
# Integrar y comparar los resultados del modelado, la evaluación
# aparente y la validación cruzada anidada de RNA-seq,
# metilación, proteómica y fosfoproteómica.
#
# Este script no ajusta nuevos modelos.
###############################################################


# 1. Directorio del proyecto

# Se define la carpeta principal para cargar los resultados
# obtenidos en los scripts anteriores.

project_directory <- "C:/Users/milen/OneDrive/Desktop/tesis/TFM_LUAD"

setwd(project_directory)


# 2. Cargar los resultados previos

# Se recuperan los resultados del modelado LASSO-Cox,
# la evaluación aparente y la validación cruzada anidada.

omic_models <- readRDS(
  "results/LUAD_omic_models.rds"
)

apparent_evaluations <- readRDS(
  "results/LUAD_apparent_evaluations.rds"
)

nested_validation_results <- readRDS(
  "results/LUAD_nested_validation.rds"
)


# 3. Resumir la selección de variables

# Se reúnen las principales etapas de reducción de dimensionalidad
# y selección de variables de las cuatro modalidades.

feature_selection_summary <- rbind(
  omic_models$rna$Summary,
  omic_models$methylation$Summary,
  omic_models$proteomics$Summary,
  omic_models$phosphoproteomics$Summary
)

rownames(feature_selection_summary) <- NULL

feature_selection_summary <- feature_selection_summary[
  ,
  c(
    "Modality",
    "Patients",
    "Events",
    "Initial_features",
    "Valid_univariate_models",
    "Features_with_p_less_0_05",
    "Features_entered_in_LASSO",
    "Features_selected_lambda_min",
    "Features_selected_lambda_1se"
  ),
  drop = FALSE
]


# 4. Extraer la evaluación aparente

# Se recuperan las medidas obtenidas al evaluar los modelos
# sobre la misma cohorte utilizada para construirlos.

apparent_summary <- apparent_evaluations$summary

apparent_summary <- apparent_summary[
  ,
  c(
    "Modality",
    "Hazard_ratio",
    "CI95_lower",
    "CI95_upper",
    "Cox_p_value",
    "Logrank_p_value",
    "Apparent_C_index",
    "Apparent_C_index_SE",
    "Evaluation_status"
  ),
  drop = FALSE
]


# 5. Extraer la validación cruzada anidada

# Se recupera el rendimiento fuera de muestra de cada modalidad.

validation_summary <- nested_validation_results$summary

validation_summary <- validation_summary[
  ,
  c(
    "Modality",
    "Mean_validated_C_index",
    "SD_validated_C_index",
    "Minimum_validated_C_index",
    "Maximum_validated_C_index",
    "Valid_folds",
    "Total_folds"
  ),
  drop = FALSE
]


# 6. Integrar los resultados

# Las tres fuentes de resultados se unen mediante el nombre
# de la modalidad ómica.

final_comparison <- merge(
  feature_selection_summary,
  apparent_summary,
  by = "Modality"
)

final_comparison <- merge(
  final_comparison,
  validation_summary,
  by = "Modality"
)


# 7. Calcular el optimismo del modelo

# El optimismo corresponde a la diferencia entre el C-index
# aparente y el C-index obtenido mediante validación anidada.
# Valores mayores indican una mayor pérdida de rendimiento
# al evaluar el modelo fuera de la cohorte de entrenamiento.

final_comparison$C_index_optimism <- round(
  final_comparison$Apparent_C_index -
    final_comparison$Mean_validated_C_index,
  4
)


# 8. Calcular el porcentaje de pliegues válidos

# Se calcula el porcentaje de las 15 evaluaciones externas
# en las que fue posible obtener un C-index.

final_comparison$Valid_fold_percentage <- round(
  final_comparison$Valid_folds /
    final_comparison$Total_folds * 100,
  2
)


# 9. Ordenar las modalidades

# Las modalidades se ordenan de mayor a menor según su
# rendimiento validado.

final_comparison <- final_comparison[
  order(
    -final_comparison$Mean_validated_C_index,
    na.last = TRUE
  ),
  ,
  drop = FALSE
]

rownames(final_comparison) <- NULL

print(final_comparison)


# 10. Crear la tabla principal para la redacción

# Se conservan únicamente los indicadores que aportan información
# directa para comparar el rendimiento pronóstico de las cuatro
# modalidades ómicas.

manuscript_summary <- data.frame(
  
  Modalidad =
    final_comparison$Modality,
  
  Variables_LASSO =
    final_comparison$Features_selected_lambda_min,
  
  HR_IC95 = paste0(
    round(
      final_comparison$Hazard_ratio,
      2
    ),
    " (",
    round(
      final_comparison$CI95_lower,
      2
    ),
    "–",
    round(
      final_comparison$CI95_upper,
      2
    ),
    ")"
  ),
  
  Valor_p_Cox =
    final_comparison$Cox_p_value,
  
  C_index_aparente =
    round(
      final_comparison$Apparent_C_index,
      3
    ),
  
  C_index_validado =
    round(
      final_comparison$Mean_validated_C_index,
      3
    ),
  
  Optimismo =
    round(
      final_comparison$C_index_optimism,
      3
    ),
  
  stringsAsFactors = FALSE
)

print(manuscript_summary)