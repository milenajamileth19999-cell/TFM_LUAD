###############################################################
# SCRIPT 4: MODELADO PRONÓSTICO MEDIANTE LASSO-COX
#
# Objetivo:
# Aplicar Cox univariado y LASSO-Cox de forma independiente
# a RNA-seq, metilación, proteómica y fosfoproteómica.
###############################################################


# 1. Directorio y paquetes

# Se define la carpeta del proyecto y se cargan los paquetes
# necesarios para el análisis de supervivencia y LASSO-Cox.

project_directory <- "C:/Users/milen/OneDrive/Desktop/tesis/TFM_LUAD"

setwd(project_directory)

library(survival)
library(glmnet)


# 2. Cargar los datos preprocesados

# Se recuperan las matrices generadas en el Script 3. Todas tienen
# pacientes en las filas, variables moleculares en las columnas
# y el mismo orden que la información de supervivencia.

processed_data <- readRDS(
  "results/LUAD_processed_data.rds"
)

survival_complete <- processed_data$survival
rna_complete <- processed_data$rna
methylation_complete <- processed_data$methylation
proteomics_complete <- processed_data$proteomics
phosphoproteomics_complete <- processed_data$phosphoproteomics


# 3. Función para ejecutar Cox univariado

# Se ajusta un modelo de Cox para cada variable molecular.
# Se conserva el coeficiente, el hazard ratio y el valor p.
# tryCatch permite continuar si una variable concreta no puede
# producir un modelo válido.

run_univariate_cox <- function(
    omic_matrix,
    survival_data
) {
  
  results <- lapply(
    seq_len(ncol(omic_matrix)),
    function(i) {
      
      model <- tryCatch(
        coxph(
          Surv(
            survival_data$OS_days,
            survival_data$OS_event
          ) ~ omic_matrix[, i],
          ties = "efron"
        ),
        error = function(e) NULL
      )
      
      if (is.null(model)) {
        return(
          c(
            Coefficient = NA_real_,
            Hazard_ratio = NA_real_,
            P_value = NA_real_
          )
        )
      }
      
      model_summary <- summary(model)
      
      c(
        Coefficient =
          model_summary$coefficients[1, "coef"],
        
        Hazard_ratio =
          model_summary$coefficients[1, "exp(coef)"],
        
        P_value =
          model_summary$coefficients[1, "Pr(>|z|)"]
      )
    }
  )
  
  results <- do.call(
    rbind,
    results
  )
  
  results <- data.frame(
    Feature_ID = colnames(omic_matrix),
    Coefficient = results[, "Coefficient"],
    Hazard_ratio = results[, "Hazard_ratio"],
    P_value = results[, "P_value"],
    stringsAsFactors = FALSE
  )
  
  results <- results[
    order(
      results$P_value,
      na.last = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  rownames(results) <- NULL
  
  return(results)
}


# 4. Función para extraer variables seleccionadas por LASSO

# Se recuperan las variables con coeficientes distintos de cero.
# El signo del coeficiente indica si contribuyen a aumentar o
# disminuir el puntaje de riesgo pronóstico.

extract_lasso_features <- function(
    lasso_model,
    lambda_value
) {
  
  coefficients <- as.matrix(
    coef(
      lasso_model,
      s = lambda_value
    )
  )
  
  selected <- which(
    coefficients[, 1] != 0
  )
  
  if (length(selected) == 0) {
    
    return(
      data.frame(
        Feature_ID = character(0),
        Coefficient = numeric(0),
        Absolute_coefficient = numeric(0),
        Association = character(0),
        stringsAsFactors = FALSE
      )
    )
  }
  
  selected_features <- data.frame(
    Feature_ID = rownames(coefficients)[selected],
    Coefficient = coefficients[selected, 1],
    stringsAsFactors = FALSE
  )
  
  selected_features$Absolute_coefficient <- abs(
    selected_features$Coefficient
  )
  
  selected_features$Association <- ifelse(
    selected_features$Coefficient > 0,
    "Mayor riesgo",
    "Menor riesgo"
  )
  
  selected_features <- selected_features[
    order(
      -selected_features$Absolute_coefficient
    ),
    ,
    drop = FALSE
  ]
  
  rownames(selected_features) <- NULL
  
  return(selected_features)
}


# 5. Función para modelar una modalidad ómica

# Primero se realiza el análisis de Cox univariado y se conservan
# las variables con p < 0,05. Las 200 con menor valor p se utilizan
# como entrada del LASSO-Cox. El parámetro lambda se selecciona
# mediante validación cruzada estratificada de cinco pliegues.

fit_omic_model <- function(
    omic_matrix,
    survival_data,
    omic_name,
    seed = 2026
) {
  
  message(
    "Modelando: ",
    omic_name
  )
  
  univariate_results <- run_univariate_cox(
    omic_matrix,
    survival_data
  )
  
  significant_results <- univariate_results[
    !is.na(univariate_results$P_value) &
      univariate_results$P_value < 0.05,
    ,
    drop = FALSE
  ]
  
  number_candidates <- min(
    200,
    nrow(significant_results)
  )
  
  if (number_candidates == 0) {
    stop(
      paste(
        "No se encontraron variables significativas en",
        omic_name
      )
    )
  }
  
  top_features <- significant_results[
    seq_len(number_candidates),
    ,
    drop = FALSE
  ]
  
  selected_feature_ids <- top_features$Feature_ID
  
  selected_matrix <- omic_matrix[
    ,
    selected_feature_ids,
    drop = FALSE
  ]
  
  
  # Crear pliegues con eventos y censuras distribuidos
  # aproximadamente de manera equilibrada.
  
  set.seed(seed)
  
  fold_id <- integer(
    nrow(survival_data)
  )
  
  event_positions <- which(
    survival_data$OS_event == 1
  )
  
  censored_positions <- which(
    survival_data$OS_event == 0
  )
  
  fold_id[event_positions] <- sample(
    rep(
      1:5,
      length.out = length(event_positions)
    )
  )
  
  fold_id[censored_positions] <- sample(
    rep(
      1:5,
      length.out = length(censored_positions)
    )
  )
  
  
  # Ajustar el modelo LASSO-Cox.
  
  lasso_model <- cv.glmnet(
    x = as.matrix(selected_matrix),
    y = Surv(
      survival_data$OS_days,
      survival_data$OS_event
    ),
    family = "cox",
    alpha = 1,
    foldid = fold_id,
    standardize = TRUE
  )
  
  
  # Extraer las variables seleccionadas con ambos criterios.
  
  selected_lambda_min <- extract_lasso_features(
    lasso_model,
    "lambda.min"
  )
  
  selected_lambda_1se <- extract_lasso_features(
    lasso_model,
    "lambda.1se"
  )
  
  
  # Crear una fila resumen de la modalidad.
  
  model_summary <- data.frame(
    Modality = omic_name,
    Patients = nrow(omic_matrix),
    Events = sum(survival_data$OS_event == 1),
    Initial_features = ncol(omic_matrix),
    Valid_univariate_models =
      sum(!is.na(univariate_results$P_value)),
    Features_with_p_less_0_05 =
      nrow(significant_results),
    Features_entered_in_LASSO =
      nrow(top_features),
    Lambda_min =
      lasso_model$lambda.min,
    Features_selected_lambda_min =
      nrow(selected_lambda_min),
    Lambda_1se =
      lasso_model$lambda.1se,
    Features_selected_lambda_1se =
      nrow(selected_lambda_1se),
    stringsAsFactors = FALSE
  )
  
  
  # Se conservan solamente los objetos utilizados en los
  # análisis posteriores.
  
  return(
    list(
      Omic_name = omic_name,
      Omic_matrix = omic_matrix,
      Survival_data = survival_data,
      Univariate_results = univariate_results,
      Significant_results = significant_results,
      Top_features = top_features,
      Selected_feature_ids = selected_feature_ids,
      Selected_matrix = selected_matrix,
      LASSO_model = lasso_model,
      Selected_features_lambda_min = selected_lambda_min,
      Selected_features_lambda_1se = selected_lambda_1se,
      Summary = model_summary
    )
  )
}


# 6. Ajustar los modelos de las cuatro modalidades

# La misma estrategia se aplica de forma independiente a cada
# modalidad, permitiendo una comparación metodológicamente justa.

rna_model <- fit_omic_model(
  rna_complete,
  survival_complete,
  "RNA-seq",
  seed = 2026
)

methylation_model <- fit_omic_model(
  methylation_complete,
  survival_complete,
  "Methylation",
  seed = 2026
)

proteomics_model <- fit_omic_model(
  proteomics_complete,
  survival_complete,
  "Proteomics",
  seed = 2026
)

phosphoproteomics_model <- fit_omic_model(
  phosphoproteomics_complete,
  survival_complete,
  "Phosphoproteomics",
  seed = 2026
)


# 7. Tabla resumen del modelado

# Esta tabla muestra la reducción de dimensionalidad obtenida
# mediante Cox univariado y LASSO-Cox en cada modalidad.

modeling_summary <- rbind(
  rna_model$Summary,
  methylation_model$Summary,
  proteomics_model$Summary,
  phosphoproteomics_model$Summary
)

rownames(modeling_summary) <- NULL

print(modeling_summary)


# 8. Tabla de variables seleccionadas

# Se reúnen las firmas obtenidas con lambda.min para conocer
# cuáles variables conforman el modelo global de cada modalidad.

selected_features_summary <- rbind(
  transform(
    rna_model$Selected_features_lambda_min,
    Modality = "RNA-seq"
  ),
  transform(
    methylation_model$Selected_features_lambda_min,
    Modality = "Methylation"
  ),
  transform(
    proteomics_model$Selected_features_lambda_min,
    Modality = "Proteomics"
  ),
  transform(
    phosphoproteomics_model$Selected_features_lambda_min,
    Modality = "Phosphoproteomics"
  )
)

selected_features_summary <- selected_features_summary[
  ,
  c(
    "Modality",
    "Feature_ID",
    "Coefficient",
    "Absolute_coefficient",
    "Association"
  ),
  drop = FALSE
]

rownames(selected_features_summary) <- NULL

print(selected_features_summary)


# 9. Guardar los resultados

# Los cuatro modelos se guardan juntos porque los siguientes
# scripts utilizarán sus coeficientes, matrices y firmas.

omic_models <- list(
  rna = rna_model,
  methylation = methylation_model,
  proteomics = proteomics_model,
  phosphoproteomics = phosphoproteomics_model,
  summary = modeling_summary
)

saveRDS(
  omic_models,
  "results/LUAD_omic_models.rds"
)

write.csv(
  modeling_summary,
  "results/modeling_summary.csv",
  row.names = FALSE
)

write.csv(
  selected_features_summary,
  "results/LASSO_selected_features.csv",
  row.names = FALSE
)
