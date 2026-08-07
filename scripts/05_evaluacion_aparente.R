###############################################################
# SCRIPT 5: EVALUACIÓN APARENTE DE LOS MODELOS ÓMICOS
#
# Objetivo:
# Evaluar el rendimiento aparente de los modelos LASSO-Cox
# construidos para RNA-seq, metilación, proteómica y
# fosfoproteómica.
#
# La evaluación se realiza sobre los mismos pacientes utilizados
# para construir cada modelo.
###############################################################


# 1. Directorio y paquetes

# Se define la carpeta del proyecto y se cargan los paquetes
# necesarios para calcular las medidas de supervivencia y obtener
# las predicciones de los modelos LASSO-Cox.

project_directory <- "C:/Users/milen/OneDrive/Desktop/tesis/TFM_LUAD"

setwd(project_directory)

library(survival)
library(glmnet)


# 2. Cargar los modelos

# Se recuperan los modelos obtenidos en el Script 4, junto con
# las matrices seleccionadas y la información de supervivencia.

omic_models <- readRDS(
  "results/LUAD_omic_models.rds"
)


# 3. Función de evaluación aparente

# Para cada modalidad se calcula un puntaje de riesgo, se divide
# a los pacientes según la mediana del puntaje y se evalúa la
# asociación con la supervivencia mediante Kaplan-Meier, log-rank,
# Cox y C-index.

evaluate_apparent_model <- function(
    model_results,
    modality
) {
  
  selected_features <- nrow(
    model_results$Selected_features_lambda_min
  )
  
  
  # Cuando LASSO no selecciona variables, no es posible calcular
  # un puntaje de riesgo ni evaluar el modelo.
  
  if (selected_features == 0) {
    
    summary_table <- data.frame(
      Modality = modality,
      Patients = nrow(model_results$Omic_matrix),
      Events = sum(model_results$Survival_data$OS_event == 1),
      Censored = sum(model_results$Survival_data$OS_event == 0),
      Low_risk_patients = NA_integer_,
      High_risk_patients = NA_integer_,
      Logrank_p_value = NA_real_,
      Hazard_ratio = NA_real_,
      CI95_lower = NA_real_,
      CI95_upper = NA_real_,
      Cox_p_value = NA_real_,
      Apparent_C_index = NA_real_,
      Apparent_C_index_SE = NA_real_,
      Evaluation_status = "LASSO no seleccionó variables",
      stringsAsFactors = FALSE
    )
    
    return(
      list(
        Risk_table = NULL,
        Kaplan_Meier_model = NULL,
        Summary = summary_table
      )
    )
  }
  
  
  # Calcular el puntaje pronóstico de cada paciente.
  
  risk_score <- as.numeric(
    predict(
      model_results$LASSO_model,
      newx = as.matrix(model_results$Selected_matrix),
      s = "lambda.min",
      type = "link"
    )
  )
  
  
  # Estandarizar el puntaje para interpretar el hazard ratio
  # como el cambio asociado a una desviación estándar.
  
  standardized_score <- as.numeric(
    scale(risk_score)
  )
  
  
  # Clasificar a los pacientes en grupos de riesgo utilizando
  # la mediana del puntaje como punto de corte.
  
  risk_group <- ifelse(
    risk_score <= median(risk_score),
    "Bajo riesgo",
    "Alto riesgo"
  )
  
  risk_group <- factor(
    risk_group,
    levels = c(
      "Bajo riesgo",
      "Alto riesgo"
    )
  )
  
  
  # Crear la tabla individual con supervivencia y riesgo.
  
  risk_table <- data.frame(
    Patient_ID = rownames(model_results$Omic_matrix),
    OS_days = model_results$Survival_data$OS_days,
    OS_event = model_results$Survival_data$OS_event,
    Risk_score = risk_score,
    Standardized_score = standardized_score,
    Risk_group = risk_group,
    stringsAsFactors = FALSE
  )
  
  
  # Estimar las curvas de supervivencia para los grupos de riesgo.
  
  kaplan_meier_model <- survfit(
    Surv(
      OS_days,
      OS_event
    ) ~ Risk_group,
    data = risk_table
  )
  
  
  # Comparar las curvas mediante la prueba log-rank.
  
  logrank_test <- survdiff(
    Surv(
      OS_days,
      OS_event
    ) ~ Risk_group,
    data = risk_table
  )
  
  logrank_p_value <- pchisq(
    logrank_test$chisq,
    df = 1,
    lower.tail = FALSE
  )
  
  
  # Evaluar el puntaje continuo mediante un modelo de Cox.
  
  risk_score_cox <- coxph(
    Surv(
      OS_days,
      OS_event
    ) ~ Standardized_score,
    data = risk_table,
    ties = "efron"
  )
  
  cox_summary <- summary(
    risk_score_cox
  )
  
  
  # Calcular la capacidad discriminativa aparente.
  
  concordance_result <- concordance(
    Surv(
      OS_days,
      OS_event
    ) ~ Standardized_score,
    data = risk_table,
    reverse = TRUE
  )
  
  
  # Crear una fila resumen de la modalidad.
  
  summary_table <- data.frame(
    Modality = modality,
    Patients = nrow(risk_table),
    Events = sum(risk_table$OS_event == 1),
    Censored = sum(risk_table$OS_event == 0),
    Low_risk_patients =
      sum(risk_table$Risk_group == "Bajo riesgo"),
    High_risk_patients =
      sum(risk_table$Risk_group == "Alto riesgo"),
    Logrank_p_value = logrank_p_value,
    Hazard_ratio =
      cox_summary$coefficients[1, "exp(coef)"],
    CI95_lower =
      cox_summary$conf.int[1, "lower .95"],
    CI95_upper =
      cox_summary$conf.int[1, "upper .95"],
    Cox_p_value =
      cox_summary$coefficients[1, "Pr(>|z|)"],
    Apparent_C_index =
      as.numeric(concordance_result$concordance),
    Apparent_C_index_SE =
      sqrt(concordance_result$var),
    Evaluation_status = "Evaluación completada",
    stringsAsFactors = FALSE
  )
  
  
  return(
    list(
      Risk_table = risk_table,
      Kaplan_Meier_model = kaplan_meier_model,
      Summary = summary_table
    )
  )
}


# 4. Evaluar las cuatro modalidades

# Se aplica exactamente el mismo procedimiento a cada modelo.

model_names <- c(
  "rna",
  "methylation",
  "proteomics",
  "phosphoproteomics"
)

modality_names <- c(
  "RNA-seq",
  "Methylation",
  "Proteomics",
  "Phosphoproteomics"
)

file_names <- c(
  "RNAseq",
  "Methylation",
  "Proteomics",
  "Phosphoproteomics"
)

apparent_evaluations <- list()

for (i in seq_along(model_names)) {
  
  apparent_evaluations[[model_names[i]]] <-
    evaluate_apparent_model(
      model_results = omic_models[[model_names[i]]],
      modality = modality_names[i]
    )
}


# 5. Tabla resumen de la evaluación aparente

# Esta tabla compara el rendimiento observado en la misma cohorte
# utilizada para construir los modelos.

apparent_evaluation_summary <- rbind(
  
  apparent_evaluations$rna$Summary,
  
  apparent_evaluations$methylation$Summary,
  
  apparent_evaluations$proteomics$Summary,
  
  apparent_evaluations$phosphoproteomics$Summary
  
)

rownames(apparent_evaluation_summary) <- NULL

print(apparent_evaluation_summary)

# 6. Crear carpetas de resultados

# Se crean carpetas separadas para organizar las figuras,
# las tablas de riesgo y los resultados generales.

figure_directory <- file.path(
  "results",
  "figures",
  "Kaplan_Meier"
)

risk_directory <- file.path(
  "results",
  "tables",
  "Risk_scores"
)

summary_directory <- file.path(
  "results",
  "tables",
  "Summary"
)

dir.create(
  figure_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  risk_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  summary_directory,
  recursive = TRUE,
  showWarnings = FALSE
)


# 7. Guardar las curvas de Kaplan-Meier

# Se genera una figura para cada modalidad que produjo
# una firma LASSO evaluable.

for (i in seq_along(model_names)) {
  
  evaluation <- apparent_evaluations[[model_names[i]]]
  
  if (!is.null(evaluation$Kaplan_Meier_model)) {
    
    figure_file <- file.path(
      figure_directory,
      paste0(
        file_names[i],
        "_Kaplan_Meier.png"
      )
    )
    
    png(
      filename = figure_file,
      width = 1800,
      height = 1400,
      res = 200
    )
    
    plot(
      evaluation$Kaplan_Meier_model,
      lty = c(1, 2),
      xlab = "Tiempo de seguimiento (días)",
      ylab = "Probabilidad de supervivencia global",
      main = paste(
        "Kaplan-Meier -",
        modality_names[i]
      ),
      mark.time = TRUE
    )
    
    legend(
      "bottomleft",
      legend = c(
        "Bajo riesgo",
        "Alto riesgo"
      ),
      lty = c(1, 2),
      bty = "n"
    )
    
    dev.off()
  }
}


# 8. Guardar las tablas de riesgo

# Se guarda una tabla por modalidad con el puntaje pronóstico
# y el grupo de riesgo asignado a cada paciente.

for (i in seq_along(model_names)) {
  
  risk_table <-
    apparent_evaluations[[model_names[i]]]$Risk_table
  
  if (!is.null(risk_table)) {
    
    risk_file <- file.path(
      risk_directory,
      paste0(
        file_names[i],
        "_risk_scores.csv"
      )
    )
    
    write.csv(
      risk_table,
      risk_file,
      row.names = FALSE
    )
  }
}


# 9. Guardar el resumen de la evaluación

# La tabla comparativa se guarda en una carpeta específica
# para los resultados generales del análisis.

summary_file <- file.path(
  summary_directory,
  "apparent_evaluation_summary.csv"
)

write.csv(
  apparent_evaluation_summary,
  summary_file,
  row.names = FALSE
)


# 10. Guardar el objeto completo

# Los resultados completos se conservan en formato RDS porque
# serán utilizados en los scripts posteriores.

apparent_evaluations$summary <-
  apparent_evaluation_summary

saveRDS(
  apparent_evaluations,
  file.path(
    "results",
    "LUAD_apparent_evaluations.rds"
  )
)
