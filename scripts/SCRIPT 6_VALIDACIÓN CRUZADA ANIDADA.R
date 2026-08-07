###############################################################
# TFM LUAD
# SCRIPT 6: VALIDACIÓN CRUZADA ANIDADA
#
# Autora: Milena Vera García
#
# Objetivo:
# Evaluar el rendimiento fuera de muestra de los modelos
# pronósticos de RNA-seq, metilación, proteómica y
# fosfoproteómica mediante validación cruzada anidada.
#
# Procedimiento:
# 1. Dividir la cohorte en cinco pliegues externos.
# 2. Realizar la selección de variables solo en entrenamiento.
# 3. Ajustar LASSO-Cox mediante validación interna.
# 4. Predecir el riesgo en el pliegue externo.
# 5. Calcular el C-index fuera de muestra.
# 6. Repetir el procedimiento con tres semillas.
###############################################################


# 1. DIRECTORIO Y PAQUETES

project_directory <- "C:/Users/milen/OneDrive/Desktop/tesis/TFM_LUAD"

setwd(project_directory)

library(survival)
library(glmnet)


# 2. CARGAR LOS DATOS PREPROCESADOS

processed_data <- readRDS(
  "results/LUAD_processed_data.rds"
)

survival_complete <- processed_data$survival

rna_complete <- processed_data$rna

methylation_complete <- processed_data$methylation

proteomics_complete <- processed_data$proteomics

phosphoproteomics_complete <-
  processed_data$phosphoproteomics


# 3. COMPROBAR QUE LAS FUNCIONES DE MODELADO ESTÉN DISPONIBLES

if (!exists("fit_omic_model")) {
  
  stop(
    paste(
      "No se encontró la función fit_omic_model().",
      "Ejecuta primero SCRIPT4_MODELADO.R",
      "en esta misma sesión de R."
    )
  )
}


# 4. CREAR PLIEGUES ESTRATIFICADOS EXTERNOS

create_outer_folds <- function(
    event,
    number_folds = 5,
    seed = 2026
) {
  
  set.seed(seed)
  
  fold_id <- integer(
    length(event)
  )
  
  for (event_value in sort(unique(event))) {
    
    event_positions <- which(
      event == event_value
    )
    
    fold_numbers <- rep(
      seq_len(number_folds),
      length.out = length(event_positions)
    )
    
    fold_id[event_positions] <- sample(
      fold_numbers
    )
  }
  
  return(fold_id)
}


# 5. VALIDACIÓN CRUZADA ANIDADA

evaluate_fold <- function(
    training_matrix,
    test_matrix,
    training_survival,
    test_survival,
    omic_name,
    seed
) {
  
  fold_model <- tryCatch(
    fit_omic_model(
      omic_matrix = training_matrix,
      survival_data = training_survival,
      omic_name = omic_name,
      number_top_features = 200,
      number_inner_folds = 5,
      seed = seed,
      show_messages = FALSE
    ),
    error = function(e) NULL
  )
  
  if (is.null(fold_model)) {
    return(
      list(
        C_index = NA_real_,
        Selected_features = character(0),
        Significant_features = NA_integer_,
        Features_entered_LASSO = NA_integer_,
        Status = "No se pudo ajustar el modelo"
      )
    )
  }
  
  selected_features <-
    fold_model$Selected_features_lambda_min$Feature_ID
  
  if (length(selected_features) == 0) {
    return(
      list(
        C_index = NA_real_,
        Selected_features = character(0),
        Significant_features =
          nrow(fold_model$Significant_results),
        Features_entered_LASSO =
          nrow(fold_model$Top_features),
        Status = "LASSO no seleccionó variables"
      )
    )
  }
  
  test_matrix_selected <- test_matrix[
    ,
    fold_model$Selected_feature_ids,
    drop = FALSE
  ]
  
  training_prediction <- as.numeric(
    predict(
      fold_model$LASSO_model,
      newx = as.matrix(
        fold_model$Selected_matrix
      ),
      s = "lambda.min",
      type = "link"
    )
  )
  
  test_prediction <- as.numeric(
    predict(
      fold_model$LASSO_model,
      newx = as.matrix(
        test_matrix_selected
      ),
      s = "lambda.min",
      type = "link"
    )
  )
  
  training_sd <- sd(
    training_prediction
  )
  
  if (!is.finite(training_sd) || training_sd == 0) {
    return(
      list(
        C_index = NA_real_,
        Selected_features = selected_features,
        Significant_features =
          nrow(fold_model$Significant_results),
        Features_entered_LASSO =
          nrow(fold_model$Top_features),
        Status = "Puntaje sin variación"
      )
    )
  }
  
  standardized_test_prediction <- (
    test_prediction -
      mean(training_prediction)
  ) / training_sd
  
  c_index <- concordance(
    Surv(
      test_survival$OS_days,
      test_survival$OS_event
    ) ~ standardized_test_prediction,
    reverse = TRUE
  )$concordance
  
  return(
    list(
      C_index = as.numeric(c_index),
      Selected_features = selected_features,
      Significant_features =
        nrow(fold_model$Significant_results),
      Features_entered_LASSO =
        nrow(fold_model$Top_features),
      Status = "Pliegue evaluado"
    )
  )
}


run_nested_validation <- function(
    omic_matrix,
    survival_data,
    omic_name,
    validation_seeds = c(
      2026,
      2030,
      2040
    ),
    number_outer_folds = 5
) {
  
  all_fold_results <- list()
  
  result_position <- 1
  
  for (current_seed in validation_seeds) {
    
    outer_fold_id <- create_outer_folds(
      event = survival_data$OS_event,
      number_folds = number_outer_folds,
      seed = current_seed
    )
    
    for (outer_fold in seq_len(number_outer_folds)) {
      
      message(
        omic_name,
        " - semilla ",
        current_seed,
        " - pliegue ",
        outer_fold
      )
      
      test_positions <- which(
        outer_fold_id == outer_fold
      )
      
      training_positions <- which(
        outer_fold_id != outer_fold
      )
      
      fold_result <- evaluate_fold(
        training_matrix =
          omic_matrix[
            training_positions,
            ,
            drop = FALSE
          ],
        
        test_matrix =
          omic_matrix[
            test_positions,
            ,
            drop = FALSE
          ],
        
        training_survival =
          survival_data[
            training_positions,
            ,
            drop = FALSE
          ],
        
        test_survival =
          survival_data[
            test_positions,
            ,
            drop = FALSE
          ],
        
        omic_name =
          paste0(
            omic_name,
            "_seed_",
            current_seed,
            "_fold_",
            outer_fold
          ),
        
        seed =
          current_seed + outer_fold
      )
      
      all_fold_results[[result_position]] <- data.frame(
        Modality = omic_name,
        Seed = current_seed,
        Fold = outer_fold,
        Training_patients =
          length(training_positions),
        Test_patients =
          length(test_positions),
        Training_events =
          sum(
            survival_data$OS_event[
              training_positions
            ] == 1
          ),
        Test_events =
          sum(
            survival_data$OS_event[
              test_positions
            ] == 1
          ),
        Significant_features =
          fold_result$Significant_features,
        Features_entered_LASSO =
          fold_result$Features_entered_LASSO,
        Selected_features =
          length(
            fold_result$Selected_features
          ),
        Test_C_index =
          fold_result$C_index,
        Fold_status =
          fold_result$Status,
        stringsAsFactors = FALSE
      )
      
      all_fold_results[[result_position]]$
        Selected_feature_IDs <- list(
          fold_result$Selected_features
        )
      
      result_position <- result_position + 1
    }
  }
  
  fold_summary <- do.call(
    rbind,
    all_fold_results
  )
  
  c_index_by_seed <- aggregate(
    Test_C_index ~ Seed,
    data = fold_summary,
    FUN = function(x) {
      mean(
        x,
        na.rm = TRUE
      )
    }
  )
  
  names(c_index_by_seed)[2] <-
    "Mean_C_index"
  
  valid_c_indices <- fold_summary$Test_C_index[
    is.finite(
      fold_summary$Test_C_index
    )
  ]
  
  general_summary <- data.frame(
    Modality = omic_name,
    Repetitions =
      length(validation_seeds),
    Outer_folds =
      number_outer_folds,
    Mean_validated_C_index =
      mean(valid_c_indices),
    SD_validated_C_index =
      sd(valid_c_indices),
    Minimum_validated_C_index =
      min(valid_c_indices),
    Maximum_validated_C_index =
      max(valid_c_indices),
    Valid_folds =
      length(valid_c_indices),
    Total_folds =
      nrow(fold_summary),
    stringsAsFactors = FALSE
  )
  
  selected_features <- unlist(
    fold_summary$Selected_feature_IDs,
    use.names = FALSE
  )
  
  feature_frequency <- as.data.frame(
    table(selected_features),
    stringsAsFactors = FALSE
  )
  
  names(feature_frequency) <- c(
    "Feature_ID",
    "Times_selected"
  )
  
  feature_frequency$Selection_percentage <- round(
    feature_frequency$Times_selected /
      nrow(fold_summary) * 100,
    2
  )
  
  feature_frequency <- feature_frequency[
    order(
      feature_frequency$Times_selected,
      decreasing = TRUE
    ),
    ,
    drop = FALSE
  ]
  
  return(
    list(
      Fold_summary =
        fold_summary,
      C_index_by_seed =
        c_index_by_seed,
      General_summary =
        general_summary,
      Feature_frequency =
        feature_frequency
    )
  )
}

# 6. VALIDAR RNA-SEQ

rna_validation <- run_nested_validation(
  omic_matrix =
    rna_complete,
  
  survival_data =
    survival_complete,
  
  omic_name =
    "RNA-seq"
)


# 7. VALIDAR METILACIÓN

methylation_validation <- run_nested_validation(
  omic_matrix =
    methylation_complete,
  
  survival_data =
    survival_complete,
  
  omic_name =
    "Methylation"
)


# 8. VALIDAR PROTEÓMICA

proteomics_validation <- run_nested_validation(
  omic_matrix =
    proteomics_complete,
  
  survival_data =
    survival_complete,
  
  omic_name =
    "Proteomics"
)


# 9. VALIDAR FOSFOPROTEÓMICA

phosphoproteomics_validation <- run_nested_validation(
  omic_matrix =
    phosphoproteomics_complete,
  
  survival_data =
    survival_complete,
  
  omic_name =
    "Phosphoproteomics"
)


# 10. CREAR LA TABLA RESUMEN DE VALIDACIÓN

nested_validation_summary <- rbind(
  
  rna_validation$General_summary,
  
  methylation_validation$General_summary,
  
  proteomics_validation$General_summary,
  
  phosphoproteomics_validation$General_summary
)

rownames(
  nested_validation_summary
) <- NULL

print(
  nested_validation_summary
)


# 11. CREAR TABLA POR REPETICIÓN

nested_repetition_summary <- rbind(
  
  rna_validation$Repetition_summary,
  
  methylation_validation$Repetition_summary,
  
  proteomics_validation$Repetition_summary,
  
  phosphoproteomics_validation$Repetition_summary
)

rownames(
  nested_repetition_summary
) <- NULL

print(
  nested_repetition_summary
)


# 12. AGRUPAR LOS RESULTADOS

nested_validation_results <- list(
  
  rna =
    rna_validation,
  
  methylation =
    methylation_validation,
  
  proteomics =
    proteomics_validation,
  
  phosphoproteomics =
    phosphoproteomics_validation,
  
  summary =
    nested_validation_summary,
  
  repetition_summary =
    nested_repetition_summary
)


# 13. GUARDAR LOS RESULTADOS GENERALES

saveRDS(
  nested_validation_results,
  file = "results/LUAD_nested_validation.rds"
)

write.csv(
  nested_validation_summary,
  file = "results/nested_validation_summary.csv",
  row.names = FALSE
)

write.csv(
  nested_repetition_summary,
  file = "results/nested_validation_by_seed.csv",
  row.names = FALSE
)


# 14. GUARDAR FRECUENCIAS DE SELECCIÓN

write.csv(
  rna_validation$Feature_frequency,
  file = "results/RNAseq_feature_frequency.csv",
  row.names = FALSE
)

write.csv(
  methylation_validation$Feature_frequency,
  file = "results/Methylation_feature_frequency.csv",
  row.names = FALSE
)

write.csv(
  proteomics_validation$Feature_frequency,
  file = "results/Proteomics_feature_frequency.csv",
  row.names = FALSE
)

write.csv(
  phosphoproteomics_validation$Feature_frequency,
  file = "results/Phosphoproteomics_feature_frequency.csv",
  row.names = FALSE
)


# 15. GUARDAR INFORMACIÓN DE LA SESIÓN

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = "results/sessionInfo_Script06.txt"
)

