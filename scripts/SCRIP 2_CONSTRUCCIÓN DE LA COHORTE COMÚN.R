###############################################################
# TFM LUAD
# SCRIPT 2: CONSTRUCCIÓN DE LA COHORTE COMÚN
#
# Autora: Milena Vera García
#
# Objetivo:
# Identificar los pacientes presentes simultáneamente en los
# datos de supervivencia, la lista de casos tumorales y las
# cuatro modalidades ómicas.
#
# Este script no filtra variables moleculares, no imputa valores
# faltantes y no ajusta modelos estadísticos.
###############################################################


#==============================================================
# 1. DIRECTORIO DEL PROYECTO
#==============================================================

project_directory <- "C:/Users/milen/OneDrive/Desktop/tesis/TFM_LUAD"

setwd(project_directory)

output_directory <- "results"


#==============================================================
# 2. CARGAR LOS DATOS GENERADOS POR EL SCRIPT 1
#==============================================================

raw_data <- readRDS(
  "results/LUAD_raw_data_complete.rds"
)

survival <- raw_data$survival
phenotype <- raw_data$phenotype
tumor_cases <- raw_data$tumor_cases

rnaseq <- raw_data$rnaseq
methylation <- raw_data$methylation
proteomics <- raw_data$proteomics
phosphoproteomics <- raw_data$phosphoproteomics


#==============================================================
# 3. EXTRAER LOS IDENTIFICADORES DE PACIENTES
#==============================================================

# En supervivencia y en la lista tumoral, los pacientes se
# encuentran en la columna case_id.

survival_ids <- as.character(
  survival$case_id
)

tumor_case_ids <- as.character(
  tumor_cases$case_id
)


# En las matrices ómicas, la primera columna es idx y las
# columnas restantes corresponden a los pacientes.

rnaseq_ids <- colnames(rnaseq)[
  colnames(rnaseq) != "idx"
]

methylation_ids <- colnames(methylation)[
  colnames(methylation) != "idx"
]

proteomics_ids <- colnames(proteomics)[
  colnames(proteomics) != "idx"
]

phosphoproteomics_ids <- colnames(phosphoproteomics)[
  colnames(phosphoproteomics) != "idx"
]


#==============================================================
# 4. RESUMIR LOS PACIENTES DISPONIBLES
#==============================================================

cohort_summary <- data.frame(
  
  Fuente = c(
    "Supervivencia",
    "Casos tumorales",
    "RNA-seq",
    "Metilación",
    "Proteómica",
    "Fosfoproteómica"
  ),
  
  Pacientes_disponibles = c(
    length(unique(survival_ids)),
    length(unique(tumor_case_ids)),
    length(unique(rnaseq_ids)),
    length(unique(methylation_ids)),
    length(unique(proteomics_ids)),
    length(unique(phosphoproteomics_ids))
  ),
  
  stringsAsFactors = FALSE
)

print(cohort_summary)


#==============================================================
# 5. IDENTIFICAR LOS PACIENTES COMUNES
#==============================================================

common_all <- Reduce(
  intersect,
  list(
    survival_ids,
    tumor_case_ids,
    rnaseq_ids,
    methylation_ids,
    proteomics_ids,
    phosphoproteomics_ids
  )
)


# Mantener el mismo orden que aparece en supervivencia.

common_all <- survival_ids[
  survival_ids %in% common_all
]

length(common_all)

head(common_all)


#==============================================================
# 6. PREPARAR LA SUPERVIVENCIA DE LA COHORTE COMÚN
#==============================================================

survival_common <- survival[
  match(
    common_all,
    survival$case_id
  ),
  ,
  drop = FALSE
]

rownames(survival_common) <- survival_common$case_id


#==============================================================
# 7. PREPARAR LA INFORMACIÓN FENOTÍPICA
#==============================================================

phenotype_common <- phenotype[
  match(
    common_all,
    phenotype$case_id
  ),
  ,
  drop = FALSE
]

rownames(phenotype_common) <- phenotype_common$case_id


#==============================================================
# 8. COMPROBAR LA COHORTE RESULTANTE
#==============================================================

cohort_result <- data.frame(
  
  Indicador = c(
    "Pacientes en supervivencia",
    "Pacientes comunes a todas las modalidades",
    "Pacientes perdidos por falta de alguna modalidad"
  ),
  
  Resultado = c(
    length(unique(survival_ids)),
    length(common_all),
    length(unique(survival_ids)) - length(common_all)
  ),
  
  stringsAsFactors = FALSE
)

print(cohort_result)


# Revisar que el orden de los pacientes sea correcto.

head(
  rownames(survival_common)
)

head(
  common_all
)

identical(
  rownames(survival_common),
  common_all
)


#==============================================================
# 9. AGRUPAR LOS RESULTADOS DE LA COHORTE
#==============================================================

common_cohort <- list(
  
  patient_ids = common_all,
  
  survival = survival_common,
  
  phenotype = phenotype_common,
  
  cohort_summary = cohort_summary,
  
  cohort_result = cohort_result
  
)


#==============================================================
# 10. GUARDAR LA COHORTE COMÚN
#==============================================================

saveRDS(
  common_cohort,
  file = file.path(
    output_directory,
    "LUAD_common_cohort.rds"
  )
)

write.csv(
  cohort_summary,
  file = file.path(
    output_directory,
    "cohort_available_patients.csv"
  ),
  row.names = FALSE
)

write.csv(
  cohort_result,
  file = file.path(
    output_directory,
    "cohort_result.csv"
  ),
  row.names = FALSE
)

write.table(
  common_all,
  file = file.path(
    output_directory,
    "common_patient_ids.txt"
  ),
  row.names = FALSE,
  col.names = FALSE,
  quote = FALSE
)


#==============================================================
# 11. INFORMACIÓN DE LA SESIÓN
#==============================================================

writeLines(
  capture.output(
    sessionInfo()
  ),
  con = file.path(
    output_directory,
    "sessionInfo_Script02.txt"
  )
)

