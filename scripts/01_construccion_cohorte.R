###############################################################
# TFM LUAD
# SCRIPT 1: CARGA DE DATOS Y CONSTRUCCIÓN DE LA COHORTE COMÚN
#
# Objetivo:
# Cargar los datos clínicos y multi-ómicos del proyecto LUAD,
# comprobar sus dimensiones e identificar los pacientes
# presentes simultáneamente en las cuatro modalidades ómicas.
###############################################################


# 1. Directorio del proyecto

# Se define la carpeta principal para localizar los archivos
# originales y guardar los resultados generados.

project_directory <- "C:/Users/milen/OneDrive/Desktop/tesis/TFM_LUAD"

setwd(project_directory)


# 2. Carga de los datos clínicos

# Se cargan los archivos de supervivencia, fenotipo y casos
# tumorales disponibles para la cohorte LUAD.

survival <- read.delim(
  "data_raw/LUAD_survival.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

phenotype <- read.delim(
  "data_raw/LUAD_phenotype.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

tumor_cases <- read.delim(
  "data_raw/LUAD_Tumor_CaseList.txt",
  header = FALSE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE
)

colnames(tumor_cases) <- "case_id"


# 3. Carga de las matrices ómicas

# Se importan las matrices de RNA-seq, metilación, proteómica
# y fosfoproteómica. La primera columna contiene el identificador
# molecular y las columnas restantes corresponden a pacientes.

rnaseq <- read.delim(
  "data_raw/LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

methylation <- read.delim(
  "data_raw/LUAD_methylation_gene_beta_value_Tumor.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

proteomics <- read.delim(
  "data_raw/LUAD_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)

phosphoproteomics <- read.delim(
  "data_raw/LUAD_phospho_site_abundance_log2_reference_intensity_normalized_Tumor.txt",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("", "NA")
)


# 4. Resumen de los datos cargados

# Se construye una tabla con el número de filas y columnas de
# cada archivo para comprobar sus dimensiones iniciales.

data_summary <- data.frame(
  Dataset = c(
    "Supervivencia",
    "Fenotipo",
    "Casos tumorales",
    "RNA-seq",
    "Metilación",
    "Proteómica",
    "Fosfoproteómica"
  ),
  Filas = c(
    nrow(survival),
    nrow(phenotype),
    nrow(tumor_cases),
    nrow(rnaseq),
    nrow(methylation),
    nrow(proteomics),
    nrow(phosphoproteomics)
  ),
  Columnas = c(
    ncol(survival),
    ncol(phenotype),
    ncol(tumor_cases),
    ncol(rnaseq),
    ncol(methylation),
    ncol(proteomics),
    ncol(phosphoproteomics)
  ),
  stringsAsFactors = FALSE
)

print(data_summary)


# 5. Extraer los identificadores de pacientes

# Los identificadores clínicos se obtienen de la columna case_id.
# En las matrices ómicas se excluye la primera columna porque
# corresponde al identificador molecular y no a un paciente.

survival_ids <- survival$case_id

tumor_case_ids <- tumor_cases$case_id

rnaseq_ids <- colnames(rnaseq)[-1]

methylation_ids <- colnames(methylation)[-1]

proteomics_ids <- colnames(proteomics)[-1]

phosphoproteomics_ids <- colnames(phosphoproteomics)[-1]


# 6. Resumen de pacientes disponibles

# Esta tabla muestra el número de pacientes disponibles en cada
# fuente antes de realizar la integración multi-ómica.

cohort_availability <- data.frame(
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

print(cohort_availability)


# 7. Identificar los pacientes comunes

# Se conservan únicamente los pacientes presentes simultáneamente
# en supervivencia, casos tumorales y las cuatro modalidades
# ómicas, garantizando una cohorte comparable.

common_patient_ids <- Reduce(
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

# Se mantiene el orden original del archivo de supervivencia.

common_patient_ids <- survival_ids[
  survival_ids %in% common_patient_ids
]


# 8. Preparar los datos de supervivencia

# Se extraen los datos clínicos de los pacientes presentes en
# todas las modalidades y se conserva el mismo orden de la cohorte.

survival_common <- survival[
  match(
    common_patient_ids,
    survival$case_id
  ),
  ,
  drop = FALSE
]

rownames(survival_common) <- survival_common$case_id


# 9. Resumen de la cohorte común

# Esta tabla documenta el número inicial de pacientes, el tamaño
# final de la cohorte integrada y las exclusiones producidas
# durante la intersección de las modalidades.

cohort_result <- data.frame(
  Indicador = c(
    "Pacientes iniciales con supervivencia",
    "Pacientes comunes a todas las modalidades",
    "Pacientes excluidos durante la integración"
  ),
  Resultado = c(
    length(unique(survival_ids)),
    length(common_patient_ids),
    length(unique(survival_ids)) -
      length(common_patient_ids)
  ),
  stringsAsFactors = FALSE
)

print(cohort_result)


# 10. Guardar los resultados

# Se guardan las tablas resumen y los dos objetos necesarios para
# continuar con el preprocesamiento de las matrices ómicas.

write.csv(
  data_summary,
  "results/data_summary.csv",
  row.names = FALSE
)

write.csv(
  cohort_availability,
  "results/cohort_available_patients.csv",
  row.names = FALSE
)

write.csv(
  cohort_result,
  "results/cohort_result.csv",
  row.names = FALSE
)

saveRDS(
  common_patient_ids,
  "results/common_patient_ids.rds"
)

saveRDS(
  survival_common,
  "results/survival_common.rds"
)
