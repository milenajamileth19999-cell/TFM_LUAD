# TFM LUAD
# Construcción preliminar de cohorte multi-ómica
# Autor: Milena Vera García

# Este script realiza la primera fase del proyecto:
# cargar los datos clínicos y ómicos de LUAD, armonizar los
# identificadores de pacientes y construir una cohorte común
# disponible simultáneamente para supervivencia, RNA-seq,
# metilación, proteómica y fosfoproteómica.

# 1. Definición de directorio de trabajo
# Cambiar esta ruta según la ubicación local del proyecto.
setwd("C:/Users/milen/OneDrive/Desktop/tesis/TFM_LUAD")


# 2. Carga de datos clínicos

# El archivo survival contiene la información clínica necesaria
# para definir el desenlace pronóstico: OS_days y OS_event.

survival <- read.delim(
  "data_raw/LUAD_survival.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# La lista de casos tumorales permite verificar qué muestras
# corresponden a tumores dentro de la cohorte LUAD.

tumor_cases <- read.delim(
  "data_raw/LUAD_Tumor_CaseList.txt",
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

# 3. Carga de matrices ómicas

# En estas matrices, la primera columna corresponde al identificador
# molecular de cada feature, por ejemplo genes, sitios CpG, proteínas
# o fosfositios. Las columnas restantes corresponden a pacientes.

rna <- read.delim(
  "data_raw/LUAD_RNAseq_gene_RSEM_coding_UQ_1500_log2_Tumor.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = TRUE
)

methylation <- read.delim(
  "data_raw/LUAD_Methylation_Tumor.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = TRUE
)

proteomics <- read.delim(
  "data_raw/LUAD_Proteomics_Tumor.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = TRUE
)

phosphoproteomics <- read.delim(
  "data_raw/LUAD_Phosphoproteomics_Tumor.txt",
  header = TRUE,
  sep = "\t",
  stringsAsFactors = FALSE,
  check.names = TRUE
)

# 4. Revición dd dimensiones iniciales

# Esta revisión permite conocer el número inicial de pacientes y
# variables por modalidad. Es una etapa importante porque los datos
# ómicos suelen presentar muchas más variables que pacientes.

dim(survival)
dim(tumor_cases)
dim(rna)
dim(methylation)
dim(proteomics)
dim(phosphoproteomics)

# 5. Armonización de identificadores de pacientes

# R modifica automáticamente los nombres de columnas que empiezan
# con números, por ejemplo "11LU013" puede convertirse en "X11LU013".
# Para evitar pérdidas artificiales de pacientes durante la integración,
# se aplica make.names() a los identificadores clínicos y tumorales.

survival$case_id_R <- make.names(survival$case_id)
tumor_cases$case_id_R <- make.names(tumor_cases$V1)

# 6. Extraer identificadores de pacientes desde matrices ómicas

# La primera columna de cada matriz no representa un paciente,
# sino una variable molecular. Por eso se excluye al extraer
# los identificadores de pacientes.

rna_ids <- names(rna)[-1]
methylation_ids <- names(methylation)[-1]
proteomics_ids <- names(proteomics)[-1]
phosphoproteomics_ids <- names(phosphoproteomics)[-1]

# 7. Construir cohorte integrada

# Se conservan únicamente los pacientes presentes simultáneamente
# en supervivencia, lista de casos tumorales y las cuatro modalidades
# ómicas. Esto garantiza que la comparación entre RNA-seq, metilación,
# proteómica y fosfoproteómica se realice sobre la misma cohorte.

common_all <- Reduce(
  intersect,
  list(
    survival$case_id_R,
    tumor_cases$case_id_R,
    rna_ids,
    methylation_ids,
    proteomics_ids,
    phosphoproteomics_ids
  )
)

# 8. Creación de resumen de disponibilidad

# El resumen documenta cuántos pacientes y variables están disponibles
# en cada fuente de información antes de los análisis posteriores.

cohort_summary <- data.frame(
  Fuente = c(
    "Survival",
    "Tumor case list",
    "RNA-seq",
    "Metilacion",
    "Proteomica",
    "Fosfoproteomica",
    "Cohorte integrada final"
  ),
  N_pacientes = c(
    nrow(survival),
    nrow(tumor_cases),
    length(rna_ids),
    length(methylation_ids),
    length(proteomics_ids),
    length(phosphoproteomics_ids),
    length(common_all)
  ),
  N_variables = c(
    5,
    1,
    nrow(rna),
    nrow(methylation),
    nrow(proteomics),
    nrow(phosphoproteomics),
    NA
  )
)

print(cohort_summary)

# 9. Filtrar datos clínicos para la cohorte integrada

# Se conserva la información de supervivencia únicamente para los
# pacientes presentes en todas las modalidades. Esta tabla clínica
# será la base para los modelos de supervivencia posteriores.

survival_common <- survival[
  survival$case_id_R %in% common_all,
]

# 10. Verificación de datos faltantes del desenlace clínico

# OS_days representa el tiempo de seguimiento en días.
# OS_event indica si ocurrió el evento de interés o si el dato
# corresponde a una observación censurada.

missing_OS_days <- sum(is.na(survival_common$OS_days))
missing_OS_event <- sum(is.na(survival_common$OS_event))

missing_OS_days
missing_OS_event

# 11. resumen de la cohorte

# El archivo cohort_summary.csv permite documentar de forma reproducible
# las dimensiones iniciales y el tamaño de la cohorte multi-ómica final.

if (!dir.exists("results")) {
  dir.create("results")
}

write.csv(
  cohort_summary,
  "results/cohort_summary.csv",
  row.names = FALSE
)
