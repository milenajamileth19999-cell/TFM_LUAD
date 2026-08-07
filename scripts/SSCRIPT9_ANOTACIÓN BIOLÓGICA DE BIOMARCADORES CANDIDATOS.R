###############################################################
# TFM LUAD
# SCRIPT 9: ANOTACIÓN BIOLÓGICA DE BIOMARCADORES CANDIDATOS
#
# Autora: Milena Vera García
#
# Objetivo:
# Anotar los biomarcadores candidatos priorizados en el Script 8
# y organizar la información necesaria para su interpretación.
###############################################################


# 1. Directorio y paquetes

project_directory <- "C:/Users/milen/OneDrive/Desktop/tesis/TFM_LUAD"

setwd(project_directory)

library(AnnotationDbi)
library(org.Hs.eg.db)


# 2. Cargar candidatos priorizados

priority_biomarkers <- read.csv(
  "results/tables/Biomarkers/priority_biomarkers.csv",
  stringsAsFactors = FALSE,
  check.names = FALSE
)


# 3. Preparar identificadores moleculares

# En RNA-seq, metilación y proteómica se elimina la versión
# del identificador Ensembl.

priority_biomarkers$Ensembl_gene <- sub(
  "\\..*$",
  "",
  priority_biomarkers$Feature_ID
)

priority_biomarkers$Ensembl_protein <- NA_character_
priority_biomarkers$Phosphosite <- NA_character_
priority_biomarkers$Peptide <- NA_character_


# En fosfoproteómica se separan gen, proteína, sitio y péptido.

phospho_positions <- which(
  priority_biomarkers$Modality == "Phosphoproteomics"
)

for (i in phospho_positions) {
  
  parts <- strsplit(
    priority_biomarkers$Feature_ID[i],
    "\\|"
  )[[1]]
  
  priority_biomarkers$Ensembl_gene[i] <- sub(
    "\\..*$",
    "",
    parts[1]
  )
  
  if (length(parts) >= 2) {
    priority_biomarkers$Ensembl_protein[i] <- sub(
      "\\..*$",
      "",
      parts[2]
    )
  }
  
  if (length(parts) >= 3) {
    priority_biomarkers$Phosphosite[i] <- parts[3]
  }
  
  if (length(parts) >= 4) {
    priority_biomarkers$Peptide[i] <- parts[4]
  }
}


# 4. Anotar genes

# Se convierten los identificadores Ensembl a símbolo génico
# y nombre oficial.

gene_annotation <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = unique(priority_biomarkers$Ensembl_gene),
  keytype = "ENSEMBL",
  columns = c(
    "ENSEMBL",
    "SYMBOL",
    "GENENAME"
  )
)

gene_annotation <- gene_annotation[
  !duplicated(gene_annotation$ENSEMBL),
]

priority_biomarkers$Gene_symbol <-
  gene_annotation$SYMBOL[
    match(
      priority_biomarkers$Ensembl_gene,
      gene_annotation$ENSEMBL
    )
  ]

priority_biomarkers$Gene_name <-
  gene_annotation$GENENAME[
    match(
      priority_biomarkers$Ensembl_gene,
      gene_annotation$ENSEMBL
    )
  ]


# 5. Definir tipo molecular

priority_biomarkers$Molecular_type <- ifelse(
  priority_biomarkers$Modality == "RNA-seq",
  "Expresión génica",
  ifelse(
    priority_biomarkers$Modality == "Methylation",
    "Metilación",
    ifelse(
      priority_biomarkers$Modality == "Proteomics",
      "Abundancia proteica",
      "Sitio de fosforilación"
    )
  )
)


# 6. Clasificar estabilidad de selección

# La prioridad refleja la frecuencia de selección en los
# 15 pliegues externos.

priority_biomarkers$Selection_priority <- cut(
  priority_biomarkers$Selection_percentage,
  breaks = c(
    -Inf,
    20,
    30,
    50,
    Inf
  ),
  labels = c(
    "Baja",
    "Moderada",
    "Alta",
    "Muy alta"
  ),
  right = FALSE
)


# 7. Ordenar y crear ranking

modality_order <- c(
  "RNA-seq",
  "Methylation",
  "Proteomics",
  "Phosphoproteomics"
)

priority_biomarkers$Modality <- factor(
  priority_biomarkers$Modality,
  levels = modality_order
)

priority_biomarkers <- priority_biomarkers[
  order(
    priority_biomarkers$Modality,
    -priority_biomarkers$Selection_percentage,
    -abs(priority_biomarkers$LASSO_coefficient)
  ),
]

priority_biomarkers$Rank <- ave(
  seq_len(nrow(priority_biomarkers)),
  priority_biomarkers$Modality,
  FUN = seq_along
)


# 8. Crear tabla maestra

master_annotation <- data.frame(
  
  Rank = priority_biomarkers$Rank,
  
  Modality = priority_biomarkers$Modality,
  
  Molecular_type = priority_biomarkers$Molecular_type,
  
  Gene_symbol = priority_biomarkers$Gene_symbol,
  
  Gene_name = priority_biomarkers$Gene_name,
  
  Ensembl_gene = priority_biomarkers$Ensembl_gene,
  
  Ensembl_protein = priority_biomarkers$Ensembl_protein,
  
  Phosphosite = priority_biomarkers$Phosphosite,
  
  Peptide = priority_biomarkers$Peptide,
  
  LASSO_beta = round(
    priority_biomarkers$LASSO_coefficient,
    3
  ),
  
  Association = priority_biomarkers$Association,
  
  Hazard_ratio = round(
    priority_biomarkers$Univariate_HR,
    3
  ),
  
  P_value = signif(
    priority_biomarkers$Univariate_p_value,
    3
  ),
  
  Times_selected = priority_biomarkers$Times_selected,
  
  Selection_percentage = round(
    priority_biomarkers$Selection_percentage,
    1
  ),
  
  Selection_priority =
    priority_biomarkers$Selection_priority,
  
  stringsAsFactors = FALSE
)

rownames(master_annotation) <- NULL

print(master_annotation)


# 9. Resumen de anotación

annotation_summary <- data.frame(
  
  Modality = modality_order,
  
  Candidates = sapply(
    modality_order,
    function(x) {
      sum(master_annotation$Modality == x)
    }
  ),
  
  Annotated_genes = sapply(
    modality_order,
    function(x) {
      sum(
        master_annotation$Modality == x &
          !is.na(master_annotation$Gene_symbol)
      )
    }
  ),
  
  stringsAsFactors = FALSE
)

print(annotation_summary)


# 10. Guardar resultado

annotation_directory <- file.path(
  "results",
  "tables",
  "Biological_Annotation"
)

dir.create(
  annotation_directory,
  recursive = TRUE,
  showWarnings = FALSE
)

write.csv(
  master_annotation,
  file.path(
    annotation_directory,
    "Master_annotation.csv"
  ),
  row.names = FALSE
)



# 11. Seleccionar biomarcadores para la discusión

# Se conservan únicamente los biomarcadores con prioridad
# "Alta" o "Muy alta", equivalentes a una frecuencia de
# selección igual o superior al 30%.

discussion_biomarkers <- subset(
  master_annotation,
  Selection_priority == "Muy alta"
)

# Ordenar por modalidad y estabilidad

discussion_biomarkers <- discussion_biomarkers[
  order(
    discussion_biomarkers$Modality,
    -discussion_biomarkers$Selection_percentage,
    discussion_biomarkers$Rank
  ),
]

rownames(discussion_biomarkers) <- NULL

print(discussion_biomarkers)

# Guardar la tabla

write.csv(
  discussion_biomarkers,
  file.path(
    annotation_directory,
    "Discussion_biomarkers.csv"
  ),
  row.names = FALSE
)