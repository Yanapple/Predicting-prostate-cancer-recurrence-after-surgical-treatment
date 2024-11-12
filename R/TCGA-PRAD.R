#TCGA-PRAD

# packages
#install.packages('BiocManager')
library(BiocManager)

library(TCGAbiolinks)
install.packages('magrittr')
library(magrittr)
#install('org.Hs.eg.db')
library("org.Hs.eg.db", character.only = TRUE)
library(SummarizedExperiment)
#install.packages("xml2")
library(xml2)

# UTF-8

# downloading expression var 1 (with filtering)
prad_query <- GDCquery(project = "TCGA-PRAD",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification"
)
prad_query$filters <- list(c("gene_type","in",c("miRNA","iRNA","lncRNA"))) 
GDCdownload(prad_query)

# check count of sample types
as.vector(as.data.frame(prad_query[[1]][[1]])$sample_type) %>%
  table()

# data preparing
prad_data <- GDCprepare(prad_query)
prad_data@rowRanges@elementMetadata@listData$gene_type
prad_expMat <- assay(prad_data)

# save the matrix
write.csv(prad_expMat, 'C:/.../prad_3RNAtype_expMat.csv')

# downloading expression var 2 (with three specific datasets: iRNA, miRNA, lncRNA)
# downloading miRNA expression
prad_miRNA_query <- GDCquery(project = "TCGA-PRAD",
                       data.category = "Transcriptome Profiling",
                       data.type = "miRNA Expression Quantification"
)
GDCdownload(prad_miRNA_query)

# data preparing
prad_miRNA_data <- GDCprepare(prad_miRNA_query)
prad_miRNA_expMat <- assay(prad_miRNA_data)

#но как загрузить iRNA, lncRNA? Как они названы в data.type?
# Attantion!!!
# luad_expMat has not normalized data!

# downloading full clinical data
prad_clinical_query <- GDCquery(project = "TCGA-PRAD",
                                data.category = "Clinical",
                                data.type = "Clinical Supplement"
)

GDCdownload(prad_clinical_query)

# data preparing
prad_clinical_data <- GDCprepare(prad_clinical_query)
prad_clinical_expMat <- assay(prad_clinical_data)

#!! prad_clinical_data has no column "bcr_patient_barcode, vital_status, disease_free_status, recurrence"
clinical_relapse_data <- prad_clinical_data %>%
  select(bcr_patient_barcode, vital_status, disease_free_status, recurrence)

# Assume list and the tibbles are in first column
tibbles <- lapply(prad_clinical_data, function(x) x[[1]])

# Bind all tibbles into one data frame
combined_df <- do.call(rbind, tibbles)

# Convert the data frame to a matrix
combined_matrix <- as.matrix(combined_df)

#extract value from xml data
extract_value <- function(xml_string) {
  xml <- read_xml(xml_string)
  xml_text(xml)
}

combined_matrix <- apply(combined_matrix, c(1, 2), extract_value)

# downloading indexed clinical data
prad_indexed_clinical_query <- GDCquery_clinic(project = "TCGA-PRAD",
                                type = "Clinical"
)

GDCdownload(prad_indexed_clinical_query)

# data preparing
prad_indexed_clinical_data <- GDCprepare(prad_indexed_clinical_query)
prad_clinical_expMat <- assay(prad_clinical_data)
# save the matrix
write.csv(combined_matrix, 'C:/.../combined_matrix.csv')
