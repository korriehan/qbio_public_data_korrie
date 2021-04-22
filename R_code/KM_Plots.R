# Loading packages
library(TCGAbiolinks)
library(SummarizedExperiment)
library(survival)
library(survminer)

# Barcodes
# barcodes_rnaseq <- c("TCGA-BH-A0DG-01A-21R-A12P-07","TCGA-A2-A0YF-01A-21R-A109-07",
#                     "TCGA-AN-A04A-01A-21R-A034-07","TCGA-AR-A1AP-01A-11R-A12P-07",
#                    "TCGA-A2-A0T3-01A-21R-A115-07", "TCGA-E2-A154-01A-11R-A115-07" )
# barcodes_clinic <- c("TCGA-BH-A0DG","TCGA-A2-A0YF","TCGA-AN-A04A","TCGA-AR-A1AP", "TCGA-A2-A0T3",
#                     "TCGA-E2-A154", "TCGA-AO-A12F", "TCGA-A2-A0YM", "TCGA-BH-A0E0", "TCGA-AN-A0FL")

# Accessing RNAseq data 
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")
# barcode = barcodes_rnaseq
# GDCdownload(query)
sum_exp <- GDCprepare(query)
patient_data <- colData(sum_exp)

# Accessing Clinical data
clin_query <- GDCquery(project = "TCGA-BRCA", data.category="Clinical", file.type = "xml")
# barcode=barcodes_clinic
# GDCdownload( clin_query ) #only need this command once. This downloads the files onto your system.
clinic <- GDCprepare_clinic(clin_query, clinical.info="patient")
names(clinic)[names(clinic) == "days_to_last_followup"] = "days_to_last_follow_up"

# Getting counts for specific genes, replacing column names with shorter barcodes (to match clinical)
htseq_counts <- assays(sum_exp)$"HTSeq - Counts"
colnames(htseq_counts) <- patient_data$patient

# # Matching clinical data sample order to RNAseq sample order
row_order <- match(colnames(htseq_counts), clinic$bcr_patient_barcode) 
# # match function takes two parameters: first is the order you wish to match something to, and the second is the dataframe
# # you wish to alter the order of.
clinic_ordered  <- clinic[row_order, ]

# # Get rid of nonmatching samples in clinical and htseq
# matching <- which(clinic_ordered$bcr_patient_barcode %in% colnames(htseq_counts))
# # which function basically only takes the values of clinic_ordered$bcr_patient_barcode that are also found in the colnames of htseq
# clinic_matched <- clinic_ordered[matching,]

# Adding age to clinical data
age_clinical = clinic$age_at_initial_pathologic_diagnosis
clinic$age_category = ifelse(age_clinical < 50, "Young", ifelse(age_clinical >= 50, "Old"))

# Accessing counts data for ESR1, categorize expression, and add to clinical data
ESR1_mask <- rowData(sum_exp)$external_gene_name == "ESR1"
ESR1_ENSG_name <- rowData(sum_exp)$ensembl_gene_id[ ESR1_mask ]
ESR1_counts <- htseq_counts[ESR1_ENSG_name, clinic$bcr_patient_barcode]
ESR1_quartiles <- quantile(ESR1_counts) # Categorizing the expression level based on quartile analysis
ESR1_expression_level <- ifelse(ESR1_counts > ESR1_quartiles[4], "High", ifelse(ESR1_counts < ESR1_quartiles[2], "Low", "Mid"))
clinic$ESR1_expression = ESR1_expression_level



# Splitting clinic_matched based on age
clinic_old <- subset(clinic, clinic$age_category=="Old")
clinic_young <- subset(clinic, clinic$age_category=="Young")

# Making KM plots for ESR1
TCGAanalyze_survival(clinic_old, "ESR1_expression", main="Kaplan-Meier Survival Curves for Old Patients with Varying ESR1 Expression", filename = "ESR1old.pdf")
TCGAanalyze_survival(clinic_young, "ESR1_expression", main="Kaplan-Meier Survival Curves for Young Patients with Varying ESR1 Expression", filename = "ESR1young.pdf")

