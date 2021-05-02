#Installs TCGAbiolinks if not already present
if (!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")

#What line of code needs to be added before using TCGAbiolinks?
library(TCGAbiolinks)

#install.packages("survival", repos = "http://cran.us.r-project.org")
#install.packages("survminer", repos = "http://cran.us.r-project.org")
#install.packages("arsenal", repos = "http://cran.us.r-project.org")
library(survival)
library(survminer)
library(arsenal)
library(SummarizedExperiment)

#creating age categories: young, mid, old
clin_query <- GDCquery(project = "TCGA-BRCA", data.category="Clinical", file.type="xml")
GDCdownload( clin_query )
clinic <- GDCprepare_clinic(clin_query, clinical.info="patient")
names(clinic)[names(clinic) == "days_to_last_followup"] = "days_to_last_follow_up"

age_clinical = clinic$age_at_initial_pathologic_diagnosis
clinic$age_category = ifelse(age_clinical < 40, "Young", ifelse(age_clinical >= 60, "Old", "Mid"))

#install.packages("tableone", repos = "http://cran.us.r-project.org")
library(tableone)
clinic_summary <- CreateTableOne(data = clinic)
subtypes <- TCGAquery_subtype(tumor = "BRCA")

age_subs <- subtypes$age_at_initial_pathologic_diagnosis
subtypes$age_category = ifelse(age_subs < 40, "Young", ifelse(age_subs >= 60, "Old", "Mid"))

names(subtypes) = make.names(names(subtypes))
table_arse <- tableby(age_category ~ (pathologic_stage)+ (mRNA.Clusters)+ (BRCA_Pathology),
                      data = subtypes, numeric.test="kwt", cat.test="chisq",
                      numeric.stats = c("Nmiss", "meansd"), total=FALSE)
df <- as.data.frame(summary(table_arse, text=TRUE, pfootnote=TRUE))
write.csv(df,'summarize_by_age.csv', row.names=FALSE)


###creating Maf Survival KM plot###
#BiocManager::install("maftools")
library(maftools)
#mutation <- GDCquery_Maf(tumor = "BRCA",save.csv=TRUE, pipeline="mutect2") #choose a pipeline
# after running ^^, navigate to the saved csv file. Open the csv file with below Code
# only need to query once. For repeating code, you can just read in the saved dataframe you created

#assigns 1 to each value
clinic$Overall_Survival_Status <- 1 

#assigns 0 when vital status is not equal to "Dead"
clinic$Overall_Survival_Status[which(clinic$vital_status != "Dead")] <- 0
clinic$time <- clinic$days_to_death

#when days_to_death is NA, that value in the time column is replaced by the corresponding value in days_to_last_follow_up
clinic$time[is.na(clinic$days_to_death)] <- clinic$days_to_last_follow_up[is.na(clinic$days_to_death)]
colnames(clinic)[1] <- "Tumor_Sample_Barcode"

#reads in as csv
my_maf <- data.table::fread("GDCdata/TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf.csv")

#reads in as maf
maf <- read.maf(maf = my_maf, clinicalData = clinic, isTCGA = TRUE)

pdf("mafSurvival.pdf")

#survival probability over time for gene, mutant vs. wild type
mafSurvival(maf = maf, genes = 'ESR1', time = 'time', Status = 'Overall_Survival_Status', isTCGA = TRUE)

dev.off()

###creating young and old patient KM plots with low, mid, and high expression level###
# Accessing RNAseq data 
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")
#GDCdownload(query)
sum_exp <- GDCprepare(query)
patient_data <- colData(sum_exp)

# Accessing Clinical data
clin_query <- GDCquery(project = "TCGA-BRCA", data.category="Clinical", file.type = "xml")
#GDCdownload( clin_query ) #only need this command once. This downloads the files onto your system.
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
matching <- which(clinic_ordered$bcr_patient_barcode %in% colnames(htseq_counts))
# # which function basically only takes the values of clinic_ordered$bcr_patient_barcode that are also found in the colnames of htseq
clinic_matched <- clinic_ordered[matching,]

# Adding age to clinical data
age_clinical = clinic_matched$age_at_initial_pathologic_diagnosis
clinic_matched$age_category = ifelse(age_clinical < 50, "Young", "Old")

# Accessing counts data for ESR1, categorize expression, and add to clinical data
ESR1_mask <- rowData(sum_exp)$external_gene_name == "ESR1"
ESR1_ENSG_name <- rowData(sum_exp)$ensembl_gene_id[ ESR1_mask ]
ESR1_counts <- htseq_counts[ESR1_ENSG_name, clinic_matched$bcr_patient_barcode]
ESR1_quartiles <- quantile(ESR1_counts) # Categorizing the expression level based on quartile analysis
ESR1_expression_level <- ifelse(ESR1_counts > ESR1_quartiles[4], "High", ifelse(ESR1_counts < ESR1_quartiles[2], "Low", "Mid"))
clinic_matched$ESR1_expression = ESR1_expression_level

# Splitting clinic_matched based on age
clinic_old <- subset(clinic_matched, clinic_matched$age_category=="Old")
clinic_young <- subset(clinic_matched, clinic_matched$age_category=="Young")

# Making KM plots for ESR1
TCGAanalyze_survival(clinic_old, "ESR1_expression", main="Kaplan-Meier Survival Curves for Old Patients with Varying ESR1 Expression", filename = "ESR1old.pdf")
TCGAanalyze_survival(clinic_young, "ESR1_expression", main="Kaplan-Meier Survival Curves for Young Patients with Varying ESR1 Expression", filename = "ESR1young.pdf")
