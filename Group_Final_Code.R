#Final Project - Ellison Data Analysis Group -  Multi-Omic Analysis of ESR1 Expression#

# We first need to download all of the required packages to access the data and generate
# the graphs created in this analysis. 

if (!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")
if(!requireNamespace("SummarizedExperiment"))BiocManager::install(c("SummarizedExperiment"))
if(!requireNamespace("maftools"))BiocManager::install(c("maftools"))
if(!requireNamespace("arsenal"))install.packages(c("arsenal"))
if(!requireNamespace("survival"))install.packages(c("survival"))
if(!requireNamespace("survminer"))install.packages(c("survminer"))
if (!require(DESeq2)) BiocManager::install("DESeq2")

# We then use the "library()" command to load all of our packages into the current 
# R environment. 

library(TCGAbiolinks) # load TCGAbiolinks library
library(maftools) # load maftools library
library(SummarizedExperiment)# load SummarizedExperiment library
library(arsenal) # load arsenal library
library(survival) # load survival library
library(survminer) # load survminer library
library(DESeq2) # load DESeq2 library

### Mutation Rates ###

# Lolipop2 Graphs #
query <- GDCquery(project = "TCGA-BRCA",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification",
                  workflow.type = "HTSeq - Counts")
GDCdownload(query) #only need this line of code ONCE to download the data
sum_exp <- GDCprepare(query)

mutation <- GDCquery_Maf(tumor = "BRCA",save.csv=TRUE, pipeline="mutect2")
maf_dataframe = read.maf(mutation)

# Here we are creating a NEW COLUMN in patient_data called "age_category"
patient_data <- colData(sum_exp)
patient_ages <- patient_data$paper_age_at_initial_pathologic_diagnosis
# NOTE: This will NOT be added to colData(sum_exp). Instead it will only be added 
# to the patient_data data table.

# Create new column with age categories for all patient samples
patient_data$age_category = ifelse(patient_ages < 50, "Young", "Old")

# Get shortened patient barcodes so that we can compare with
short_maf <- substr(maf_dataframe@clinical.data$Tumor_Sample_Barcode, 1,12)

# Create a new column in maf_dataframe
maf_dataframe@clinical.data$short_barcodes <- short_maf

# Age vector has either young or old based on age category at that barcode
maf_ages <- patient_data[short_maf, "age_category"]

# Add new column to the maf dataframe containing age info
maf_dataframe@clinical.data$Ages <- maf_ages

# Extract codes for each age group
young_codes <- maf_dataframe@clinical.data[maf_dataframe@clinical.data$Ages =="Young",]
old_codes <- maf_dataframe@clinical.data[maf_dataframe@clinical.data$Ages =="Old",]

# Create maf subsets for each age group
young_maf <- subsetMaf(maf_dataframe, tsb = young_codes$Tumor_Sample_Barcode)
old_maf <- subsetMaf(maf_dataframe, tsb = old_codes$Tumor_Sample_Barcode)

# Create the Lollipop Plot
jpeg("LolliPlot2OldYoung.jpeg")
lollipopPlot2(m1=young_maf, m2=old_maf, m1_name = "Young", m2_name = "Old", gene="ESR1")

### RNA Expression ###

patient_data <- colData(sum_exp)  # save colData in patient_data variable
patient_ages <- patient_data$paper_age_at_initial_pathologic_diagnosis  # set patient_ages as age of initial diagnosis
patient_data$age_category = ifelse(patient_ages < 50, "Young", "Old")  # create new age_category column
htseq_counts <- assays(sum_exp)$"HTSeq - Counts"  # get gene expression counts

patient_data$ESR1_counts = htseq_counts["ENSG00000091831",]  # create column of ESR1 counts in patient_data

patient_data$ESR1_counts_log = sapply(htseq_counts["ENSG00000091831",], log10) # create column of log(ESR1) counts in patient_data

all_young <- sum(patient_data$age_category == "Young", na.rm = TRUE)  # get # of young patients
all_old <- sum(patient_data$age_category == "Old", na.rm=TRUE)  # get # of old patients

patient_data <- subset(patient_data, patient_data$paper_BRCA_Subtype_PAM50 == "LumB" | patient_data$paper_BRCA_Subtype_PAM50 == "LumA")  # create subset of data with only ER-positive
ER_young <- sum(patient_data$age_category == "Young", na.rm = TRUE)  # get # of er-positive young patients
ER_old <- sum(patient_data$age_category == "Old", na.rm=TRUE) # get # of er-positive old patients

ER_young/all_young  # get percentage of young patients that are ER-positive
ER_old/all_old  # get percentage of old patients that are ER-positive

require(stats)  # load stats
jpeg("log(ESR1)_RNA_scatterplot")  # create jpeg
reg_log <- lm(ESR1_counts_log ~ paper_age_at_initial_pathologic_diagnosis, data = patient_data)  # perform regression on log(ESR1)
coeff_log <- coefficients(reg_log)  # get coefficients to get equation
coeff_log  # print coefficients
plot(x = patient_data$paper_age_at_initial_pathologic_diagnosis, y = patient_data$ESR1_counts_log, main = "y = 0.01x + 3.87")  # plot scatterplot
abline(reg_log)  # plot regression line
dev.off()  # plot complete

jpeg("ESR1_RNA_scatterplot")  # create jpeg
reg <- lm(ESR1_counts ~ paper_age_at_initial_pathologic_diagnosis, data = patient_data)  # regression on ESR1
coeff <- coefficients(reg)  # get coefficients to get equation
coeff  # print coefficients
plot(x = patient_data$paper_age_at_initial_pathologic_diagnosis, y = patient_data$ESR1_counts, main = "y = 1015.7x - 14201.4")  # plot scatterplot
abline(reg)  # plot regression line 
dev.off()  # plot complete

jpeg("log(ESR1)_RNA_boxplot")  # create jpeg
boxplot(ESR1_counts_log~factor(age_category, levels=c("Young", "Old")), data = patient_data, main = "Boxplot of log(HTSeq - Counts) for ESR1 by Age Category")  # plot boxplot of log(ESR1)
dev.off()  # plot complete

library(ggpubr)  # load ggpubr library
patient_data$age_category_factor <- as.factor(patient_data$age_category)  # make age_category a factor
patient_data <- as.data.frame(patient_data) # make sure patient_data is a dataframe
compare_means(ESR1_counts ~ age_category, data=patient_data)  # calculates p-value

###creating Maf Survival KM plot###

# From here, we create young and old patient KM plots with low, mid, and high expression level#

# Accessing Clinical data
clin_query <- GDCquery(project = "TCGA-BRCA", data.category="Clinical", file.type = "xml")
#GDCdownload( clin_query ) #only need this command once. This downloads the files onto your system.
clinic <- GDCprepare_clinic(clin_query, clinical.info="patient")
names(clinic)[names(clinic) == "days_to_last_followup"] = "days_to_last_follow_up"

# Getting counts for specific genes, replacing column names with shorter barcodes (to match clinical)
htseq_counts <- assays(sum_exp)$"HTSeq - Counts"
colnames(htseq_counts) <- patient_data$patient

# Matching clinical data sample order to RNAseq sample order
row_order <- match(colnames(htseq_counts), clinic$bcr_patient_barcode) 

# match function takes two parameters: first is the order you wish to match something to, and the second 
# is the dataframe you wish to alter the order of.
clinic_ordered  <- clinic[row_order, ]

# Get rid of nonmatching samples in clinical and htseq
matching <- which(clinic_ordered$bcr_patient_barcode %in% colnames(htseq_counts))

# which function basically only takes the values of clinic_ordered$bcr_patient_barcode that are also found in the colnames of htseq
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

# From here, we create a Suvival analysis Plot for ERS1 comparing mutant and wild type genes. 

# Assigns 1 to each value
clinic$Overall_Survival_Status <- 1 

# Assigns 0 when vital status is not equal to "Dead"
clinic$Overall_Survival_Status[which(clinic$vital_status != "Dead")] <- 0
clinic$time <- clinic$days_to_death

# When days_to_death is NA, that value in the time column is replaced by the corresponding value in days_to_last_follow_up
clinic$time[is.na(clinic$days_to_death)] <- clinic$days_to_last_follow_up[is.na(clinic$days_to_death)]
colnames(clinic)[1] <- "Tumor_Sample_Barcode"

# Reads in as csv
my_maf <- data.table::fread("GDCdata/TCGA.BRCA.mutect.995c0111-d90b-4140-bee7-3845436c3b42.DR-10.0.somatic.maf.csv")

# Reads in as maf
maf <- read.maf(maf = my_maf, clinicalData = clinic, isTCGA = TRUE)

#survival probability over time for gene, mutant vs. wild type
pdf("mafSurvival.pdf")
mafSurvival(maf = maf, genes = 'ESR1', time = 'time', Status = 'Overall_Survival_Status', isTCGA = TRUE)
dev.off()

### DESeq2 ###

#access counts
counts <- assays(sum_exp)$"HTSeq - Counts"
# remove patients with NA in age, PAM50, pathology
patients_no_NA_mask <- ( !is.na(colData(sum_exp)$paper_age_at_initial_pathologic_diagnosis)
                         & !is.na(colData(sum_exp)$paper_BRCA_Subtype_PAM50)
                         & !is.na(colData(sum_exp)$paper_BRCA_Pathology)
                         & !colData(sum_exp)$paper_BRCA_Pathology == "NA" )
#access the patient_data from coldata and apply NA mask
patient_data <- colData(sum_exp)[ patients_no_NA_mask, ]

# apply NA mask to counts
counts <- counts[ , patients_no_NA_mask]

# create age_category column
patient_data$age_category <- ifelse(patient_data$paper_age_at_initial_pathologic_diagnosis < 50, "Young", "Old")


# make categorical columns factors 
patient_data$age_category <- factor(patient_data$age_category, levels=c("Young", "Old"))
patient_data$paper_BRCA_Subtype_PAM50 <- factor( patient_data$paper_BRCA_Subtype_PAM50, levels=c("Her2","LumA","LumB","Basal","Normal") )
patient_data$paper_BRCA_Pathology <- factor( patient_data$paper_BRCA_Pathology, levels=c("IDC","Other","Mixed","ILC") )

# prepare DESeq2 data matrix
dds_with_adjustment <- DESeqDataSetFromMatrix(countData = counts, colData = patient_data, design = ~paper_BRCA_Pathology+ paper_BRCA_Subtype_PAM50 +age_category)
# perform DESeq2 analysis
dds_obj <- DESeq(dds_with_adjustment)
# load results 
results <- results(dds_obj, contrast=c("age_category", "Young",'Old'))
# create mask of NA values
results_NA <- !is.na(results[,"padj"])
# apply mask to results
results <- results[results_NA,]
# create new column of fold change
results$FoldChange <- 2^results$log2FoldChange
# save results csv file
write.csv(results, "/PATHWAY/DESeq2.csv")
# separate significant genes from insignificant
results_significant_adjp <- results[results$padj < padj_threshold,]  #separate significant genes from insignificant ones
# separate up regulated in young
results_sig_up_regulated <- results_significant_adjp[results_significant_adjp$log2FoldChange > log2FC_threshold,]  
# separate down regulated in young
results_sig_down_regulated <- results_significant_adjp[results_significant_adjp$log2FoldChange < -log2FC_threshold,]

gene_information <- rowData(sum_exp)  #access gene information from sum_exp rowData
matching_rows_up <- gene_information[gene_information$ensembl_gene_id %in% rownames(results_sig_up_regulated),]  # create dataframe with only gene information for genes that were up regulated
results_sig_up_regulated$gene_name <- matching_rows_up$external_gene_name # create new column in up regulated df for common gene name
matching_rows_down <- gene_information[gene_information$ensembl_gene_id %in% rownames(results_sig_down_regulated),] # create dataframe with only gene information for genes that were down regulated
results_sig_down_regulated$gene_name <- matching_rows_down$external_gene_name  # create new column in down regulated df for common gene names

#save as .csv files
write.csv(results_sig_up_regulated, "/PATHWAY/high_in_young.csv")
write.csv(results_sig_down_regulated, "/PATHWAY/low_in_young.csv")


### GSEA .rnk File ###

# create variable of gene names of genes in results
matching_rows_up <- gene_information[gene_information$ensembl_gene_id %in% rownames(results),]  
# create new gene_name column in results 
results$gene_name <- matching_rows_up$external_gene_name 
# create gsea_rnk dataframe with gene names column
gsea_rnk <- subset(results, select=c("gene_name"))
# add log2FC column
gsea_rnk$log2FC <- results$log2FoldChange
# write dataframe to .rnk file
write.table(gsea_rnk, "/PATHWAY/file.rnk", sep="\t", row.names=FALSE, col.names = FALSE, quote=FALSE)