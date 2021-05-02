#Installs TCGAbiolinks if not already present
if (!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")

#What line of code needs to be added before using TCGAbiolinks?
library(TCGAbiolinks)

barcodes_clinic <- c("TCGA-BH-A0DG","TCGA-A2-A0YF","TCGA-AN-A04A","TCGA-AR-A1AP", "TCGA-A2-A0T3",
                     "TCGA-E2-A154", "TCGA-AO-A12F", "TCGA-A2-A0YM", "TCGA-BH-A0E0", "TCGA-AN-A0FL")

#install.packages("survival", repos = "http://cran.us.r-project.org")
#install.packages("survminer", repos = "http://cran.us.r-project.org")
#install.packages("arsenal", repos = "http://cran.us.r-project.org")
library(survival)
library(survminer)
library(arsenal)

#remove.packages("ggplot2")

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


###MAF###
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


#survival probability over time for geneset, mutant vs. wild type
mafSurvGroup(maf = maf, geneSet = c("TP53", "PIK3CA"), time = "time", Status = "Overall_Survival_Status")

mafSurvGroup(maf = maf, geneSet = c("TP53", "MUC16"), time = "time", Status = "Overall_Survival_Status")

mafSurvGroup(maf = maf, geneSet = c("MUC16", "PIK3CA"), time = "time", Status = "Overall_Survival_Status")


#survival probability over time for gene, mutant vs. wild type
mafSurvival(maf = maf, genes = 'TP53', time = 'time', Status = 'Overall_Survival_Status', isTCGA = TRUE)

mafSurvival(maf = maf, genes = 'PIK3CA', time = 'time', Status = 'Overall_Survival_Status', isTCGA = TRUE)

mafSurvival(maf = maf, genes = 'MUC16', time = 'time', Status = 'Overall_Survival_Status', isTCGA = TRUE)

mafSurvival(maf = maf, genes = 'ESR1', time = 'time', Status = 'Overall_Survival_Status', isTCGA = TRUE)

dev.off()
