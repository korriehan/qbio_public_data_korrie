#Installs TCGAbiolinks if not already present
if (!require(TCGAbiolinks)) BiocManager::install("TCGAbiolinks")

#What line of code needs to be added before using TCGAbiolinks?
library(TCGAbiolinks)

barcodes_clinic <- c("TCGA-BH-A0DG","TCGA-A2-A0YF","TCGA-AN-A04A","TCGA-AR-A1AP", "TCGA-A2-A0T3",
                     "TCGA-E2-A154", "TCGA-AO-A12F", "TCGA-A2-A0YM", "TCGA-BH-A0E0", "TCGA-AN-A0FL")

install.packages("survival")
install.packages("survminer")
install.packages("arsenal")
library(survival)
library(survminer)
library(arsenal)

clin_query <- GDCquery(project = "TCGA-BRCA", data.category="Clinical", file.type="xml")
GDCdownload( clin_query )
clinic <- GDCprepare_clinic(clin_query, clinical.info="patient")
names(clinic)[names(clinic) == "days_to_last_followup"] = "days_to_last_follow_up"

age_clinical = clinic$age_at_initial_pathologic_diagnosis
clinic$age_category = ifelse(age_clinical < 40, "Young", ifelse(age_clinical >= 60, "Old", "Mid"))

install.packages("tableone")
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
write.csv(df,"data/tableby.csv", row.names=FALSE)
