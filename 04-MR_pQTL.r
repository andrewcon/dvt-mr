# Set the working directory------
setwd("#")

# Load libraries------
library(TwoSampleMR)
library(MRInstruments)
library(ggplot2)
library(openxlsx)
library(ieugwasr)
library(data.table)

# pQTL data prep------
# Load the summary statistics for pQTL associated with BMI from Lucy's analysis
pqtl_files <- list.files("pqtl_sum_stats", pattern = "*.tsv", full.names = TRUE)
pqtl_bmi <- read.xlsx("BMI_protein_supplementary_Tables_Int_J_O.xlsx", sheet = 9, startRow = 2)
pqtl_bmi <- pqtl_bmi[pqtl_bmi$P_val < 1.4e-5,]
pqtl_bmi <- pqtl_bmi[, c("SomaID", "Target.full.name", "UniProtID", "Beta_coefficient", "SE", "P_val", "N")]
pqtl_bmi2 <- read.xlsx("41467_2021_21542_MOESM8_ESM.xlsx", sheet = 1, startRow = 3)
pqtl_bmi2 <- pqtl_bmi2[pqtl_bmi2$P_value_protein < 5.43e-5,]
pqtl_bmi2 <- pqtl_bmi2[, c("SOMA_SeqID", "TargetFullName", "UniProt", 
                           "Beta_protein", "SE_protein", "P_value_protein", "N")]
setnames(pqtl_bmi2, 
         c("SOMA_SeqID", "TargetFullName", "UniProt", "Beta_protein", "SE_protein", "P_value_protein"), 
         c("SomaID", "Target.full.name", "UniProtID", "Beta_coefficient", "SE", "P_val")) 
pqtl_bmi <- rbind(pqtl_bmi, pqtl_bmi2)

# Load the MR-Base GWAS catalog
ao <- available_outcomes()

# Get the names of pQTL to extract from MR-Base
pqtl_names <- c("Leptin", "Receptor-type tyrosine-protein phosphatase delta", "PILR alpha-associated neural protein", 
                "Inhibin beta C chain", "Fumarylacetoacetase", "C5", "fatty acid binding protein 4", 
                "SHBG", "Insulin-like growth factor-binding protein 1", "Insulin-like growth factor-binding protein 2", 
                "Plasminogen activator inhibitor 1", 
                "WAP, Kazal, immunoglobulin, Kunitz and NTR domain-containing protein 2", 
                "Dickkopf-related protein 3", "galectin 3", "Growth hormone receptor", 
                "GDF2", "Netrin receptor UNC5D", 
                "Neurogenic locus notch homolog protein 1", "Hepatocyte growth factor receptor", 
                "Antithrombin-III", "C-Reactive protein level", 
                "C-reactive protein", "Neural cell adhesion molecule 1, 120 kDa isoform", 
                "Protein jagged-1", "Cystatin-M", 
                "Endothelial cell-specific molecule 1", 
                
                "Inhibin beta A chain", "Inhibin beta A chain:Inhibin beta B chain heterodimer", 
                "Inhibin beta A chain:Inhibin beta B chain heterodimer"
)
ao_pqtl <- subset(ao[ao$trait %in% pqtl_names,], select=c(trait, id, sample_size, group_name, author))
#ao_pqtl <- ao_pqtl[!(startsWith(ao_pqtl$id, "ukb") | 
#                       startsWith(ao_pqtl$id, "ieu") | 
#                       grepl("prot-c-2575_5_5", ao_pqtl$id)),]
pqtl_exp_dat <- extract_instruments(outcomes=ao_pqtl$id, clump = TRUE)
pqtl_exp_dat <- as.data.table(pqtl_exp_dat)
#pqtl_exp_dat[, temp_EA := ifelse(beta.exposure >= 0, effect_allele.exposure, other_allele.exposure)]
#pqtl_exp_dat[, temp_NEA := ifelse(beta.exposure < 0, effect_allele.exposure, other_allele.exposure)]
#pqtl_exp_dat[, eaf.exposure := ifelse(beta.exposure >= 0, eaf.exposure, 1 - eaf.exposure)]
#setnames(pqtl_exp_dat, c("temp_EA", "temp_NEA"), c("effect_allele.exposure", "other_allele.exposure"))
#pqtl_exp_dat <- pqtl_exp_dat[, -c(9, 10)]
#pqtl_exp_dat[, beta.exposure := ifelse(beta.exposure >= 0, beta.exposure, -1*beta.exposure)]

#Perform MR 
id_outcome <- "ukb-a-65"
outcome_dat <- extract_outcome_data(pqtl_exp_dat$SNP, id_outcome)
dat <- harmonise_data(pqtl_exp_dat, outcome_dat)
res <- mr(dat, method_list = c("mr_ivw", "mr_weighted_median", 
                               "mr_weighted_mode", "mr_simple_mode", 
                               "mr_wald_ratio"))
res <- as.data.table(res)
res <- split_exposure(res)
res[res$exposure == "fatty acid binding protein 4", "exposure"] <- "Fatty acid binding protein 4"
pqtl_bmi[pqtl_bmi$Target.full.name == "Fatty acid-binding protein", "Target.full.name"] <- "Fatty acid binding protein 4"
res <- merge(res, pqtl_bmi[,c("Target.full.name", "Beta_coefficient", "SE")], by.x="exposure", by.y="Target.full.name")
setnames(res, c("Beta_coefficient", "SE"), c("BETA_BMI_pQTL", "SE_BMI_pQTL"))

data <- res
data$id.outcome <- gsub('\\:', '-', data$id.outcome)
data$id.outcome <- tolower(data$id.outcome)

#find number of cases and controls
ao$id.outcome <- ao$id
#ao <- as.data.table(ao)

df <- merge(x=data, y=ao, by='id.outcome', all.x=T)
df <- subset(df, select = c("id.outcome", "exposure", "outcome","id.exposure", "method", 
                            "sample_size","b","se","pval", 
                            "nsnp.x","ncase","ncontrol", "BETA_BMI_pQTL", 
                            "SE_BMI_pQTL"))
df <- merge(x=df, y=ao[, c("id", "author", "consortium", "population")], by.x='id.exposure', by.y = 'id', all.x=T)
df <- df[population == "European",]
df <- df[df$id.exposure != "ukb-d-30710_raw",]

res <- df
res$outcome <- "Deep vein thrombosis"
res$beta_BMI_indirect_DVT <- res$lnrr*res$BETA_BMI_pQTL
res$se_BMI_indirect_DVT <- res$lnrr_se*res$SE_BMI_pQTL
res$lci_BMI_indirect_DVT <- res$beta_BMI_indirect_DVT - (1.96*res$se_BMI_indirect_DVT)
res$uci_BMI_indirect_DVT <- res$beta_BMI_indirect_DVT + (1.96*res$se_BMI_indirect_DVT)
res <- res[!is.nan(res$rr_p),]
#res$or_BMI_indirect_DVT <- exp(res$beta_BMI_indirect_DVT)
#res$lci_BMI_indirect_DVT <- res$or_BMI_indirect_DVT - (1.96*res$se_BMI_indirect_DVT)
#res$uci_BMI_indirect_DVT <- res$or_BMI_indirect_DVT + (1.96*res$se_BMI_indirect_DVT)

# Write all results
fwrite(res, file = "mr_results.csv", sep = ",", quote = F, row.names = F)
write.xlsx(res, "mr_results.xlsx", overwrite = T)


