library(ggplot2)
library(TwoSampleMR)
library(dplyr)
library(ggrepel)
library(readxl)
library(openxlsx)
library(data.table)
library(eQTpLot)
library(ggpubr)

setwd("#")

#################
#Hypothesis-free#
#################
myhf <- read.xlsx("MR_HF.xlsx")
myhf <- myhf[!is.na(myhf$rr_p),]
myhf <- myhf[myhf$rr_p < 5.43e-5,]
myhf$rr_p <- as.numeric(myhf$rr_p)
myhf <- arrange(myhf, desc(lnrr), rr_p)
myhf$rr_p <- format.pval(pv=myhf$rr_p, digits = 3, eps = 1e-200)
colnames(myhf)[12] <- "SNP"
myhf$rr_p[myhf$rr_p > 0.05 & myhf$rr_p <= 1.1] <- lapply(myhf$rr_p[myhf$rr_p > 0.05 & myhf$rr_p <= 1.1], function(x){signif(as.numeric(x), 3)})
myhf$rr_p <- as.character(myhf$rr_p)
myhf<-split_exposure(myhf)
myhf[myhf$exposure == "Treatment/medication code: warfarin",1] <- "Warfarin"
myhf[myhf$exposure == "Non-cancer illness code  self-reported: mania/bipolar disorder/manic depression",1] <- "Mania/bipolar/manic depression"
myhf[myhf$exposure == "Treatment/medication code: carbimazole",1] <- "Carbimazole"
myhf[myhf$exposure == "Non-cancer illness code  self-reported: hyperthyroidism/thyrotoxicosis",1] <- "Hyperthyroidism/thyrotoxicosis"
myhf[myhf$exposure == "Diagnoses - main ICD10: I83 Varicose veins of lower extremities",1] <- "Varicose veins of lower extremities"
myhf[myhf$exposure == "Non-cancer illness code  self-reported: varicose veins",1] <- "Varicose veins"
myhf[myhf$exposure == "Long-standing illness  disability or infirmity",1] <- "Disability or infirmity"
myhf[myhf$exposure == "Non-cancer illness code  self-reported: chronic obstructive airways disease/copd",1] <- "Chronic obstructive pulmonary disease"
myhf <- myhf[myhf$exposure != "Adrenic acid (22:4n6)",]

hf_plot <- forest_plot_1_to_many(myhf, b="lnrr", se="lnrr_se",
                                   ao_slc=FALSE, lo=-1.1, up=5.5, col1_title = "Trait",
                                   TraitM="exposure",col1_width=3.0,by=NULL,
                                   addcols=c("SNP","rr_p"), addcol_widths=c(0.5,0.8),
                                   addcol_titles=c("SNP","P-val"),
                                   col_text_size = 3, subheading_size = 11, shape_points = 1,
                                   xlab = "Log Risk Ratio for DVT per unit increase in trait")

ggsave(
  "hypothesisFreeMR/hf_mr.tiff",
  plot = hf_plot,
  device = "tiff",
  width = 173,
  height = 210,
  units = c("mm"),
  dpi = 300,
)

#################
#   pQTL BMI    #
#################
# Load Lucy and Zaghlool MR data
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

mypqtl <- read.xlsx("s_tables.xlsx", sheet = 3, rows = 2:27)
mypqtl <- as.data.table(mypqtl)
setnames(mypqtl, c("P-value*", "Log.Risk.Ratio*"), c("RR_P", "logRR"))
mypqtl <- mypqtl[, -c(7, 8)]
mypqtl <- mypqtl[!is.na(mypqtl$RR_P),]
mypqtl$RR_P <- as.numeric(mypqtl$RR_P)
mypqtl <- arrange(mypqtl, desc(logRR), RR_P)
mypqtl$RR_P <- format.pval(pv=mypqtl$RR_P, digits = 3, eps = 1e-200)
colnames(mypqtl)[4] <- "SNP"
colnames(mypqtl)[7] <- "logRR.SE"
mypqtl$RR_P[mypqtl$RR_P > 0.05 & mypqtl$RR_P <= 1.1] <- lapply(mypqtl$RR_P[mypqtl$RR_P > 0.05 & mypqtl$RR_P <= 1.1], function(x){signif(as.numeric(x), 3)})
mypqtl$RR_P <- as.character(mypqtl$RR_P)
mypqtl<-split_exposure(mypqtl)
mypqtl <- mypqtl[mypqtl$GWAS.Author != "Neale lab", ]
test <- mypqtl[, ! names(mypqtl) %in% c("beta", "se", "P_val")]

pqtl_plot <- forest_plot_1_to_many(mypqtl, b="logRR", se="logRR.SE",
                                 ao_slc=F,lo=-2.0, up=1.3, col1_title = "Protein",
                                 TraitM="Gene.symbol", col1_width=0.8, by=NULL,
                                 addcols=c("GWAS.Author", "MR_method", "SNP","RR_P"), addcol_widths=c(0.8,1.7,0.5,0.7),
                                 addcol_titles=c("GWAS Author", "MR method", "No.SNP", "P-val"),
                                 col_text_size = 5, subheading_size = 12, 
                                 xlab = "Log Risk Ratio for DVT per unit increase in circulating protein level")

ggsave(
  "s_pqtl.tiff",
  plot = pqtl_plot,
  device = "tiff",
  width = 14,
  height = 8,
  units = c("in"),
  dpi = 300,
)

# Just the important pQTLs
mypqtl$RR_P <- as.numeric(mypqtl$RR_P)
mypqtl <- mypqtl[RR_P < 0.003,]
pqtl_bmi <- pqtl_bmi[pqtl_bmi$Target.full.name %in% mypqtl$Exposure, ]
setnames(pqtl_bmi, c("Target.full.name", "Beta_coefficient", "SE", "P_val"), 
         c("Exposure", "BETA_BMI", "SE_BMI", "P_BMI"))
mypqtl <- merge(mypqtl, pqtl_bmi[, c("Exposure", "BETA_BMI", "SE_BMI", "P_BMI")])
mypqtl <- melt(mypqtl, 
     id.vars = c("Exposure", "Gene.symbol"), 
     measure.vars = list(c("BETA_BMI", "logRR", "Proportion.mediated.(%).by.protein"), 
                         c("SE_BMI", "logRR.SE"),
                         c("P_BMI", "RR_P")),
     variable.name = c("Type"),
     value.name = c("Effect", "SE", "P")
)
mypqtl$Exposure <- "Body mass index"
mypqtl[, Type := fcase(
  mypqtl$Type == 1, "BMI->Protein",
  mypqtl$Type == 2, "Protein->DVT",
  mypqtl$Type == 3, "BMI->Protein->DVT"
)]
mypqtl[, Analysis := fcase(
  mypqtl$Type == "BMI->Protein", "One-sample MR",
  mypqtl$Type == "Protein->DVT", "Two-sample MR",
  mypqtl$Type == "BMI->Protein->DVT", "Mediation MR"
)]
mypqtl$P <- format.pval(pv=mypqtl$P, digits = 3, eps = 1e-200, na.form = "N/A")
mypqtl[9, "P"] <- "N/A"

pqtl_dvt <- forest_plot_1_to_many(mypqtl[-c(7:9), ], b="Effect", se="SE",
                                   ao_slc=F,lo=-2.0, up=1.0, col1_title = "Relationship type",
                                   TraitM="Type", col1_width=1.2, by=c("Gene.symbol"),
                                   addcols=c("Analysis", "P"), addcol_widths=c(1.0, 0.7),
                                   addcol_titles=c("Analysis type", "P-val"),
                                   col_text_size = 5, subheading_size = 12, shape_points = 1,
                                   xlab = "Log Risk Ratio on DVT per unit increase in trait")

ggsave(
  "pqtl.tiff",
  plot = pqtl_dvt,
  device = "tiff",
  width = 14,
  height = 8,
  units = c("in"),
  dpi = 300,
)
