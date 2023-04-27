library(devtools)
install_github("MRCIEU/TwoSampleMR")
library(ggplot2)
#link to the test version
library(TwoSampleMR)
#library(MRInstruments)
library("readxl")
rm(list=ls(all=TRUE)) 

##read in files
setwd("#")

data <- read.table("MR_results.txt",
                   header=T,sep="\t",quote="")
data$id.outcome <- gsub('\\:', '-', data$id.outcome)
data$id.outcome <- tolower(data$id.outcome)

#find number of cases and controls
ao<-available_outcomes()
ao$id.outcome<-ao$id

df <- merge(x=data, y=ao, by='id.outcome', all.x=T)
df <- subset(df, select = c("id.outcome", "exposure", "outcome","id.exposure",
                            "sample_size","b","se","pval",
                            "nsnp.x","ncase","ncontrol", "method"))

#work out log risk ratio
#b -- absolute change in risk due to exposure in uk biobank from linear regression model of case control status (0 or 1) on SNP

#risk.u = ncases/(ncases+ncontrols) #absolute risk in unexposed group
df$risk.u <- df$ncase/(df$ncase+df$ncontrol)

#risk.e = risk.u+b # absolute risk in exposed group
df$risk.e <- df$risk.u + df$b

#rr<-risk.e/risk.u #risk ratio 
df$rr <- df$risk.e/df$risk.u

#lnrr<-log(rr) # log risk ratio
df$lnrr <- log(df$rr)

#work out standard error for log risk ratio
df$b_lci = df$b-(1.96*df$se) #lower 95% CI for b
df$b_uci <- df$b+(1.96*df$se) #upper 95% CI for b
df$risk.e_lci <- df$risk.u+df$b_lci # lower CI for absolute risk in exposed group
df$risk.e_uci <- df$risk.u+df$b_uci # upper CI for absolute risk in exposed group
df$lnrr_lci <- log(df$risk.e_lci/df$risk.u) #lower CI for risk ratio
df$lnrr_uci <- log(df$risk.e_uci/df$risk.u) #upper CI for risk ratio
df$lnrr_se <- (df$lnrr_uci-df$lnrr_lci)/(1.96*2) # standard error for log risk ratio

######re-estimate p values#####
df$rr_z <- df$lnrr/df$lnrr_se

df$rr_p <- 2*pnorm(-abs(df$rr_z))

df <- subset(df, select = c("exposure", "outcome","id.exposure","id.outcome", 
                            "method", "sample_size", "pval", "rr",
                            "lnrr","lnrr_se","rr_p","nsnp.x", 
                            "lnrr_lci", "lnrr_uci"))
df$bonferroni <- p.adjust(df$rr_p, method = "bonferroni")

#result_file0 <- paste0("results/combined/summary/MR.single.ukb.v3.txt")
write.xlsx(df, "MR_results.xlsx")
