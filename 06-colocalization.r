##################################################
#### Set up R workspace 
##################################################
# Assign path to user-installed packages

library(devtools)
library(ggplot2)
library(TwoSampleMR)
library(readxl)
library(data.table)
library(snpStats)
library(truncnorm)
library(rrcov)
library(coloc)

rm(list=ls(all=TRUE))

setwd("#")

##################################################
#### Import table of all files which contain list
#### of SNPs for coloc analysis. 
####
#### Modified to also read outcome data from
#### opt_all
##################################################

opt_all <- read.table("3-colocalization_list.txt", header = TRUE, 
                      stringsAsFactors=F, sep = ",")
opt_all$id <- NA

for (i in 1:nrow(opt_all)) {
  # Load instrument file
  instrument_file <- as.vector(opt_all$instrument_file[i])
  exposure_dat <- read.delim(instrument_file, sep="\t", 
                             header=TRUE, stringsAsFactors=F)
  if(c("A2.x") %in% colnames(exposure_dat)) {
    setnames(exposure_dat, "A2.x", "A2")
  }
  
  # Flip eaf if > 0.5 as coloc requires maf
  exposure_dat$minor <- ifelse(exposure_dat$AAF > 0.5, exposure_dat$A2, exposure_dat$A1)
  exposure_dat$major <- ifelse(exposure_dat$AAF > 0.5, exposure_dat$A1, exposure_dat$A2)
  exposure_dat$maf <- ifelse(exposure_dat$AAF > 0.5, 1 - exposure_dat$AAF, exposure_dat$AAF)
  
  # Format data
  exposure_dat <- format_data(exposure_dat, type = "exposure", header = TRUE,
                              snp_col = "SNP", beta_col = "BETA",
                              se_col = "SE", eaf_col = "maf", effect_allele_col = "minor",
                              other_allele_col = "major", pval_col = "P", chr_col = "CHR", 
                              pos_col = "BP")
  
  # Load the outcome file
  outcome_file <- as.vector(opt_all$outcome_file[i])
  outcome_dat <- read.delim(outcome_file, header=TRUE, sep="\t", stringsAsFactors=F)
  if(c("A2.x") %in% colnames(outcome_dat)) {
    setnames(outcome_dat, "A2.x", "A2")
  }
  #colnames(outcome_dat)[8] <- "n"
  #outcome_dat$allele1 <- as.character(outcome_dat$allele1)
  #outcome_dat$allele0 <- as.character(outcome_dat$allele0)
  
  # Flip eaf if > 0.5 as coloc requires maf
  outcome_dat$minor <- ifelse(outcome_dat$AAF > 0.5, outcome_dat$A2, outcome_dat$A1)
  outcome_dat$major <- ifelse(outcome_dat$AAF > 0.5, outcome_dat$A1, outcome_dat$A2)
  outcome_dat$maf <- ifelse(outcome_dat$AAF > 0.5, 1 - outcome_dat$AAF, outcome_dat$AAF)

  outcome_dat <- format_data(outcome_dat, type = "outcome", header = TRUE,
                             snp_col = "SNP", beta_col = "BETA", se_col = "SE",
                             eaf_col = "maf", effect_allele_col = "minor", other_allele_col = "major", 
                             pval_col = "P", samplesize_col = "n", 
                             chr_col = "CHR", pos_col = "BP")
  outcome.name <- as.vector(opt_all$outcome[i])
  
  # Harmonise!
  dat <- harmonise_data(exposure_dat, outcome_dat)
  dat <- dat[dat$remove == "FALSE", ]
  cols = c("SNP", "beta.exposure", "beta.outcome", "se.exposure", "se.outcome", "eaf.exposure", "eaf.outcome")
  
  dat2 <- dat[,cols]
  trait_name <- tools::file_path_sans_ext(opt_all[i, "output"])
  trait_name <- sub(".*_", "", trait_name)
  opt_all[i, c("id")] <- trait_name
  if (exists("dat2") == TRUE) {
    write.table(dat2, file=paste0("cond_coloc/", opt_all[i, "output"]), 
                sep = "\t", col.names = T, row.names = F, quote = F)
  }
  
}

write.table(opt_all, "coloc_list.txt", row.names=F, quote=F, sep = ",")

##################################################
#### Imports and function neccessary for coloc
#### NB: MAF is MINOR allele frequency!
##################################################
#biocLite("snpStats")
#biocLite("truncnorm")
#biocLite("coloc")
library(snpStats)
library(truncnorm)
library(rrcov)
library(coloc)

coloc.analysis <- function(beta1, beta2, se1, se2, MAF1, MAF2, N2, s) {
	# Convert the inputs in order to run in coloc function.
	# type: quant (quantitative) for pQTL study, cc (case control) for CHD
	dataset1 <- list(beta = beta1, varbeta = se1^2, MAF = MAF1, type = "quant", N = 3300)
	
	# N for CHD cases is 60801, proportion which are cases is 60801/184305 = 0.329893383
	dataset2 <- list(beta = beta2, varbeta = se2^2, MAF = MAF2, type = "cc", s = s, N = N2)
	
	# Run the coloc analysis, setting the prior probabilities for association with each trait
	# (p1, p2) and both traits together (p12) as 0.002
	result <- coloc.abf(dataset1, dataset2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
	
	# Format the data for saving
	df <- data.frame(matrix(unlist(result$summary), nrow = 1, byrow = T))
	
	# Label the columns in df
	names(df) <- c("nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H1.abf", "PP.H3.abf", "PP.H4.abf")
	
	return(df)
}

coloc.analysis.quant <- function(beta1, beta2, se1, se2, MAF1, MAF2, N2) {
	# Conver the inputs in order to run in coloc function
	# type, quant for pQTL study
	dataset1 <- list(beta = beta1, varbeta = se1^2, MAF = MAF1, type = "quant", N = 3300)
	
	# type quant for quantitative trait study
	dataset2 <- list(beta = beta2, varbeta = se2^2, MAF = MAF2, type = "quant", N = N2)
	
	# Run the coloc analysis, setting the prior probabilities for association with each trait
	# (p1, p2) and both traits together (p12) as 0.002
	result <- coloc.abf(dataset1, dataset2, p1 = 1e-4, p2 = 1e-4, p12 = 1e-5)
	
	# Format the data for saving
	df <- data.frame(matrix(unlist(result$summary), nrow = 1, byrow = T))
	
	# Label the columns in df
	names(df) <- c("nsnps", "PP.H0.abf", "PP.H1.abf", "PP.H2.abf", "PP.H3.abf", "PP.H4.abf")
	
	return(df)
}

LD <- ieugwasr::ld_matrix(df$SNP, df$SNP, with_alleles = F)
LD_snps <- dimnames(LD)[[1]]
df <- df[df$SNP %in% LD_snps, ]
d1 <- list(beta = df$beta.exposure, varbeta = df$se.exposure^2, 
           LD = LD, snp = df$SNP, type = "quant", MAF = df$MAF1, N = 3300)
d2 <- list(snp = df$SNP, beta = df$beta.outcome, varbeta = df$se.outcome^2,
           LD = LD, snp = df$SNP, type = "cc", MAF = df$MAF2, N = 462933)

dataset1 <- runsusie(d1, suffix = 1)
dataset2 <- runsusie(d2, suffix = 2)
r_susie <- coloc.susie(dataset1 = dataset1, dataset2 = dataset2)

##################################################
#### Coloc using real data
##################################################
#setwd("")

opt_all <- read.table("cond_coloc/coloc_list.txt", header = T, stringsAsFactors = F, sep = ",")
results_df <- data.frame()

for (i in 1:nrow(opt_all)) {
  input <- as.vector(opt_all$output[i])
	input_file <- as.vector(opt_all$output[i]) # harmonised file
	
	try(df <- read.table(paste0("cond_coloc/", input_file), header = TRUE, stringsAsFactors = F))
	df <- df[complete.cases(df), ]
	
	if (nrow(df) != 0) {
		df$MAF1 <- NULL
		df$MAF2 <- NULL
		
		for (j in 1:nrow(df)) {
			if (df$eaf.exposure[j] > 0.5) {
				df$MAF1[j] = 1 - df$eaf.exposure[j]
			}
			else {
				df$MAF1[j] = df$eaf.exposure[j]
			}
			if (df$eaf.outcome[j] > 0.5) {
				df$MAF2[j] = 1 - df$eaf.outcome[j]
			}
			else {
				df$MAF2[j] = df$eaf.outcome[j]
			}
		}
		
		ids <- as.vector(opt_all$id[i])
		
		#if (ao[ao$id == ids, ]$category == "Disease") {
		if (0) {
			# sample size from local file
		  N1 <- as.numeric(opt_all$n_case[i])
			N2 <- as.numeric(opt_all$n_total[i])
			
			# sample size from MR-Base traits
			#N1 <- as.numeric(ao[ao$id == ids, ]$ncase)
			#N2 <- as.numeric(ao[ao$id == ids, ]$sample_size)
			s <- N1/N2
			
			# andrei added
			N1 <- 462933
			N2 <- 3300
			
			
			if (is.na(N1) == FALSE) {
				result <- coloc.analysis(df$beta.exposure, df$beta.outcome, df$se.exposure, df$se.outcome, df$MAF1, df$MAF2, N2, s)
			}
			else {
				result <- cbind(input, nrow(df), "NA", "NA", "NA", "NA", "NA")
			}
		}
		else {
			# sample size from local file
			N2 <- as.numeric(opt_all$n_total[i])
			
			# sample sie from MR-Base traits
			#N2 <- as.numeric(ao[ao$id == ids, ]$sample_size)
			result <- coloc.analysis.quant(df$beta.exposure, df$beta.outcome, df$se.exposure, df$se.outcome, df$MAF1, df$MAF2, N2)
		}
		
		result <- cbind(input, result)
	}
	else {
	  result <- cbind(input, 0, "NA", "NA", "NA", "NA", "NA")
	}
	
	results_df <- rbind(results_df, result)
	
	result_file <- paste0("cond_coloc/", "all.txt")
	if (exists("results_df")) {
		write.table(results_df, file = result_file, sep = "\t", col.names = T, row.names = F, quote = F)
	}
}

openxlsx::write.xlsx(results_df, "all.xlsx", overwrite = TRUE)

