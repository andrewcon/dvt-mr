#Hypothesis Free analysis of all instruments in MRBase to DVT
#Can be modified to use any trait as an outcome

#Functions
#************************************************************************************
setWorkDir <- function(){
  n <- tolower(readline(prompt="Would you like to set the working directory? (Y / N)"))
  if(!grepl("[yn]",n) || nchar(n) > 1)
  {
    print("Please type 'Y' or 'N'")
    return(instPkg())
  }else{
    if(n == 'y'){
      wD <- readline(prompt="Set working directory: ")
      y <- tolower(readline(prompt="Make sure you have no spelling errors. Proceed? (Y / N)"))
      if(!grepl("y",y) || nchar(y) > 1)
      {
        print("Going back...")
        return(setWorkDir())
      }
      setwd(wD)
    }
  }
}


instPkg <- function(){ 
  n <- tolower(readline(prompt="Would you like to install the packages neeeded for the analysis? (Y / N)"))
  if(!grepl("[yn]",n) || nchar(n) > 1)
  {
    print("Please type 'Y' or 'N'")
    return(instPkg())
  }else{
    
    if(n == 'y'){
      install.packages("devtools")
      install.packages("biomaRt")
      devtools::install_github("MRCIEU/TwoSampleMR")
      install.packages("plyr")
      install.packages("svMisc")
      library(devtools)
      devtools::install_github("MRCIEU/MRInstruments")
    }
    library(MRInstruments)
    library(TwoSampleMR)
    library(plyr)
    library(dplyr)
    library(openxlsx)
    library(svMisc)
    return('Libraries loaded.')
  }
}

curateAo <- function(dataframe){
  tts <- readline(prompt="Enter the trait you wish to use as an exposure (e.g Body Mass Index) NOT CASE SENSITIVE: ")
  n <- tolower(readline(prompt="Make sure you have no spelling errors. Proceed? (Y / N)"))
  if(!grepl("y",n) || nchar(n) > 1)
  {
    print("Going back...")
    return(curateAo(dataframe))
  }else{
    #Curate ao to remove non-European populations
    return(dataframe %>% 
             filter(substr(population,1,4) == "Euro") %>%
             filter(!grepl(tts, trait, ignore.case = TRUE)) %>%
             group_by(trait) %>%
             slice(which.max(sample_size)) %>%
             arrange(trait))
  }
}

exposureID <- function(){
  idnr <- readline(prompt="Enter the id for the exposure (e.g UKB-a:382 is a Waist circumference study) CASE SENSITIVE: ")
  n <- tolower(readline(prompt="Make sure you have no spelling errors. Proceed? (Y / N)"))
  if(!grepl("y",n) || nchar(n) > 1)
  {
    print("Going back...")
    return(exposureID())
  }else{
    return(idnr)
  }
}

partitionNrows <- function(dataframe){
  nCl <- ceiling(nrow(aoC) / 10)
  divs <- c()
  divs[1] <- 0; divs[11] <- nrow(aoC)
  for(x in 2:10){
    divs[x] <- nCl*x
  }
  return(divs)
}

#************************************************************************************
setWorkDir()
instPkg()

print("Loading available outcomes from MRBase...")
ao <- available_outcomes()
aoC <- ao
aoC <- aoC %>% filter(substr(population,1,4) == "Euro") %>%
  filter(!grepl('dvt', trait, ignore.case = TRUE)) %>%
  group_by(trait) %>%
  slice(which.max(sample_size)) %>%
  arrange(trait)

#Get vector with partitioned aoC rows
aR <- partitionNrows(aoC)

#Switch to test server
toggle_api("test")

#Chose outcome (by id)
id_nr <- exposureID()

#Initialize dat
dat <- data.frame()

#Extracting exposure data for the best-fitting DVT study
exposure_dat <- extract_instruments(outcomes="UKB-a:65")
clumped_exp_dat <- clump_data(exposure_dat)


for(x in 1:10){
  print(paste0("*** Current outcome chunk: ", x, " ***"))
  outChk <- paste0('outcome_chunk',x)
  assign(outChk, data.frame(stringsAsFactors = FALSE))
  datVar <- paste0('dat',x)
  
  print("Extracting outcome data...")
  print(outChk)
  i <- 0
  for(y in aoC$id[c((aR[x]+1):aR[x+1])]){
    print(paste0('Current study ID: ', y))
    #Progress bar
    i <- i + 1
    progress(i, max.value = length(aoC$id[c((aR[x]+1):aR[x+1])]), progress.bar = TRUE, char = "|",
             init = TRUE, console = TRUE, gui = FALSE)
    cat("\n")
    assign(outChk, rbind.fill(get(outChk), extract_outcome_data(clumped_exp_dat$SNP, c(y))))
    
  }
  
  print("Harmonising data...")
  assign(paste0('dat',x), harmonise_data(clumped_exp_dat, get(outChk))) #Harmonise each chunk
  
  dat <- rbind.fill(dat, get(datVar)) #Binding the dataframes into one large dataframe
  
  write.csv(get(datVar), file=paste(datVar, 'csv', sep='.')) #Save each chunk as a csv file
}

dat_names <- list.files("traitsToDVT/", pattern = "^dat")
dat <- rbindlist(lapply(dat_names, fread), fill = T)
exp_names <- read.xlsx("bi-directional.xlsx", sheet = 1)
dat <- dat[dat$outcome %in% exp_names$outcome,]

#Using the power.prune function (see docs)
print("Using power.prune function")
dat <- power_prune(dat, method = 1)

# Load the downloaded RData object. This loads the rf object
#load("rf.rdata")

#Perform MR and sort the dataframe by pval
print("Preforming MR - one with Wald Ratio/IVW,
      and one with everything")
mr_results_WrIvw <- mr(dat, method_list=c("mr_wald_ratio", "mr_ivw"))
mr_results_WrIvw <- mr_results_WrIvw %>% arrange(pval)
het <- mr_heterogeneity(dat, method_list=c("mr_egger_regression"))
ple <- mr_pleiotropy_test(dat)
res <- merge(mr_results_WrIvw, het, by = "outcome")
res <- merge(res, ple, by = "outcome")
res <- split_exposure(res)
res <- split_outcome(res)
res <- res[,c("outcome", "nsnp", "method.x", "b", "se.x", "pval.x", "Q_pval", "pval.y")]

#Export as xlsx table
write.xlsx(res,"bi-directional.xlsx", overwrite = T)

