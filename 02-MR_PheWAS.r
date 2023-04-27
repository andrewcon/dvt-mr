#Hypothesis Free analysis of all instruments in MRBase to DVT
#This script can be modified for any desired outcome

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
      install.packages("dplyr")
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
  tts <- readline(prompt="Enter the trait you wish to study (e.g Body Mass Index) NOT CASE SENSITIVE: ")
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

outcomeID <- function(){
  idnr <- readline(prompt="Enter the id number for the outcome (e.g UKB-a:382 is a Waist circumference study) CASE SENSITIVE: ")
  n <- tolower(readline(prompt="Make sure you have no spelling errors. Proceed? (Y / N)"))
  if(!grepl("y",n) || nchar(n) > 1)
  {
    print("Going back...")
    return(outcomeID())
  }else{
    return(idnr)
  }
}

partitionNrows <- function(dataframe){
  nCl <- ceiling(nrow(aoC) / 6)
  divs <- c()
  divs[1] <- 0; divs[7] <- nrow(aoC)
  for(x in 2:6){
    divs[x] <- nCl*x
  }
  return(divs)
}

#************************************************************************************
setWorkDir()

instPkg()

#Switch to test server
toggle_api("test")

print("Loading available outcomes from MRBase...")
ao <- available_outcomes()
aoC <- curateAo(ao)

#Get vector with partitioned aoC rows
aR <- partitionNrows(aoC)

#Chose outcome (by id)
id_nr <- outcomeID()

#Initialize dat
dat <- data.frame()

for(x in 1:6){
  print(paste0("*** Current chunk: ", x, " ***"))
  chkVar <- paste0('exposure_chunk',x)
  datVar <- paste0('dat',x)
  
  progress(x, max.value = NULL, progress.bar = TRUE, char = "|",
           init = TRUE, console = TRUE, gui = FALSE)
  cat("\n")
  
  print("Extracting instruments...") #Read in exposure data
  exposure_dat <- extract_instruments(aoC$id[c((aR[x]+1):aR[x+1])])
  
  print("Clumping exposure data...") #Prune SNPs for LD 
  assign(chkVar, clump_data(exposure_dat)) #Create each chunk
  
  print("Harmonising data...")
  outcome_dat <- extract_outcome_data(get(chkVar)$SNP, c(id_nr))
  assign(paste0('dat',x), harmonise_data(get(chkVar), outcome_dat)) #Harmonise each chunk
  dat <- rbind.fill(dat, get(datVar)) #Binding the dataframes into one large dataframe
  write.csv(get(datVar), file=paste(datVar, 'csv', sep='.')) #Save each chunk as a csv file
}

#Using the power.prune function (see docs)
print("Using power.prune function")
dat <- power.prune(dat, method.size=T)

#Perform MR and sort the dataframe by pval
print("Preforming MR - one with Wald Ratio/IVW,
      and one with everything")
mr_results_WrIvw <- mr(dat, method_list=c("mr_wald_ratio", "mr_ivw"))
mr_results_other <- mr(dat)
mr_results_other <- mr_results_other %>% arrange(pval) %>% group_by(method)
mr_results_WrIvw <- mr_results_WrIvw %>% arrange(pval)

#Export as xlsx table
write.xlsx(mr_results_other,"./HypothesisFreeMR_DVT_all.xlsx")
write.xlsx(mr_results_WrIvw,"./HypothesisFreeMR_DVT_WrIvw.xlsx")
