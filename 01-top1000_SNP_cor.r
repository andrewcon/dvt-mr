#Creates a correlation matrix to be analysed using PhenoSpD

#Set the working directory
setwd("#")

#Read dataframe which contains a correlation table between traits using the top 1000 SNPs available
read_df <- read.delim("top_1000_gwas_hits_limit_10000.txt.cor.flat", header = FALSE)
colnames(read_df) <- c("trait.1", "trait.2", "cor")
read_df_bkup <- read_df

#mr is a dataframe with the MR results 
mr_exposures <- strsplit(gsub("[|]{1} id:", "", mr$exposure), "\\s{1}\\|")
mr_exposures <- sapply(mr_exposures,function(x) paste(x[2], x[1], sep=':'))

read_df_curated <- data.frame()
i <- 0
for(x in mr_exposures){
  read_df_curated <- rbind.fill(read_df_curated, read_df %>% filter(tolower(x) == tolower(trait.1) | 
                                                                      tolower(x) == tolower(trait.2)))
  progress(i, max.value = length(mr_exposures), progress.bar = TRUE, char = "|",
            console = TRUE, gui = FALSE)
  i <- i + 1
}

read_df_curated <- read_df_curated %>% filter(cor >= 0.0) %>% 
                                        distinct(trait.1, trait.2, .keep_all = TRUE)


read_df_curated <- arrange(read_df_curated, gsub('.*[:]','', read_df_curated$trait.1), 
                           gsub('.*[:]','', read_df_curated$trait.2))

mr_exposures_df <- data.frame(mr_exposures)
mr_exposures_df$mr_exposures <- as.character(mr_exposures_df$mr_exposures)


missing_val <- !mr_exposures_df$mr_exposures %in% read_df_curated$trait.1 & 
  !mr_exposures_df$mr_exposures %in% read_df_curated$trait.2

missing_val <- mr_exposures_df[missing_val,]
missing_val <- data.frame(missing_val)
colnames(missing_val) <- c("trait.1")
missing_val$trait.2 <- c(rep("", length(missing_val$trait.1)))
missing_val$cor <- c(rep(0, length(missing_val$trait.1)))
missing_val[] <- lapply(missing_val, as.character)


links <- read_df_curated
#links <- rbind.fill(links, missing_val)

#For PhenoSpd
interdf <- read_df_curated[, c("trait.2", "trait.1", "cor")]
colnames(interdf) <- c("trait.1", "trait.2", "cor")
interdf <- rbind(read_df_curated, interdf)

interdf2 <- interdf
interdf2[, "trait.2"] <- interdf2[, "trait.1"]
interdf2$cor <- 1.0
missing_val <- unlist(missing_val)
interdf22 <- data.frame(missing_val, missing_val, cor = 1, stringsAsFactors = FALSE)
colnames(interdf22) <- colnames(interdf2)
interdf2 <- rbind(interdf2, interdf22)
interdf2 <- unique(interdf2)

interdf3 <- rbind(interdf, interdf2)

cor_matrix <- acast(interdf3, trait.1~trait.2, value.var="cor")

write.csv(cor_matrix, file = "cor_matrix.csv")

