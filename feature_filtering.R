##load packages
library(tidyverse)
library(ALDEx2)


##feature filtering based on wilcox test,
##only check features that were detected in more than 10% of total samples

###filename is the name of your taxonomy/phylotype file
###metadata is metadata file
###outcome: 'was_preterm', 'was_early_preterm'
###outname: output file name



feature_filtering <- function(filename,metadata,outcome,outname, prev = 0.1, seed = 1000){
  
  meta_data <- read.csv(metadata, header = T)
  outcome_vec <- meta_data[,outcome]
  
  raw <- read.csv(filename,header = T,row.names = 1)[meta_data$specimen,]
  
  
  reads <- as.data.frame(t(raw))
  
  
  reads_clean <- reads[rowSums(reads != 0) > prev*ncol(reads),]
  
  rm(raw)
  rm(reads)
  
  set.seed(seed)
  
  x <- aldex.clr(reads_clean, outcome_vec,
                 mc.samples = 128,
                 denom = "all", verbose = F)
  
  
  rm(reads_clean)
  
  x_tt <- aldex.ttest(x, paired.test = F, verbose = F)
  
  x_effect <- aldex.effect(x, CI = T, verbose = F)
  
  aldex_out <- data.frame(x_tt,x_effect)
  
  
  test <- rownames_to_column(aldex_out, "feature_selected") %>%
    filter(wi.eBH <= 0.05)  %>% 
    dplyr::select(feature_selected, we.eBH, wi.eBH, c(effect:overlap)) %>%
    dplyr::arrange(desc(abs(effect)))
  
  write.csv(test, outname, row.names = F, na = "")
  
}

feature_filtering("raw_data/training_data_2022-07-21/taxonomy/taxonomy_nreads.family.csv",
                  "metadata_imputed.csv",
                  "was_preterm",
                  "tax_fam_preterm.csv")
feature_filtering("raw_data/training_data_2022-07-21/taxonomy/taxonomy_nreads.genus.csv",
                  "metadata_imputed.csv",
                  "was_preterm",
                  "tax_genus_preterm.csv")
feature_filtering("raw_data/training_data_2022-07-21/taxonomy/taxonomy_nreads.species.csv",
                  "metadata_imputed.csv",
                  "was_preterm",
                  "tax_species_preterm.csv")

feature_filtering("raw_data/training_data_2022-07-21/phylotypes/phylotype_nreads.1e_1.csv",
                  "metadata_imputed.csv",
                  "was_preterm",
                  "phylo_.1_preterm.csv")
feature_filtering("raw_data/training_data_2022-07-21/phylotypes/phylotype_nreads.1e0.csv",
                  "metadata_imputed.csv",
                  "was_preterm",
                  "phylo_1_preterm.csv")
feature_filtering("raw_data/training_data_2022-07-21/phylotypes/phylotype_nreads.5e_1.csv",
                  "metadata_imputed.csv",
                  "was_preterm",
                  "phylo_.5_preterm.csv")


feature_filtering("raw_data/training_data_2022-07-21/taxonomy/taxonomy_nreads.family.csv",
                  "metadata_imputed.csv",
                  "was_early_preterm",
                  "tax_fam_early_preterm.csv")
feature_filtering("raw_data/training_data_2022-07-21/taxonomy/taxonomy_nreads.genus.csv",
                  "metadata_imputed.csv",
                  "was_early_preterm",
                  "tax_genus_early_preterm.csv")
feature_filtering("raw_data/training_data_2022-07-21/taxonomy/taxonomy_nreads.species.csv",
                  "metadata_imputed.csv",
                  "was_early_preterm",
                  "tax_species_early_preterm.csv")

feature_filtering("raw_data/training_data_2022-07-21/phylotypes/phylotype_nreads.1e_1.csv",
                  "metadata_imputed.csv",
                  "was_early_preterm",
                  "phylo_.1_early_preterm.csv")
feature_filtering("raw_data/training_data_2022-07-21/phylotypes/phylotype_nreads.1e0.csv",
                  "metadata_imputed.csv",
                  "was_early_preterm",
                  "phylo_1_early_preterm.csv")
feature_filtering("raw_data/training_data_2022-07-21/phylotypes/phylotype_nreads.5e_1.csv",
                  "metadata_imputed.csv",
                  "was_early_preterm",
                  "phylo_.5_early_preterm.csv")


