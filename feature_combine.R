##load packages
library(tidyverse)


feature_combine <- function(filename,metadata,CST,outfile){
  
  meta_data <- read.csv(metadata, header = T)
  
  
  multi_data <- meta_data %>% group_by(participant_id, collect_wk) %>% 
    filter(n() > 1) %>% ungroup() %>% dplyr::select(participant_id,specimen,collect_wk)
  
  
  CST_data <- read.csv(CST, header = T)
  
  CST_multi <- merge(multi_data, CST_data,by = "specimen", all.x = T) %>% 
    arrange(participant_id, collect_wk) %>% 
    mutate(uid = paste(participant_id,collect_wk,sep = "_"))
  
  
  CST_multi_sum <- CST_multi %>% group_by(uid) %>% 
    summarise(same_CST = n_distinct(CST) == 1) %>% ungroup()
  
  inconsistent_uid <- CST_multi_sum %>% filter(!same_CST) %>% pull(uid)
  
  
  raw <- read.csv(filename,header = T)
  
  combo_data_filter_inconsistent <- merge(meta_data,raw,by = "specimen",all.x = T) %>%
    mutate(uid = paste(participant_id,collect_wk, sep = "_")) %>% 
    filter(!uid %in% inconsistent_uid)
  
  meta_clean <- combo_data_filter_inconsistent %>% dplyr::select(c(colnames(meta_data),"uid")) %>% 
    arrange(specimen) %>% group_by(uid) %>% filter(row_number() == 1) %>% ungroup() %>% 
    arrange(uid) %>% select(uid,specimen)
  
  raw_clean <- combo_data_filter_inconsistent %>% dplyr::select(c(colnames(raw),"uid")) %>%
    arrange(specimen) %>% column_to_rownames("specimen") %>% 
    group_by(uid) %>% summarise_all("mean") %>% ungroup() %>% arrange(uid)
  
  
  output_file <- merge(raw_clean, meta_clean, by = "uid", all.x = T) %>%
                  select(-uid) %>% relocate(specimen, .before = 1)
  
  write.csv(output_file,outfile,row.names = F, na = "")

}


feature_combine("data/training_data_2022-07-21/taxonomy/taxonomy_relabd.family.csv",
                "data/metadata_imputed.csv",
                "data/training_data_2022-07-21/community_state_types/cst_valencia.csv",
                "data/raw_combine/comb_taxonomy_relabd.family.csv")

feature_combine("data/training_data_2022-07-21/taxonomy/taxonomy_relabd.genus.csv",
                "data/metadata_imputed.csv",
                "data/training_data_2022-07-21/community_state_types/cst_valencia.csv",
                "data/raw_combine/comb_taxonomy_relabd.genus.csv")

feature_combine("data/training_data_2022-07-21/taxonomy/taxonomy_relabd.species.csv",
               "data/metadata_imputed.csv",
               "data/training_data_2022-07-21/community_state_types/cst_valencia.csv",
               "data/raw_combine/comb_taxonomy_relabd.species.csv")

feature_combine("data/training_data_2022-07-21/phylotypes/phylotype_relabd.1e_1.csv",
                "data/metadata_imputed.csv",
                "data/training_data_2022-07-21/community_state_types/cst_valencia.csv",
                "data/raw_combine/comb_phylotype_relabd.1e_1.csv")

feature_combine("data/training_data_2022-07-21/phylotypes/phylotype_relabd.1e0.csv",
                "data/metadata_imputed.csv",
                "data/training_data_2022-07-21/community_state_types/cst_valencia.csv",
                "data/raw_combine/comb_phylotype_relabd.1e0.csv")

feature_combine("data/training_data_2022-07-21/phylotypes/phylotype_relabd.5e_1.csv",
                "data/metadata_imputed.csv",
                "data/training_data_2022-07-21/community_state_types/cst_valencia.csv",
                "data/raw_combine/comb_phylotype_relabd.5e_1.csv")






##
meta_alpha_combine <- function(metadata,CST,alpha,outfile){
  
  meta_data <- read.csv(metadata, header = T)
  
  
  multi_data <- meta_data %>% group_by(participant_id, collect_wk) %>% 
    filter(n() > 1) %>% ungroup() %>% dplyr::select(participant_id,specimen,collect_wk)
  
  
  CST_data <- read.csv(CST, header = T)
  
  CST_multi <- merge(multi_data, CST_data,by = "specimen", all.x = T) %>% 
    arrange(participant_id, collect_wk) %>% 
    mutate(uid = paste(participant_id,collect_wk,sep = "_"))
  
  
  CST_multi_sum <- CST_multi %>% group_by(uid) %>% 
    summarise(same_CST = n_distinct(CST) == 1) %>% ungroup()
  
  inconsistent_uid <- CST_multi_sum %>% filter(!same_CST) %>% pull(uid)
  
  
  alpha <- read.csv(alpha,header = T)
  
  combo_data_filter_inconsistent <- merge(meta_data,alpha,by = "specimen",all.x = T) %>%
    mutate(uid = paste(participant_id,collect_wk, sep = "_")) %>% 
    filter(!uid %in% inconsistent_uid)
  
  meta_clean <- merge(combo_data_filter_inconsistent %>% dplyr::select(c(colnames(meta_data),"uid")) %>% 
    arrange(specimen) %>% group_by(uid) %>% filter(row_number() == 1) %>% ungroup() %>% 
    arrange(uid),CST_data %>% select(CST,specimen),by = "specimen",all.x = T)
  
  raw_clean <- combo_data_filter_inconsistent %>% dplyr::select(c(colnames(alpha),"uid")) %>%
    arrange(specimen) %>% column_to_rownames("specimen") %>% 
    group_by(uid) %>% summarise_all("mean") %>% ungroup() %>% arrange(uid)
  
  
  output_file <- merge(raw_clean, meta_clean, by = "uid", all.x = T) %>%
    select(specimen,participant_id,project,was_preterm,was_early_preterm,
           delivery_wk,collect_wk,CST,shannon,inv_simpson,bwpd,phylo_entropy,quadratic,
           unrooted_pd,rooted_pd)
           
  
  write.csv(output_file,outfile,row.names = F, na = "")
  
}

meta_alpha_combine("data/metadata_imputed.csv",
                "data/training_data_2022-07-21/community_state_types/cst_valencia.csv",
                "data/training_data_2022-07-21/alpha_diversity/alpha_diversity.csv",
                "data/raw_combine/comb_meta_alpha.csv")
