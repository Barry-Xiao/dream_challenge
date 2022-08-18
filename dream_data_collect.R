
#load the packages
library(tidyverse)



## read metadata


meta_data <- read.csv("data/metadata_imputed.csv", header = T) %>%
                select(participant_id,specimen,collect_wk,project,
                       was_preterm, was_early_preterm,age_imp,
                       delivery_wk,collect_tri)

multi_data <- meta_data %>% group_by(participant_id, collect_wk) %>% 
                filter(n() > 1) %>% ungroup() %>% select(participant_id,specimen,collect_wk)



alpha_data <- read.csv("data/raw_data/training_data_2022-07-21/alpha_diversity/alpha_diversity.csv", header = T)

CST_data <- read.csv("data/raw_data/training_data_2022-07-21/community_state_types/cst_valencia.csv", header = T)






taxo_genus <- read.csv("data/raw_data/training_data_2022-07-21/taxonomy/taxonomy_relabd.genus.csv", header = T) %>%
                select(specimen, Lactobacillus)


## aggregate multiple tests

CST_multi <- merge(multi_data, CST_data,by = "specimen", all.x = T) %>% 
                arrange(participant_id, collect_wk) %>% 
                  mutate(uid = paste(participant_id,collect_wk, sep = "_"))


CST_multi_sum <- CST_multi %>% group_by(uid) %>% 
  summarise(same_CST = n_distinct(CST) == 1) %>% ungroup()

inconsistent_uid <- CST_multi_sum %>% filter(!same_CST) %>% pull(uid)




combo_data <- Reduce(function(x,y){merge(x,y,by = "specimen")},
                     list(meta_data,alpha_data,CST_data %>% select(specimen,CST),taxo_genus)) %>%
                  mutate(uid = paste(participant_id,collect_wk, sep = "_"))


combo_clean_data <- combo_data %>% filter(!uid %in% inconsistent_uid) %>%
                    group_by(uid) %>% mutate(shannon = mean(shannon),
                                             bwpd = mean(bwpd),
                                             Lactobacillus = mean(Lactobacillus)) %>% 
                    filter(row_number() == 1) %>% ungroup() %>% select(-uid)
                    

combo_inconsistent_data <- combo_data %>% filter(uid %in% inconsistent_uid)
                


write.csv(combo_clean_data, "combo_clean_data.csv", row.names = F, na = "")


##check multiple measurement per week


 