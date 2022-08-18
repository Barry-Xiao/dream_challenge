
#load the packages
library(tidyverse)

args = commandArgs(trailingOnly = TRUE)


##
outcome = args[1]

data_path = args[2]

meta_path = args[3]

profile_path = args[4]





##


data_filtered = read.csv(data_path, header = T)



## read metadata

meta_data <- read.csv("data/metadata_imputed.csv", header = T) %>%
                select(participant_id,specimen,collect_wk,project,
                       was_preterm, was_early_preterm,delivery_wk,
                       age_imp)

multi_data <- meta_data %>% group_by(participant_id, collect_wk) %>% 
                filter(n() > 1) %>% ungroup() %>% select(participant_id,specimen,collect_wk)



alpha_data <- read.csv("data/training_data_2022-07-21/alpha_diversity/alpha_diversity.csv", header = T) 

CST_data <- read.csv("data/training_data_2022-07-21/community_state_types/cst_valencia.csv", header = T) 


## aggregate multiple tests

CST_multi <- merge(multi_data, CST_data,by = "specimen", all.x = T) %>% 
                arrange(participant_id, collect_wk) %>% 
                  mutate(uid = paste(participant_id,collect_wk, sep = "_"))


CST_multi_sum <- CST_multi %>% group_by(uid) %>% 
  summarise(same_CST = n_distinct(CST) == 1) %>% ungroup()

inconsistent_uid <- CST_multi_sum %>% filter(!same_CST) %>% pull(uid)


combo_data <- Reduce(function(x,y){merge(x,y,by = "specimen")},
                     list(meta_data,alpha_data,CST_data %>% select(specimen,CST),data_filtered)) %>%
                  mutate(uid = paste(participant_id,collect_wk, sep = "_"))


combo_data_fix <- combo_data %>% select(participant_id,specimen,project,
                                        was_preterm,was_early_preterm,delivery_wk,
                                        age_imp,CST,uid) %>% 
                        filter(!uid %in% inconsistent_uid) %>% group_by(uid) %>%
                        filter(row_number() == 1) %>% ungroup()


combo_data_change <- combo_data %>% select(-c(participant_id,specimen,project,
                                              was_preterm,was_early_preterm,age_imp,
                                              delivery_wk,CST)) %>%
                    filter(!uid %in% inconsistent_uid) %>%
                    group_by(uid) %>% summarise_all("mean") %>% ungroup()
  
combo_data_clean <- merge(combo_data_fix,combo_data_change,by = "uid") %>% select(-uid)




if (outcome == "was_preterm") {
  
  combo_data_clean <- combo_data_clean %>% filter(collect_wk <= 32)
  
} else {
  
  combo_data_clean <- combo_data_clean %>% filter(collect_wk <= 28)
}


combo_data_meta <- combo_data_clean %>% select(c(participant_id:rooted_pd))


combo_data_profile <- combo_data_clean %>% select(-c(participant_id:rooted_pd))


if (!file.exists(meta_path)){
  write.csv(combo_data_meta, meta_path, row.names = F, na = "")
}


write.table(combo_data_profile, profile_path, sep = ',', row.names = F, col.names = F)




##check multiple measurement per week


 