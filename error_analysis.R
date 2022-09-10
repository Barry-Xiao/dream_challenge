

library(tidyverse)

metadata <- read.csv("data/combo_clean_data.csv", header = T)

metadata_preterm <- metadata %>% filter(collect_wk <= 32)

metadata_early_preterm <- metadata %>% filter(collect_wk <= 28)


##read error data

# error_data <- rbind(read.csv("test_error_analysis.csv", header = T)[,-1] %>% mutate(dataset = "test"),
#                     read.csv("train_val_error_analysis.csv", header = T)[,-1]%>% mutate(dataset = "train_val")) %>%
#               rename(specimen = participant_id)


error_data <- read.csv("test_error_analysis.csv", header = T)[,-1] %>% mutate(correct = labels == predicted_labels)

error_data_analysis_preterm <- merge(error_data, metadata_preterm %>% group_by(participant_id) %>% summarise(n = n()) %>%
                                       ungroup(),
                             by = "participant_id", all.x = T)

trupos <- error_data_analysis_preterm %>% filter(!correct & labels == 1)

