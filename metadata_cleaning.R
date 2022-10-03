
#load the packages
library(tidyverse)
library(mice)

#unzip the downloaded data to raw_data folder

untar("training_data_2022-08-25.tgz", exdir = "raw_data")


#metadata cleaning

##import raw metadata

metadata_sample_raw <- read.csv("raw_data/training_data_2022-07-21/metadata/metadata.csv", header = T) 

##initial clean metadata
###generate `was_preterm` and `was_early_preterm`
###generate `collect_tri`: trimester of sample collection
###remove .0 in age

metadata_sample <- metadata_sample_raw  %>%
  mutate(was_preterm = as.logical(was_preterm),
         was_early_preterm = as.logical(was_early_preterm),
         collect_tri = factor(case_when(collect_wk <= 13 ~ 1,
                                        collect_wk > 13 &
                                          collect_wk <= 26 ~ 2,
                                        collect_wk > 26 ~ 3)),
         NIH.Racial.Category = 
           factor(case_when(NIH.Racial.Category == "American Indian or Alaska Native" ~ 1,
                            NIH.Racial.Category == "Asian" ~ 2,
                            NIH.Racial.Category == "Black or African American" ~ 3,
                            NIH.Racial.Category == "Native Hawaiian or Other Pacific Islander" ~ 4,
                            NIH.Racial.Category == "White" ~ 5), 
                  labels = c("American Indian or Alaska Native",
                             "Asian",
                             "Black or African American",
                             "Native Hawaiian or Other Pacific Islander",
                             "White")),
         age_tmp = str_replace_all(age, "\\..*","")) %>%
  arrange(project,participant_id,specimen)


##take demo/outcome data from first row as meta_ind

metadata_ind <- metadata_sample %>% group_by(participant_id) %>%
  filter(row_number() == 1) %>% ungroup() %>% 
  select(-c(specimen,collect_wk,collect_tri,race,NIH.Ethnicity.Category)) %>%
  mutate(age_imp = as.numeric(age_tmp))


##impute age


###set seed

set.seed(41072202)

###impute categorical

age_imp_Below_18 <- metadata_ind %>% filter((is.na(age_imp) & age == "Below_18")|
                                         (!is.na(age_imp) & age < 18)) %>% 
              select(age_imp,participant_id)
  
imp_Below_18 <- mice(age_imp_Below_18, m = 3, method = "sample", printFlag = F)

age_imp_Below_18$age_imp <- round(rowMeans(cbind(complete(imp_Below_18,1)$age_imp,
                             complete(imp_Below_18,2)$age_imp,
                             complete(imp_Below_18,3)$age_imp)))

age_imp_18_to_28 <- metadata_ind %>% filter((is.na(age_imp) & age == "18_to_28")|
                                              (!is.na(age_imp) & age >= 18 & age <=28)) %>% 
  select(age_imp,participant_id)

imp_18_to_28 <- mice(age_imp_18_to_28, m = 3, method = "sample", printFlag = F)

age_imp_18_to_28$age_imp <- round(rowMeans(cbind(complete(imp_18_to_28,1)$age_imp,
                                                 complete(imp_18_to_28,2)$age_imp,
                                                 complete(imp_18_to_28,3)$age_imp)))


age_imp_29_38 <- metadata_ind %>% filter((is.na(age_imp) & age == "29-38")|
                                              (!is.na(age_imp) & age >= 29 & age <=38)) %>% 
  select(age_imp,participant_id)

imp_29_38 <- mice(age_imp_29_38, m = 3, method = "sample", printFlag = F)

age_imp_29_38$age_imp <- round(rowMeans(cbind(complete(imp_29_38,1)$age_imp,
                                                 complete(imp_29_38,2)$age_imp,
                                                 complete(imp_29_38,3)$age_imp)))

age_imp_Above_38 <- metadata_ind %>% filter((is.na(age_imp) & age == "Above_38")|
                                           (!is.na(age_imp) & age > 38)) %>% 
  select(age_imp,participant_id)

imp_Above_38 <- mice(age_imp_Above_38, m = 3, method = "sample", printFlag = F)

age_imp_Above_38$age_imp <- round(rowMeans(cbind(complete(imp_Above_38,1)$age_imp,
                                              complete(imp_Above_38,2)$age_imp,
                                              complete(imp_Above_38,3)$age_imp)))



age_imp_cat <- rbind(age_imp_Below_18,age_imp_18_to_28,age_imp_29_38,age_imp_Above_38) %>%
                  rename(age_imp2 = age_imp)

##impute rest of unknow age

age_imp_all <- merge(metadata_ind %>% 
                       select(age_imp,participant_id),
                     age_imp_cat, by = "participant_id", all.x = T) %>%
               mutate(age_imp3 = ifelse(is.na(age_imp),age_imp2,age_imp)) %>%
                select(participant_id,age_imp3)



imp_all <- mice(age_imp_all, m = 3, method = "sample", printFlag = F)

age_imp_all$age_imp3<- round(rowMeans(cbind(complete(imp_all,1)$age_imp3,
                                                 complete(imp_all,2)$age_imp3,
                                                 complete(imp_all,3)$age_imp3)))




##impute race

metadata_race_imp <- metadata_ind %>% select(participant_id,NIH.Racial.Category) %>% 
                                mutate(race_imp = as.numeric(NIH.Racial.Category)) %>%
                                select(-NIH.Racial.Category)


imp_all_race <-  mice(metadata_race_imp, m = 3, method = "sample", printFlag = F)

metadata_race_imp$race_imp<- round(rowMeans(cbind(complete(imp_all_race,1)$race_imp,
                                            complete(imp_all_race ,2)$race_imp,
                                            complete(imp_all_race ,3)$race_imp)))


metadata_sample_imp <- reduce(list(metadata_sample,age_imp_all,
                                   metadata_race_imp), function(x,y){merge(x,y,by = "participant_id",all.x = T)}) %>%
                        rename(age_imp = age_imp3) %>%
                        mutate(race_imp = factor(race_imp, 
                                                 labels = c("American Indian or Alaska Native",
                                                            "Asian","Black or African American",
                                                             "Native Hawaiian or Other Pacific Islander",
                                                             "White")),
                               age_imp_cat = factor(case_when(age_imp < 18 ~ 1,
                                                              age_imp >= 18 & 
                                                                age_imp <= 28 ~ 2,
                                                              age_imp >= 29 &
                                                                age_imp <= 38~ 3,
                                                              age_imp > 38 ~ 4),
                                                    labels = c("below_18","from_18_to_28","from_29_to_38","above_38"))) %>%
                          select(-age_tmp)


write.csv(metadata_sample_imp,"metadata_imputed.csv", row.names = F, na = "")




library(reshape2)

long_sv <- read.csv("raw_data/training_data_2022-07-21/sv_counts/sp_sv_long.csv", header = T)


wide_sv_nreads <- dcast(long_sv %>% select(specimen,sv,nreads),
                        specimen ~ sv, value.var = "nreads")

write.csv(wide_sv_nreads, "sp_sv_wide_nreads.csv",
          row.names = F, na = "0")

wide_sv_fract <- dcast(long_sv %>% select(specimen,sv,fract),
                       specimen ~ sv, value.var = "fract")

write.csv(wide_sv_fract, "sp_sv_wide_fract.csv",
          row.names = F,na = "0")






