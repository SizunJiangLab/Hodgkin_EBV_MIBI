library(tidyverse)


# Load data
hodgkin_annotated <- read_csv('~/hodgkinebvmibi/clustering_data/corrected_annotation_032423.csv') 

# Create EBV class
hodgkin_annotated <- hodgkin_annotated %>% 
  mutate(ebv_status = case_when(patientID %in% c('2_rLN_ctrl', '1_rLN_ctrl') ~ 'control',
                                patientID %in% unique(hodgkin_annotated$patientID)[c(4,6,7,8,9,10,11,12,13,14,15,16,17,18)] ~ 'Negative',
                                patientID %in% unique(hodgkin_annotated$patientID)[c(1,2,3,5,21,22)] ~ 'Positive'))

hodgkin_annotated <- hodgkin_annotated %>% 
  mutate(ebv_status = factor(ebv_status, levels = c('Positive', 'Negative', 'control')))

# Filter out Other
hodgkin_annotated <- hodgkin_annotated %>% 
  dplyr::filter(Annotation != 'Other') 

# Filter out control

hodgkin_annotated <- hodgkin_annotated %>% 
  dplyr::filter(ebv_status != 'control')

tumor_effect_score_df <- read_csv('/mnt/nfs/home/huayingqiu/hodgkinebvmibi/anchor_data/tumor_effect_score_100um.csv')

hodgkin_tumor_score_100um <- tumor_effect_score_df %>% 
  left_join(hodgkin_annotated, by = c('pointNum', 'cellLabel', 'centroidX', 'centroidY'))


hodgkin_tumor_score_100um_no_0 <- hodgkin_tumor_score_100um %>% 
  filter(score != 0)

density_est <- density(hodgkin_tumor_score_100um_no_0$score)


deriva <- diff(density_est$y)/diff(density_est$x)


bin_deriva_threshold <- density_est$x[which.min(deriva)]

hodgkin_tumor_score_100um <- hodgkin_tumor_score_100um %>% 
  mutate(bin_deriva = ifelse(score == 0, 'Far', ifelse(score <= bin_deriva_threshold, 'Tumor Sparse', 'Tumor Dense')))

bin_score_df <- hodgkin_tumor_score_100um %>% 
  dplyr::select(pointNum, cellLabel, centroidX, centroidY, score, bin_deriva)


hodgkin_annotated <- hodgkin_annotated %>% 
  left_join(bin_score_df, by = c('pointNum', 'cellLabel', 'centroidX', 'centroidY'))


hodgkin_annotated <- hodgkin_annotated %>% 
  mutate(bin_deriva = ifelse(is.na(bin_deriva), 'Tumor', bin_deriva))

hodgkin_annotated <- hodgkin_annotated %>% 
  dplyr::select(-patientID, -tissue_block)

write_csv(hodgkin_annotated, '~/Hodgkin_github/data/hodgkin_DFCI_noID.csv')


####################################################################################


df_topic <- read_csv('~/hodgkinebvmibi/clustering_data/final_df_with_topic_032123.csv')

df_topic <- df_topic %>% 
  dplyr::select(-patientID, -tissue_block)

write_csv(df_topic, '~/Hodgkin_github/data/df_toipic_noID.csv')



