hodgkin_annotated <- read_csv('~/Hodgkin_github/data/hodgkin_DFCI_noID.csv')

hodgkin_annotated <- hodgkin_annotated %>% 
  mutate(ebv_status = factor(ebv_status, levels = c('Positive', 'Negative')))

###############################################################################################


test_df <- hodgkin_annotated %>% 
  select(pointNum, Annotation, ebv_status) %>% 
  group_by(pointNum, ebv_status, Annotation) %>% 
  dplyr::count() %>% 
  group_by(pointNum, ebv_status) %>% 
  mutate(total_n = sum(n)) %>% 
  mutate(prop = n/total_n) %>% 
  ungroup() %>% 
  select(pointNum, Annotation, ebv_status, prop)

test_result_df <- data.frame(Annotation = unique(test_df$Annotation),
                             Test = rep('EBV+ vs EBV-, Wilcoxon Rank Sum, Two-sided', 13),
                             p.value = rep(NA, 13))

for (type in test_result_df$Annotation){
  sub_df <- test_df %>% 
    filter(Annotation == type)
  test_result <- wilcox.test(sub_df$prop[sub_df$ebv_status == 'Positive'],
                                sub_df$prop[sub_df$ebv_status == 'Negative'])
  test_result_df$p.value[test_result_df$Annotation == type] <- test_result$p.value
}


test_result_df <- test_result_df %>% 
  mutate(rank = rank(p.value)) %>% 
  mutate(p.adj = (rank/13) * 0.05) %>% 
  mutate(reject = ifelse(p.value < 0.05, 1, 0),
         reject.adj = ifelse(p.value <= p.adj, 1, 0))

write_csv(test_result_df, '~/hodgkinebvmibi/paper_figures/figure2/figure2c_test.csv')
  




#################################################################################################

test_df <- hodgkin_annotated %>% 
  dplyr::select(pointNum, Annotation, ebv_status, all_of(marker_list)) %>% 
  pivot_longer(cols = -c('pointNum', 'Annotation', 'ebv_status'), names_to = 'marker', values_to = 'value') %>% 
  group_by(pointNum, Annotation, ebv_status, marker) %>% 
  summarise(marker_fov_avg = mean(value)) %>% 
  ungroup()


test_result_df <- data.frame(marker = rep(NA, 100),
                             celltype = rep(NA, 100),
                             p.value = rep(NA, 100))


k <- 0

for (i in unique(test_df$marker)){
  df <- test_df %>% 
    dplyr::filter(marker == i)
  if (i %in% c('B2-Microglobulin', 'HLA1')){
    for (j in unique(df$Annotation)){
      sub_df <- df %>% 
        filter(Annotation == j)
      test_result <- wilcox.test(sub_df$marker_fov_avg[sub_df$ebv_status == 'Positive'],
                                 sub_df$marker_fov_avg[sub_df$ebv_status == 'Negative'],
                                 alternative = 'greater')
      k <- k + 1
      
      test_result_df$marker[k] <- i
      test_result_df$celltype[k] <- j
      test_result_df$p.value[k] <- test_result$p.value
      
    }
  } else if (i == 'HLA-DR'){
    df <- df %>% 
      dplyr::filter(Annotation %in% c('B', 'DC', 'M1', 'M2', 'Tumor'))
    
    for (j in unique(df$Annotation)){
      sub_df <- df %>% 
        filter(Annotation == j)
      test_result <- wilcox.test(sub_df$marker_fov_avg[sub_df$ebv_status == 'Positive'],
                                 sub_df$marker_fov_avg[sub_df$ebv_status == 'Negative'],
                                 alternative = 'greater')
      k <- k + 1
      
      test_result_df$marker[k] <- i
      test_result_df$celltype[k] <- j
      test_result_df$p.value[k] <- test_result$p.value
    
    } 
  } else if (i %in% c('CD45RO', 'Tox', 'Lag3', 'PD-1')){
    df <- df %>% 
      dplyr::filter(Annotation %in% c('CD4', 'CD8', 'Cytotoxic CD4', 'Cytotoxic CD8', 'Treg'))
    
      for (j in unique(df$Annotation)){
        sub_df <- df %>% 
          filter(Annotation == j)
        test_result <- wilcox.test(sub_df$marker_fov_avg[sub_df$ebv_status == 'Positive'],
                                   sub_df$marker_fov_avg[sub_df$ebv_status == 'Negative'],
                                   alternative = 'greater')
        k <- k + 1
        
        test_result_df$marker[k] <- i
        test_result_df$celltype[k] <- j
        test_result_df$p.value[k] <- test_result$p.value
    } 
  } else if (i == 'PD-L1'){
    df <- df %>% 
      dplyr::filter(Annotation %in% c('DC', 'M1', 'M2', 'Tumor'))
    
      for (j in unique(df$Annotation)){
        sub_df <- df %>% 
          filter(Annotation == j)
        test_result <- wilcox.test(sub_df$marker_fov_avg[sub_df$ebv_status == 'Positive'],
                                   sub_df$marker_fov_avg[sub_df$ebv_status == 'Negative'],
                                   alternative = 'greater')
        k <- k + 1
        
        test_result_df$marker[k] <- i
        test_result_df$celltype[k] <- j
        test_result_df$p.value[k] <- test_result$p.value
      }
    
  }
}

test_result_df <- test_result_df %>% 
  drop_na() %>% 
  arrange(p.value) 

test_result_df <- test_result_df %>% 
  mutate(rank = rank(p.value)) %>% 
  mutate(p.adj = (rank/55) * 0.05) %>% 
  mutate(reject = ifelse(p.value < 0.05, 1, 0),
         reject.adj = ifelse(p.value <= p.adj, 1, 0))

write_csv(test_result_df, '~/hodgkinebvmibi/paper_figures/figure2/test.csv')



#################################################################################


CN_prop_test_df <- df_topic %>% 
  select(pointNum, topic, ebv_status) %>% 
  group_by(pointNum, ebv_status, topic) %>% 
  dplyr::count() %>% 
  group_by(pointNum, ebv_status) %>% 
  mutate(total_n = sum(n)) %>% 
  mutate(prop = n/total_n) %>% 
  ungroup() %>% 
  select(pointNum, topic, ebv_status, prop)

test_result_df <- data.frame(topic = unique(CN_prop_test_df$topic),
                             Test = rep('EBV+ vs EBV-, Wilcoxon Rank Sum, Two-sided', 8),
                             p.value = rep(NA, 8))

for (type in CN_prop_test_df$topic){
  sub_df <- CN_prop_test_df %>% 
    filter(topic == type)
  
  if (sum(sub_df$ebv_status == 'Positive') == 0 | sum(sub_df$ebv_status == 'Negative') == 0) {
    next
  }
  
  test_result <- wilcox.test(sub_df$prop[sub_df$ebv_status == 'Positive'],
                             sub_df$prop[sub_df$ebv_status == 'Negative'])
  test_result_df$p.value[test_result_df$topic == type] <- test_result$p.value
}


test_result_df <- test_result_df %>% 
  mutate(rank = rank(p.value)) %>% 
  mutate(p.adj = (rank/13) * 0.05) %>% 
  mutate(reject = ifelse(p.value < 0.05, 1, 0),
         reject.adj = ifelse(p.value <= p.adj, 1, 0))

write_csv(test_result_df, '~/hodgkinebvmibi/paper_figures/figure3/figure3c_test.csv')





###################################################################################





CN_test_df <- df_topic %>% 
  dplyr::select(pointNum, Annotation, ebv_status, all_of(marker_list), topic) %>% 
  pivot_longer(cols = -c('pointNum', 'Annotation', 'ebv_status', 'topic'), names_to = 'marker', values_to = 'value') %>% 
  group_by(pointNum, Annotation, ebv_status, marker, topic) %>% 
  summarise(marker_fov_avg = mean(value)) %>% 
  ungroup()


lm_test_df <- CN_test_df %>% 
  pivot_wider(names_from = marker, values_from = marker_fov_avg) %>% 
  filter(topic %in% c('Topic-0', 'Topic-1'))

fit1 <- lm(Tox ~ Annotation + ebv_status + topic + ebv_status * topic, data = lm_test_df)


summary(fit1)

wilcox.test(CN_test_df$marker_fov_avg[CN_test_df$marker == 'B2-Microglobulin' & CN_test_df$Annotation == 'CD4' & 
                                        CN_test_df$topic == 'Topic-0' & CN_test_df$ebv_status == 'Negative'],
            CN_test_df$marker_fov_avg[CN_test_df$marker == 'B2-Microglobulin' & CN_test_df$Annotation == 'CD4' & 
                                        CN_test_df$topic == 'Topic-1' & CN_test_df$ebv_status == 'Negative'])

CN_test_result_df <- data.frame(marker = rep(NA, 200),
                                celltype = rep(NA, 200),
                             test = rep(NA, 200),
                             p.value = rep(NA, 200))


k <- 0


for (i in unique(CN_test_df$marker)){
  df <- CN_test_df %>% 
    dplyr::filter(marker == i)
  if (i %in% c('CD45RO', 'Tox', 'Lag3', 'PD-1')){
    for (j in c('CD4', 'CD8', 'Cytotoxic CD4', 'Cytotoxic CD8', 'Treg')){
      sub_df <- df %>% 
        filter(Annotation == j)
      
      k <- k + 1
      
      test_result1 <- wilcox.test(sub_df$marker_fov_avg[sub_df$topic == 'Topic-0' & sub_df$ebv_status == 'Positive'],
                                 sub_df$marker_fov_avg[sub_df$topic == 'Topic-1' & sub_df$ebv_status == 'Positive'])
      
      CN_test_result_df$marker[k] <- i
      CN_test_result_df$celltype[k] <- j
      CN_test_result_df$test[k] <- paste0('Positive, Topic-0 vs. Topic-1')
      CN_test_result_df$p.value[k] <- test_result1$p.value
      
      k <- k + 1
      
      test_result2 <- wilcox.test(sub_df$marker_fov_avg[sub_df$topic == 'Topic-0' & sub_df$ebv_status == 'Negative'],
                                  sub_df$marker_fov_avg[sub_df$topic == 'Topic-1' & sub_df$ebv_status == 'Negative'])
      
      CN_test_result_df$marker[k] <- i
      CN_test_result_df$celltype[k] <- j
      CN_test_result_df$test[k] <- paste0('Negative, Topic-0 vs. Topic-1')
      CN_test_result_df$p.value[k] <- test_result2$p.value
      
    }
  } else if (i %in% c('HLA-DR', 'PD-L1')){
    df <- CN_test_df %>% 
      dplyr::filter(Annotation %in% c('DC', 'M1', 'M2', 'Tumor'))
    
    for (j in unique(df$Annotation)){
      sub_df <- df %>% 
        filter(Annotation == j)
      
      k <- k + 1
      
      test_result1 <- wilcox.test(sub_df$marker_fov_avg[sub_df$topic == 'Topic-0' & sub_df$ebv_status == 'Positive'],
                                  sub_df$marker_fov_avg[sub_df$topic == 'Topic-1' & sub_df$ebv_status == 'Positive'])
      
      CN_test_result_df$marker[k] <- i
      CN_test_result_df$celltype[k] <- j
      CN_test_result_df$test[k] <- paste0('Positive, Topic-0 vs. Topic-1')
      CN_test_result_df$p.value[k] <- test_result1$p.value
      
      k <- k + 1
      
      test_result2 <- wilcox.test(sub_df$marker_fov_avg[sub_df$topic == 'Topic-0' & sub_df$ebv_status == 'Negative'],
                                  sub_df$marker_fov_avg[sub_df$topic == 'Topic-1' & sub_df$ebv_status == 'Negative'])
      
      CN_test_result_df$marker[k] <- i
      CN_test_result_df$celltype[k] <- j
      CN_test_result_df$test[k] <- paste0('Negative, Topic-0 vs. Topic-1')
      CN_test_result_df$p.value[k] <- test_result2$p.value
      
    } 
  } else if (i %in% c('B2-Microglobulin', 'HLA1')){
      for (j in c('CD4', 'CD8', 'Cytotoxic CD4', 'Cytotoxic CD8', 'Treg', 'DC', 'M1', 'M2', 'Tumor')){
        sub_df <- df %>% 
          filter(Annotation == j)
        
        k <- k + 1
        
        test_result1 <- wilcox.test(sub_df$marker_fov_avg[sub_df$topic == 'Topic-0' & sub_df$ebv_status == 'Positive'],
                                    sub_df$marker_fov_avg[sub_df$topic == 'Topic-1' & sub_df$ebv_status == 'Positive'])
        
        CN_test_result_df$marker[k] <- i
        CN_test_result_df$celltype[k] <- j
        CN_test_result_df$test[k] <- paste0('Positive, Topic-0 vs. Topic-1')
        CN_test_result_df$p.value[k] <- test_result1$p.value
        
        k <- k + 1
        
        test_result2 <- wilcox.test(sub_df$marker_fov_avg[sub_df$topic == 'Topic-0' & sub_df$ebv_status == 'Negative'],
                                    sub_df$marker_fov_avg[sub_df$topic == 'Topic-1' & sub_df$ebv_status == 'Negative'])
        
        CN_test_result_df$marker[k] <- i
        CN_test_result_df$celltype[k] <- j
        CN_test_result_df$test[k] <- paste0('Negative, Topic-0 vs. Topic-1')
        CN_test_result_df$p.value[k] <- test_result2$p.value
        
      }
  }
}





CN_test_result_df <- CN_test_result_df %>% 
  drop_na() %>% 
  mutate(rank = rank(p.value, ties.method = "max")) %>% 
  mutate(stat = (rank/124) * 0.05) %>% 
  mutate(reject = ifelse(p.value < 0.05, 1, 0),
         reject.adj = ifelse(p.value <= stat, 1, 0))


write_csv(CN_test_result_df, '~/hodgkinebvmibi/paper_figures/figure3/fig3e_test.csv')



##########################################################################################################


CN_cellcount_df <- df_topic %>% 
  group_by(pointNum, topic, Annotation, ebv_status) %>% 
  dplyr::count() %>% 
  group_by(pointNum, topic, ebv_status) %>% 
  mutate(total_n = sum(n),
         prop = n/sum(n))

CN_count_test_result <- data.frame(celltype = rep(NA, 100),
                                   test = rep(NA, 100),
                                   p.value = rep(NA, 100))

k <- 0
for (i in c('Tumor', 'CD4', 'CD8', 'Cytotoxic CD4', 'Cytotoxic CD8', 'DC', 'M1', 'M2', 'Treg', 'B',
            'Endothelial', 'Neutrophil', 'NK')){
  
  sub_df <- CN_cellcount_df %>% 
    filter(Annotation == i) %>% 
    filter(topic %in% c('Topic-0', 'Topic-1'))
  
  test_result1 <- wilcox.test(sub_df$prop[sub_df$ebv_status == 'Positive' & sub_df$topic == 'Topic-0'],
                              sub_df$prop[sub_df$ebv_status == 'Negative' & sub_df$topic == 'Topic-0'])
  
  k <- k + 1
  
  CN_count_test_result$celltype[k] <- i
  CN_count_test_result$test[k] <- paste0('Topic-0, Positive vs. Negative')
  CN_count_test_result$p.value[k] <- test_result1$p.value
  
  test_result2 <- wilcox.test(sub_df$prop[sub_df$ebv_status == 'Positive' & sub_df$topic == 'Topic-1'],
                              sub_df$prop[sub_df$ebv_status == 'Negative' & sub_df$topic == 'Topic-1'])
  
  k <- k + 1
  
  CN_count_test_result$celltype[k] <- i
  CN_count_test_result$test[k] <- paste0('Topic-1, Positive vs. Negative')
  CN_count_test_result$p.value[k] <- test_result2$p.value
  
  
  test_result3 <- wilcox.test(sub_df$prop[sub_df$ebv_status == 'Positive' & sub_df$topic == 'Topic-0'],
                              sub_df$prop[sub_df$ebv_status == 'Positive' & sub_df$topic == 'Topic-1'])
  
  k <- k + 1
  
  CN_count_test_result$celltype[k] <- i
  CN_count_test_result$test[k] <- paste0('Positive, Topic-0 vs Topic-1')
  CN_count_test_result$p.value[k] <- test_result3$p.value
  
  test_result4 <- wilcox.test(sub_df$prop[sub_df$ebv_status == 'Negative' & sub_df$topic == 'Topic-0'],
                              sub_df$prop[sub_df$ebv_status == 'Negative' & sub_df$topic == 'Topic-1'])
  
  k <- k + 1
  
  CN_count_test_result$celltype[k] <- i
  CN_count_test_result$test[k] <- paste0('Negative, Topic-0 vs Topic-1')
  CN_count_test_result$p.value[k] <- test_result4$p.value
}

CN_count_test_result <- CN_count_test_result %>% 
  drop_na() %>% 
  mutate(rank = rank(p.value, ties.method = 'max')) %>% 
  mutate(stat = (rank/32) * 0.05) %>% 
  mutate(reject = ifelse(p.value < 0.05, 1, 0),
         reject.adj = ifelse(p.value <= stat, 1, 0))


write_csv(CN_count_test_result, '~/hodgkinebvmibi/paper_figures/figure3/fig3d_test.csv')



#########################################################################



CN_test_result_df_ebv <- data.frame(marker = rep(NA, 200),
                                celltype = rep(NA, 200),
                                test = rep(NA, 200),
                                p.value = rep(NA, 200))


k <- 0


for (i in unique(CN_test_df$marker)){
  df <- CN_test_df %>% 
    dplyr::filter(marker == i)
  if (i %in% c('CD45RO', 'Tox', 'Lag3', 'PD-1')){
    for (j in c('CD4', 'CD8', 'Cytotoxic CD4', 'Cytotoxic CD8', 'Treg')){
      sub_df <- df %>% 
        filter(Annotation == j)
      
      k <- k + 1
      
      test_result1 <- wilcox.test(sub_df$marker_fov_avg[sub_df$topic == 'Topic-0' & sub_df$ebv_status == 'Positive'],
                                  sub_df$marker_fov_avg[sub_df$topic == 'Topic-0' & sub_df$ebv_status == 'Negative'],
                                  alternative = 'greater')
      
      CN_test_result_df_ebv$marker[k] <- i
      CN_test_result_df_ebv$celltype[k] <- j
      CN_test_result_df_ebv$test[k] <- paste0('Topic-0, Positive vs. Negative')
      CN_test_result_df_ebv$p.value[k] <- test_result1$p.value
      
      k <- k + 1
      
      test_result2 <- wilcox.test(sub_df$marker_fov_avg[sub_df$topic == 'Topic-1' & sub_df$ebv_status == 'Positive'],
                                  sub_df$marker_fov_avg[sub_df$topic == 'Topic-1' & sub_df$ebv_status == 'Negative'],
                                  alternative = 'greater')
      
      CN_test_result_df_ebv$marker[k] <- i
      CN_test_result_df_ebv$celltype[k] <- j
      CN_test_result_df_ebv$test[k] <- paste0('Topic-1, Positive vs. Negative')
      CN_test_result_df_ebv$p.value[k] <- test_result2$p.value
      
    }
  } else if (i %in% c('HLA-DR', 'PD-L1')){
    df <- CN_test_df %>% 
      dplyr::filter(Annotation %in% c('DC', 'M1', 'M2', 'Tumor'))
    
    for (j in unique(df$Annotation)){
      sub_df <- df %>% 
        filter(Annotation == j)
      
      k <- k + 1
      
      test_result1 <- wilcox.test(sub_df$marker_fov_avg[sub_df$topic == 'Topic-0' & sub_df$ebv_status == 'Positive'],
                                  sub_df$marker_fov_avg[sub_df$topic == 'Topic-0' & sub_df$ebv_status == 'Negative'],
                                  alternative = 'greater')
      
      CN_test_result_df_ebv$marker[k] <- i
      CN_test_result_df_ebv$celltype[k] <- j
      CN_test_result_df_ebv$test[k] <- paste0('Topic-0, Positive vs. Negative')
      CN_test_result_df_ebv$p.value[k] <- test_result1$p.value
      
      k <- k + 1
      
      test_result2 <- wilcox.test(sub_df$marker_fov_avg[sub_df$topic == 'Topic-1' & sub_df$ebv_status == 'Positive'],
                                  sub_df$marker_fov_avg[sub_df$topic == 'Topic-1' & sub_df$ebv_status == 'Negative'],
                                  alternative = 'greater')
      
      CN_test_result_df_ebv$marker[k] <- i
      CN_test_result_df_ebv$celltype[k] <- j
      CN_test_result_df_ebv$test[k] <- paste0('Topic-1, Positive vs. Negative')
      CN_test_result_df_ebv$p.value[k] <- test_result2$p.value
      
    } 
  } else if (i %in% c('B2-Microglobulin', 'HLA1')){
    for (j in c('CD4', 'CD8', 'Cytotoxic CD4', 'Cytotoxic CD8', 'Treg', 'DC', 'M1', 'M2', 'Tumor')){
      sub_df <- df %>% 
        filter(Annotation == j)
      
      k <- k + 1
      
      test_result1 <- wilcox.test(sub_df$marker_fov_avg[sub_df$topic == 'Topic-0' & sub_df$ebv_status == 'Positive'],
                                  sub_df$marker_fov_avg[sub_df$topic == 'Topic-0' & sub_df$ebv_status == 'Negative'],
                                  alternative = 'greater')
      
      CN_test_result_df_ebv$marker[k] <- i
      CN_test_result_df_ebv$celltype[k] <- j
      CN_test_result_df_ebv$test[k] <- paste0('Topic-0, Positive vs. Negative')
      CN_test_result_df_ebv$p.value[k] <- test_result1$p.value
      
      k <- k + 1
      
      test_result2 <- wilcox.test(sub_df$marker_fov_avg[sub_df$topic == 'Topic-1' & sub_df$ebv_status == 'Positive'],
                                  sub_df$marker_fov_avg[sub_df$topic == 'Topic-1' & sub_df$ebv_status == 'Negative'],
                                  alternative = 'greater')
      
      CN_test_result_df_ebv$marker[k] <- i
      CN_test_result_df_ebv$celltype[k] <- j
      CN_test_result_df_ebv$test[k] <- paste0('Topic-1, Positive vs. Negative')
      CN_test_result_df_ebv$p.value[k] <- test_result2$p.value
      
    }
  }
}





CN_test_result_df_ebv <- CN_test_result_df_ebv %>% 
  drop_na() %>% 
  mutate(rank = rank(p.value, ties.method = "max")) %>% 
  mutate(stat = (rank/124) * 0.05) %>% 
  mutate(reject = ifelse(p.value < 0.05, 1, 0),
         reject.adj = ifelse(p.value <= stat, 1, 0))


write_csv(CN_test_result_df_ebv, '~/hodgkinebvmibi/paper_figures/figure3/fig3e_test_ebv.csv')





######################################################################################################


exhaustion_test_result_df <- data.frame(celltype = rep(NA, 200),
                                        test = rep(NA, 200),
                                        p.value = rep(NA, 200))


k <- 0


for (i in unique(exhaustion_df$Annotation)){
  df <- exhaustion_df %>% 
    dplyr::filter(Annotation == i) %>% 
    dplyr::filter(bin_deriva != 'Far')
  
  summary_data <- df %>% 
    dplyr::select(pointNum, bin_deriva, ebv_status, exhaustion_score_v3) %>% 
    group_by(pointNum, bin_deriva, ebv_status) %>% 
    summarise(score_avg = mean(exhaustion_score_v3)) %>% 
    ungroup() %>% 
    group_by(pointNum, ebv_status) %>% 
    mutate(diff = lag(score_avg) - score_avg) %>% 
    dplyr::select(pointNum, diff) %>% 
    drop_na()
  
  
  k <- k + 1
  
  test_result1 <- t.test(summary_data$diff[summary_data$ebv_status == 'Positive'], alternative = 'greater')
  
  exhaustion_test_result_df$celltype[k] <- i
  exhaustion_test_result_df$test[k] <- paste0(i, ', EBV+')
  exhaustion_test_result_df$p.value[k] <- test_result1$p.value
  
  k <- k + 1
  
  test_result2 <- t.test(summary_data$diff[summary_data$ebv_status == 'Negative'], alternative = 'greater')
  
  
  exhaustion_test_result_df$celltype[k] <- i
  exhaustion_test_result_df$test[k] <- paste0(i, ', EBV-')
  exhaustion_test_result_df$p.value[k] <- test_result2$p.value
}

exhaustion_test_result_df <- exhaustion_test_result_df %>% 
  drop_na() %>% 
  mutate(rank = rank(p.value, ties.method = 'max')) %>% 
  mutate(stat = (rank/124) * 0.05) %>% 
  mutate(reject = ifelse(p.value < 0.05, 1, 0),
         reject.adj = ifelse(p.value <= stat, 1, 0))



write_csv(exhaustion_test_result_df, '~/hodgkinebvmibi/paper_figures/figure4/fig4c_test.csv')












  
