#### Calculate Metrics ####
cl <- multidplyr::new_cluster(workers)
# Yeast
roc_yeast_mix_baldur_results <- create_roc('err', yeast_mix_baldur_results, 'YEAST', cl) %>%
  mutate(
    method = 'LGMR-Baldur EB'
  )
roc_yeast_sin_baldur_results <- create_roc('err', yeast_sin_baldur_results, 'YEAST', cl) %>%
  mutate(
    method = 'GR-Baldur EB'
  )

roc_yeast_mix_baldur_wi_results <- create_roc('err', yeast_mix_baldur_wi_results, 'YEAST', cl) %>%
  mutate(
    method = 'LGMR-Baldur WI'
  )
roc_yeast_sin_baldur_wi_results <- create_roc('err', yeast_sin_baldur_wi_results, 'YEAST', cl) %>%
  mutate(
    method = 'GR-Baldur WI'
  )


roc_yeast_trend <- create_roc('p_val', yeast_trend, 'YEAST', cl) %>%
  mutate(
    method = 'Limma-Trend',
    comparison = 'ng50 vs ng100'
  )
roc_yeast_ttest <- create_roc('p_val', yeast_ttest, 'YEAST', cl) %>%
  mutate(
    method = 't-test',
    comparison = 'ng50 vs ng100'
  )

# UPS
roc_ups_mix_baldur_results <- create_roc('err', ups_mix_baldur_results, 'UPS', cl) %>%
  mutate(
    method = 'LGMR-Baldur EB'
  )
roc_ups_sin_baldur_results <- create_roc('err', ups_sin_baldur_results, 'UPS', cl) %>%
  mutate(
    method = 'GR-Baldur EB'
  )

roc_ups_mix_baldur_wi_results <- create_roc('err', ups_mix_baldur_wi_results, 'UPS', cl) %>%
  mutate(
    method = 'LGMR-Baldur WI'
  )
roc_ups_sin_baldur_wi_results <- create_roc('err', ups_sin_baldur_wi_results, 'UPS', cl) %>%
  mutate(
    method = 'GR-Baldur WI'
  )


roc_ups_trend <- create_roc('p_val', ups_trend, 'UPS', cl) %>%
  mutate(
    method = 'Limma-Trend'
  )
roc_ups_ttest <- create_roc('p_val', ups_ttest, 'UPS', cl) %>%
  mutate(
    method = 't-test'
  )

# Ramus
roc_ramus_mix_baldur_results <- create_roc('err', ramus_mix_baldur_results, 'UPS', cl) %>%
  mutate(
    method = 'LGMR-Baldur EB'
  )
roc_ramus_sin_baldur_results <- create_roc('err', ramus_sin_baldur_results, 'UPS', cl) %>%
  mutate(
    method = 'GR-Baldur EB'
  )

roc_ramus_mix_baldur_wi_results <- create_roc('err', ramus_mix_baldur_wi_results, 'UPS', cl) %>%
  mutate(
    method = 'LGMR-Baldur WI'
  )
roc_ramus_sin_baldur_wi_results <- create_roc('err', ramus_sin_baldur_wi_results, 'UPS', cl) %>%
  mutate(
    method = 'GR-Baldur WI'
  )

roc_ramus_trend <- create_roc('p_val', ramus_trend, 'UPS', cl) %>%
  mutate(
    method = 'Limma-Trend'
  )
roc_ramus_ttest <- create_roc('p_val', ramus_ttest, 'UPS', cl) %>%
  mutate(
    method = 't-test'
  )

# Human
roc_human_mix_baldur_results <- create_roc('err', human_mix_baldur_results, 'ECOLI', cl) %>%
  mutate(
    method = 'LGMR-Baldur EB'
  )
roc_human_sin_baldur_results <- create_roc('err', human_sin_baldur_results, 'ECOLI', cl) %>%
  mutate(
    method = 'GR-Baldur EB'
  )

roc_human_mix_wi_baldur_wi_results <- create_roc('err', human_mix_baldur_wi_results, 'ECOLI', cl) %>%
  mutate(
    method = 'LGMR-Baldur WI'
  )
roc_human_sin_wi_baldur_wi_results <- create_roc('err', human_sin_baldur_wi_results, 'ECOLI', cl) %>%
  mutate(
    method = 'GR-Baldur WI'
  )


roc_human_trend <- create_roc('p_val', human_trend, 'ECOLI', cl) %>%
  mutate(
    method = 'Limma-Trend'
  )
roc_human_ttest <- create_roc('p_val', human_ttest, 'ECOLI', cl) %>%
  mutate(
    method = 't-test'
  )

# Bruderer
roc_bruder_mix_baldur_results <- create_roc('err', bruder_results$bruder_mix_baldur_results, 'TRUE', cl) %>%
  mutate(
    method = 'LGMR-Baldur EB'
  )
roc_bruder_sin_baldur_results <- create_roc('err', bruder_results$bruder_sin_baldur_results, 'TRUE', cl) %>%
  mutate(
    method = 'GR-Baldur EB'
  )

roc_bruder_mix_baldur_wi_results <- create_roc('err', bruder_results$bruder_mix_baldur_wi_results, 'TRUE', cl) %>%
  mutate(
    method = 'LGMR-Baldur WI'
  )
roc_bruder_sin_baldur_wi_results <- create_roc('err', bruder_results$bruder_sin_baldur_wi_results, 'TRUE', cl) %>%
  mutate(
    method = 'GR-Baldur WI'
  )

roc_bruder_trend <- create_roc('p_val', bruder_results$bruder_trend, 'TRUE', cl) %>%
  mutate(
    method = 'Limma-Trend'
  )
roc_bruder_ttest <- create_roc('p_val', bruder_results$bruder_ttest, 'TRUE', cl) %>%
  mutate(
    method = 't-test'
  )

# Navarro
roc_nav_mix_baldur_results <- create_roc('err', nav_mix_baldur_results, '^((?!HUMAN).)*$', cl) %>%
  mutate(
    method = 'LGMR-Baldur EB'
  )
roc_nav_sin_baldur_results <- create_roc('err', nav_sin_baldur_results, '^((?!HUMAN).)*$', cl) %>%
  mutate(
    method = 'GR-Baldur EB'
  )

roc_nav_mix_baldur_wi_results <- create_roc('err', nav_mix_baldur_wi_results, '^((?!HUMAN).)*$', cl) %>%
  mutate(
    method = 'LGMR-Baldur WI'
  )
roc_nav_sin_baldur_wi_results <- create_roc('err', nav_sin_baldur_wi_results, '^((?!HUMAN).)*$', cl) %>%
  mutate(
    method = 'GR-Baldur WI'
  )

roc_nav_trend <- create_roc('p_val', nav_trend, '^((?!HUMAN).)*$', cl) %>%
  mutate(
    method = 'Limma-Trend',
    comparison = 'condi1 vs condi2'
  )
roc_nav_ttest <- create_roc('p_val', nav_ttest, '^((?!HUMAN).)*$', cl) %>%
  mutate(
    method = 't-test',
    comparison = 'condi1 vs condi2'
  )

rm(cl)
gc()


 #### Integrate Curves ####
# yeast
yeast_roc <- ls(pattern = 'roc_yeast.*') %>%
  map(get) %>%
  bind_rows()

yeast_auroc <- yeast_roc %>%
  group_by(comparison, method) %>%
  summarise(
    auROC = integrater(FPR, TPR)
  ) %>%
  arrange(desc(auROC)) %>%
  mutate(
    FPR = seq(0.1, .98, length.out = n()),
    TPR = 0
  )

# UPS
ups_roc <- ls(pattern = 'roc_ups.*') %>%
  map(get) %>%
  bind_rows() %>%
  mutate(
    comparison = str_replace(comparison, 'fmol100 and fmol25', 'fmol25 and fmol100')
  )

ups_auroc <- ups_roc %>%
  group_by(comparison, method) %>%
  summarise(
    auROC = integrater(FPR, TPR)
  ) %>%
  arrange(desc(auROC)) %>%
  mutate(
    FPR = seq(0.1, .98, length.out = n()),
    TPR = 0
  )


# Ramus
ramus_roc <- ls(pattern = 'roc_ramus.*') %>%
  map(get) %>%
  bind_rows() %>%
  mutate(
    comparison = str_replace_all(comparison, ramus_str_replace)
  )

ramus_auroc <- ramus_roc %>%
  group_by(comparison, method) %>%
  summarise(
    auROC = integrater(FPR, TPR)
  ) %>%
  arrange(desc(auROC)) %>%
  mutate(
    FPR = seq(0.1, .98, length.out = n()),
    TPR = 0
  )

# Human
human_roc <- ls(pattern = 'roc_human.*') %>%
  map(get) %>%
  bind_rows() %>%
  mutate(
    comparison = str_replace_all(comparison, human_str_replace)
  )

human_auroc <- human_roc %>%
  group_by(comparison, method) %>%
  summarise(
    auROC = integrater(FPR, TPR)
  ) %>%
  arrange(desc(auROC)) %>%
  mutate(
    FPR = seq(0.1, .98, length.out = n()),
    TPR = 0
  )

# Bruderer
bruder_roc <- ls(pattern = 'roc_bruder.*') %>%
  map(get) %>%
  bind_rows() %>%
  mutate(
    comparison = "All"
  )

bruder_auroc <- bruder_roc %>%
  group_by(method) %>%
  summarise(
    auROC = integrater(FPR, TPR)
  ) %>%
  arrange(desc(auROC)) %>%
  mutate(
    FPR = seq(0.1, .98, length.out = n()),
    TPR = 0
  )

# Navarro
nav_roc <- ls(pattern = 'roc_nav.*') %>%
  map(get) %>%
  bind_rows()

nav_auroc <- nav_roc %>%
  group_by(comparison, method) %>%
  summarise(
    auROC = integrater(FPR, TPR)
  ) %>%
  arrange(desc(auROC)) %>%
  mutate(
    FPR = seq(0.1, .98, length.out = n()),
    TPR = 0
  )
