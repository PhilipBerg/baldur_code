#### Yeast ####
yeast_design <- model.matrix(~0+factor(rep(1:2, each = 3)))
colnames(yeast_design) <- paste0('ng', c(50, 100))

yeast_prnn <- yeast_prnn %>%
  trend_partitioning(yeast_design, h = .1)

#### UPS ####
ups_design <- model.matrix(~0+factor(rep(1:3, each = 4)))
colnames(ups_design) <- paste0('fmol', c(25, 50, 100))

ups_prnn <- ups_prnn %>%
  trend_partitioning(ups_design, h = .1)

#### Ramus ####
ramus_design <- model.matrix(~0+factor(rep(1:9, each = 3)))
colnames(ramus_design) <- paste0('condi', 1:9)

ramus_prnn <- ramus_prnn %>%
  trend_partitioning(ramus_design, h = .1)

#### Human ####
human_design <- model.matrix(~0+factor(rep(1:3, each = 23)))
colnames(human_design) <- paste0('spike_prop_', c(25, 12, 6))

human_prnn <- human_prnn %>%
  trend_partitioning(human_design, h = .1)

#### Calculate McFadden’s Pseudo R2 ####
# Single trend
r2sin <- ls(pattern = '_prnn') %>%
  setNames(., .) %>%
  map(get) %>%
  map(fit_gamma_regression, sd ~ mean) %>%
  map_dbl(calc_mcfadden) %>%
  round(2)

r2mix <- ls(pattern = '_prnn') %>%
  setNames(., .) %>%
  map(get) %>%
  map(fit_gamma_regression) %>%
  map_dbl(calc_mcfadden) %>%
  round(2)

#### Clustering table ####
ls(pattern = '_prnn') %>%
  setNames(., .) %>%
  map(get) %>%
  map2(c('ECOLI', 'UPS', 'UPS', 'YEAST'), ~.x %$%
         table(c, str_detect(identifier, .y)) %>%
         as.data.frame.matrix() %>%
         as_tibble(rownames = 'cluster')
  ) %>%
  imap(
    ~ .x %>%
      mutate(
        Dataset = .y
      )
  ) %>%
  bind_rows() %>%
  rename(`False Positive` = 2, `True Positive` = 3) %>%
  mutate(
    Dataset = str_remove(Dataset, '_prnn') %>%
      str_replace_all(c('ups' = 'UPS', 'ramus' = 'Ramus'))
  ) %>%
  unite(tmp, 3:2, sep = '/') %>%
  spread(cluster, tmp) %>%
  rename_with(
    ~ paste0(., ' (TP/TN)') %>%
      str_replace_all(c('L' = 'Lower', 'U' = 'Upper')),
    -1
  )

ls(pattern = '_prnn') %>%
  setNames(., .) %>%
  map(get) %>%
  map2(c('ECOLI', 'UPS', 'UPS', 'YEAST'), ~.x %$%
         table(c, str_detect(identifier, .y)) %>%
         as.data.frame.matrix() %>%
         as_tibble(rownames = 'cluster')
  ) %>%
  imap(
    ~ .x %>%
      mutate(
        Dataset = .y,
        `TRUE` = round(`TRUE`/sum(`TRUE`), 3),
        `FALSE` = round(`FALSE`/sum(`FALSE`), 3)
      )
  ) %>%
  bind_rows() %>%
  rename(`False Positive` = 2, `True Positive` = 3) %>%
  mutate(
    Dataset = str_remove(Dataset, '_prnn') %>%
      str_replace_all(c('ups' = 'UPS', 'ramus' = 'Ramus'))
  ) %>%
  unite(tmp, 3:2, sep = '/') %>%
  spread(cluster, tmp) %>%
  rename_with(
    ~ paste0(., ' (TP/TN)') %>%
      str_replace_all(c('L' = 'Lower', 'U' = 'Upper')),
    -1
  )