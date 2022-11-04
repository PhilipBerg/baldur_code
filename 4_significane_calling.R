#### YEAST ####
#### Gamma regressions ####
# Single trend
yeast_sin_gam <- yeast_prnn %>%
  fit_gamma_regression(sd ~ mean)
# Mixed trend
yeast_mix_gam <- yeast_prnn %>%
  fit_gamma_regression(sd ~ mean + c)

#### Limma & t-test ####
yeast_contrast <- limma::makeContrasts(contrasts = 'ng50-ng100', levels = yeast_design)

yeast_trend <- yeast_prnn %>%
  select(-mean:-c) %>%
  limma_wrapper(yeast_design, yeast_contrast)

# yeast_mv_weights <- yeast_prnn %>%
#   estimate_uncertainty('identifier', yeast_design, yeast_mix_gam) %>%
#   magrittr::raise_to_power(-2)
#
# yeast_limma_mix <- yeast_prnn %>%
#   select(-mean:-c) %>%
#   limma_wrapper(yeast_design, yeast_contrast, yeast_mv_weights)
#
# yeast_weights <- yeast_prnn %>%
#   estimate_uncertainty('identifier', yeast_design, yeast_sin_gam) %>%
#   magrittr::raise_to_power(-2)
#
# yeast_limma_gamma <- yeast_prnn %>%
#   select(-mean:-c) %>%
#   limma_wrapper(yeast_design, yeast_contrast, yeast_weights)

yeast_ttest <- yeast_prnn %>%
  drop_na() %>%
  group_by(identifier) %>%
  mutate(
    lfc = mean(c(ng50_1, ng50_2, ng50_3)) - mean(c(ng100_1, ng100_2, ng100_3)),
    p_val = t.test(c(ng50_1, ng50_2, ng50_3), c(ng100_1, ng100_2, ng100_3), var.equal = T)$p.value,
    p_val = p.adjust(p_val, 'fdr'),
    comparison = 'ng100 vs ng50'
  )

#### Baldur ####
yeast_cont <- matrix(1:2, 1)

yeast_sin_baldur_results <- baldur_wrapper(
  yeast_prnn,
  yeast_design,
  yeast_cont,
  yeast_sin_gam,
  baldur:::stanmodels$empirical_bayes,
  workers
)

yeast_mix_baldur_results <- baldur_wrapper(
  yeast_prnn,
  yeast_design,
  yeast_cont,
  yeast_mix_gam,
  baldur:::stanmodels$empirical_bayes,
  workers
)

#### UPS ####
#### Gamma regressions ####
# Single trend
ups_sin_gam <- ups_prnn %>%
  fit_gamma_regression(sd ~ mean)
# Mixed trend
ups_mix_gam <- ups_prnn %>%
  fit_gamma_regression(sd ~ mean + c)

#### Limma & t-test ####
ups_contrast <- combn(colnames(ups_design), 2) %>%
  t() %>%
  apply(1, str_flatten, '-') %>%
  limma::makeContrasts(contrasts = ., levels = ups_design)

ups_trend <- ups_prnn %>%
  select(-mean:-c) %>%
  limma_wrapper(ups_design, ups_contrast)

# ups_mv_weights <- ups_prnn %>%
#   estimate_uncertainty('identifier', yeast_design, yeast_mix_gam) %>%
#   magrittr::raise_to_power(-2)
#
# ups_limma_mix <- ups_prnn %>%
#   limma_wrapper(ups_design, ups_contrast, ups_mv_weights)

# ups_weights <- ups_prnn %>%
#   fit_gamma_weights(ups_design, 'identifier') %>%
#   pair::calc_weights(ups_prnn, .) %>%
#   select(-1) %>%
#   as.matrix()
#
# ups_limma_gamma <- ups_prnn %>%
#   limma_wrapper(ups_design, ups_contrast, ups_weights)

ups_ttest <- ups_prnn %>%
  group_by(identifier) %>%
  transmute(
    p_val_25_50 = t.test(c(fmol25_1, fmol25_2, fmol25_3, fmol25_4), c(fmol50_1, fmol50_2, fmol50_3, fmol50_4), var.equal = T)$p.value,
    p_val_25_100 = t.test(c(fmol25_1, fmol25_2, fmol25_3, fmol25_4), c(fmol100_1, fmol100_2, fmol100_3, fmol100_4), var.equal = T)$p.value,
    p_val_50_100 = t.test(c(fmol50_1, fmol50_2, fmol50_3, fmol50_4), c(fmol100_1, fmol100_2, fmol100_3, fmol100_4), var.equal = T)$p.value,
    across(contains('p_val'), p.adjust, 'fdr')
  ) %>%
  ungroup() %>%
  pivot_longer(where(is.numeric), values_to = 'p_val', names_to = 'comparison') %>%
  mutate(
    comparison = str_replace(comparison, '^.*_([0-9]*)_(.*)$', 'fmol\\1 vs fmol\\2')
  )



#### Baldur ####
ups_cont <- combn(1:3, 2) %>%
  t()

ups_sin_baldur_results <- baldur_wrapper(
  ups_prnn,
  ups_design,
  ups_cont,
  ups_sin_gam,
  baldur:::stanmodels$empirical_bayes,
  workers
)

ups_mix_baldur_results <- baldur_wrapper(
  ups_prnn,
  ups_design,
  ups_cont,
  ups_mix_gam,
  baldur:::stanmodels$empirical_bayes,
  workers
)


#### Ramus ####
ramus_sin_gam <- ramus_prnn %>%
  fit_gamma_regression(sd ~ mean)
# Mixed trend
ramus_mix_gam <- ramus_prnn %>%
  fit_gamma_regression(sd ~ mean + c)

#### Limma & t-test ####
ramus_contrast <- combn(colnames(ramus_design), 2) %>%
  t() %>%
  apply(1, str_flatten, '-') %>%
  limma::makeContrasts(contrasts = ., levels = ramus_design)

ramus_trend <- ramus_prnn %>%
  select(-mean:-c) %>%
  limma_wrapper(ramus_design, ramus_contrast)

ramus_ttest <- map(colnames(ramus_contrast), ttest_wrapper, ramus_prnn) %>%
  bind_rows()

#### Baldur ####
ramus_cont <- combn(1:9, 2) %>%
  t()

ramus_sin_baldur_results <- baldur_wrapper(
  ramus_prnn,
  ramus_design,
  ramus_cont,
  ramus_sin_gam,
  baldur:::stanmodels$empirical_bayes,
  workers
)

ramus_mix_baldur_results <- baldur_wrapper(
  ramus_prnn,
  ramus_design,
  ramus_cont,
  ramus_mix_gam,
  baldur:::stanmodels$empirical_bayes,
  workers
)


#### Human ####
#### Gamma regressions ####
# Single trend
human_sin_gam <- human_prnn %>%
  fit_gamma_regression(sd ~ mean)
# Mixed trend
human_mix_gam <- human_prnn %>%
  fit_gamma_regression(sd ~ mean + c)

#### Limma & t-test ####
human_contrast <- combn(colnames(human_design), 2) %>%
  t() %>%
  apply(1, str_flatten, '-') %>%
  limma::makeContrasts(contrasts = ., levels = ups_design)

human_trend <- human_prnn %>%
  select(-mean:-c) %>%
  limma_wrapper(human_design, human_contrast)

human_ttest <- map(colnames(human_contrast), ttest_wrapper, human_prnn) %>%
  bind_rows()


#### Baldur ####
human_cont <- combn(1:3, 2) %>%
  t()

human_sin_baldur_results <- baldur_wrapper(
  human_prnn,
  human_design,
  human_cont,
  human_sin_gam,
  baldur:::stanmodels$empirical_bayes,
  workers
)

human_mix_baldur_results <- baldur_wrapper(
  human_prnn,
  human_design,
  human_cont,
  human_mix_gam,
  baldur:::stanmodels$empirical_bayes,
  workers
)
