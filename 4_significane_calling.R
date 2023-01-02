#### YEAST ####
#### Gamma regressions ####
# Single trend
yeast_sin_gam <- yeast_prnn %>%
  fit_gamma_regression(sd ~ mean)

#### Limma & t-test ####
yeast_contrast <- limma::makeContrasts(contrasts = 'ng50-ng100', levels = yeast_design)

yeast_trend <- load_data('yeast') %>%
  limma_wrapper(yeast_design, yeast_contrast)

yeast_ttest <- load_data('yeast') %>%
  group_by(identifier) %>%
  mutate(
    lfc = mean(c(ng50_1, ng50_2, ng50_3)) - mean(c(ng100_1, ng100_2, ng100_3)),
    p_val = t.test(c(ng50_1, ng50_2, ng50_3), c(ng100_1, ng100_2, ng100_3), var.equal = T)$p.value,
    p_val = p.adjust(p_val, 'fdr'),
    comparison = 'ng100 vs ng50'
  )

#### Baldur ####
yeast_cont <- matrix(1:2, 1)

tictoc::tic('LGMR-Baldur Yeast')
yeast_baldur_results <- sample_posterior(
  data = yeast_prnn,
  id_col_name = 'identifier',
  design_matrix = yeast_design,
  contrast_matrix = yeast_cont,
  uncertainty_matrix = yeast_unc,
  bayesian_model = baldur:::stanmodels$empirical_bayes,
  clusters = workers
)
tictoc::toc()

tictoc::tic('GR-Baldur Yeast')
yeast_baldur_single_results <- baldur_wrapper(
  yeast_prnn,
  yeast_design,
  yeast_cont,
  yeast_sin_gam,
  baldur:::stanmodels$empirical_bayes,
  workers
)
tictoc::toc()

#### UPS ####
#### Gamma regressions ####
# Single trend
ups_sin_gam <- ups_prnn %>%
  fit_gamma_regression(sd ~ mean)

#### Limma & t-test ####
ups_contrast <- combn(colnames(ups_design), 2) %>%
  t() %>%
  apply(1, str_flatten, '-') %>%
  limma::makeContrasts(contrasts = ., levels = ups_design)

ups_trend <- load_data('ups') %>%
  limma_wrapper(ups_design, ups_contrast)

ups_ttest <- load_data('ups') %>%
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

ups_mix_baldur_results <- sample_posterior(
  data = ups_prnn,
  id_col_name = 'identifier',
  design_matrix = ups_design,
  contrast_matrix = ups_cont,
  uncertainty_matrix = ups_unc,
  bayesian_model = baldur:::stanmodels$empirical_bayes,
  clusters = workers
)

ups_baldur_single_results <- baldur_wrapper(
  ups_prnn,
  ups_design,
  ups_cont,
  ups_sin_gam,
  baldur:::stanmodels$empirical_bayes,
  workers
)

#### Ramus ####
ramus_sin_gam <- ramus_prnn %>%
  fit_gamma_regression(sd ~ mean)

#### Limma & t-test ####
ramus_contrast <- combn(colnames(ramus_design), 2) %>%
  t() %>%
  apply(1, str_flatten, '-') %>%
  limma::makeContrasts(contrasts = ., levels = ramus_design)

ramus_trend <- load_data('ramus') %>%
  limma_wrapper(ramus_design, ramus_contrast)

ramus_ttest <- map(colnames(ramus_contrast), ttest_wrapper, ramus_prnn) %>%
  bind_rows()

#### Baldur ####
ramus_cont <- combn(1:9, 2) %>%
  t()

tictoc::tic('LGMR-Baldur Ramus')
ramus_mix_baldur_results <- sample_posterior(
  data = ramus_prnn,
  id_col_name = 'identifier',
  design_matrix = ramus_design,
  contrast_matrix = ramus_cont,
  uncertainty_matrix = ramus_unc,
  bayesian_model = baldur:::stanmodels$empirical_bayes,
  clusters = workers
)
tictoc::toc()

tictoc::tic('GR-Baldur Ramus')
ramus_baldur_single_results <- baldur_wrapper(
  ramus_prnn,
  ramus_design,
  ramus_cont,
  ramus_sin_gam,
  baldur:::stanmodels$empirical_bayes,
  workers
)
tictoc::toc()

#### Human ####
#### Gamma regressions ####
# Single trend
human_sin_gam <- human_prnn %>%
  fit_gamma_regression(sd ~ mean)

#### Limma & t-test ####
human_contrast <- combn(colnames(human_design), 2) %>%
  t() %>%
  apply(1, str_flatten, '-') %>%
  limma::makeContrasts(contrasts = ., levels = human_design)

human_trend <- load_data('human') %>%
  limma_wrapper(human_design, human_contrast)

human_ttest <- map(colnames(human_contrast), ttest_wrapper, human_prnn) %>%
  bind_rows()


#### Baldur ####
human_cont <- combn(1:3, 2) %>%
  t()

human_mix_baldur_results <- sample_posterior(
  data = human_prnn,
  id_col_name = 'identifier',
  design_matrix = human_design,
  contrast_matrix = human_cont,
  uncertainty_matrix = human_unc,
  bayesian_model = baldur:::stanmodels$empirical_bayes,
  clusters = workers
)

human_baldur_single_results <- baldur_wrapper(
  human_prnn,
  human_design,
  human_cont,
  human_sin_gam,
  baldur:::stanmodels$empirical_bayes,
  workers
)
