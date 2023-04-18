#### YEAST ####
#### Limma & t-test ####
yeast_contrast <- limma::makeContrasts(contrasts = 'ng50-ng100', levels = yeast_design)

yeast_trend <- load_data('yeast') %>%
  limma_wrapper(yeast_design, yeast_contrast)

yeast_ttest <- map(colnames(yeast_contrast), ttest_wrapper, yeast_prnn) %>%
  bind_rows()

#### Baldur ####
# Contrast #
yeast_cont <- matrix(c(1, -1), 2)

# LGMR #
tictoc::tic('LGMR-Baldur Yeast')
yeast_mix_baldur_results <- baldur_wrapper(
  yeast_prnn,
  yeast_design,
  yeast_cont,
  samp,
  empirical_bayes,
  workers
)
tictoc::toc()

tictoc::tic('LGMR-Baldur Yeast (WI)')
yeast_mix_baldur_wi_results <- baldur_wrapper(
  yeast_prnn,
  yeast_design,
  yeast_cont,
  samp,
  weakly_informative,
  workers
)
tictoc::toc()

# GR #
yeast_sin_gam <- yeast_prnn %>%
  fit_gamma_regression(sd ~ mean)

tictoc::tic('GR-Baldur Yeast')
yeast_sin_baldur_results <- baldur_wrapper(
  yeast_prnn,
  yeast_design,
  yeast_cont,
  yeast_sin_gam,
  empirical_bayes,
  workers
)
tictoc::toc()

tictoc::tic('GR-Baldur Yeast (WI)')
yeast_sin_baldur_wi_results <- baldur_wrapper(
  yeast_prnn,
  yeast_design,
  yeast_cont,
  yeast_sin_gam,
  weakly_informative,
  workers
)
tictoc::toc()

#### UPS ####
#### Limma & t-test ####
complex_designs <- map(0:2, ~ (1:3 + .x) %% 3 + 1) %>%
  map( ~ colnames(ups_design)[.x]) %>%
  map_chr( ~ sprintf('%s-(%s + %s)/2', .x[1], .x[2], .x[3]))
ups_contrast <- combn(colnames(ups_design), 2) %>%
  t() %>%
  apply(1, str_flatten, '-') %>%
  c(., complex_designs) %>%
  limma::makeContrasts(
    contrasts = .,
    levels = ups_design
  )

ups_trend <- load_data('ups') %>%
  limma_wrapper(ups_design, ups_contrast) %>%
  mutate(
    comparison = str_remove(comparison, '\\('),
    comparison = str_remove(comparison, '\\)/2'),
    comparison = str_replace(comparison, '\\+', 'and')
  )

ups_ttest <- map(colnames(ups_contrast)[1:3], ttest_wrapper, ups_prnn) %>%
  bind_rows()

#### Baldur ####
ups_cont <- matrix(
  c(
    1,   -1,    0,
    1,    0,   -1,
    0,    1,   -1,
    1,   -0.5, -0.5,
   -0.5,  1,   -0.5,
   -0.5, -0.5,  1
  ),, 3, byrow = T
) %>%
  t()

# LGMR #

tictoc::tic('LGMR-Baldur UPS')
ups_mix_baldur_results <- baldur_wrapper(
  ups_prnn,
  ups_design,
  ups_cont,
  samp_ups,
  empirical_bayes,
  workers
)
tictoc::toc()

tictoc::tic('LGMR-Baldur UPS (WI)')
ups_mix_baldur_wi_results <- baldur_wrapper(
  ups_prnn,
  ups_design,
  ups_cont,
  samp_ups,
  weakly_informative,
  workers
)
tictoc::toc()

# GR #
ups_sin_gam <- ups_prnn %>%
  fit_gamma_regression(sd ~ mean)

tictoc::tic('GR-Baldur UPS')
ups_sin_baldur_results <- baldur_wrapper(
  ups_prnn,
  ups_design,
  ups_cont,
  ups_sin_gam,
  empirical_bayes,
  workers
)
tictoc::toc()

tictoc::tic('GR-Baldur UPS (WI)')
ups_sin_baldur_wi_results <- baldur_wrapper(
  ups_prnn,
  ups_design,
  ups_cont,
  ups_sin_gam,
  weakly_informative,
  workers
)
tictoc::toc()

#### Ramus ####
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
  apply(2, make_pairwise_contrast, 9)

# LGMR #
tictoc::tic('LGMR-Baldur Ramus')
ramus_mix_baldur_results <- baldur_wrapper(
  ramus_prnn,
  ramus_design,
  ramus_cont,
  samp_ramus,
  empirical_bayes,
  workers
)
tictoc::toc()

tictoc::tic('LGMR-Baldur Ramus (WI)')
ramus_mix_baldur_wi_results <- baldur_wrapper(
  ramus_prnn,
  ramus_design,
  ramus_cont,
  samp_ramus,
  weakly_informative,
  workers
)
tictoc::toc()

# GR #
ramus_sin_gam <- ramus_prnn %>%
  fit_gamma_regression(sd ~ mean)

tictoc::tic('GR-Baldur Ramus')
ramus_sin_baldur_results <- baldur_wrapper(
  ramus_prnn,
  ramus_design,
  ramus_cont,
  ramus_sin_gam,
  empirical_bayes,
  workers
)
tictoc::toc()

tictoc::tic('GR-Baldur Ramus (WI)')
ramus_sin_baldur_wi_results <- baldur_wrapper(
  ramus_prnn,
  ramus_design,
  ramus_cont,
  ramus_sin_gam,
  weakly_informative,
  workers
)
tictoc::toc()

#### Human ####
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
  apply(2, make_pairwise_contrast, 3)

# LGMR #
tictoc::tic('LGMR-Baldur Human')
human_mix_baldur_results <- baldur_wrapper(
  human_prnn,
  human_design,
  human_cont,
  samp_human,
  empirical_bayes,
  workers
)
tictoc::toc()

tictoc::tic('LGMR-Baldur Human (WI)')
human_mix_baldur_wi_results <- baldur_wrapper(
  human_prnn,
  human_design,
  human_cont,
  samp_human,
  weakly_informative,
  workers
)
tictoc::toc()

# GR #
human_sin_gam <- human_prnn %>%
  fit_gamma_regression(sd ~ mean)

tictoc::tic('GR-Baldur Human')
human_sin_baldur_results <- baldur_wrapper(
  human_prnn,
  human_design,
  human_cont,
  human_sin_gam,
  empirical_bayes,
  workers
)
tictoc::toc()

tictoc::tic('GR-Baldur Human (WI)')
human_sin_baldur_wi_results <- baldur_wrapper(
  human_prnn,
  human_design,
  human_cont,
  human_sin_gam,
  weakly_informative,
  workers
)
tictoc::toc()
