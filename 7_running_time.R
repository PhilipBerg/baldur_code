#### Clocking data and decision model ####
## Yeast-DS
yeast_cont <- matrix(c(1, -1), 2)
yeast_sin_gam <- yeast_prnn %>%
  fit_gamma_regression(sd ~ mean)

print('yeast')
yeast <- microbenchmark::microbenchmark(
  yeast2 = baldur_wrapper(
    yeast_prnn,
    yeast_design,
    yeast_cont,
    yeast_sin_gam,
    empirical_bayes,
    2
  ),
  yeast4 = baldur_wrapper(
    yeast_prnn,
    yeast_design,
    yeast_cont,
    yeast_sin_gam,
    empirical_bayes,
    4
  ),
  yeast6 = baldur_wrapper(
    yeast_prnn,
    yeast_design,
    yeast_cont,
    yeast_sin_gam,
    empirical_bayes,
    6
  ),
  yeast8 = baldur_wrapper(
    yeast_prnn,
    yeast_design,
    yeast_cont,
    yeast_sin_gam,
    empirical_bayes,
    8
  ),
  yeast10 = baldur_wrapper(
    yeast_prnn,
    yeast_design,
    yeast_cont,
    yeast_sin_gam,
    empirical_bayes,
    10
  ),
  yeast16 = baldur_wrapper(
    yeast_prnn,
    yeast_design,
    yeast_cont,
    yeast_sin_gam,
    empirical_bayes,
    16
  ),
  yeast32 = baldur_wrapper(
    yeast_prnn,
    yeast_design,
    yeast_cont,
    yeast_sin_gam,
    empirical_bayes,
    32
  ),
  yeast64 = baldur_wrapper(
    yeast_prnn,
    yeast_design,
    yeast_cont,
    yeast_sin_gam,
    empirical_bayes,
    64
  ),
  times = 10
)

## UPS-DS
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

ups_sin_gam <- ups_prnn %>%
  fit_gamma_regression(sd ~ mean)

print('ups')
ups <- microbenchmark::microbenchmark(
  ups2 = baldur_wrapper(
    ups_prnn,
    ups_design,
    ups_cont,
    ups_sin_gam,
    empirical_bayes,
    2
  ),
  ups4 = baldur_wrapper(
    ups_prnn,
    ups_design,
    ups_cont,
    ups_sin_gam,
    empirical_bayes,
    4
  ),
  ups6 = baldur_wrapper(
    ups_prnn,
    ups_design,
    ups_cont,
    ups_sin_gam,
    empirical_bayes,
    6
  ),
  ups8 = baldur_wrapper(
    ups_prnn,
    ups_design,
    ups_cont,
    ups_sin_gam,
    empirical_bayes,
    8
  ),
  ups10 = baldur_wrapper(
    ups_prnn,
    ups_design,
    ups_cont,
    ups_sin_gam,
    empirical_bayes,
    10
  ),
  ups16 = baldur_wrapper(
    ups_prnn,
    ups_design,
    ups_cont,
    ups_sin_gam,
    empirical_bayes,
    16
  ),
  ups32 = baldur_wrapper(
    ups_prnn,
    ups_design,
    ups_cont,
    ups_sin_gam,
    empirical_bayes,
    32
  ),
  ups64 = baldur_wrapper(
    ups_prnn,
    ups_design,
    ups_cont,
    ups_sin_gam,
    empirical_bayes,
    64
  ),
  times = 10
)

## Ramus-DS
ramus_cont <- combn(1:9, 2) %>%
  apply(2, make_pairwise_contrast, 9)

ramus_sin_gam <- ramus_prnn %>%
  fit_gamma_regression(sd ~ mean)

print('ramus')
ramus <- microbenchmark::microbenchmark(
  ramus2 = baldur_wrapper(
    ramus_prnn,
    ramus_design,
    ramus_cont,
    ramus_sin_gam,
    empirical_bayes,
    2
  ),
  ramus4 = baldur_wrapper(
    ramus_prnn,
    ramus_design,
    ramus_cont,
    ramus_sin_gam,
    empirical_bayes,
    4
  ),
  ramus6 = baldur_wrapper(
    ramus_prnn,
    ramus_design,
    ramus_cont,
    ramus_sin_gam,
    empirical_bayes,
    6
  ),
  ramus8 = baldur_wrapper(
    ramus_prnn,
    ramus_design,
    ramus_cont,
    ramus_sin_gam,
    empirical_bayes,
    8
  ),
  ramus10 = baldur_wrapper(
    ramus_prnn,
    ramus_design,
    ramus_cont,
    ramus_sin_gam,
    empirical_bayes,
    10
  ),
  ramus16 = baldur_wrapper(
    ramus_prnn,
    ramus_design,
    ramus_cont,
    ramus_sin_gam,
    empirical_bayes,
    16
  ),
  ramus32 = baldur_wrapper(
    ramus_prnn,
    ramus_design,
    ramus_cont,
    ramus_sin_gam,
    empirical_bayes,
    32
  ),
  ramus64 = baldur_wrapper(
    ramus_prnn,
    ramus_design,
    ramus_cont,
    ramus_sin_gam,
    empirical_bayes,
    64
  ),
  times = 10
)

## Human-DS
human_cont <- combn(1:3, 2) %>%
  apply(2, make_pairwise_contrast, 3)

human_sin_gam <- human_prnn %>%
  fit_gamma_regression(sd ~ mean)

print('human')
human <- microbenchmark::microbenchmark(
  human2 = baldur_wrapper(
    human_prnn,
    human_design,
    human_cont,
    human_sin_gam,
    empirical_bayes,
    2
  ),
  human4 = baldur_wrapper(
    human_prnn,
    human_design,
    human_cont,
    human_sin_gam,
    empirical_bayes,
    4
  ),
  human6 = baldur_wrapper(
    human_prnn,
    human_design,
    human_cont,
    human_sin_gam,
    empirical_bayes,
    6
  ),
  human8 = baldur_wrapper(
    human_prnn,
    human_design,
    human_cont,
    human_sin_gam,
    empirical_bayes,
    8
  ),
  human10 = baldur_wrapper(
    human_prnn,
    human_design,
    human_cont,
    human_sin_gam,
    empirical_bayes,
    10
  ),
  human16 = baldur_wrapper(
    human_prnn,
    human_design,
    human_cont,
    human_sin_gam,
    empirical_bayes,
    16
  ),
  human32 = baldur_wrapper(
    human_prnn,
    human_design,
    human_cont,
    human_sin_gam,
    empirical_bayes,
    32
  ),
  human64 = baldur_wrapper(
    human_prnn,
    human_design,
    human_cont,
    human_sin_gam,
    empirical_bayes,
    64
  ),
  times = 10
)

print('Bruderer')
bruder_cont <- combn(1:8, 2) %>%
  apply(2, make_pairwise_contrast, 8)

bruder_sin_gam <- bruder_prnn %>%
  fit_gamma_regression(sd ~ mean)

bruder <- microbenchmark::microbenchmark(
  bruder2 = baldur_wrapper(
    bruder_prnn,
    bruder_design,
    bruder_cont,
    bruder_sin_gam,
    empirical_bayes,
    2
  ),
  bruder4 = baldur_wrapper(
    bruder_prnn,
    bruder_design,
    bruder_cont,
    bruder_sin_gam,
    empirical_bayes,
    4
  ),
  bruder6 = baldur_wrapper(
    bruder_prnn,
    bruder_design,
    bruder_cont,
    bruder_sin_gam,
    empirical_bayes,
    6
  ),
  bruder8 = baldur_wrapper(
    bruder_prnn,
    bruder_design,
    bruder_cont,
    bruder_sin_gam,
    empirical_bayes,
    8
  ),
  bruder10 = baldur_wrapper(
    bruder_prnn,
    bruder_design,
    bruder_cont,
    bruder_sin_gam,
    empirical_bayes,
    10
  ),
  bruder16 = baldur_wrapper(
    bruder_prnn,
    bruder_design,
    bruder_cont,
    bruder_sin_gam,
    empirical_bayes,
    16
  ),
  bruder32 = baldur_wrapper(
    bruder_prnn,
    bruder_design,
    bruder_cont,
    bruder_sin_gam,
    empirical_bayes,
    32
  ),
  bruder64 = baldur_wrapper(
    bruder_prnn,
    bruder_design,
    bruder_cont,
    bruder_sin_gam,
    empirical_bayes,
    64
  ),
  times = 10
)


time_data_dec_model <- mget(c('yeast', 'ups', 'ramus', 'human', 'bruder')) %>%
  map(
    ~ tibble(
      tmp = .x$expr,
      time = .x$time*1e-9
    )
  ) %>%
  bind_rows() %>%
  tidyr::extract(tmp, c('dataset', 'workers'), '([a-z]*)([0-9]*)', convert = T) %>%
  mutate(
    dataset = str_to_title(dataset) %>%
      str_replace('Ups', "UPS")
  )


#### LGMR model ####
lgmr_time <- microbenchmark::microbenchmark(
  yeast  = fit_lgmr(yeast_prnn,  'identifier', simplify = TRUE, cores = 5),
  ups    = fit_lgmr(ups_prnn,    'identifier', simplify = TRUE, cores = 5),
  ramus  = fit_lgmr(ramus_prnn,  'identifier', simplify = TRUE, cores = 5),
  human  = fit_lgmr(human_prnn,  'identifier', simplify = TRUE, cores = 5),
  bruder = fit_lgmr(bruder_prnn, 'identifier', simplify = TRUE, cores = 5),
  times = 10
) %>%
  as_tibble() %>%
  rename(dataset = expr) %>%
  mutate(
    dataset = str_to_title(dataset) %>%
      str_replace('Ups', "UPS"),
    time = time*1e-9/60,
    p = case_when(
      dataset == 'Yeast' ~ nrow(yeast_prnn),
      dataset == 'UPS'   ~ nrow(ups_prnn),
      dataset == 'Ramus' ~ nrow(ramus_prnn),
      dataset == 'Human' ~ nrow(human_prnn),
      T                  ~ nrow(bruder_prnn),
    )
  )
