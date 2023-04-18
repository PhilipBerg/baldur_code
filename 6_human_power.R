#!/usr/bin/env Rscript

options(warn = 1)
source('1_packages_and_function.R')
source('2_read_data.R')

produce_sets <- function(n, repeats) {
  sets <- map(1:repeats, ~ sample(1:23, n)) %>%
    map(sort)
  unique_sets(sets, n, repeats)
}

unique_sets <- function(sets, n, repeats) {
  sets <- unique(sets)
  sets2 <- map(seq_len(repeats - length(sets)), ~ sample(1:23, n)) %>%
    map(sort)
  sets <- c(sets, sets2)
  if(length(unique(sets)) < repeats) {
    unique_sets(sets, n, repeats)
  }
  return(sets)
}

run_performance <- function(dataset, replicate, human_design) {
  workers   <- floor(parallel::detectCores()/2)
  human_str_replace <- c(
    'spike_prop_6' = '1:6',
    'spike_prop_12' = '1:12',
    'spike_prop_25' = '1:25'
  )


  human_prnn <- dataset %>%
    calculate_mean_sd_trends(human_design)

	tictoc::tic('human')
	samp_human <- fit_lgmr(human_prnn, 'identifier', simplify = TRUE, cores = min(20, workers), chains = 20, warmup = 1000, iter = 2500)
	tictoc::toc()

	human_prnn <- samp_human %>%
  	estimate_gamma_hyperparameters(human_prnn, "identifier")

	human_unc <- samp_human %>%
  		estimate_uncertainty(human_prnn, 'identifier', human_design)

  # Single trend
  human_sin_gam <- human_prnn %>%
    fit_gamma_regression(sd ~ mean)

  #### Limma & t-test ####
  human_contrast <- combn(colnames(human_design), 2) %>%
    t() %>%
    apply(1, str_flatten, '-') %>%
    limma::makeContrasts(contrasts = ., levels = human_design)

  cols_to_select <- str_flatten(colnames(human_design), '|')
  human_trend <- human_prnn %>%
    select(identifier, matches(cols_to_select)) %>%
    limma_wrapper(human_design, human_contrast)

  human_ttest <- map(colnames(human_contrast), ttest_wrapper, human_prnn) %>%
    bind_rows()

  #### Baldur ####
  human_cont <- combn(1:3, 2) %>%
  apply(2, make_pairwise_contrast, 3)

	tictoc::tic('LGMR-Baldur Human')
	human_mix_baldur_results <- infer_data_and_decision_model(
  		data = human_prnn,
  		id_col = 'identifier',
  		design_matrix = human_design,
  		contrast_matrix = human_cont,
  		uncertainty_matrix = human_unc,
  		stan_model = empirical_bayes,
  		clusters = workers
	)
	tictoc::toc()

	tictoc::tic('GR-Baldur Human')
	human_sin_baldur_results <- baldur_wrapper(
  		human_prnn,
  		human_design,
  		human_cont,
  		human_sin_gam,
  		baldur:::stanmodels$empirical_bayes,
  		workers
	)
	tictoc::toc()

	tictoc::tic('LGMR-Baldur Human (WI)')
	human_mix_baldur_wi_results <- infer_data_and_decision_model(
  		data = human_prnn,
  		id_col = 'identifier',
  		design_matrix = human_design,
  		contrast_matrix = human_cont,
  		uncertainty_matrix = human_unc,
  		stan_model = weakly_informative,
  		clusters = workers
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

  cl <- multidplyr::new_cluster(workers)
  roc_human_mix_baldur_results <- create_roc('err', human_mix_baldur_results, 'ECOLI', cl) %>%
    mutate(
      method = 'LGMR-Baldur'
    )
  roc_human_sin_baldur_results <- create_roc('err', human_sin_baldur_results, 'ECOLI', cl) %>%
    mutate(
      method = 'GR-Baldur'
    )

  roc_human_mix_wi_baldur_results <- create_roc('err', human_mix_baldur_wi_results, 'ECOLI', cl) %>%
    mutate(
      method = 'LGMR-Baldur WI'
    )
  roc_human_sin_wi_baldur_results <- create_roc('err', human_sin_baldur_wi_results, 'ECOLI', cl) %>%
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
  rm(cl)
  gc()

  print(paste0('Finished: ', replicate))

  env <- rlang::caller_env(n = 0)
  ls(pattern = 'roc_human.*') %>%
    mget(envir = env) %>%
    bind_rows() %>%
    mutate(
      comparison = str_replace_all(comparison, human_str_replace)
    ) %>%
    mutate(
      method = str_replace(method, '(Baldur$)', '\\1 EB'),
      rep = as.integer(word(replicate, 2, sep = '_'))
    )
}

n <- as.integer(commandArgs(TRUE)[1])
repeats <- as.integer(commandArgs(TRUE)[2])

set.seed(1)
sets <- produce_sets(n, repeats)

columns_sets <- sets %>%
  map(
    ~ paste0(
      rep(c("spike_prop_25", "spike_prop_12", "spike_prop_6"), each = n), '_', .x
    )
  )

datasets <- map(columns_sets, ~ select(human_prnn, identifier, all_of(.x))) %>%
  set_names(paste0('replicate_', 1:repeats))

# Human
human_design <- model.matrix(~0+factor(rep(1:3, each = n)))
colnames(human_design) <- paste0('spike_prop_', c(25, 12, 6))

results <- imap(datasets,
               run_performance,
               human_design
) %>%
  bind_rows() %>%
  mutate(
    columns = n
  )
rm(workers)
image_name <- paste0('human_', n, '_columns.RData')
save.image(image_name)
warnings()
