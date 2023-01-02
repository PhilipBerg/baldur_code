if (!("pacman" %in% .packages(all.available = T))) {
  install.packages("pacman")
  library("pacman")
} else if(!("pacman" %in% (.packages()))) {
  library("pacman")
}
if (!("devtools" %in% .packages(all.available = T))) {
  install.packages("devtools")
}

if (!("BiocManager" %in% .packages(all.available = T))) {
  install.packages("BiocManager")
  BiocManager::install(version = "3.12")
}
if (!("baldur" %in% .packages(all.available = T))) {
  devtools::install_github("PhilipBerg/baldur")
}

p_load(fitdistrplus, vroom, magrittr, tibble, plyr, dplyr,
       purrr, tidyr, stringr, readxl, multidplyr,
       ggplot2, future, furrr, janitor, baldur, rstan)
if (!("limma" %in% .packages(all.available = T))) {
  BiocManager::install("limma")
  library("limma")
}else{
  library("limma")
}

calc_mcfadden <- function(model){
  1 - model$deviance / model$null.deviance
}

load_data <- function(data) {
  if (data == 'yeast') {
    if (file.exists('yeast_data.RData')) {
      load('yeast_data.RData')
      return(yeast_prnn)
    } else {
      yeast_prnn <- yeast %>%
        psrn(load_info = F, id_col = 'identifier') %>%
        mf_wrapper()
      save(yeast_prnn, file = 'yeast_data.RData')
      return(yeast_prnn)
    }
  } else if (data == 'ups') {
    if (file.exists('ups_data.RData')) {
      load('ups_data.RData')
      return(ups_prnn)
    } else {
      ups_prnn <- ups %>%
        psrn('identifier') %>%
        mf_wrapper()
      save(ups_prnn, file = 'ups_data.RData')
      return(ups_prnn)
    }
  } else if (data == 'ramus') {
    if (file.exists('ramus_data.RData')) {
      load('ramus_data.RData')
      return(ramus_prnn)
    } else {
      ramus_prnn <- readr::read_csv('ramus_clean.csv') %>%
        rename_with(~paste0('condi', rep(1:9, each = 3), '_', rep(1:3, lenght.out = 9*3)), where(is.numeric)) %>%
        psrn(load_info = F, id_col = 'identifier') %>%
        mf_wrapper()
      save(ramus_prnn, file = 'ramus_data.RData')
      return(ramus_prnn)
    }
  } else if (data == 'human') {
    if (file.exists('human_data.RData')) {
      load('human_data.RData')
      return(human_prnn)
    } else {
      human_prnn <- readxl::read_excel('diaWorkflowResults_allDilutions.xlsx', na = "NA") %>% #, sheet = 'OSW_DIANN_AI_GPF') %>%
        dplyr::select(-2:-24) %>%
        rename_with(
          ~str_replace_all(.,
                           c(
                             '\\.\\.\\.1' = 'identifier',
                             '1-' = 'spike_prop_'
                           )
          ),
          everything()
        ) %>%
        rename_with(
          ~str_remove(., '[0-9]*$') %>%
            paste0(., 1:23),
          where(is.numeric)
        ) %>%
        filter(
          rowSums(across(-1, is.na)) != 69
        ) %>%
        mf_wrapper()
      save(human_prnn, file = 'human_data.RData')
      return(human_prnn)
    }
  } else {
    stop('Dataset does not exist')
  }
}

fit_lgmr <- function(data, model, mixed = TRUE, p_beta = c(1, 1), prop = .01, bound = .5) {
  if (!'sd' %in% names(data)) {
    stop('sd is not a column in the data\n Did you forget to calculate the Mean-Variance trend?')
  } else if(is.unsorted(data$sd) & is.unsorted(data$sd, strictly = TRUE)) {
    stop('Data is not sorted from smallest to largest sd.\n Did you forget to sort the data?')
  }
  n <- nrow(data)
  input <- list(
    N = n, U = round(n*prop), M = mixed,
    y = data$sd, x = data$mean,
    p_beta = p_beta, lb_bound = 1 - bound, ub_bound = bound
  )
  rstan::sampling(
    model, data = input,
    cores = 5, chains = 5,
    iter = 2500, warmup = 500,
    control = list(adapt_delta = .9)
  )
}

reg_fun <- function(p, reg_pars, y_bar, mixed = TRUE){
  q <- 1 - p
  out <- exp(q*(reg_pars['I_U'] - reg_pars['S_U']*y_bar))
  if (mixed) {
    exp(p*(reg_pars['I_L'] - reg_pars['S_L']*y_bar))*out
  } else {
    exp(reg_pars['I_L'] - reg_pars['S_L']*y_bar)*out
  }
}

plot_regression <- function(data, reg_pars, mixed = TRUE) {
  anno_pars <- round(reg_pars, 2)
  if(mixed) {
    eq_part1 <- paste0(
      'bold(\u03bc) == exp(bold(p)%.%(',
      sprintf(
        '%1$s - %2$s*bold(bar(y)))',
        anno_pars['I_L'], anno_pars['S_L']
      )
    )
  } else {
    eq_part1 <- paste0(
      'bold(\u03bc) == exp(',
      sprintf(
        '%1$s - %2$s*bold(bar(y))',
        anno_pars['I_L'], anno_pars['S_L']
      )
    )
  }
  eq <- paste0(
    eq_part1,
    ' + bold(q)%.%',
    sprintf(
      '(%1$s - %2$s*bold(bar(y))))',
      anno_pars['I_U'], anno_pars['S_U']
    )
  )
  upper <- ~ reg_fun(0,   reg_pars, .x, mixed)
  middl <- ~ reg_fun(0.5, reg_pars, .x, mixed)
  lower <- ~ reg_fun(1,   reg_pars, .x, mixed)
  data %>%
    ggplot(aes(mean, sd, color = p)) +
    geom_point(size = .1) +
    theme_classic() +
    scale_color_viridis_c(end = .9) +
    theme(
      legend.position = c(.75, .75)
    ) +
    stat_function(aes(color = 1),   fun = lower, linewidth = .9, n = 10000) +
    stat_function(aes(color = .5),  fun = middl, linewidth = .9, n = 10000) +
    stat_function(aes(color = 0),   fun = upper, linewidth = .9, n = 10000) +
    labs(
      x = expression(bold(bar(y))),
      y = expression(bold(s))
    ) +
    annotate('text', Inf, Inf, hjust = 1, vjust = 1, label = eq, parse = T)
}

#### RF Imputation ####
mf_wrapper <- function(data) {
  auxilary <- data %>%
    select(-where(is.numeric))
  cl <- parallel::makeCluster(
    min(ncol(data) - 1, round(parallel::detectCores()/2))
  )
  doParallel::registerDoParallel(cl)
  out <- data %>%
    select(where(is.numeric)) %>%
    as.data.frame() %>%
    missForest::missForest(maxiter = 20, parallelize = 'forests', verbose = TRUE) %>%
    use_series(ximp) %>%
    as_tibble() %>%
    bind_cols(auxilary, .)
  rm(cl)
  gc(F)
  return(out)
}

#### Imputation ####
trend_imputation <- function (data, design) {
  conditions <- get_conditions(design)
  LOQ <- data %>%
    dplyr::select(dplyr::matches(conditions)) %>%
    unlist(T, F) %>%
    {
      stats::quantile(., 0.25, na.rm = T) - 1.5 * stats::IQR(., na.rm = T)
    } %>%
    unname()
  order <- data %>%
    colnames()
  gam_reg <- estimate_trend(data, design)
  order <- data %>%
    colnames()
  tmp_cols <- data %>%
    select(-matches(conditions)) %>%
    colnames()
  means <- data %>%
    select(matches(conditions)) %>%
    split.default(str_extract(colnames(.), conditions)) %>%
    map(rowMeans, na.rm = T)
  data %>%
    bind_cols(means) %>%
    select(matches(conditions)) %>%
    split.default(str_extract(colnames(.), conditions)) %>%
    map(rename, mean = last_col()) %>%
    map(
      mutate,
      mean = replace_na(mean, LOQ),
      mean = if_else(mean > LOQ, mean, LOQ),
      sd_pred = predict.glm(gam_reg, newdata = data.frame(mean = mean), type = 'response'),
      imp = apply(across(everything()), 1, impute_row, LOQ)
    ) %>%
    map(select, imp) %>%
    map(unnest, imp) %>%
    bind_cols(data[tmp_cols], .) %>%
    select(all_of(order))
}

get_conditions <- function(design) {
  colnames(design) %>%
    str_flatten('|')
}

pooled_sd <- function(data, design) {
  condi <- get_conditions(design)
  n <- colSums(design)
  data %>%
    select(matches(condi)) %>%
    split.default(str_extract(colnames(.), condi)) %>%
    map2(n,
         ~ apply(.x, 1, sd)*(.y - 1)
    ) %>%
    as_tibble() %>%
    summarise(
      sd = apply(across(everything()), 1, sum)/(sum(n) - length(n))
    ) %>%
    bind_cols(select(data, -sd), .)
}

estimate_trend <- function(data, design) {
  data %>%
    drop_na() %>%
    calculate_mean_sd_trends(design) %>%
    #pooled_sd(design) %>%
    fit_gamma_regression(sd ~ mean)
}

impute_row <- function (data, LOQ) {
  if (anyNA(data)) {
    missing <- is.na(data)
    data[missing] <- stats::rnorm(
      n = sum(missing),
      mean = max(data['mean'], LOQ), sd = data['sd_pred']
    )
  }
  as_tibble(as.list(data[seq_len(length(data)-2)]))
}

#### Data related functions ####
get_human_data <- function() {
  download.file('https://zenodo.org/record/6379087/files/Source%20Data.zip?download=1', 'human_data.zip', method = 'curl')
  unzip('human_data.zip')
  file.rename('Source Data/MinimumDataset/diaWorkflowResults_allDilutions.xlsx', 'diaWorkflowResults_allDilutions.xlsx')
  file.remove('human_data.zip')
  unlink('__MACOSX', TRUE)
  unlink('Source Data', TRUE)
}
get_ramus_data <- function() {
  readr::read_csv('https://figshare.com/ndownloader/files/35592290?private_link=28e837bfe865e8f13479', show_col_types = FALSE) %>%
    janitor::clean_names() %>%
    mutate(
      across(where(is.numeric), ~na_if(.x, 0))
    ) %>%
    readr::write_csv('ramus_clean.csv')
}
get_data <- function() {
  if(!file.exists('diaWorkflowResults_allDilutions.xlsx')) {
    get_human_data()
  }
  if(!file.exists('ramus_clean.csv')) {
    get_ramus_data()
  }
}

#### Performance ####
calc_tpr_fpr <- function(alpha, hits, p, n, p_val_col){
  p_val <- pull(hits, p_val_col)
  tp <- hits$tp
  tn <- hits$tn
  tibble(
    TP = sum((p_val<=alpha) & tp),
    TN = sum((p_val>=alpha) & tn),
    FP = sum((p_val<=alpha) & tn),
    FN = sum((p_val>=alpha) & tp)
  ) %>%
    mutate(
      TPR = TP/p,
      FPR = FP/n,
      NPV = TN/(TN+FN),
      precision = TP/(TP + FP),
      MCC = sqrt(precision)*sqrt(TPR)*sqrt(1-FPR)*sqrt(NPV) - sqrt(1-precision)*sqrt(1-TPR)*sqrt(FPR)*sqrt(1-NPV)
    )
}

create_roc <- function(p_val_col, data, regex, cl = NULL){
  data <- data %>%
    split.data.frame(.$comparison)
  out <- map(data,
             pull, p_val_col
  ) %>%
    map(unique) %>%
    map(~c(-0.1, .)) %>%
    map(enframe, value = 'alpha', name = NULL)
  data <- data %>%
    map(
      mutate,
      tp = str_detect(identifier, regex),
      tn = !tp
    )
  p <- data %>%
    map(
      use_series, tp
    ) %>%
    map(sum)
  n <- data %>%
    map(
      use_series, tn
    ) %>%
    map(sum)
  pmap(list(out, data, p, n), create_roc_helper, cl, p_val_col) %>%
    bind_rows()
}

create_roc_helper <- function(out, data, p, n, cl, p_val_col) {
  if(!is.null(cl)){
    multidplyr::cluster_library(cl,
                                c("dplyr",
                                  "stringr",
                                  "tidyr",
                                  "purrr",
                                  "tibble"
                                )
    )
    multidplyr::cluster_copy(cl,
                             c(
                               "calc_tpr_fpr",
                               "data",
                               "regex",
                               "p",
                               "n"
                             )
    )
    out <- out %>%
      multidplyr::partition(cl)
  }
  out <- out %>%
    mutate(
      results = map(alpha, calc_tpr_fpr, data, p, n, p_val_col),
      comparison = data$comparison[[1]]
    )
  if(!is.null(cl)){
    out <- out %>%
      collect()
  }
  out %>%
    unnest(results)
}

integrater <- function(x, y) {
  x_na <- is.na(x)
  y_na <- is.na(y)
  x <- x[!x_na&!y_na]
  y <- y[!x_na&!y_na]
  x_order <- order(x)
  x <- x[x_order]
  y <- y[x_order]
  dx <- diff(x)
  end <- length(y)
  my <- (y[1:(end - 1)] + y[2:end]) / 2
  sum(dx * my)
}

#### Testing ####
limma_wrapper <- function(data, design, contrast, weights = NULL){
  data <- data %>%
    dplyr::select(-1) %>%
    as.data.frame() %>%
    magrittr::set_rownames(
      data %>%
        magrittr::use_series(identifier)
    )
  hit <- limma::lmFit(data, design, weights = weights) %>%
    limma::contrasts.fit(contrast) %>%
    limma::eBayes(robust = T, trend = is.null(weights))
  colnames(contrast) %>%
    stats::setNames(colnames(contrast)) %>%
    purrr::map(
      ~limma::topTable(hit, .x, number = Inf, adjust.method = "fdr") %>%
        tibble::rownames_to_column('identifier') %>%
        tibble::as_tibble()
    ) %>%
    purrr::imap(dplyr::mutate) %>%
    purrr::map(dplyr::rename, comparison = dplyr::last_col()) %>%
    dplyr::bind_rows() %>%
    dplyr::select(
      all_of('identifier'), lfc = logFC, mean = AveExpr,
      p_val = P.Value, comparison
    ) %>%
    dplyr::mutate(comparison = stringr::str_replace(comparison, "-", " vs "))
}

calc_weights <- function(df, gm){
  id <- df$identifier
  df %>%
    mutate(
      dplyr::across(where(is.numeric),
                    ~stats::predict.glm(
                      gm,
                      newdata = data.frame(mean = ., c = c),
                      type = "response"
                    ),
      )
    ) %>%
    select(-c,-1,-identifier) %>%
    as.matrix() %>%
    set_rownames(id)
}
ttest_wrapper <- function(contrast, data) {
  comp <- str_split(contrast, '-')[[1]]
  data %>%
    group_by(identifier) %>%
    mutate(
      comparison = str_flatten(comp, ' vs '),
      p_val = t.test(across(contains(comp[1])), across(contains(comp[2])), var.equal = T)$p.value,
      p_val = p.adjust(p_val, 'fdr')
    ) %>%
    select(identifier, comparison, p_val)
}
baldur_wrapper <- function(data, design, contrast, gam_model, bald_model, workers) {
  uncertainty <- data %>%
    estimate_uncertainty('identifier', design, gam_model)
  data %>%
    estimate_gamma_priors(design, gam_model) %>%
    sample_posterior(
      id_col_name = 'identifier',
      design_matrix = design,
      contrast_matrix = contrast,
      uncertainty_matrix = uncertainty,
      bayesian_model = bald_model,
      clusters = workers
    )
}
#### Figure saving wrapper ####
ggsave_wrapper <- function(file_name, width, height = NA){
  dir.create('figures', F)
  ggsave(paste0('figures/', file_name, '.tiff'), width = width, height = height, units = 'mm', dpi = 320, compression = 'lzw')
  ggsave(paste0('figures/', file_name, '.pdf'), width = width, height = height, units = 'mm', dpi = 320)
}
