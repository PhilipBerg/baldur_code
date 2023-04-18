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
       ggplot2, janitor, baldur, rstan)
if (!("limma" %in% .packages(all.available = T))) {
  BiocManager::install("limma")
  library("limma")
} else {
  library("limma")
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


#### Data related functions ####
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
      human_prnn <- readxl::read_excel('diaWorkflowResults_allDilutions.xlsx', na = "NA") %>%
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
        drop_na() %>%
        mutate(
          across(where(is.numeric), ~ raise_to_power(2, .x))
        ) %>%
        psrn('identifier')

      save(human_prnn, file = 'human_data.RData')
      return(human_prnn)
    }
  } else {
    stop('Dataset does not exist')
  }
}
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
  }
  pmap(list(out, data, p, n), create_roc_helper, cl, p_val_col) %>%
    bind_rows()
}

create_roc_helper <- function(out, data, p, n, cl, p_val_col) {
  if(!is.null(cl)){
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
      p_val = adj.P.Val, comparison
    ) %>%
    dplyr::mutate(comparison = stringr::str_replace(comparison, "-", " vs "))
}

ttest_wrapper <- function(contrast, data) {
  comp <- str_split(contrast, '-')[[1]]
  data %>%
    group_by(identifier) %>%
    mutate(
      comparison = str_flatten(comp, ' vs '),
      p_val = t.test(across(contains(comp[1])), across(contains(comp[2])), var.equal = T)$p.value
    ) %>%
    group_by(comparison) %>%
    mutate(
      p_val = p.adjust(p_val, 'fdr')
    ) %>%
    select(identifier, comparison, p_val)
}

baldur_wrapper <- function(data, design, contrast, gam_model, bald_model, workers) {
  uncertainty <- gam_model %>%
    estimate_uncertainty(data, 'identifier', design)
  gam_model %>%
    estimate_gamma_hyperparameters(data, id_col = 'identifier') %>%
    infer_data_and_decision_model(
      id_col = 'identifier',
      design_matrix = design,
      contrast_matrix = contrast,
      uncertainty_matrix = uncertainty,
      stan_model = bald_model,
      clusters = workers
    )
}

make_pairwise_contrast <- function(x, conditions) {
  m <- matrix(0, conditions, 1)
  m[x] <- c(1, -1)
  return(m)
}

#### Plotting ####
ggsave_wrapper <- function(file_name, width, height = NA, dpi = 1200){
  dir.create('figures', F)
  ggsave(paste0('figures/', file_name, '.tiff'), width = width, height = height, units = 'mm', dpi = dpi, compression = 'lzw')
  ggsave(paste0('figures/', file_name, '.pdf'), width = width, height = height, units = 'mm', dpi = dpi)
}
set_shape <- function(plt) {
  plt$layers[[1]]$aes_params$shape <- 20
  return(plt)
}

#### Table calculations ####
calc_nrmse <- function(model) {
  residuals.glm(model, type = 'response') %>%
    raise_to_power(2) %>%
    mean() %>%
    divide_by(var(model$data$sd)) %>%
    sqrt()
}
