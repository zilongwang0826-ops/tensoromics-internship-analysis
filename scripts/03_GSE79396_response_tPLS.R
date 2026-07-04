options(stringsAsFactors = FALSE)
analysis_seed <- 2026
set.seed(analysis_seed)

# GSE79396 response analysis with tensorOmics tPLS.
# Early PBMC transcriptomic changes are compared with IgG and Tfh responses.

top_n_features <- 5000
ncomp_tpls <- 2
n_perm <- 199
k_folds <- 5

post_days <- c(1, 3, 7)
primary_outcome <- "D30/0_IgG"
secondary_outcome <- "D7/0_TFH.CXCR3pos"

root_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
supp_dir <- file.path(
  root_dir,
  "data",
  "GSE79396_Cell2017_supplement",
  "SupplementalPackage-MMRN"
)

if (!dir.exists(supp_dir)) {
  stop("Could not find SupplementalPackage-MMRN. Check the extracted Cell supplemental package path.")
}

source_data_dir <- file.path(supp_dir, "sourcedata")
btm_member_dir <- file.path(supp_dir, "BTM.members", "SupplementaryData_TutorialPackage")
tensoromics_src_dir <- file.path(root_dir, "external", "tensorOmics", "R")

result_dir <- file.path(root_dir, "results", "GSE79396_response")
figure_dir <- file.path(result_dir, "figures")
table_dir <- file.path(result_dir, "tables")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

install_if_missing <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
  }
}

install_if_missing(c(
  "data.table", "dplyr", "tidyr", "ggplot2", "ggrepel",
  "readr", "tibble", "matrixStats", "scales", "gsignal", "pls"
))

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(readr)
  library(tibble)
  library(matrixStats)
  library(scales)
})

ensure_tensoromics <- function() {
  if (
    requireNamespace("tensorOmics", quietly = TRUE) &&
      all(c("tpls", "tplsda", "dctii_m_transforms") %in%
            getNamespaceExports("tensorOmics"))
  ) {
    message("Using installed tensorOmics package.")
    return(list(
      tpls = tensorOmics::tpls,
      tplsda = tensorOmics::tplsda,
      dctii_m_transforms = tensorOmics::dctii_m_transforms
    ))
  }

  local_files <- file.path(
    tensoromics_src_dir,
    c(
      "vendor.R",
      "names.R",
      "tens.mproduct.R",
      "tens.tsvdm.R",
      "tens.tpls.R",
      "tens.tplsda.R"
    )
  )

  if (all(file.exists(local_files))) {
    message("Sourcing tensorOmics functions from local GitHub clone.")
    for (f in local_files) source(f)
  } else {
    message("Local tensorOmics source not found. Sourcing from GitHub.")
    raw_base <- "https://raw.githubusercontent.com/brendanlu/tensorOmics/main/R"
    remote_files <- paste0(
      raw_base, "/",
      c(
        "vendor.R",
        "names.R",
        "tens.mproduct.R",
        "tens.tsvdm.R",
        "tens.tpls.R",
        "tens.tplsda.R"
      )
    )
    for (u in remote_files) {
      message("Sourcing: ", u)
      source(u)
    }
  }

  if (!exists("tpls") || !exists("tplsda") || !exists("dctii_m_transforms")) {
    stop("Could not load tensorOmics tpls/tplsda functions.")
  }

  list(
    tpls = get("tpls"),
    tplsda = get("tplsda"),
    dctii_m_transforms = get("dctii_m_transforms")
  )
}

tensor_tools <- ensure_tensoromics()
tpls_fun <- tensor_tools$tpls
tplsda_fun <- tensor_tools$tplsda
dctii_fun <- tensor_tools$dctii_m_transforms

###############################################################################
# Helper functions
###############################################################################

normalize_zv_id <- function(x) {
  out <- as.character(x)
  out <- gsub("^ZV([0-9]+)-0*([0-9]+)$", "ZV\\1\\2", out)
  out
}

transpose_feature_table <- function(path, feature_col_name = NULL) {
  dt <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  if (!is.null(feature_col_name) && !(feature_col_name %in% names(dt))) {
    names(dt)[1] <- feature_col_name
  }
  feature_names <- dt[[1]]
  mat <- as.matrix(dt[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- make.unique(feature_names)
  out <- as.data.frame(t(mat), check.names = FALSE)
  out$subject_id <- rownames(out)
  rownames(out) <- NULL
  out
}

safe_cor_test <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4 || sd(x[ok]) == 0 || sd(y[ok]) == 0) {
    return(tibble(n = sum(ok), r = NA_real_, p_value = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = method))
  tibble(n = sum(ok), r = unname(ct$estimate), p_value = ct$p.value)
}

scale_tensor_with_reference <- function(X, ref_idx = seq_len(dim(X)[1])) {
  out <- X
  for (k in seq_len(dim(X)[3])) {
    ref_mat <- X[ref_idx, , k, drop = FALSE][, , 1]
    mu <- colMeans(ref_mat, na.rm = TRUE)
    sig <- apply(ref_mat, 2, sd, na.rm = TRUE)
    sig[!is.finite(sig) | sig == 0] <- 1

    mat <- X[, , k, drop = FALSE][, , 1]
    out[, , k] <- sweep(sweep(mat, 2, mu, "-"), 2, sig, "/")
  }
  out[!is.finite(out)] <- 0
  out
}

make_response_tensor <- function(y_scaled, subject_ids, outcome_name, time_names) {
  # tensorOmics::tpls() can drop dimensions internally when Y has only one
  # response column. The second column contains the same response with opposite
  # sign, so no new biological information is added; it only keeps Y matrix-like.
  y_mat <- cbind(y_scaled, -y_scaled)
  Y <- array(
    rep(y_mat, times = length(time_names)),
    dim = c(length(y_scaled), 2, length(time_names)),
    dimnames = list(
      subject_ids,
      c(outcome_name, paste0(outcome_name, "_opposite")),
      time_names
    )
  )
  Y
}

orient_tpls_fit <- function(fit, y_scaled) {
  for (i in seq_len(ncol(fit$x_projected))) {
    r <- suppressWarnings(cor(fit$x_projected[, i], y_scaled, use = "complete.obs"))
    if (is.finite(r) && r < 0) {
      fit$x_projected[, i] <- -fit$x_projected[, i]
      fit$y_projected[, i] <- -fit$y_projected[, i]
      fit$x_loadings[, i] <- -fit$x_loadings[, i]
      fit$y_loadings[, i] <- -fit$y_loadings[, i]
    }
  }
  fit
}

project_tpls_component1 <- function(fit, X_new) {
  transforms <- dctii_fun(dim(X_new)[3])
  Xhat <- transforms$m(X_new)
  face <- fit$faces[1]
  loading <- fit$x_loadings[, 1]
  mat <- matrix(Xhat[, , face], nrow = dim(Xhat)[1], ncol = dim(Xhat)[2])
  as.numeric(mat %*% loading)
}

cv_tpls_component1 <- function(X_raw, y, nfold = 5) {
  n <- length(y)
  if (n < nfold * 2) {
    return(tibble(
      cv_n = n,
      cv_r = NA_real_,
      cv_r2 = NA_real_,
      cv_rmse = NA_real_
    ))
  }

  fold_id <- sample(rep(seq_len(nfold), length.out = n))
  pred <- rep(NA_real_, n)

  for (fold in seq_len(nfold)) {
    test_idx <- which(fold_id == fold)
    train_idx <- setdiff(seq_len(n), test_idx)

    X_scaled <- scale_tensor_with_reference(X_raw, ref_idx = train_idx)
    X_train <- X_scaled[train_idx, , , drop = FALSE]
    X_test <- X_scaled[test_idx, , , drop = FALSE]

    y_train <- y[train_idx]
    y_mean <- mean(y_train, na.rm = TRUE)
    y_sd <- sd(y_train, na.rm = TRUE)
    if (!is.finite(y_sd) || y_sd == 0) y_sd <- 1
    y_train_scaled <- (y_train - y_mean) / y_sd

    Y_train <- make_response_tensor(
      y_train_scaled,
      dimnames(X_train)[[1]],
      "response",
      dimnames(X_train)[[3]]
    )

    fit <- try(
      tpls_fun(
        x = X_train,
        y = Y_train,
        ncomp = 1,
        mode = "regression",
        center = FALSE
      ),
      silent = TRUE
    )

    if (inherits(fit, "try-error")) next

    train_scores <- fit$x_projected[, 1]
    test_scores <- project_tpls_component1(fit, X_test)

    score_model <- lm(y_train_scaled ~ train_scores)
    pred_scaled <- predict(
      score_model,
      newdata = data.frame(train_scores = test_scores)
    )
    pred[test_idx] <- as.numeric(pred_scaled) * y_sd + y_mean
  }

  ok <- is.finite(pred) & is.finite(y)
  cv_r <- if (sum(ok) >= 4) cor(pred[ok], y[ok]) else NA_real_
  cv_r2 <- 1 - sum((y[ok] - pred[ok])^2) / sum((y[ok] - mean(y[ok]))^2)
  cv_rmse <- sqrt(mean((y[ok] - pred[ok])^2))

  tibble(
    cv_n = sum(ok),
    cv_r = cv_r,
    cv_r2 = cv_r2,
    cv_rmse = cv_rmse
  )
}

permutation_tpls_test <- function(X_scaled, y_scaled, ncomp, n_perm) {
  Y <- make_response_tensor(
    y_scaled,
    dimnames(X_scaled)[[1]],
    "response",
    dimnames(X_scaled)[[3]]
  )

  fit_obs <- tpls_fun(
    x = X_scaled,
    y = Y,
    ncomp = ncomp,
    mode = "regression",
    center = FALSE
  )

  stat_fun <- function(scores, yy) {
    cors <- apply(scores, 2, function(z) {
      suppressWarnings(cor(z, yy, use = "complete.obs"))
    })
    max(abs(cors), na.rm = TRUE)
  }

  observed_stat <- stat_fun(fit_obs$x_projected[, seq_len(ncomp), drop = FALSE], y_scaled)

  perm_stats <- rep(NA_real_, n_perm)
  if (n_perm > 0) {
    for (b in seq_len(n_perm)) {
      y_perm <- sample(y_scaled)
      Y_perm <- make_response_tensor(
        y_perm,
        dimnames(X_scaled)[[1]],
        "response",
        dimnames(X_scaled)[[3]]
      )
      fit_perm <- try(
        tpls_fun(
          x = X_scaled,
          y = Y_perm,
          ncomp = ncomp,
          mode = "regression",
          center = FALSE
        ),
        silent = TRUE
      )
      if (!inherits(fit_perm, "try-error")) {
        perm_stats[b] <- stat_fun(
          fit_perm$x_projected[, seq_len(ncomp), drop = FALSE],
          y_perm
        )
      }
    }
  }

  p_perm <- (1 + sum(perm_stats >= observed_stat, na.rm = TRUE)) /
    (1 + sum(is.finite(perm_stats)))

  list(
    observed_stat = observed_stat,
    p_perm = p_perm,
    perm_stats = perm_stats,
    fit_obs = fit_obs
  )
}

read_gmt_genes <- function(gmt_path, pattern) {
  if (!file.exists(gmt_path)) return(character(0))
  lines <- readLines(gmt_path, warn = FALSE)
  hit <- grep(pattern, lines, ignore.case = TRUE, value = TRUE)
  if (length(hit) == 0) return(character(0))
  parts <- strsplit(hit[1], "\t", fixed = TRUE)[[1]]
  unique(parts[-c(1, 2)])
}

###############################################################################
# Read original Cell supplementary data
###############################################################################

expr_path <- file.path(source_data_dir, "genetable_vax05_rma_ordered.txt")
igg_path <- file.path(source_data_dir, "formated2_IgG_NLSrenalyzed.txt")
tfh_path <- file.path(source_data_dir, "formatted2_TFH.txt")
btm_path <- file.path(source_data_dir, "formatted_BTMs.txt")
age_path <- file.path(source_data_dir, "VZV_age_gender2.txt")
gmt_path <- file.path(btm_member_dir, "BTM_for_GSEA_20131008.gmt")

message("Reading transcriptome table...")
expr_dt <- data.table::fread(expr_path, data.table = FALSE, check.names = FALSE)
gene_id <- make.unique(expr_dt[[1]])
expr_mat <- as.matrix(expr_dt[, -1, drop = FALSE])
storage.mode(expr_mat) <- "numeric"
rownames(expr_mat) <- gene_id
rm(expr_dt)
gc()

sample_info <- tibble(sample_col = colnames(expr_mat)) %>%
  mutate(
    subject_id = sub("_D[0-9]+$", "", sample_col),
    day = as.integer(sub("^.*_D([0-9]+)$", "\\1", sample_col))
  )

complete_subjects <- sample_info %>%
  filter(day %in% c(0, post_days)) %>%
  dplyr::count(subject_id) %>%
  filter(n == length(c(0, post_days))) %>%
  pull(subject_id) %>%
  sort()

message("Complete transcriptome subjects: ", length(complete_subjects))

igg_wide <- transpose_feature_table(igg_path, "subj")
tfh_wide <- transpose_feature_table(tfh_path, "Subj")
btm_wide <- transpose_feature_table(btm_path, "BTM")

age_df <- data.table::fread(age_path, data.table = FALSE, check.names = FALSE) %>%
  mutate(
    subject_id = normalize_zv_id(ZV),
    age = as.numeric(AGE),
    gender = factor(Gender),
    age_group = factor(
      ifelse(`Age Group` == "Y", "young", "elderly"),
      levels = c("young", "elderly")
    )
  ) %>%
  select(subject_id, age, gender, age_group)

srebf1_module_genes <- read_gmt_genes(gmt_path, "M178|SREBF1")
sterol_marker_genes <- c(
  "SREBF1", "SREBF2", "HMGCR", "HMGCS1", "MVK", "MVD",
  "IDI1", "SQLE", "LDLR", "INSIG1", "INSIG2", "DHCR7",
  "DHCR24", "FASN", "SC5D", "MSMO1", "ACAT2"
)
focus_genes <- unique(c(srebf1_module_genes, sterol_marker_genes))

###############################################################################
# Build transcriptomic tensor: subjects x genes x post-vaccination changes
###############################################################################

complete_cols <- sample_info %>%
  filter(subject_id %in% complete_subjects, day %in% c(0, post_days)) %>%
  pull(sample_col)

gene_var <- matrixStats::rowVars(expr_mat[, complete_cols, drop = FALSE], na.rm = TRUE)
names(gene_var) <- rownames(expr_mat)

top_variable_genes <- names(sort(gene_var, decreasing = TRUE))[
  seq_len(min(top_n_features, length(gene_var)))
]

features_keep <- unique(c(
  top_variable_genes,
  intersect(focus_genes, rownames(expr_mat))
))

message("Features in tensor: ", length(features_keep))
message("SREBF1/M178/focus genes retained: ", sum(features_keep %in% focus_genes))

build_delta_tensor <- function(subjects, features, post_days) {
  X <- array(
    NA_real_,
    dim = c(length(subjects), length(features), length(post_days)),
    dimnames = list(subjects, features, paste0("D", post_days, "-D0"))
  )

  for (i in seq_along(subjects)) {
    sid <- subjects[i]
    d0_col <- paste0(sid, "_D0")
    for (k in seq_along(post_days)) {
      day_col <- paste0(sid, "_D", post_days[k])
      X[i, , k] <- expr_mat[features, day_col] - expr_mat[features, d0_col]
    }
  }

  X
}

X_delta_all <- build_delta_tensor(complete_subjects, features_keep, post_days)

###############################################################################
# Main analysis function
###############################################################################

run_response_analysis <- function(outcome_label, outcome_table, outcome_col) {
  message("Running tPLS analysis for: ", outcome_label)

  if (!(outcome_col %in% names(outcome_table))) {
    stop("Outcome column not found: ", outcome_col)
  }

  outcome_df <- outcome_table %>%
    select(subject_id, all_of(outcome_col))
  names(outcome_df)[names(outcome_df) == outcome_col] <- "outcome"
  outcome_df <- outcome_df %>%
    mutate(outcome = as.numeric(outcome)) %>%
    filter(is.finite(outcome))

  subjects <- complete_subjects[complete_subjects %in% outcome_df$subject_id]
  y <- outcome_df$outcome[match(subjects, outcome_df$subject_id)]

  X_raw <- X_delta_all[subjects, , , drop = FALSE]
  X_scaled <- scale_tensor_with_reference(X_raw)
  y_scaled <- as.numeric(scale(y))
  names(y_scaled) <- subjects

  Y_tensor <- make_response_tensor(
    y_scaled,
    subjects,
    outcome_label,
    dimnames(X_scaled)[[3]]
  )

  tpls_res <- tpls_fun(
    x = X_scaled,
    y = Y_tensor,
    ncomp = ncomp_tpls,
    mode = "regression",
    center = FALSE
  )
  tpls_res <- orient_tpls_fit(tpls_res, y_scaled)

  score_df <- tibble(
    subject_id = subjects,
    outcome = y,
    outcome_z = y_scaled,
    Comp1 = tpls_res$x_projected[, 1],
    Comp2 = tpls_res$x_projected[, 2]
  ) %>%
    left_join(age_df, by = "subject_id") %>%
    mutate(
      response_group = factor(
        ifelse(outcome >= median(outcome, na.rm = TRUE), "high", "low"),
        levels = c("low", "high")
      )
    )

  readr::write_csv(
    score_df,
    file.path(table_dir, paste0(outcome_label, "_tpls_scores.csv"))
  )

  comp_cor <- bind_rows(
    safe_cor_test(score_df$Comp1, score_df$outcome) %>% mutate(component = "Comp1"),
    safe_cor_test(score_df$Comp2, score_df$outcome) %>% mutate(component = "Comp2")
  ) %>%
    select(component, n, r, p_value)

  lm_basic <- summary(lm(outcome ~ Comp1 + Comp2, data = score_df))
  lm_age <- try(summary(lm(outcome ~ Comp1 + age + gender, data = score_df)), silent = TRUE)

  lm_stats <- tibble(
    model = c("outcome ~ Comp1 + Comp2", "outcome ~ Comp1 + age + gender"),
    n = c(nrow(score_df), nrow(score_df %>% filter(!is.na(age), !is.na(gender)))),
    r_squared = c(lm_basic$r.squared, if (!inherits(lm_age, "try-error")) lm_age$r.squared else NA_real_),
    adj_r_squared = c(lm_basic$adj.r.squared, if (!inherits(lm_age, "try-error")) lm_age$adj.r.squared else NA_real_)
  )

  readr::write_csv(
    comp_cor,
    file.path(table_dir, paste0(outcome_label, "_component_correlations.csv"))
  )
  readr::write_csv(
    lm_stats,
    file.path(table_dir, paste0(outcome_label, "_linear_model_summary.csv"))
  )

  perm <- permutation_tpls_test(X_scaled, y_scaled, ncomp_tpls, n_perm)
  perm_summary <- tibble(
    outcome = outcome_label,
    n_subjects = length(subjects),
    n_features = dim(X_scaled)[2],
    n_time_points = dim(X_scaled)[3],
    observed_max_abs_cor = perm$observed_stat,
    permutation_p = perm$p_perm,
    n_perm = n_perm
  )
  readr::write_csv(
    perm_summary,
    file.path(table_dir, paste0(outcome_label, "_permutation_summary.csv"))
  )
  readr::write_csv(
    tibble(perm_stat = perm$perm_stats),
    file.path(table_dir, paste0(outcome_label, "_permutation_stats.csv"))
  )

  cv_summary <- cv_tpls_component1(X_raw, y, nfold = k_folds) %>%
    mutate(outcome = outcome_label, .before = 1)
  readr::write_csv(
    cv_summary,
    file.path(table_dir, paste0(outcome_label, "_tpls_comp1_cv_summary.csv"))
  )

  loading_df <- tibble(
    gene = rownames(tpls_res$x_loadings),
    Comp1_loading = tpls_res$x_loadings[, 1],
    Comp2_loading = tpls_res$x_loadings[, 2],
    abs_Comp1 = abs(Comp1_loading),
    abs_Comp2 = abs(Comp2_loading),
    in_M178_SREBF1_BTM = gene %in% srebf1_module_genes,
    sterol_focus_gene = gene %in% sterol_marker_genes
  ) %>%
    arrange(desc(abs_Comp1))

  readr::write_csv(
    loading_df,
    file.path(table_dir, paste0(outcome_label, "_tpls_gene_loadings.csv"))
  )

  top_loading_comp1 <- loading_df %>%
    arrange(desc(abs_Comp1)) %>%
    slice_head(n = 30)

  top_loading_comp2 <- loading_df %>%
    arrange(desc(abs_Comp2)) %>%
    slice_head(n = 30)

  readr::write_csv(
    top_loading_comp1,
    file.path(table_dir, paste0(outcome_label, "_top30_loadings_comp1.csv"))
  )
  readr::write_csv(
    top_loading_comp2,
    file.path(table_dir, paste0(outcome_label, "_top30_loadings_comp2.csv"))
  )

  # Matrix PLS baseline on the flattened subject x (gene-time) matrix.
  X_flat <- do.call(cbind, lapply(seq_len(dim(X_scaled)[3]), function(k) {
    mat <- X_scaled[, , k, drop = FALSE][, , 1]
    colnames(mat) <- paste0(colnames(mat), "_", dimnames(X_scaled)[[3]][k])
    mat
  }))
  colnames(X_flat) <- make.names(colnames(X_flat), unique = TRUE)
  pls_df <- data.frame(y = y, X_flat, check.names = FALSE)
  matrix_pls <- try(
    pls::plsr(y ~ ., data = pls_df, ncomp = 2, validation = "LOO", scale = FALSE),
    silent = TRUE
  )
  if (!inherits(matrix_pls, "try-error")) {
    pred <- drop(matrix_pls$validation$pred[, , 1])
    matrix_pls_summary <- tibble(
      outcome = outcome_label,
      model = "matrix PLSR flattened gene-time matrix",
      cv_n = sum(is.finite(pred)),
      cv_r = cor(pred, y, use = "complete.obs"),
      cv_r2 = 1 - sum((y - pred)^2) / sum((y - mean(y))^2),
      cv_rmse = sqrt(mean((y - pred)^2))
    )
  } else {
    matrix_pls_summary <- tibble(
      outcome = outcome_label,
      model = "matrix PLSR flattened gene-time matrix",
      cv_n = NA_integer_,
      cv_r = NA_real_,
      cv_r2 = NA_real_,
      cv_rmse = NA_real_
    )
  }
  readr::write_csv(
    matrix_pls_summary,
    file.path(table_dir, paste0(outcome_label, "_matrix_plsr_cv_summary.csv"))
  )

  # BTM module correlations, including the SREBF1 target module M178.
  btm_sub <- btm_wide %>%
    filter(subject_id %in% subjects) %>%
    arrange(match(subject_id, subjects))

  btm_features <- setdiff(names(btm_sub), "subject_id")

  btm_cor_outcome <- bind_rows(lapply(btm_features, function(feature_name) {
    z <- as.numeric(btm_sub[[feature_name]])
    safe_cor_test(z, y) %>%
      mutate(BTM = feature_name, target = "outcome")
  })) %>%
    select(BTM, target, n, r, p_value) %>%
    arrange(p_value)

  btm_cor_comp1 <- bind_rows(lapply(btm_features, function(feature_name) {
    z <- as.numeric(btm_sub[[feature_name]])
    safe_cor_test(z, score_df$Comp1) %>%
      mutate(BTM = feature_name, target = "tPLS_Comp1")
  })) %>%
    select(BTM, target, n, r, p_value) %>%
    arrange(p_value)

  btm_cor_all <- bind_rows(btm_cor_outcome, btm_cor_comp1)
  readr::write_csv(
    btm_cor_all,
    file.path(table_dir, paste0(outcome_label, "_BTM_correlations.csv"))
  )

  srebf_btm_cor <- btm_cor_all %>%
    filter(grepl("M178|SREBF1|sterol|cholesterol|lipid", BTM, ignore.case = TRUE)) %>%
    mutate(
      module_time = ifelse(
        grepl("^D[0-9]+/0", BTM),
        sub("^(D[0-9]+/0).*", "\\1", BTM),
        "D0"
      )
    ) %>%
    arrange(BTM, target)

  readr::write_csv(
    srebf_btm_cor,
    file.path(table_dir, paste0(outcome_label, "_SREBF1_lipid_BTM_correlations.csv"))
  )

  # tPLS-DA with median-split response, useful as an intuitive supervised plot.
  da_res <- try(
    tplsda_fun(
      x = X_scaled,
      y = score_df$response_group,
      ncomp = 2,
      center = FALSE
    ),
    silent = TRUE
  )

  if (!inherits(da_res, "try-error")) {
    da_score_df <- score_df %>%
      mutate(
        DA1 = da_res$x_projected[, 1],
        DA2 = da_res$x_projected[, 2]
      )
    readr::write_csv(
      da_score_df,
      file.path(table_dir, paste0(outcome_label, "_tplsda_high_low_scores.csv"))
    )
  }

  # Figures --------------------------------------------------------------------

  p_scores <- ggplot(score_df, aes(x = Comp1, y = Comp2, color = outcome)) +
    geom_point(aes(shape = age_group), size = 3, alpha = 0.9) +
    ggrepel::geom_text_repel(
      data = score_df %>%
        arrange(desc(abs(outcome - median(outcome, na.rm = TRUE)))) %>%
        slice_head(n = 8),
      aes(label = subject_id),
      size = 3,
      max.overlaps = Inf,
      show.legend = FALSE
    ) +
    scale_color_gradient(low = "#2b8cbe", high = "#d95f0e") +
    theme_bw(base_size = 12) +
    labs(
      title = paste0("tensorOmics tPLS: ", outcome_label),
      subtitle = "Transcriptome tensor uses D1-D0, D3-D0 and D7-D0 changes",
      x = "tPLS Component 1",
      y = "tPLS Component 2",
      color = outcome_label,
      shape = "Age group"
    )

  ggsave(
    file.path(figure_dir, paste0(outcome_label, "_01_tpls_scores.png")),
    p_scores,
    width = 8,
    height = 6,
    dpi = 320
  )

  p_corr <- ggplot(score_df, aes(x = Comp1, y = outcome)) +
    geom_point(aes(color = age_group), size = 3, alpha = 0.9) +
    geom_smooth(method = "lm", se = TRUE, color = "black") +
    theme_bw(base_size = 12) +
    labs(
      title = paste0("Association between tPLS Component 1 and ", outcome_label),
      subtitle = paste0(
        "r = ", round(comp_cor$r[comp_cor$component == "Comp1"], 3),
        ", p = ", signif(comp_cor$p_value[comp_cor$component == "Comp1"], 3),
        "; permutation p = ", signif(perm_summary$permutation_p, 3)
      ),
      x = "tPLS Component 1",
      y = outcome_label,
      color = "Age group"
    )

  ggsave(
    file.path(figure_dir, paste0(outcome_label, "_02_comp1_outcome_correlation.png")),
    p_corr,
    width = 7,
    height = 5,
    dpi = 320
  )

  if (!inherits(da_res, "try-error")) {
    p_da <- ggplot(da_score_df, aes(x = DA1, y = DA2, color = response_group)) +
      geom_point(aes(shape = age_group), size = 3, alpha = 0.9) +
      theme_bw(base_size = 12) +
      scale_color_manual(values = c(low = "#2b8cbe", high = "#d95f0e")) +
      labs(
        title = paste0("tPLS-DA median split: ", outcome_label),
        subtitle = "This plot is only an intuitive high/low response view",
        x = "tPLS-DA Component 1",
        y = "tPLS-DA Component 2",
        color = "Response group",
        shape = "Age group"
      )

    ggsave(
      file.path(figure_dir, paste0(outcome_label, "_03_tplsda_high_low.png")),
      p_da,
      width = 7,
      height = 5,
      dpi = 320
    )
  }

  srebf_m178_plot <- srebf_btm_cor %>%
    filter(grepl("M178|SREBF1", BTM, ignore.case = TRUE)) %>%
    mutate(
      module_time = factor(module_time, levels = c("D1/0", "D3/0", "D7/0")),
      target = factor(target, levels = c("outcome", "tPLS_Comp1"))
    )

  if (nrow(srebf_m178_plot) > 0) {
    p_srebf <- ggplot(
      srebf_m178_plot,
      aes(x = module_time, y = r, fill = target)
    ) +
      geom_col(position = position_dodge(width = 0.75), width = 0.65) +
      geom_hline(yintercept = 0, linewidth = 0.4) +
      theme_bw(base_size = 12) +
      scale_fill_manual(values = c(outcome = "#7b3294", tPLS_Comp1 = "#008837")) +
      labs(
        title = paste0("SREBF1 target BTM M178 correlation: ", outcome_label),
        subtitle = "Comparison with the SREBF1/sterol-metabolism signal from Li et al. 2017",
        x = "BTM time contrast",
        y = "Pearson correlation",
        fill = "Target"
      )

    ggsave(
      file.path(figure_dir, paste0(outcome_label, "_04_SREBF1_M178_BTM_correlations.png")),
      p_srebf,
      width = 7,
      height = 5,
      dpi = 320
    )
  }

  top_genes_for_plot <- loading_df %>%
    arrange(desc(abs_Comp1)) %>%
    slice_head(n = 8) %>%
    pull(gene)

  trajectory_df <- bind_rows(lapply(top_genes_for_plot, function(gene_name) {
    vals <- X_raw[, gene_name, , drop = FALSE]
    tibble(
      subject_id = rep(subjects, times = length(post_days)),
      gene = gene_name,
      time = rep(dimnames(X_raw)[[3]], each = length(subjects)),
      delta_expression = as.numeric(vals)
    )
  })) %>%
    left_join(score_df %>% select(subject_id, response_group), by = "subject_id") %>%
    group_by(gene, time, response_group) %>%
    summarise(
      mean_delta = mean(delta_expression, na.rm = TRUE),
      se_delta = sd(delta_expression, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(time = factor(time, levels = paste0("D", post_days, "-D0")))

  p_traj <- ggplot(
    trajectory_df,
    aes(x = time, y = mean_delta, group = response_group, color = response_group)
  ) +
    geom_line(linewidth = 0.8) +
    geom_point(size = 2) +
    geom_errorbar(
      aes(ymin = mean_delta - se_delta, ymax = mean_delta + se_delta),
      width = 0.12,
      linewidth = 0.4
    ) +
    facet_wrap(~ gene, scales = "free_y", ncol = 4) +
    theme_bw(base_size = 11) +
    scale_color_manual(values = c(low = "#2b8cbe", high = "#d95f0e")) +
    labs(
      title = paste0("Trajectories of top tPLS Component 1 loading genes: ", outcome_label),
      x = "Post-vaccination contrast",
      y = "Mean log2 expression change",
      color = "Response group"
    )

  ggsave(
    file.path(figure_dir, paste0(outcome_label, "_05_top_loading_gene_trajectories.png")),
    p_traj,
    width = 10,
    height = 6,
    dpi = 320
  )

  summary_row <- tibble(
    outcome = outcome_label,
    outcome_column = outcome_col,
    n_subjects = length(subjects),
    n_features = dim(X_scaled)[2],
    n_time_points = dim(X_scaled)[3],
    comp1_r = comp_cor$r[comp_cor$component == "Comp1"],
    comp1_p = comp_cor$p_value[comp_cor$component == "Comp1"],
    comp2_r = comp_cor$r[comp_cor$component == "Comp2"],
    comp2_p = comp_cor$p_value[comp_cor$component == "Comp2"],
    permutation_p = perm_summary$permutation_p,
    tpls_cv_r = cv_summary$cv_r,
    tpls_cv_r2 = cv_summary$cv_r2,
    matrix_pls_cv_r = matrix_pls_summary$cv_r,
    matrix_pls_cv_r2 = matrix_pls_summary$cv_r2,
    top_comp1_gene_1 = loading_df$gene[1],
    top_comp1_gene_2 = loading_df$gene[2],
    top_comp1_gene_3 = loading_df$gene[3]
  )

  summary_row
}

###############################################################################
# Run IgG and Tfh response analyses
###############################################################################

summary_primary <- run_response_analysis(
  outcome_label = "IgG_D30_0",
  outcome_table = igg_wide,
  outcome_col = primary_outcome
)

summary_tfh <- run_response_analysis(
  outcome_label = "TFH_CXCR3pos_D7_0",
  outcome_table = tfh_wide,
  outcome_col = secondary_outcome
)

all_summary <- bind_rows(summary_primary, summary_tfh)
readr::write_csv(all_summary, file.path(table_dir, "response_analysis_summary.csv"))

###############################################################################
# Dataset summary
###############################################################################

summary_txt <- file.path(table_dir, "dataset_summary_response_analysis.txt")
sink(summary_txt)
cat("GSE79396 tPLS response analysis\n")
cat("================================\n\n")
cat("Supplemental package directory:\n", supp_dir, "\n\n", sep = "")
cat("Transcriptome table:", expr_path, "\n")
cat("Expression matrix dimensions:", nrow(expr_mat), "genes x", ncol(expr_mat), "samples\n")
cat("Complete transcriptome subjects:", length(complete_subjects), "\n")
cat("Tensor time contrasts:", paste(dimnames(X_delta_all)[[3]], collapse = ", "), "\n")
cat("Features retained:", length(features_keep), "\n")
cat("Focus SREBF1/M178/sterol genes retained:", sum(features_keep %in% focus_genes), "\n")
cat("IgG subjects overlapping complete transcriptome:",
    sum(complete_subjects %in% igg_wide$subject_id), "\n")
cat("Tfh subjects overlapping complete transcriptome:",
    sum(complete_subjects %in% tfh_wide$subject_id), "\n\n")
cat("Main outputs are in:\n")
cat("Tables:", table_dir, "\n")
cat("Figures:", figure_dir, "\n\n")
print(all_summary)
sink()

message("Done.")
message("Tables: ", table_dir)
message("Figures: ", figure_dir)
message("Dataset summary: ", summary_txt)
