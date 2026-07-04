

options(stringsAsFactors = FALSE)
set.seed(1234)

###############################################################################
# User settings
###############################################################################

gse_id <- "GSE79396"
geo_matrix_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE79nnn/GSE79396/matrix/GSE79396_series_matrix.txt.gz"

days_keep <- c(0, 1, 3, 7)
top_n_features <- 2000
ncomp_tpca <- 5
n_perm <- 999

base_dir <- getwd()
data_dir <- file.path(base_dir, "data", "GSE79396")
result_dir <- file.path(base_dir, "results", "GSE79396_tPCA")
figure_dir <- file.path(result_dir, "figures")
table_dir <- file.path(result_dir, "tables")

dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

matrix_path <- file.path(data_dir, "GSE79396_series_matrix.txt.gz")

# Package installation/loading

install_if_missing <- function(pkgs) {
  missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    install.packages(missing_pkgs, repos = "https://cloud.r-project.org")
  }
}

install_if_missing(c(
  "data.table", "dplyr", "tidyr", "ggplot2", "ggrepel",
  "patchwork", "scales", "matrixStats", "tibble", "readr", "remotes"
))

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(scales)
  library(matrixStats)
  library(tibble)
  library(readr)
})

ensure_biocparallel <- function() {
  if (!requireNamespace("BiocParallel", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org")
    }
    BiocManager::install("BiocParallel", ask = FALSE, update = FALSE)
  }
}

ensure_tensoromics <- function() {
  if (requireNamespace("tensorOmics", quietly = TRUE)) {
    message("Using installed tensorOmics package.")
    return(tensorOmics::tpca)
  }

  message("tensorOmics package not found. Trying remotes::install_github().")
  ensure_biocparallel()
  install_if_missing(c("gsignal", "rARPACK", "MASS"))

  ok <- FALSE
  try({
    remotes::install_github(
      "brendanlu/tensorOmics",
      dependencies = TRUE,
      upgrade = "never"
    )
    ok <- requireNamespace("tensorOmics", quietly = TRUE)
  }, silent = TRUE)

  if (ok) {
    message("Using tensorOmics installed from GitHub.")
    return(tensorOmics::tpca)
  }

  message("GitHub installation failed or R version is too old.")
  message("Fallback: sourcing the author's R source files directly from GitHub.")
  ensure_biocparallel()
  install_if_missing(c("gsignal", "rARPACK", "MASS"))

  source_urls <- c(
    "https://raw.githubusercontent.com/brendanlu/tensorOmics/main/R/vendor.R",
    "https://raw.githubusercontent.com/brendanlu/tensorOmics/main/R/names.R",
    "https://raw.githubusercontent.com/brendanlu/tensorOmics/main/R/tens.mproduct.R",
    "https://raw.githubusercontent.com/brendanlu/tensorOmics/main/R/tens.tsvdm.R",
    "https://raw.githubusercontent.com/brendanlu/tensorOmics/main/R/tens.tpca.R"
  )

  for (u in source_urls) {
    message("Sourcing: ", u)
    source(u)
  }

  if (!exists("tpca")) {
    stop("Could not load tpca() from tensorOmics source files.")
  }
  tpca
}

tpca_fun <- ensure_tensoromics()

# Download and parse GEO series matrix

if (!file.exists(matrix_path)) {
  message("Downloading ", gse_id, " series matrix from GEO...")
  download.file(geo_matrix_url, destfile = matrix_path, mode = "wb")
} else {
  message("Using existing file: ", matrix_path)
}

clean_geo_value <- function(x) {
  x <- sub("^\"", "", x)
  x <- sub("\"$", "", x)
  x
}

read_geo_metadata <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con))

  meta <- list()
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0) break
    if (startsWith(line, "!series_matrix_table_begin")) break

    if (
      startsWith(line, "!Sample_title") ||
        startsWith(line, "!Sample_geo_accession") ||
        startsWith(line, "!Sample_source_name_ch1") ||
        startsWith(line, "!Sample_characteristics_ch1") ||
        startsWith(line, "!Series_title") ||
        startsWith(line, "!Series_summary") ||
        startsWith(line, "!Series_overall_design")
    ) {
      parts <- strsplit(line, "\t", fixed = TRUE)[[1]]
      key <- parts[1]
      values <- clean_geo_value(parts[-1])
      meta[[key]] <- append(meta[[key]], list(values))
    }
  }
  meta
}

meta <- read_geo_metadata(matrix_path)

sample_title <- meta[["!Sample_title"]][[1]]
geo_accession <- meta[["!Sample_geo_accession"]][[1]]
source_name <- meta[["!Sample_source_name_ch1"]][[1]]
characteristics <- meta[["!Sample_characteristics_ch1"]]

extract_characteristic <- function(prefix) {
  out <- rep(NA_character_, length(sample_title))
  for (row in characteristics) {
    hit <- startsWith(row, prefix)
    out[hit] <- trimws(sub(prefix, "", row[hit], fixed = TRUE))
  }
  out
}

sample_info <- tibble(
  geo_accession = geo_accession,
  sample_title = sample_title,
  source_name = source_name,
  subject_id = extract_characteristic("subject id:"),
  vaccine = extract_characteristic("vaccine:"),
  cohort = extract_characteristic("cohort:"),
  visit = extract_characteristic("visit:"),
  sex = extract_characteristic("Sex:"),
  age = as.numeric(extract_characteristic("age:")),
  age_group = extract_characteristic("age group:"),
  tissue = extract_characteristic("tissue/cell type:")
) %>%
  mutate(
    day = as.integer(sub(".*_D([0-9]+)_.*", "\\1", sample_title)),
    age_group = factor(age_group, levels = c("young", "elderly")),
    sex = factor(sex),
    cohort = factor(cohort)
  )

if (any(is.na(sample_info$subject_id)) || any(is.na(sample_info$day))) {
  stop("Could not parse subject_id or day from GEO metadata.")
}

message("Samples in GEO matrix metadata: ", nrow(sample_info))
message("Unique subjects: ", length(unique(sample_info$subject_id)))
message("Time points: ", paste(sort(unique(sample_info$day)), collapse = ", "))

# Read expression table

message("Reading expression matrix. This can take a few minutes...")
expr_dt <- data.table::fread(
  matrix_path,
  skip = "!series_matrix_table_begin",
  data.table = TRUE,
  check.names = FALSE,
  showProgress = TRUE
)

if (!("ID_REF" %in% names(expr_dt))) {
  setnames(expr_dt, old = names(expr_dt)[1], new = "ID_REF")
}

expr_dt <- expr_dt[!grepl("^!", ID_REF)]
sample_cols <- intersect(sample_info$geo_accession, names(expr_dt))
if (length(sample_cols) != nrow(sample_info)) {
  warning("Not all metadata samples were found as expression columns.")
}

probe_ids <- expr_dt$ID_REF
expr_mat <- as.matrix(expr_dt[, ..sample_cols])
mode(expr_mat) <- "numeric"
rownames(expr_mat) <- probe_ids

rm(expr_dt)
gc()

# Select complete subjects and construct subject x feature x time tensor
sample_info_complete_time <- sample_info %>%
  filter(day %in% days_keep, geo_accession %in% sample_cols)

complete_subjects <- sample_info_complete_time %>%
  dplyr::group_by(subject_id) %>%
  dplyr::summarise(n_time = dplyr::n_distinct(day), .groups = "drop") %>%
  dplyr::filter(n_time == length(days_keep)) %>%
  dplyr::pull(subject_id)

sample_info_complete <- sample_info_complete_time %>%
  filter(subject_id %in% complete_subjects) %>%
  arrange(subject_id, day)

subject_meta <- sample_info_complete %>%
  group_by(subject_id) %>%
  summarise(
    age_group = first(age_group),
    age = first(age),
    sex = first(sex),
    cohort = first(cohort),
    .groups = "drop"
  ) %>%
  arrange(subject_id)

subject_order <- subject_meta$subject_id

message("Complete subjects for days ", paste(days_keep, collapse = ", "), ": ",
        length(subject_order))
message("Age groups among complete subjects:")
print(table(subject_meta$age_group))

# Feature filtering and scaling
message("Filtering probes with missing values and selecting top variable probes...")

expr_subset <- expr_mat[, sample_info_complete$geo_accession, drop = FALSE]
finite_keep <- matrixStats::rowAlls(is.finite(expr_subset))
expr_subset <- expr_subset[finite_keep, , drop = FALSE]

probe_vars <- matrixStats::rowVars(expr_subset)
probe_vars[!is.finite(probe_vars)] <- 0
names(probe_vars) <- rownames(expr_subset)

top_n <- min(top_n_features, nrow(expr_subset))
top_probes <- names(sort(probe_vars, decreasing = TRUE))[seq_len(top_n)]
expr_top <- expr_subset[top_probes, , drop = FALSE]

probe_mean <- rowMeans(expr_top)
probe_sd <- matrixStats::rowSds(expr_top)
probe_sd[probe_sd == 0 | !is.finite(probe_sd)] <- 1
expr_scaled <- sweep(expr_top, 1, probe_mean, "-")
expr_scaled <- sweep(expr_scaled, 1, probe_sd, "/")
expr_scaled[!is.finite(expr_scaled)] <- 0

message("Selected probes: ", nrow(expr_scaled))

# Construct tensor: subjects x probes x tim
X_tensor <- array(
  NA_real_,
  dim = c(length(subject_order), nrow(expr_scaled), length(days_keep)),
  dimnames = list(
    subject_order,
    rownames(expr_scaled),
    paste0("D", days_keep)
  )
)

for (k in seq_along(days_keep)) {
  d <- days_keep[k]
  sample_order_day <- sample_info_complete %>%
    filter(day == d) %>%
    arrange(match(subject_id, subject_order))

  if (!identical(sample_order_day$subject_id, subject_order)) {
    stop("Subject order mismatch while constructing tensor.")
  }

  X_tensor[, , k] <- t(expr_scaled[, sample_order_day$geo_accession, drop = FALSE])
}

# Matrix PCA baseline
message("Running classical matrix PCA baseline...")

matrix_sample_info <- sample_info_complete %>%
  arrange(day, subject_id)

X_matrix_samples <- t(expr_scaled[, matrix_sample_info$geo_accession, drop = FALSE])

pca_matrix <- prcomp(
  X_matrix_samples,
  center = FALSE,
  scale. = FALSE,
  rank. = 5
)

matrix_scores <- as_tibble(pca_matrix$x[, 1:5, drop = FALSE]) %>%
  setNames(paste0("PC", seq_len(5))) %>%
  bind_cols(matrix_sample_info, .)

matrix_subject_scores <- matrix_scores %>%
  group_by(subject_id, age_group, age, sex, cohort) %>%
  summarise(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    PC3 = mean(PC3),
    .groups = "drop"
  )

# Tensor PCA using author's tensorOmics tpca interface
message("Running tensorOmics tPCA...")

tpca_res <- tpca_fun(
  X_tensor,
  ncomp = ncomp_tpca,
  center = TRUE,
  matrix_output = TRUE
)

tpca_scores <- as.data.frame(tpca_res$variates)
tpca_scores <- tpca_scores[, seq_len(min(ncomp_tpca, ncol(tpca_scores))), drop = FALSE]
colnames(tpca_scores) <- paste0("Comp", seq_len(ncol(tpca_scores)))

tpca_scores <- tpca_scores %>%
  rownames_to_column("subject_id") %>%
  left_join(subject_meta, by = "subject_id")

tpca_loadings <- as.data.frame(tpca_res$loadings)
tpca_loadings <- tpca_loadings[, seq_len(min(ncomp_tpca, ncol(tpca_loadings))), drop = FALSE]
colnames(tpca_loadings) <- paste0("Comp", seq_len(ncol(tpca_loadings)))
tpca_loadings <- tpca_loadings %>%
  rownames_to_column("probe_id")

# Statistical tests
extract_manova_pillai <- function(fit) {
  s <- summary(fit, test = "Pillai")
  stat <- s$stats[1, "Pillai"]
  pval <- s$stats[1, "Pr(>F)"]
  c(statistic = unname(stat), p_value = unname(pval))
}

permute_manova_pillai <- function(df, score_cols, group_col, n_perm = 999) {
  form <- as.formula(
    paste0("cbind(", paste(score_cols, collapse = ", "), ") ~ ", group_col)
  )
  observed <- extract_manova_pillai(manova(form, data = df))["statistic"]

  perm_stats <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    df_perm <- df
    df_perm[[group_col]] <- sample(df_perm[[group_col]])
    perm_stats[i] <- extract_manova_pillai(manova(form, data = df_perm))["statistic"]
  }

  p_perm <- (sum(perm_stats >= observed) + 1) / (n_perm + 1)
  list(statistic = observed, p_value = p_perm)
}

message("Running association tests...")

tpca_manova <- extract_manova_pillai(
  manova(cbind(Comp1, Comp2) ~ age_group, data = tpca_scores)
)
tpca_perm <- permute_manova_pillai(
  tpca_scores,
  score_cols = c("Comp1", "Comp2"),
  group_col = "age_group",
  n_perm = n_perm
)

matrix_manova <- extract_manova_pillai(
  manova(cbind(PC1, PC2) ~ age_group, data = matrix_subject_scores)
)
matrix_perm <- permute_manova_pillai(
  matrix_subject_scores,
  score_cols = c("PC1", "PC2"),
  group_col = "age_group",
  n_perm = n_perm
)

lm_tpca_comp1 <- summary(lm(Comp1 ~ age_group, data = tpca_scores))
lm_tpca_comp2 <- summary(lm(Comp2 ~ age_group, data = tpca_scores))

stat_tests <- tibble(
  analysis = c(
    "Tensor PCA scores Comp1-Comp2 ~ age_group",
    "Tensor PCA scores Comp1-Comp2 ~ age_group (permutation)",
    "Matrix PCA subject-averaged PC1-PC2 ~ age_group",
    "Matrix PCA subject-averaged PC1-PC2 ~ age_group (permutation)",
    "Tensor PCA Comp1 ~ age_group",
    "Tensor PCA Comp2 ~ age_group"
  ),
  test = c(
    "MANOVA Pillai",
    "Permutation MANOVA Pillai",
    "MANOVA Pillai",
    "Permutation MANOVA Pillai",
    "Linear model F-test",
    "Linear model F-test"
  ),
  statistic = c(
    tpca_manova["statistic"],
    tpca_perm$statistic,
    matrix_manova["statistic"],
    matrix_perm$statistic,
    lm_tpca_comp1$fstatistic[1],
    lm_tpca_comp2$fstatistic[1]
  ),
  p_value = c(
    tpca_manova["p_value"],
    tpca_perm$p_value,
    matrix_manova["p_value"],
    matrix_perm$p_value,
    pf(
      lm_tpca_comp1$fstatistic[1],
      lm_tpca_comp1$fstatistic[2],
      lm_tpca_comp1$fstatistic[3],
      lower.tail = FALSE
    ),
    pf(
      lm_tpca_comp2$fstatistic[1],
      lm_tpca_comp2$fstatistic[2],
      lm_tpca_comp2$fstatistic[3],
      lower.tail = FALSE
    )
  )
)

# Loadings and trajectory summaries
top_loadings <- function(loadings_df, comp = "Comp1", n = 30) {
  loadings_df %>%
    transmute(
      probe_id,
      loading = .data[[comp]],
      abs_loading = abs(.data[[comp]])
    ) %>%
    arrange(desc(abs_loading)) %>%
    slice_head(n = n)
}

top_pc1 <- top_loadings(tpca_loadings, "Comp1", n = 30)
top_pc2 <- top_loadings(tpca_loadings, "Comp2", n = 30)

make_feature_long <- function(feature_ids) {
  keep_samples <- sample_info_complete %>%
    arrange(day, subject_id)

  values <- t(expr_scaled[feature_ids, keep_samples$geo_accession, drop = FALSE])
  values_df <- as_tibble(values) %>%
    mutate(
      subject_id = keep_samples$subject_id,
      day = keep_samples$day,
      age_group = keep_samples$age_group
    )

  values_df %>%
    pivot_longer(
      cols = all_of(feature_ids),
      names_to = "probe_id",
      values_to = "scaled_expression"
    ) %>%
    group_by(age_group, day, probe_id) %>%
    summarise(
      mean_expression = mean(scaled_expression, na.rm = TRUE),
      .groups = "drop"
    )
}

traj_pc1 <- make_feature_long(top_pc1$probe_id[1:8])
traj_pc2 <- make_feature_long(top_pc2$probe_id[1:8])

# Figures
age_cols <- c("young" = "#00A9A5", "elderly" = "#F8766D")

p_matrix <- ggplot(matrix_scores, aes(x = PC1, y = PC2, color = age_group, shape = factor(day))) +
  geom_point(size = 2.4, alpha = 0.85) +
  scale_color_manual(values = age_cols, drop = FALSE) +
  labs(
    title = "Classical matrix PCA",
    subtitle = "Each point is one sample-time observation",
    x = "PC1",
    y = "PC2",
    color = "Age group",
    shape = "Day"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

p_tpca <- ggplot(tpca_scores, aes(x = Comp1, y = Comp2, color = age_group, shape = cohort)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(aes(label = subject_id), size = 2.7, max.overlaps = 20) +
  scale_color_manual(values = age_cols, drop = FALSE) +
  labs(
    title = "tensorOmics tPCA",
    subtitle = "Each point is one subject represented by the full D0-D1-D3-D7 trajectory",
    x = "Component 1",
    y = "Component 2",
    color = "Age group",
    shape = "Cohort"
  ) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(face = "bold"))

plot_loading_bar <- function(top_df, comp_name) {
  top_df %>%
    slice_head(n = 20) %>%
    mutate(probe_id = factor(probe_id, levels = rev(probe_id))) %>%
    ggplot(aes(x = probe_id, y = loading, fill = loading > 0)) +
    geom_col(width = 0.75) +
    coord_flip() +
    scale_fill_manual(values = c("TRUE" = "#e39894", "FALSE" = "#8fafd9"), guide = "none") +
    labs(
      title = paste0("Top tPCA loadings: ", comp_name),
      x = "Probe ID",
      y = "Loading"
    ) +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))
}

p_load_pc1 <- plot_loading_bar(top_pc1, "Component 1")
p_load_pc2 <- plot_loading_bar(top_pc2, "Component 2")

plot_trajectories <- function(traj_df, title) {
  ggplot(traj_df, aes(x = day, y = mean_expression, color = age_group, group = age_group)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~ probe_id, scales = "free_y", ncol = 4) +
    scale_color_manual(values = age_cols, drop = FALSE) +
    scale_x_continuous(breaks = days_keep) +
    labs(
      title = title,
      x = "Day post-vaccination",
      y = "Mean scaled expression",
      color = "Age group"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      strip.text = element_text(size = 8)
    )
}

p_traj_pc1 <- plot_trajectories(
  traj_pc1,
  "Mean trajectories of top Component 1 probes"
)
p_traj_pc2 <- plot_trajectories(
  traj_pc2,
  "Mean trajectories of top Component 2 probes"
)

ggsave(file.path(figure_dir, "01_matrix_pca_samples.png"), p_matrix, width = 8, height = 6, dpi = 300)
ggsave(file.path(figure_dir, "02_tensor_pca_subjects.png"), p_tpca, width = 8, height = 6, dpi = 300)
ggsave(file.path(figure_dir, "03_tensor_pc1_loadings.png"), p_load_pc1, width = 8, height = 6, dpi = 300)
ggsave(file.path(figure_dir, "04_tensor_pc2_loadings.png"), p_load_pc2, width = 8, height = 6, dpi = 300)
ggsave(file.path(figure_dir, "05_tensor_pc1_top_feature_trajectories.png"), p_traj_pc1, width = 10, height = 7, dpi = 300)
ggsave(file.path(figure_dir, "06_tensor_pc2_top_feature_trajectories.png"), p_traj_pc2, width = 10, height = 7, dpi = 300)

# Save tables and summary
write_csv(sample_info_complete, file.path(table_dir, "sample_metadata_complete.csv"))
write_csv(tpca_scores, file.path(table_dir, "tensor_pca_scores.csv"))
write_csv(matrix_scores, file.path(table_dir, "matrix_pca_scores.csv"))
write_csv(top_pc1, file.path(table_dir, "top_loadings_pc1.csv"))
write_csv(top_pc2, file.path(table_dir, "top_loadings_pc2.csv"))
write_csv(stat_tests, file.path(table_dir, "statistical_tests.csv"))

summary_lines <- c(
  paste0("Dataset: ", gse_id),
  "Title: Integrated transcriptomics and metabolomics profiling delineates early molecular correlates of immunity to herpes zoster vaccination in humans",
  paste0("Samples in metadata: ", nrow(sample_info)),
  paste0("Complete subjects used: ", length(subject_order)),
  paste0("Time points used: ", paste(days_keep, collapse = ", ")),
  paste0("Selected top variable probes: ", nrow(expr_scaled)),
  "",
  "Complete subject age-group counts:",
  capture.output(print(table(subject_meta$age_group))),
  "",
  "Statistical tests:",
  capture.output(print(stat_tests)),
  "",
  "Output folders:",
  paste0("Figures: ", figure_dir),
  paste0("Tables: ", table_dir)
)

writeLines(summary_lines, con = file.path(table_dir, "dataset_summary.txt"))

message("Analysis complete.")
message("Figures saved to: ", figure_dir)
message("Tables saved to: ", table_dir)
message("Please send the figures and statistical_tests.csv when you want to write the report.")
