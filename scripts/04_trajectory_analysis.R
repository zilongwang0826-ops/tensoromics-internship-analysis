options(stringsAsFactors = FALSE)
analysis_seed <- 2026
set.seed(analysis_seed)

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(matrixStats)
  library(readr)
  library(tibble)
  library(scales)
})

# Independent trajectory analyses from the local GEO files.

root_dir <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out_dir <- file.path(root_dir, "results", "trajectory_analysis")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_bw(base_size = 11))

clean_geo_value <- function(x) {
  x <- sub("^\"", "", x)
  sub("\"$", "", x)
}

read_geo_metadata <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con))

  meta <- list()
  repeat {
    line <- readLines(con, n = 1, warn = FALSE)
    if (length(line) == 0 || startsWith(line, "!series_matrix_table_begin")) break

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

extract_characteristic <- function(characteristics, prefix, n) {
  out <- rep(NA_character_, n)
  for (row in characteristics) {
    hit <- startsWith(row, prefix)
    out[hit] <- trimws(sub(prefix, "", row[hit], fixed = TRUE))
  }
  out
}

safe_cor <- function(x, y, method = "pearson") {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 4 || sd(x[ok]) == 0 || sd(y[ok]) == 0) {
    return(tibble(n = sum(ok), r = NA_real_, p_value = NA_real_))
  }
  ct <- suppressWarnings(cor.test(x[ok], y[ok], method = method))
  tibble(n = sum(ok), r = unname(ct$estimate), p_value = ct$p.value)
}

cohens_d <- function(x, g, positive_group) {
  g <- as.character(g)
  x1 <- x[g == positive_group]
  x0 <- x[g != positive_group]
  sp <- sqrt(((length(x1) - 1) * var(x1) + (length(x0) - 1) * var(x0)) /
    (length(x1) + length(x0) - 2))
  (mean(x1) - mean(x0)) / sp
}

pls1_components <- function(X, y, ncomp = 2) {
  X <- scale(X)
  y <- as.numeric(scale(y))
  X[!is.finite(X)] <- 0
  y[!is.finite(y)] <- 0

  X_work <- X
  y_work <- y
  scores <- matrix(NA_real_, nrow(X), ncomp)
  loadings <- matrix(NA_real_, ncol(X), ncomp)
  rownames(scores) <- rownames(X)
  rownames(loadings) <- colnames(X)

  for (h in seq_len(ncomp)) {
    w <- as.numeric(crossprod(X_work, y_work))
    if (!any(is.finite(w)) || sqrt(sum(w^2)) == 0) {
      stop("PLS component could not be estimated; check input variance.")
    }
    w <- w / sqrt(sum(w^2))
    t_score <- as.numeric(X_work %*% w)
    denom <- sum(t_score^2)
    p <- as.numeric(crossprod(X_work, t_score) / denom)
    q <- sum(y_work * t_score) / denom

    scores[, h] <- t_score
    loadings[, h] <- w
    X_work <- X_work - tcrossprod(t_score, p)
    y_work <- y_work - q * t_score
  }

  colnames(scores) <- paste0("Comp", seq_len(ncomp))
  colnames(loadings) <- paste0("Comp", seq_len(ncomp))
  list(scores = scores, loadings = loadings)
}

extract_pillai <- function(df, score_cols, group_col) {
  form <- as.formula(paste0("cbind(", paste(score_cols, collapse = ", "), ") ~ ", group_col))
  st <- summary(manova(form, data = df), test = "Pillai")$stats
  tibble(statistic = unname(st[1, "Pillai"]), p_value = unname(st[1, "Pr(>F)"]))
}

permute_pillai <- function(df, score_cols, group_col, n_perm = 499) {
  obs <- extract_pillai(df, score_cols, group_col)$statistic
  stats <- numeric(n_perm)
  for (i in seq_len(n_perm)) {
    tmp <- df
    tmp[[group_col]] <- sample(tmp[[group_col]])
    stats[i] <- extract_pillai(tmp, score_cols, group_col)$statistic
  }
  tibble(
    statistic = obs,
    p_value = (sum(stats >= obs) + 1) / (n_perm + 1),
    n_perm = n_perm
  )
}

download_if_missing <- function(url, path) {
  if (!file.exists(path)) {
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    download.file(url, destfile = path, mode = "wb", quiet = FALSE)
  }
}

read_platform_annotation <- function(path) {
  annot <- data.table::fread(
    path,
    skip = 27,
    sep = "\t",
    data.table = FALSE,
    fill = TRUE,
    quote = "",
    check.names = FALSE
  )
  names(annot)[1] <- "probe_id"
  annot %>%
    filter(!startsWith(probe_id, "!")) %>%
    transmute(
      probe_id,
      gene_title = .data[["Gene title"]],
      gene_symbol = .data[["Gene symbol"]],
      gene_id = .data[["Gene ID"]],
      go_process = .data[["GO:Process"]],
      go_function = .data[["GO:Function"]]
    )
}

normalize_zv_id <- function(x) {
  gsub("^ZV([0-9]+)-0*([0-9]+)$", "ZV\\1\\2", as.character(x))
}

transpose_feature_table <- function(path, first_col_name) {
  dt <- data.table::fread(path, data.table = FALSE, check.names = FALSE)
  names(dt)[1] <- first_col_name
  feature_names <- make.unique(as.character(dt[[1]]))
  mat <- as.matrix(dt[, -1, drop = FALSE])
  storage.mode(mat) <- "numeric"
  rownames(mat) <- feature_names
  out <- as.data.frame(t(mat), check.names = FALSE)
  out$subject_id <- rownames(out)
  rownames(out) <- NULL
  out
}

message("Fish aging trajectory analysis")

fish_counts_path <- file.path(root_dir, "data", "GSE150318", "GSE150318_counts.csv.gz")
fish_meta_path <- file.path(root_dir, "data", "GSE150318", "GSE150318_series_matrix.txt.gz")
download_if_missing(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE150nnn/GSE150318/matrix/GSE150318_series_matrix.txt.gz",
  fish_meta_path
)

fish_meta_geo <- read_geo_metadata(fish_meta_path)
fish_counts_dt <- data.table::fread(fish_counts_path, data.table = FALSE, check.names = FALSE)
names(fish_counts_dt)[1] <- "gene_id"
fish_counts <- as.matrix(fish_counts_dt[, -1, drop = FALSE])
mode(fish_counts) <- "numeric"
rownames(fish_counts) <- fish_counts_dt$gene_id
rm(fish_counts_dt)
gc()

fish_sample_title <- fish_meta_geo[["!Sample_title"]][[1]]
fish_geo_accession <- fish_meta_geo[["!Sample_geo_accession"]][[1]]
fish_characteristics <- fish_meta_geo[["!Sample_characteristics_ch1"]]

fish_sample_meta <- tibble(
  sample_title = fish_sample_title,
  geo_accession = fish_geo_accession,
  fish_id = extract_characteristic(fish_characteristics, "animal id:", length(fish_sample_title)),
  sex = extract_characteristic(fish_characteristics, "Sex:", length(fish_sample_title)),
  tissue = extract_characteristic(fish_characteristics, "tissue:", length(fish_sample_title)),
  time_week = as.numeric(extract_characteristic(
    fish_characteristics,
    "age at fin clip (weeks):",
    length(fish_sample_title)
  )),
  age_death_days = suppressWarnings(as.numeric(extract_characteristic(
    fish_characteristics,
    "age at death (days):",
    length(fish_sample_title)
  )))
) %>%
  mutate(sample_name = paste0(fish_id, "_", time_week, "w")) %>%
  relocate(sample_name, .before = sample_title)

if (!all(fish_sample_meta$sample_name %in% colnames(fish_counts))) {
  stop("Fish sample names from GEO metadata do not match the count matrix.")
}

fish_info <- fish_sample_meta %>%
  filter(is.finite(age_death_days)) %>%
  distinct(fish_id, age_death_days) %>%
  arrange(fish_id)

median_death <- median(fish_info$age_death_days)
fish_info <- fish_info %>%
  mutate(
    lifespan_group = ifelse(age_death_days <= median_death, "short_lived", "long_lived"),
    lifespan_group = factor(lifespan_group, levels = c("short_lived", "long_lived"))
  )

fish_complete <- fish_sample_meta %>%
  filter(fish_id %in% fish_info$fish_id, time_week %in% c(10, 20)) %>%
  group_by(fish_id) %>%
  filter(n_distinct(time_week) == 2) %>%
  ungroup() %>%
  arrange(fish_id, time_week) %>%
  left_join(fish_info, by = "fish_id")

counts_valid <- fish_counts[, fish_complete$sample_name, drop = FALSE]
keep_genes <- rowSums(counts_valid) > 10
counts_valid <- counts_valid[keep_genes, , drop = FALSE]

cpm <- t(t(counts_valid) / colSums(counts_valid)) * 1e6
logcpm <- log2(cpm + 1)
gene_var <- matrixStats::rowVars(logcpm)
names(gene_var) <- rownames(logcpm)
top_fish_genes <- names(sort(gene_var, decreasing = TRUE))[seq_len(1000)]
logcpm_top <- logcpm[top_fish_genes, , drop = FALSE]

fish_order <- fish_info$fish_id
time_order <- c(10, 20)
X_fish_10 <- matrix(NA_real_, length(fish_order), length(top_fish_genes),
  dimnames = list(fish_order, top_fish_genes)
)
X_fish_20 <- X_fish_10

for (i in seq_along(fish_order)) {
  sid <- fish_order[i]
  s10 <- fish_complete$sample_name[fish_complete$fish_id == sid & fish_complete$time_week == 10]
  s20 <- fish_complete$sample_name[fish_complete$fish_id == sid & fish_complete$time_week == 20]
  X_fish_10[i, ] <- logcpm_top[, s10]
  X_fish_20[i, ] <- logcpm_top[, s20]
}

X_fish_10_raw <- X_fish_10
X_fish_20_raw <- X_fish_20
X_fish_delta_raw <- X_fish_20 - X_fish_10
X_fish_delta <- X_fish_delta_raw
colnames(X_fish_10) <- paste0(colnames(X_fish_10), "|w10")
colnames(X_fish_20) <- paste0(colnames(X_fish_20), "|w20")
colnames(X_fish_delta) <- paste0(colnames(X_fish_delta), "|delta20minus10")
X_fish_traj <- cbind(X_fish_10, X_fish_20, X_fish_delta)

group_numeric <- ifelse(fish_info$lifespan_group == "long_lived", 1, 0)
fish_pls <- pls1_components(X_fish_traj, group_numeric, ncomp = 2)
fish_scores <- as_tibble(fish_pls$scores, rownames = "fish_id") %>%
  left_join(fish_info, by = "fish_id")

# Orient components toward longer lifespan for plotting and interpretation.
if (cor(fish_scores$Comp1, fish_scores$age_death_days) < 0) {
  fish_scores$Comp1 <- -fish_scores$Comp1
  fish_pls$loadings[, "Comp1"] <- -fish_pls$loadings[, "Comp1"]
}
if (cor(fish_scores$Comp2, fish_scores$age_death_days) < 0) {
  fish_scores$Comp2 <- -fish_scores$Comp2
  fish_pls$loadings[, "Comp2"] <- -fish_pls$loadings[, "Comp2"]
}

fish_pillai <- extract_pillai(fish_scores, c("Comp1", "Comp2"), "lifespan_group") %>%
  mutate(analysis = "trajectory PLS Comp1-Comp2 ~ lifespan group", .before = 1)
fish_perm <- permute_pillai(fish_scores, c("Comp1", "Comp2"), "lifespan_group", n_perm = 499) %>%
  mutate(analysis = "trajectory PLS label permutation", .before = 1)
fish_lm <- summary(lm(age_death_days ~ Comp1 + Comp2, data = fish_scores))
fish_lm_row <- tibble(
  analysis = "age at death ~ trajectory PLS Comp1 + Comp2",
  statistic = fish_lm$r.squared,
  p_value = pf(
    fish_lm$fstatistic[1],
    fish_lm$fstatistic[2],
    fish_lm$fstatistic[3],
    lower.tail = FALSE
  ),
  n_perm = NA_integer_
)

fish_component_cor <- bind_rows(
  safe_cor(fish_scores$Comp1, fish_scores$age_death_days) %>% mutate(component = "Comp1"),
  safe_cor(fish_scores$Comp2, fish_scores$age_death_days) %>% mutate(component = "Comp2")
) %>%
  select(component, everything())

fish_feature_loadings <- as_tibble(fish_pls$loadings, rownames = "feature") %>%
  mutate(
    gene_id = sub("\\|.*$", "", feature),
    time_feature = sub("^.*\\|", "", feature)
  )

fish_gene_loadings <- fish_feature_loadings %>%
  group_by(gene_id) %>%
  summarise(
    loading_comp1_rss = sqrt(sum(Comp1^2)),
    loading_comp2_rss = sqrt(sum(Comp2^2)),
    loading_max = max(abs(c(Comp1, Comp2))),
    top_time_feature = time_feature[which.max(abs(Comp1))],
    .groups = "drop"
  ) %>%
  arrange(desc(loading_comp1_rss))

fish_delta_tests <- tibble(gene_id = top_fish_genes) %>%
  rowwise() %>%
  mutate(
    mean_delta_short = mean(X_fish_delta_raw[fish_info$lifespan_group == "short_lived", gene_id]),
    mean_delta_long = mean(X_fish_delta_raw[fish_info$lifespan_group == "long_lived", gene_id]),
    difference_long_minus_short = mean_delta_long - mean_delta_short,
    p_value = t.test(
      X_fish_delta_raw[fish_info$lifespan_group == "long_lived", gene_id],
      X_fish_delta_raw[fish_info$lifespan_group == "short_lived", gene_id]
    )$p.value
  ) %>%
  ungroup() %>%
  mutate(fdr = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

top_loading50 <- fish_gene_loadings$gene_id[seq_len(50)]
top_delta50 <- fish_delta_tests$gene_id[seq_len(50)]
fish_overlap <- tibble(
  top_loading_n = 50,
  top_delta_n = 50,
  universe_n = length(top_fish_genes),
  overlap_n = length(intersect(top_loading50, top_delta50)),
  expected_overlap = 50 * 50 / length(top_fish_genes),
  p_value = phyper(
    length(intersect(top_loading50, top_delta50)) - 1,
    50,
    length(top_fish_genes) - 50,
    50,
    lower.tail = FALSE
  )
)

boot_n <- 150
boot_top <- vector("list", boot_n)
for (b in seq_len(boot_n)) {
  idx <- sample(seq_len(nrow(X_fish_traj)), replace = TRUE)
  boot_fit <- pls1_components(X_fish_traj[idx, , drop = FALSE], group_numeric[idx], ncomp = 1)
  boot_load <- as_tibble(boot_fit$loadings, rownames = "feature") %>%
    mutate(gene_id = sub("\\|.*$", "", feature)) %>%
    group_by(gene_id) %>%
    summarise(score = sqrt(sum(Comp1^2)), .groups = "drop") %>%
    arrange(desc(score)) %>%
    slice_head(n = 50)
  boot_top[[b]] <- boot_load$gene_id
}

fish_boot_stability <- tibble(gene_id = unlist(boot_top)) %>%
  dplyr::count(gene_id, name = "top50_count") %>%
  mutate(top50_frequency = top50_count / boot_n) %>%
  arrange(desc(top50_frequency)) %>%
  left_join(fish_gene_loadings, by = "gene_id")

write_csv(fish_sample_meta, file.path(tab_dir, "fish_geo_sample_metadata.csv"))
write_csv(fish_scores, file.path(tab_dir, "fish_trajectory_pls_scores.csv"))
write_csv(bind_rows(fish_pillai, fish_perm, fish_lm_row), file.path(tab_dir, "fish_trajectory_statistical_tests.csv"))
write_csv(fish_component_cor, file.path(tab_dir, "fish_component_lifespan_correlations.csv"))
write_csv(fish_gene_loadings, file.path(tab_dir, "fish_trajectory_pls_gene_loadings.csv"))
write_csv(fish_delta_tests, file.path(tab_dir, "fish_delta_group_tests.csv"))
write_csv(fish_overlap, file.path(tab_dir, "fish_loading_delta_overlap.csv"))
write_csv(fish_boot_stability, file.path(tab_dir, "fish_bootstrap_loading_stability.csv"))

fish_age_plot <- fish_scores %>%
  pivot_longer(c(Comp1, Comp2), names_to = "component", values_to = "score") %>%
  ggplot(aes(score, age_death_days, color = lifespan_group)) +
  geom_point(size = 2.4, alpha = 0.9) +
  geom_smooth(method = "lm", se = TRUE, color = "grey25", linewidth = 0.7) +
  facet_wrap(~component, scales = "free_x") +
  scale_color_manual(values = c(short_lived = "#D95F02", long_lived = "#1B9E77")) +
  labs(
    title = "Trajectory PLS of fish lifespan signal",
    x = "Component score",
    y = "Age at death (days)",
    color = "Lifespan group"
  )

ggsave(file.path(fig_dir, "fish_lifespan_component_correlations.png"), fish_age_plot,
  width = 8.5, height = 4.3, dpi = 300
)

fish_heatmap_genes <- fish_delta_tests$gene_id[seq_len(24)]
fish_heatmap <- bind_rows(
  as_tibble(X_fish_10_raw[, fish_heatmap_genes, drop = FALSE], rownames = "fish_id") %>%
    mutate(time_week = "10 weeks"),
  as_tibble(X_fish_20_raw[, fish_heatmap_genes, drop = FALSE], rownames = "fish_id") %>%
    mutate(time_week = "20 weeks")
) %>%
  left_join(fish_info %>% select(fish_id, lifespan_group), by = "fish_id") %>%
  pivot_longer(all_of(fish_heatmap_genes), names_to = "gene_id", values_to = "expression") %>%
  group_by(gene_id, lifespan_group, time_week) %>%
  summarise(mean_expression = mean(expression), .groups = "drop") %>%
  group_by(gene_id) %>%
  mutate(z_expression = as.numeric(scale(mean_expression))) %>%
  ungroup() %>%
  mutate(
    gene_id = factor(gene_id, levels = rev(fish_heatmap_genes)),
    group_time = factor(
      paste(lifespan_group, time_week, sep = "\n"),
      levels = c("short_lived\n10 weeks", "short_lived\n20 weeks", "long_lived\n10 weeks", "long_lived\n20 weeks")
    )
  ) %>%
  ggplot(aes(group_time, gene_id, fill = z_expression)) +
  geom_tile(color = "white", linewidth = 0.25) +
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", midpoint = 0) +
  labs(
    title = "Top lifespan-associated 10-to-20-week expression changes",
    x = NULL,
    y = "Gene ID",
    fill = "Mean\nz-score"
  ) +
  theme(axis.text.x = element_text(size = 8))

ggsave(file.path(fig_dir, "fish_delta_group_heatmap.png"), fish_heatmap,
  width = 7.4, height = 7.8, dpi = 300
)

fish_loading_delta <- fish_gene_loadings %>%
  left_join(fish_delta_tests, by = "gene_id") %>%
  mutate(label = ifelse(gene_id %in% fish_boot_stability$gene_id[seq_len(8)], gene_id, NA_character_))

loading_delta_cor <- safe_cor(
  fish_loading_delta$loading_comp1_rss,
  abs(fish_loading_delta$difference_long_minus_short),
  method = "spearman"
)

fish_loading_delta_plot <- ggplot(
  fish_loading_delta,
  aes(loading_comp1_rss, abs(difference_long_minus_short))
) +
  geom_point(alpha = 0.45, color = "#3B6EA8") +
  geom_text_repel(aes(label = label), size = 2.6, max.overlaps = 20, na.rm = TRUE) +
  labs(
    title = "Trajectory loadings track group-specific 10-to-20-week expression changes",
    subtitle = paste0(
      "Spearman r = ", round(loading_delta_cor$r, 3),
      ", p ", ifelse(loading_delta_cor$p_value == 0, "< 2.2e-16", paste0("= ", signif(loading_delta_cor$p_value, 3))),
      "; top-50 overlap = ", fish_overlap$overlap_n, "/50"
    ),
    x = "Aggregated Component 1 loading",
    y = "|Long-lived delta minus short-lived delta|"
  )

ggsave(file.path(fig_dir, "fish_loading_delta_relation.png"), fish_loading_delta_plot,
  width = 7.2, height = 5.2, dpi = 300
)

fish_boot_plot <- fish_boot_stability %>%
  slice_head(n = 20) %>%
  mutate(gene_id = factor(gene_id, levels = rev(gene_id))) %>%
  ggplot(aes(top50_frequency, gene_id)) +
  geom_col(fill = "#6A3D9A", width = 0.75) +
  scale_x_continuous(labels = percent_format()) +
  labs(
    title = "Bootstrap stability of fish trajectory PLS loading genes",
    x = "Frequency in bootstrap top 50",
    y = "Gene ID"
  )

ggsave(file.path(fig_dir, "fish_bootstrap_loading_stability.png"), fish_boot_plot,
  width = 7.2, height = 5.8, dpi = 300
)

message("Vaccine trajectory analysis")

gse_path <- file.path(root_dir, "data", "GSE79396", "GSE79396_series_matrix.txt.gz")
gpl_path <- file.path(root_dir, "data", "GSE79396", "GPL13158.annot.gz")
days_keep <- c(0, 1, 3, 7)
top_n_features <- 2000

vax_meta_geo <- read_geo_metadata(gse_path)
sample_title <- vax_meta_geo[["!Sample_title"]][[1]]
geo_accession <- vax_meta_geo[["!Sample_geo_accession"]][[1]]
source_name <- vax_meta_geo[["!Sample_source_name_ch1"]][[1]]
characteristics <- vax_meta_geo[["!Sample_characteristics_ch1"]]

vax_sample_info <- tibble(
  geo_accession = geo_accession,
  sample_title = sample_title,
  source_name = source_name,
  subject_id = extract_characteristic(characteristics, "subject id:", length(sample_title)),
  vaccine = extract_characteristic(characteristics, "vaccine:", length(sample_title)),
  cohort = extract_characteristic(characteristics, "cohort:", length(sample_title)),
  visit = extract_characteristic(characteristics, "visit:", length(sample_title)),
  sex = extract_characteristic(characteristics, "Sex:", length(sample_title)),
  age = as.numeric(extract_characteristic(characteristics, "age:", length(sample_title))),
  age_group = extract_characteristic(characteristics, "age group:", length(sample_title)),
  tissue = extract_characteristic(characteristics, "tissue/cell type:", length(sample_title))
) %>%
  mutate(
    day = as.integer(sub(".*_D([0-9]+)_.*", "\\1", sample_title)),
    age_group = factor(age_group, levels = c("young", "elderly")),
    cohort = factor(cohort),
    sex = factor(sex)
  )

expr_dt <- data.table::fread(
  gse_path,
  skip = "!series_matrix_table_begin",
  data.table = TRUE,
  check.names = FALSE,
  showProgress = FALSE
)
if (!("ID_REF" %in% names(expr_dt))) {
  data.table::setnames(expr_dt, names(expr_dt)[1], "ID_REF")
}
expr_dt <- expr_dt[!grepl("^!", ID_REF)]
sample_cols <- intersect(vax_sample_info$geo_accession, names(expr_dt))
expr_mat <- as.matrix(expr_dt[, ..sample_cols])
mode(expr_mat) <- "numeric"
rownames(expr_mat) <- expr_dt$ID_REF
rm(expr_dt)
gc()

vax_complete_time <- vax_sample_info %>%
  filter(day %in% days_keep, geo_accession %in% sample_cols)

complete_subjects <- vax_complete_time %>%
  group_by(subject_id) %>%
  summarise(n_time = n_distinct(day), .groups = "drop") %>%
  filter(n_time == length(days_keep)) %>%
  pull(subject_id)

vax_sample_complete <- vax_complete_time %>%
  filter(subject_id %in% complete_subjects) %>%
  arrange(subject_id, day)

vax_subject_meta <- vax_sample_complete %>%
  group_by(subject_id) %>%
  summarise(
    age_group = first(age_group),
    age = first(age),
    sex = first(sex),
    cohort = first(cohort),
    .groups = "drop"
  ) %>%
  arrange(subject_id)

expr_subset <- expr_mat[, vax_sample_complete$geo_accession, drop = FALSE]
expr_subset <- expr_subset[matrixStats::rowAlls(is.finite(expr_subset)), , drop = FALSE]
probe_vars <- matrixStats::rowVars(expr_subset)
names(probe_vars) <- rownames(expr_subset)
top_probes <- names(sort(probe_vars, decreasing = TRUE))[seq_len(top_n_features)]
expr_top <- expr_subset[top_probes, , drop = FALSE]

probe_mean <- rowMeans(expr_top)
probe_sd <- matrixStats::rowSds(expr_top)
probe_sd[!is.finite(probe_sd) | probe_sd == 0] <- 1
expr_scaled <- sweep(sweep(expr_top, 1, probe_mean, "-"), 1, probe_sd, "/")
expr_scaled[!is.finite(expr_scaled)] <- 0

subject_order <- vax_subject_meta$subject_id
X_vax <- array(
  NA_real_,
  dim = c(length(subject_order), length(top_probes), length(days_keep)),
  dimnames = list(subject_order, top_probes, paste0("D", days_keep))
)

for (k in seq_along(days_keep)) {
  d <- days_keep[k]
  day_samples <- vax_sample_complete %>%
    filter(day == d) %>%
    arrange(match(subject_id, subject_order))
  if (!identical(day_samples$subject_id, subject_order)) {
    stop("Subject order mismatch in vaccine tensor construction.")
  }
  X_vax[, , k] <- t(expr_scaled[, day_samples$geo_accession, drop = FALSE])
}

X_vax_flat <- do.call(cbind, lapply(seq_along(days_keep), function(k) {
  mat <- X_vax[, , k, drop = FALSE][, , 1]
  colnames(mat) <- paste0(colnames(mat), "|D", days_keep[k])
  mat
}))

vax_pca <- prcomp(X_vax_flat, center = TRUE, scale. = FALSE, rank. = 5)
vax_scores <- as_tibble(vax_pca$x[, 1:5, drop = FALSE], rownames = "subject_id") %>%
  setNames(c("subject_id", paste0("PC", 1:5))) %>%
  left_join(vax_subject_meta, by = "subject_id")

if (
  mean(vax_scores$PC2[vax_scores$age_group == "elderly"]) <
    mean(vax_scores$PC2[vax_scores$age_group == "young"])
) {
  vax_scores$PC2 <- -vax_scores$PC2
  vax_pca$rotation[, "PC2"] <- -vax_pca$rotation[, "PC2"]
}

vax_pillai <- extract_pillai(vax_scores, c("PC1", "PC2"), "age_group") %>%
  mutate(analysis = "trajectory PCA PC1-PC2 ~ age group", .before = 1)
vax_perm <- permute_pillai(vax_scores, c("PC1", "PC2"), "age_group", n_perm = 999) %>%
  mutate(analysis = "trajectory PCA age-group permutation", .before = 1)

component_effects <- bind_rows(lapply(paste0("PC", 1:5), function(pc) {
  f_age <- as.formula(paste0(pc, " ~ age_group"))
  f_adj <- as.formula(paste0(pc, " ~ age_group + cohort + sex"))
  fit_age <- summary(lm(f_age, data = vax_scores))
  fit_adj <- anova(lm(f_adj, data = vax_scores))
  tibble(
    component = pc,
    mean_young = mean(vax_scores[[pc]][vax_scores$age_group == "young"]),
    mean_elderly = mean(vax_scores[[pc]][vax_scores$age_group == "elderly"]),
    elderly_minus_young = mean_elderly - mean_young,
    cohen_d_elderly_vs_young = cohens_d(vax_scores[[pc]], vax_scores$age_group, "elderly"),
    age_group_p = pf(
      fit_age$fstatistic[1],
      fit_age$fstatistic[2],
      fit_age$fstatistic[3],
      lower.tail = FALSE
    ),
    adjusted_age_group_p = fit_adj["age_group", "Pr(>F)"],
    cohort_p_in_adjusted_model = fit_adj["cohort", "Pr(>F)"]
  )
})) %>%
  mutate(
    age_group_fdr = p.adjust(age_group_p, method = "BH"),
    adjusted_age_group_fdr = p.adjust(adjusted_age_group_p, method = "BH")
  )

annot <- read_platform_annotation(gpl_path)

vax_loadings <- as_tibble(vax_pca$rotation[, 1:5, drop = FALSE], rownames = "feature") %>%
  mutate(
    probe_id = sub("\\|.*$", "", feature),
    day = sub("^.*\\|", "", feature)
  )

vax_probe_loadings <- vax_loadings %>%
  group_by(probe_id) %>%
  summarise(
    pc2_loading_rss = sqrt(sum(PC2^2)),
    pc2_max_abs_loading = max(abs(PC2)),
    pc2_dominant_day = day[which.max(abs(PC2))],
    .groups = "drop"
  ) %>%
  arrange(desc(pc2_loading_rss)) %>%
  left_join(annot, by = "probe_id")

keyword_enrichment <- function(top_ids, background_ids, annot_df, pattern, label) {
  bg <- annot_df %>%
    filter(probe_id %in% background_ids) %>%
    mutate(hit = grepl(pattern, paste(go_process, gene_title, gene_symbol), ignore.case = TRUE))
  top <- bg$probe_id %in% top_ids
  mat <- table(top = top, hit = bg$hit)
  ft <- fisher.test(mat)
  tibble(
    category = label,
    top_hits = sum(bg$hit[top]),
    top_n = sum(top),
    background_hits = sum(bg$hit),
    background_n = nrow(bg),
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value
  )
}

vax_keyword_tests <- bind_rows(
  keyword_enrichment(
    vax_probe_loadings$probe_id[1:100], top_probes, annot,
    "inflamm|cytokine|NF-kappa|leukocyte|neutrophil|immune|lipopolysaccharide|defense response|interferon",
    "inflammatory/innate immune terms"
  ),
  keyword_enrichment(
    vax_probe_loadings$probe_id[1:100], top_probes, annot,
    "oxidative|superoxide|reactive oxygen",
    "oxidative-stress terms"
  ),
  keyword_enrichment(
    vax_probe_loadings$probe_id[1:100], top_probes, annot,
    "cholesterol|sterol|lipid|fatty acid|isoprenoid",
    "lipid/sterol terms"
  ),
  keyword_enrichment(
    vax_probe_loadings$probe_id[1:100], top_probes, annot,
    "antigen processing|antigen presentation|MHC",
    "antigen-presentation terms"
  )
) %>%
  mutate(fdr = p.adjust(p_value, method = "BH"))

old_tpca_pc2_path <- file.path(root_dir, "results", "GSE79396_tPCA", "tables", "top_loadings_pc2.csv")
if (file.exists(old_tpca_pc2_path)) {
  old_tpca_pc2_annot <- read_csv(old_tpca_pc2_path, show_col_types = FALSE) %>%
    left_join(annot, by = "probe_id")
  write_csv(old_tpca_pc2_annot, file.path(tab_dir, "vaccine_tensoromics_pc2_top_loadings_annotated.csv"))
}

top_pc2_plot_probes <- vax_probe_loadings %>%
  filter(!is.na(gene_symbol), gene_symbol != "") %>%
  distinct(gene_symbol, .keep_all = TRUE) %>%
  slice_head(n = 8) %>%
  pull(probe_id)

vax_traj_df <- lapply(seq_along(days_keep), function(k) {
  mat <- X_vax[, top_pc2_plot_probes, k, drop = FALSE][, , 1]
  as_tibble(mat, rownames = "subject_id") %>%
    mutate(day = days_keep[k])
}) %>%
  bind_rows() %>%
  left_join(vax_subject_meta, by = "subject_id") %>%
  pivot_longer(all_of(top_pc2_plot_probes), names_to = "probe_id", values_to = "scaled_expression") %>%
  left_join(annot %>% select(probe_id, gene_symbol), by = "probe_id") %>%
  group_by(age_group, day, probe_id, gene_symbol) %>%
  summarise(
    mean_expression = mean(scaled_expression),
    se = sd(scaled_expression) / sqrt(n()),
    .groups = "drop"
  ) %>%
  mutate(feature_label = ifelse(gene_symbol == "" | is.na(gene_symbol), probe_id, gene_symbol))

supp_dir <- file.path(
  root_dir,
  "data",
  "GSE79396_Cell2017_supplement",
  "SupplementalPackage-MMRN",
  "sourcedata"
)
age_table <- read.delim(file.path(supp_dir, "VZV_age_gender2.txt"), check.names = FALSE) %>%
  mutate(
    subject_id = normalize_zv_id(ZV),
    age_group_supp = factor(ifelse(`Age Group` == "E", "elderly", "young"), levels = c("young", "elderly"))
  )

btm_raw <- data.table::fread(file.path(supp_dir, "formatted_BTMs.txt"), data.table = FALSE, check.names = FALSE)
btm_long <- btm_raw %>%
  pivot_longer(-BTM, names_to = "subject_id", values_to = "module_score") %>%
  mutate(subject_id = normalize_zv_id(subject_id)) %>%
  left_join(age_table %>% select(subject_id, age_group_supp, AGE, Gender), by = "subject_id") %>%
  mutate(module_time = sub("_.*$", "", BTM), module_name = sub("^[^_]+_", "", BTM))

btm_age_tests <- btm_long %>%
  filter(!is.na(age_group_supp)) %>%
  group_by(BTM, module_time, module_name) %>%
  summarise(
    n = sum(is.finite(module_score)),
    mean_young = mean(module_score[age_group_supp == "young"], na.rm = TRUE),
    mean_elderly = mean(module_score[age_group_supp == "elderly"], na.rm = TRUE),
    elderly_minus_young = mean_elderly - mean_young,
    cohen_d = cohens_d(module_score, age_group_supp, "elderly"),
    p_value = t.test(module_score ~ age_group_supp)$p.value,
    .groups = "drop"
  ) %>%
  mutate(fdr = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

pc2_scores_for_btm <- vax_scores %>%
  select(subject_id, trajectory_PC2 = PC2)

btm_pc2_cor <- btm_long %>%
  inner_join(pc2_scores_for_btm, by = "subject_id") %>%
  group_by(BTM, module_time, module_name) %>%
  summarise(
    n = sum(is.finite(module_score) & is.finite(trajectory_PC2)),
    r = safe_cor(module_score, trajectory_PC2)$r,
    p_value = safe_cor(module_score, trajectory_PC2)$p_value,
    .groups = "drop"
  ) %>%
  mutate(fdr = p.adjust(p_value, method = "BH")) %>%
  arrange(p_value)

igg <- transpose_feature_table(file.path(supp_dir, "formated2_IgG_NLSrenalyzed.txt"), "outcome")
tfh <- transpose_feature_table(file.path(supp_dir, "formatted2_TFH.txt"), "outcome")
outcome_df <- vax_scores %>%
  select(subject_id, PC1, PC2, age_group, cohort) %>%
  left_join(igg %>% select(subject_id, `D30/0_IgG`), by = "subject_id") %>%
  left_join(tfh %>% select(subject_id, `D7/0_TFH.CXCR3pos`), by = "subject_id")

vax_outcome_cor <- bind_rows(
  safe_cor(outcome_df$PC2, outcome_df$`D30/0_IgG`) %>%
    mutate(outcome = "D30/0 IgG", component = "trajectory PC2"),
  safe_cor(outcome_df$PC2, outcome_df$`D7/0_TFH.CXCR3pos`) %>%
    mutate(outcome = "D7/0 CXCR3+ Tfh", component = "trajectory PC2")
) %>%
  select(outcome, component, everything())

write_csv(vax_sample_complete, file.path(tab_dir, "vaccine_geo_sample_metadata_complete.csv"))
write_csv(vax_scores, file.path(tab_dir, "vaccine_trajectory_pca_scores.csv"))
write_csv(bind_rows(vax_pillai, vax_perm), file.path(tab_dir, "vaccine_trajectory_statistical_tests.csv"))
write_csv(component_effects, file.path(tab_dir, "vaccine_component_age_effects.csv"))
write_csv(vax_probe_loadings, file.path(tab_dir, "vaccine_trajectory_pc2_probe_loadings_annotated.csv"))
write_csv(vax_keyword_tests, file.path(tab_dir, "vaccine_pc2_keyword_enrichment.csv"))
write_csv(btm_age_tests, file.path(tab_dir, "vaccine_btm_age_group_tests.csv"))
write_csv(btm_pc2_cor, file.path(tab_dir, "vaccine_btm_pc2_correlations.csv"))
write_csv(vax_outcome_cor, file.path(tab_dir, "vaccine_trajectory_pc2_response_correlations.csv"))

vax_score_plot <- ggplot(vax_scores, aes(PC1, PC2, color = age_group, shape = cohort)) +
  geom_point(size = 3, alpha = 0.9) +
  geom_text_repel(aes(label = subject_id), size = 2.4, max.overlaps = 18) +
  scale_color_manual(values = c(young = "#1B9E77", elderly = "#D95F02")) +
  labs(
    title = "Trajectory PCA of GSE79396",
    subtitle = "Each point is one subject represented by D0, D1, D3 and D7 together",
    x = "Trajectory PC1",
    y = "Trajectory PC2",
    color = "Age group",
    shape = "Cohort"
  )

ggsave(file.path(fig_dir, "vaccine_trajectory_pca_scores.png"), vax_score_plot,
  width = 8.2, height = 6.2, dpi = 300
)

component_plot <- component_effects %>%
  mutate(component = factor(component, levels = paste0("PC", 1:5))) %>%
  ggplot(aes(component, cohen_d_elderly_vs_young, fill = adjusted_age_group_fdr < 0.05)) +
  geom_hline(yintercept = 0, color = "grey45") +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0("adj. p=", signif(adjusted_age_group_p, 2))),
    vjust = ifelse(component_effects$cohen_d_elderly_vs_young >= 0, -0.4, 1.2),
    size = 3
  ) +
  scale_fill_manual(values = c("TRUE" = "#D95F02", "FALSE" = "#A6A6A6"), guide = "none") +
  labs(
    title = "Age-group effect by trajectory PCA component",
    x = NULL,
    y = "Cohen's d, elderly minus young"
  )

ggsave(file.path(fig_dir, "vaccine_age_effects_by_component.png"), component_plot,
  width = 7.4, height = 4.7, dpi = 300
)

vax_traj_plot <- ggplot(vax_traj_df, aes(day, mean_expression, color = age_group, group = age_group)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean_expression - se, ymax = mean_expression + se),
    width = 0.18, linewidth = 0.35
  ) +
  facet_wrap(~feature_label, scales = "free_y", ncol = 4) +
  scale_color_manual(values = c(young = "#1B9E77", elderly = "#D95F02")) +
  scale_x_continuous(breaks = days_keep) +
  labs(
    title = "Top annotated trajectory PC2 probes",
    x = "Day post-vaccination",
    y = "Mean scaled expression",
    color = "Age group"
  ) +
  theme(strip.text = element_text(size = 8))

ggsave(file.path(fig_dir, "vaccine_pc2_top_gene_trajectories_annotated.png"), vax_traj_plot,
  width = 10.5, height = 6.6, dpi = 300
)

btm_age_plot <- btm_age_tests %>%
  filter(is.finite(cohen_d)) %>%
  slice_min(p_value, n = 18) %>%
  mutate(
    label = paste0(module_time, " ", sub(" \\([^)]*\\)$", "", module_name)),
    label = make.unique(label),
    label = factor(label, levels = rev(label))
  ) %>%
  ggplot(aes(cohen_d, label, fill = cohen_d > 0)) +
  geom_vline(xintercept = 0, color = "grey45") +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c("TRUE" = "#D95F02", "FALSE" = "#1B9E77"), guide = "none") +
  labs(
    title = "BTM modules with the strongest young-versus-elderly differences",
    x = "Cohen's d, elderly minus young",
    y = NULL
  )

ggsave(file.path(fig_dir, "vaccine_btm_age_effects.png"), btm_age_plot,
  width = 9.2, height = 6.4, dpi = 300
)

btm_pc2_plot <- btm_pc2_cor %>%
  filter(is.finite(r)) %>%
  arrange(p_value) %>%
  slice_head(n = 18) %>%
  mutate(
    label = paste0(module_time, " ", sub(" \\([^)]*\\)$", "", module_name)),
    label = make.unique(label),
    label = factor(label, levels = rev(label))
  ) %>%
  ggplot(aes(r, label, fill = r > 0)) +
  geom_vline(xintercept = 0, color = "grey45") +
  geom_col(width = 0.72) +
  scale_fill_manual(values = c("TRUE" = "#D95F02", "FALSE" = "#1B9E77"), guide = "none") +
  labs(
    title = "BTM modules correlated with trajectory PC2",
    x = "Pearson correlation with trajectory PC2",
    y = NULL
  )

ggsave(file.path(fig_dir, "vaccine_pc2_btm_correlations.png"), btm_pc2_plot,
  width = 9.2, height = 6.4, dpi = 300
)

response_summary_path <- file.path(
  root_dir,
  "results",
  "GSE79396_response",
  "tables",
  "response_analysis_summary.csv"
)

if (file.exists(response_summary_path)) {
  response_summary <- read_csv(response_summary_path, show_col_types = FALSE) %>%
    transmute(
      outcome = recode(outcome, IgG_D30_0 = "IgG D30/0", TFH_CXCR3pos_D7_0 = "CXCR3+ Tfh D7/0"),
      in_sample_r2 = comp1_r^2,
      cv_r2 = tpls_cv_r2,
      permutation_p = permutation_p
    ) %>%
    pivot_longer(c(in_sample_r2, cv_r2), names_to = "estimate", values_to = "r2") %>%
    mutate(estimate = recode(estimate, in_sample_r2 = "In-sample r^2", cv_r2 = "Cross-validated R^2"))

  response_plot <- ggplot(response_summary, aes(outcome, r2, fill = estimate)) +
    geom_hline(yintercept = 0, color = "grey45") +
    geom_col(position = position_dodge(width = 0.75), width = 0.65) +
    scale_fill_manual(values = c("In-sample r^2" = "#377EB8", "Cross-validated R^2" = "#E41A1C")) +
    labs(
      title = "Response-level tPLS associations did not validate out of sample",
      x = NULL,
      y = "Variance explained",
      fill = NULL
    )

  ggsave(file.path(fig_dir, "vaccine_response_validation_summary.png"), response_plot,
    width = 7.0, height = 4.6, dpi = 300
  )
}

summary_lines <- c(
  "Trajectory analysis summary",
  "",
  "Fish aging:",
  paste0("  Complete fish: ", nrow(fish_info)),
  paste0("  Median age at death: ", median_death, " days"),
  paste0("  Top variable genes: ", length(top_fish_genes)),
  paste0("  Trajectory PLS MANOVA p: ", signif(fish_pillai$p_value, 4)),
  paste0("  Trajectory PLS permutation p: ", signif(fish_perm$p_value, 4)),
  paste0("  Age-at-death model R2: ", signif(fish_lm$r.squared, 4)),
  paste0("  Top-loading/top-delta overlap: ", fish_overlap$overlap_n, "/50"),
  "",
  "Vaccine response:",
  paste0("  Complete subjects: ", nrow(vax_subject_meta)),
  paste0("  Young/elderly: ", paste(capture.output(print(table(vax_subject_meta$age_group))), collapse = " ")),
  paste0("  Top variable probes: ", length(top_probes)),
  paste0("  Trajectory PCA MANOVA p: ", signif(vax_pillai$p_value, 4)),
  paste0("  Trajectory PCA permutation p: ", signif(vax_perm$p_value, 4)),
  paste0(
    "  Strongest age component: ",
    component_effects$component[which.min(component_effects$adjusted_age_group_p)],
    " (adjusted p = ",
    signif(min(component_effects$adjusted_age_group_p), 4),
    ")"
  ),
  "",
  paste0("Figures: ", fig_dir),
  paste0("Tables: ", tab_dir)
)

writeLines(summary_lines, file.path(tab_dir, "trajectory_analysis_summary.txt"))
message(paste(summary_lines, collapse = "\n"))
