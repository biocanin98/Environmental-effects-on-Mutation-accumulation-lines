# Interactive comparison of 2+ oCelloscope measurements of ONE treatment.
# CHOOSE the CSV files by number when the script starts; no filename edits needed.
# Plate layout: Layout_parents_ancestors.ods (4 wells/line/measurement).
# Each line has R1 and R2, with two wells for each replicate label.
# Each comparison creates its own time-stamped folder inside
# MA_measurement_comparison_results/, preserving earlier results.
#
# REQUIRED INPUTS in the same folder as this script:
#   Layout_parents_ancestors.ods
#   oCelloscope_functions.R
#   MA_lines.zip (automatically extracted), OR an extracted MA_lines/ folder
#
# REQUIRED R PACKAGES (install once if absent):
#   install.packages(c("readODS", "dplyr", "tidyr", "ggplot2"))

# ============================ CONFIGURATION ===============================
# Choose 2 or more files interactively after the script finds available CSVs.
# Only select files from the SAME treatment/experimental condition.

r2_cutoff <- 0.90               # Original QC script: STRICT R2 > 0.90
large_deviation_pct <- 20      # Flag |between-measurement change| >= 20%
min_window_points <- 8L        # Require adequate linear-regression observations

# Same 15,000-s window width, 5,000-s step as oCelloscope_functions.R.
limitmatrix <- matrix(
  c(0, 15000, 5000, 20000, 10000, 25000, 15000, 30000,
    20000, 35000, 25000, 40000, 30000, 45000, 35000, 50000,
    40000, 55000, 45000, 60000), ncol = 2, byrow = TRUE
)
# ========================================================================

required <- c("readODS", "dplyr", "tidyr", "ggplot2")
missing_packages <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages)) {
  stop("Missing R packages: ", paste(missing_packages, collapse = ", "),
       "\nInstall with: install.packages(c(",
       paste(sprintf('"%s"', missing_packages), collapse = ", "), "))")
}

# Original convert.ocelloscope.data.long() uses unqualified pivot_longer().
library(tidyr)

# Resolve inputs relative to this script even if it was sourced elsewhere.
script_path <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
base_dir <- if (!is.null(script_path) && nzchar(script_path)) {
  dirname(normalizePath(script_path))
} else {
  getwd()
}

layout_path <- file.path(base_dir, "Layout_parents_ancestors.ods")
function_path <- file.path(base_dir, "oCelloscope_functions.R")
zip_path <- file.path(base_dir, "MA_lines.zip")
data_dir <- file.path(base_dir, "MA_lines")
if (!file.exists(layout_path)) stop("Missing layout: ", layout_path)
if (!file.exists(function_path)) stop("Missing growth function file: ", function_path)
if (!dir.exists(data_dir) && file.exists(zip_path)) {
  message("Extracting MA_lines.zip...")
  utils::unzip(zip_path, exdir = base_dir)
}
# ================= INTERACTIVE MEASUREMENT SELECTION ======================
# Show all measurement CSVs found in MA_lines/. Enter any two or more file
# numbers, e.g. "1, 3, 8". Each selected file becomes M1, M2, M3, ...
if (!dir.exists(data_dir)) {
  stop("Cannot find MA_lines/ or MA_lines.zip beside this script: ", base_dir)
}
available_files <- sort(list.files(data_dir, pattern = "\\.csv$",
                                   recursive = TRUE, full.names = TRUE,
                                   ignore.case = TRUE))
if (length(available_files) < 2L) {
  stop("At least two measurement CSVs are required in: ", data_dir)
}

cat("\n================ AVAILABLE OCELLOSCOPE MEASUREMENTS ================\n")
for (i in seq_along(available_files)) {
  cat(sprintf("%3d. %s\n", i, basename(available_files[i])))
}
cat("====================================================================\n")
cat("Choose files from the SAME treatment. You can select 2, 3, 4, ... files.\n")
if (!interactive()) {
  stop("Interactive selection needs the RStudio/R console. Run source() in RStudio, not Rscript.")
}
repeat {
  answer <- trimws(readline(prompt =
    "Enter file numbers separated by commas (e.g. 1, 3, 5): "))
  if (!nzchar(answer)) stop("No files selected. Run the script again to choose files.")
  tokens <- strsplit(answer, "[,;[:space:]]+")[[1]]
  numbers <- suppressWarnings(as.integer(tokens))
  valid <- length(numbers) >= 2L &&
    all(!is.na(numbers)) && all(numbers >= 1L) &&
    all(numbers <= length(available_files)) &&
    !anyDuplicated(numbers)
  if (valid) break
  cat("Please enter at least TWO DIFFERENT valid numbers from the list.\n")
}
measurement_paths <- normalizePath(available_files[numbers])

# User-defined label, deliberately not inferred from potentially ambiguous
# filename conventions (e.g. "salinity_1" versus "salt_1").
treatment_name <- trimws(readline(
  prompt = "Treatment label for results folder (e.g. pH_8, salt_1): "))
if (!nzchar(treatment_name)) treatment_name <- "selected_treatment"
safe_label <- gsub("[^A-Za-z0-9_.-]+", "_", treatment_name)
output_dir <- file.path(base_dir, "MA_measurement_comparison_results",
                        paste0(safe_label, "_", format(Sys.time(), "%Y%m%d_%H%M%S")))
cat("\nSelected measurements for: ", treatment_name, "\n", sep = "")
for (i in seq_along(measurement_paths)) {
  cat(sprintf("  M%d = %s\n", i, basename(measurement_paths[i])))
}
cat("\n")
# ========================================================================

# Retain the supplied oCelloscope CSV loader/long-format conversion.
source(function_path, local = TRUE)

# Read the ODS plate layout. First row contains column numbers, first column
# contains plate letters. Ignore empty and edge wells. This gives 41 x 4 = 164
# occupied wells per measurement, with 2 technical wells per R label.
layout_raw <- as.data.frame(readODS::read_ods(layout_path, sheet = 1,
                                               col_names = FALSE),
                            stringsAsFactors = FALSE)
plate_col <- suppressWarnings(as.integer(as.character(unlist(layout_raw[1, -1]))))
plate_rows <- lapply(2:nrow(layout_raw), function(i) {
  row_letter <- trimws(as.character(layout_raw[[1]][i]))
  lapply(2:ncol(layout_raw), function(j) {
    value <- trimws(as.character(layout_raw[[j]][i]))
    if (is.na(value) || !nzchar(value) || is.na(plate_col[j-1]) ||
        is.na(row_letter) || !nzchar(row_letter)) return(NULL)
    data.frame(well = paste0(row_letter, plate_col[j-1]),
               label = value, line = sub("_R[0-9]+$", "", value),
               replicate = sub("^.*_(R[0-9]+)$", "\\1", value),
               stringsAsFactors = FALSE)
  })
})
plate_map <- dplyr::bind_rows(unlist(plate_rows, recursive = FALSE))
if (nrow(plate_map) == 0L) stop("No occupied wells found in ODS layout.")
if (anyDuplicated(plate_map$well)) stop("Duplicate well IDs in plate layout.")
if (any(!grepl("_R[0-9]+$", plate_map$label))) {
  stop("Some layout entries lack an _R# suffix; inspect the layout.")
}
planned_per_line <- table(plate_map$line)
if (any(planned_per_line != 4L)) {
  stop("Expected 4 occupied wells per line; found: ",
       paste(names(planned_per_line)[planned_per_line != 4L], collapse = ", "))
}
line_order <- c("anc_A", "anc_a", paste0("L", seq_len(40)))
line_order <- line_order[line_order %in% plate_map$line]

# Count the repetitions in the UniExplorer metadata (60 or 120, etc.).
read_repetitions <- function(csv_path) {
  con <- file(csv_path, open = "rt", encoding = "UTF-16")
  on.exit(close(con))
  header <- readLines(con, n = 60L, warn = FALSE)
  hit <- grep("^\\s*Repetitions\\s*:\\s*[0-9]+", header, value = TRUE)
  if (!length(hit)) stop("Cannot detect number of repetitions in: ", csv_path)
  as.integer(sub("^.*:\\s*([0-9]+).*$", "\\1", hit[1]))
}

# Fit every potential time window using the previous script's rules:
#   hyphal area ~ time in hours (growth = fitted slope)
#   at least one tip must exist at the START OF THE SELECTED WINDOW
#   discard negative slopes
#   select best adjusted R2; if >=2 windows have adjusted R2 > 0.995,
#     select the one with the highest positive growth slope.
# IMPORTANT: Original script records tips.at.limit1, but does NOT correct
# growth for tips. We retain the original slope AND additionally compute
# growth_per_initial_tip = slope / tips.at.limit1.
fit_one_well <- function(x) {
  x <- x[order(x$Time), , drop = FALSE]
  valid_windows <- list()
  max_time <- suppressWarnings(max(x$Time, na.rm = TRUE))
  if (!is.finite(max_time)) max_time <- -Inf
  for (j in seq_len(nrow(limitmatrix))) {
    from <- limitmatrix[j, 1]
    to <- limitmatrix[j, 2]
    # Do not fit truncated windows beyond the acquired measurement duration.
    if (max_time < to - 500) next
    w <- x[is.finite(x$Time) & x$Time >= from & x$Time <= to,
           , drop = FALSE]
    if (nrow(w) < min_window_points) next
    initial_tips <- w$no.tips[1]
    if (!is.finite(initial_tips) || initial_tips <= 0) next
    w <- w[is.finite(w$area) & is.finite(w$timestamp), , drop = FALSE]
    if (nrow(w) < min_window_points) next
    w$hours <- w$timestamp / 3600
    if (length(unique(w$hours)) < 3L) next
    model <- tryCatch(stats::lm(area ~ hours, data = w), error = function(e) NULL)
    if (is.null(model)) next
    growth <- unname(stats::coef(model)[2])
    r2 <- unname(summary(model)$adj.r.squared)
    if (!is.finite(growth) || !is.finite(r2) || growth <= 0) next
    valid_windows[[length(valid_windows) + 1L]] <- data.frame(
      gr = growth, R2 = r2, tips.at.limit1 = initial_tips,
      window = paste(from, to, sep = ":"), stringsAsFactors = FALSE
    )
  }
  if (!length(valid_windows)) {
    return(data.frame(gr = NA_real_, R2 = NA_real_,
                      tips.at.limit1 = NA_real_, window = NA_character_,
                      growth_per_initial_tip = NA_real_))
  }
  candidates <- dplyr::bind_rows(valid_windows)
  very_linear <- which(candidates$R2 > 0.995)
  best <- if (length(very_linear) >= 2L) {
    very_linear[which.max(candidates$gr[very_linear])]
  } else {
    which.max(candidates$R2)
  }
  answer <- candidates[best, , drop = FALSE]
  answer$growth_per_initial_tip <- answer$gr / answer$tips.at.limit1
  rownames(answer) <- NULL
  answer
}

# Analyze all mapped wells from all measurements; retain failed/missing wells
# instead of silently dropping them. QC: adjusted R2 STRICTLY > 0.90.
run_one_measurement <- function(csv_path, id) {
  message("Processing ", id, ": ", basename(csv_path))
  repetitions <- read_repetitions(csv_path)
  input <- load.ocelloscope.data2(csv_path, repetitions = repetitions,
                                   sep = ",", dec = ".")
  long <- convert.ocelloscope.data.long(input)
  long <- long[long$sample %in% plate_map$well, , drop = FALSE]
  by_well <- split(long, long$sample)
  fitted <- lapply(names(by_well), function(well) {
    data.frame(well = well, fit_one_well(by_well[[well]]),
               stringsAsFactors = FALSE)
  })
  fitted <- dplyr::bind_rows(fitted)
  result <- dplyr::left_join(plate_map, fitted, by = "well")
  result$measurement <- id
  result$file <- basename(csv_path)
  result$treatment <- treatment_name
  result$qc_pass <- is.finite(result$R2) & result$R2 > r2_cutoff &
    is.finite(result$growth_per_initial_tip) &
    result$growth_per_initial_tip > 0 &
    is.finite(result$tips.at.limit1) & result$tips.at.limit1 > 0
  result$qc_reason <- ifelse(result$qc_pass, "PASS",
                             ifelse(is.na(result$R2), "NO_VALID_FIT_OR_TIPS",
                                    "R2_AT_OR_BELOW_CUTOFF"))
  result
}

ids <- paste0("M", seq_along(measurement_paths))
measurement_key <- data.frame(measurement = ids,
                              filename = basename(measurement_paths),
                              treatment = treatment_name,
                              source_path = measurement_paths)
well_results <- dplyr::bind_rows(Map(run_one_measurement,
                                     measurement_paths, ids))
well_results$line <- factor(well_results$line, levels = line_order)
well_results$measurement <- factor(well_results$measurement, levels = ids)

# QC counts PER LINE and PER MEASUREMENT. Each line should have 4 wells;
# for two measurements it has 8 wells in total.
qc_by_line <- well_results |>
  dplyr::group_by(line, measurement) |>
  dplyr::summarise(
    n_planned = dplyr::n(),
    n_pass = sum(qc_pass),
    n_failed = sum(!qc_pass),
    n_R1_pass = sum(qc_pass & replicate == "R1"),
    n_R2_pass = sum(qc_pass & replicate == "R2"),
    n_without_fit = sum(is.na(R2)),
    n_R2_low = sum(is.finite(R2) & R2 <= r2_cutoff),
    .groups = "drop"
  )
qc_total <- qc_by_line |>
  dplyr::group_by(line) |>
  dplyr::summarise(n_pass_total = sum(n_pass),
                   n_planned_total = sum(n_planned),
                   n_failed_total = sum(n_failed),
                   .groups = "drop")
qc_wide <- qc_by_line |>
  dplyr::select(line, measurement, n_pass, n_failed, n_R1_pass, n_R2_pass) |>
  tidyr::pivot_wider(names_from = measurement,
                     values_from = c(n_pass, n_failed, n_R1_pass, n_R2_pass),
                     names_glue = "{.value}_{measurement}") |>
  dplyr::left_join(qc_total, by = "line") |>
  dplyr::arrange(line)

# R1 and R2 are each placed in TWO wells. Average QC-passing technical wells
# within each R-label before an exploratory significance test to avoid
# counting the repeated placement of an R-label as independent evidence.
valid_wells <- well_results[well_results$qc_pass, , drop = FALSE]
replicate_means <- valid_wells |>
  dplyr::group_by(line, measurement, replicate) |>
  dplyr::summarise(n_passing_wells = dplyr::n(),
                   mean_growth_per_tip = mean(growth_per_initial_tip),
                   mean_raw_growth = mean(gr),
                   .groups = "drop")

mean_or_na <- function(x) if (length(x)) mean(x) else NA_real_
sd_or_na <- function(x) if (length(x) >= 2L) stats::sd(x) else NA_real_
cv_or_na <- function(x) {
  m <- mean_or_na(x)
  s <- sd_or_na(x)
  if (is.finite(m) && m != 0 && is.finite(s)) 100 * s / abs(m) else NA_real_
}

# One row per line per measurement; keep all lines even if every well failed.
means_by_line <- replicate_means |>
  dplyr::group_by(line, measurement) |>
  dplyr::summarise(n_Rlabels_passing = dplyr::n(),
                   mean_growth_per_tip = mean(mean_growth_per_tip),
                   sd_across_Rlabels = sd_or_na(mean_growth_per_tip),
                   cv_across_Rlabels_pct = cv_or_na(mean_growth_per_tip),
                   mean_raw_growth = mean(mean_raw_growth),
                   .groups = "drop")
means_by_line <- qc_by_line |>
  dplyr::left_join(means_by_line, by = c("line", "measurement")) |>
  dplyr::mutate(n_Rlabels_passing = tidyr::replace_na(n_Rlabels_passing, 0L)) |>
  dplyr::arrange(line, measurement)

safe_p <- function(x, y) {
  if (length(x) < 2L || length(y) < 2L) return(NA_real_)
  if (all(x == x[1]) && all(y == y[1]) && x[1] == y[1]) return(1)
  suppressWarnings(tryCatch(stats::t.test(x, y, var.equal = FALSE)$p.value,
                            error = function(e) NA_real_))
}

# Pairwise comparisons for all measurements (M1 vs M2, M1 vs M3, ...).
# Difference is B - A. Percent change = 100 * (B / A - 1).
# p-values are exploratory: only two distinct R-labels per measurement.
pairs <- utils::combn(ids, 2, simplify = FALSE)
pair_rows <- lapply(line_order, function(this_line) {
  lapply(pairs, function(pair) {
    a <- pair[1]
    b <- pair[2]
    ax <- replicate_means$mean_growth_per_tip[
      as.character(replicate_means$line) == this_line &
        as.character(replicate_means$measurement) == a]
    bx <- replicate_means$mean_growth_per_tip[
      as.character(replicate_means$line) == this_line &
        as.character(replicate_means$measurement) == b]
    ma <- mean_or_na(ax)
    mb <- mean_or_na(bx)
    difference <- if (is.finite(ma) && is.finite(mb)) mb - ma else NA_real_
    change_pct <- if (is.finite(ma) && is.finite(mb) && ma != 0) {
      100 * (mb / ma - 1)
    } else NA_real_
    data.frame(line = this_line, measurement_A = a, measurement_B = b,
               n_Rlabels_A = length(ax), n_Rlabels_B = length(bx),
               mean_A_per_tip = ma, mean_B_per_tip = mb,
               difference_B_minus_A = difference,
               percent_change_B_vs_A = change_pct,
               absolute_percent_change = abs(change_pct),
               large_deviation = is.finite(change_pct) &
                 abs(change_pct) >= large_deviation_pct,
               welch_p = safe_p(ax, bx), stringsAsFactors = FALSE)
  })
})
pairwise <- dplyr::bind_rows(unlist(pair_rows, recursive = FALSE))
pairwise$welch_p_BH <- stats::p.adjust(pairwise$welch_p, method = "BH")
pairwise <- pairwise |>
  dplyr::left_join(measurement_key |>
                     dplyr::select(measurement, filename),
                   by = c("measurement_A" = "measurement")) |>
  dplyr::rename(file_A = filename) |>
  dplyr::left_join(measurement_key |>
                     dplyr::select(measurement, filename),
                   by = c("measurement_B" = "measurement")) |>
  dplyr::rename(file_B = filename)

# Summary across all selected measurements: QC totals and spread in mean
# per-initial-tip growth (relative to the mean of the observed measurements).
line_deviations <- means_by_line |>
  dplyr::group_by(line) |>
  dplyr::summarise(
    n_measurements_with_QC = sum(is.finite(mean_growth_per_tip)),
    minimum_mean_per_tip = if (all(is.na(mean_growth_per_tip))) NA_real_
                           else min(mean_growth_per_tip, na.rm = TRUE),
    maximum_mean_per_tip = if (all(is.na(mean_growth_per_tip))) NA_real_
                           else max(mean_growth_per_tip, na.rm = TRUE),
    overall_mean_per_tip = if (all(is.na(mean_growth_per_tip))) NA_real_
                           else mean(mean_growth_per_tip, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    range_pct_of_overall_mean = dplyr::if_else(
      n_measurements_with_QC >= 2 & overall_mean_per_tip > 0,
      100 * (maximum_mean_per_tip - minimum_mean_per_tip) / overall_mean_per_tip,
      NA_real_)) |>
  dplyr::left_join(qc_total, by = "line") |>
  dplyr::arrange(dplyr::desc(range_pct_of_overall_mean))

# ODS lines are always included, including those with 0 successful replicates.
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write.csv(measurement_key, file.path(output_dir, "00_measurement_key.csv"), row.names = FALSE)
write.csv(plate_map, file.path(output_dir, "01_plate_mapping.csv"), row.names = FALSE)
write.csv(well_results, file.path(output_dir, "02_well_growth_and_QC.csv"), row.names = FALSE)
write.csv(qc_by_line, file.path(output_dir, "03_QC_counts_by_line_and_measurement.csv"), row.names = FALSE)
write.csv(qc_wide, file.path(output_dir, "04_QC_counts_by_line_total.csv"), row.names = FALSE)
write.csv(replicate_means, file.path(output_dir, "05_R1_R2_means.csv"), row.names = FALSE)
write.csv(means_by_line, file.path(output_dir, "06_growth_summary_by_line_and_measurement.csv"), row.names = FALSE)
write.csv(pairwise, file.path(output_dir, "07_pairwise_measurement_deviations.csv"), row.names = FALSE)
write.csv(line_deviations, file.path(output_dir, "08_line_deviation_ranking.csv"), row.names = FALSE)

# Diagnostic: R2 distribution (the dashed threshold line is QC = 0.90).
p_r2 <- ggplot2::ggplot(well_results, ggplot2::aes(x = R2)) +
  ggplot2::geom_histogram(bins = 35, fill = "grey70", color = "grey20") +
  ggplot2::geom_vline(xintercept = r2_cutoff, linetype = "dashed", color = "red") +
  ggplot2::facet_wrap(~measurement, nrow = 1) +
  ggplot2::labs(title = paste(treatment_name, "— adjusted R² quality control"),
                x = "Adjusted R²", y = "Number of wells") +
  ggplot2::theme_bw()
ggplot2::ggsave(file.path(output_dir, "09_R2_diagnostics.png"),
                p_r2, width = 9, height = 4.5, dpi = 300)

# Measurement-to-measurement profiles of tip-normalized growth.
p_profile <- ggplot2::ggplot(
  means_by_line,
  ggplot2::aes(x = measurement, y = mean_growth_per_tip, group = line)
) +
  ggplot2::geom_line(color = "grey60", linewidth = 0.35, na.rm = TRUE) +
  ggplot2::geom_point(ggplot2::aes(color = n_pass), size = 1.7, na.rm = TRUE) +
  ggplot2::scale_color_gradient(low = "red3", high = "darkgreen", limits = c(0, 4)) +
  ggplot2::labs(title = paste(treatment_name, "— replicate-means across measurements"),
                x = "Measurement", y = "Hyphal area growth per hour per initial tip",
                color = "QC-passing wells") +
  ggplot2::theme_bw()
ggplot2::ggsave(file.path(output_dir, "10_measurement_profiles.png"),
                p_profile, width = 9, height = 6, dpi = 300)

# QC counts heatmap (red=0 pass, green=4 pass).
p_qc <- ggplot2::ggplot(qc_by_line,
                        ggplot2::aes(x = measurement, y = line, fill = n_pass)) +
  ggplot2::geom_tile(color = "white") +
  ggplot2::geom_text(ggplot2::aes(label = paste0(n_pass, "/", n_planned)), size = 2.7) +
  ggplot2::scale_fill_gradient(low = "red3", high = "darkgreen", limits = c(0, 4)) +
  ggplot2::labs(title = paste(treatment_name, "— QC-passing wells by line"),
                x = "Measurement", y = "Line", fill = "Passing wells") +
  ggplot2::theme_minimal()
ggplot2::ggsave(file.path(output_dir, "11_QC_heatmap.png"),
                p_qc, width = max(6, 2 + 1.5 * length(ids)), height = 11, dpi = 300)

message("\nComparison complete: ", treatment_name)
message("Total occupied wells across measurements: ", nrow(well_results))
message("QC-passing wells: ", sum(well_results$qc_pass), "/", nrow(well_results))
message("Lines flagged >= ", large_deviation_pct, "% difference in at least one pair: ",
        length(unique(as.character(pairwise$line[pairwise$large_deviation]))))
message("Files saved to: ", output_dir)
print(qc_wide, row.names = FALSE)
