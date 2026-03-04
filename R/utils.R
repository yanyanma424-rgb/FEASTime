# utils.R

#' Summarise a FEASTime result
#'
#' Returns per-stage mean, SD, and sample counts as a tidy data frame.
#'
#' @param feastime_result Output from run_feastime().
#' @return A data.frame with columns Stage, Source, Mean, SD, n_valid, n_total.
#' @export
summarise_feastime <- function(feastime_result) {
  results <- do.call(rbind, lapply(names(feastime_result$per_stage_results), function(st) {
    stage_list <- feastime_result$per_stage_results[[st]]
    valid <- Filter(Negate(is.null), stage_list)
    n_total <- length(stage_list)
    n_valid <- length(valid)
    if (n_valid == 0) return(NULL)
    mat   <- do.call(rbind, valid)
    means <- colMeans(mat, na.rm = TRUE)
    sds   <- apply(mat, 2, stats::sd, na.rm = TRUE)
    data.frame(Stage = st, Source = names(means),
               Mean = as.numeric(means), SD = as.numeric(sds),
               n_valid = n_valid, n_total = n_total,
               stringsAsFactors = FALSE)
  }))
  rownames(results) <- NULL
  results
}

#' Export FEASTime results to Excel (.xlsx)
#'
#' Requires the writexl package.
#'
#' @param feastime_result Output from run_feastime().
#' @param path File path for the .xlsx file.
#' @return Invisibly returns path.
#' @export
export_xlsx <- function(feastime_result, path) {
  if (!requireNamespace("writexl", quietly = TRUE))
    stop("writexl is required. Install with: install.packages(writexl)")
  sheets <- list()
  cumul <- feastime_result$cumulative_contributions
  sheets$cumulative <- do.call(rbind, lapply(names(cumul), function(st) {
    v <- cumul[[st]]
    data.frame(Stage = st, t(v), check.names = FALSE, stringsAsFactors = FALSE)
  }))
  for (st in names(feastime_result$per_stage_results)) {
    valid <- Filter(Negate(is.null), feastime_result$per_stage_results[[st]])
    if (length(valid) == 0) next
    df <- do.call(rbind, lapply(names(valid), function(nm) {
      data.frame(Sample = nm, t(valid[[nm]]),
                 check.names = FALSE, stringsAsFactors = FALSE)
    }))
    sheets[[paste0("stage_", st)]] <- df
  }
  writexl::write_xlsx(sheets, path)
  invisible(path)
}
