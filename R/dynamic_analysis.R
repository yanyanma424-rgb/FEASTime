# ============================================================
# dynamic_analysis.R  --  FEASTime v2
#
# Changes vs v1:
#   [A] allow_continuous_sources parameter added to:
#       run_feastime(), .analyse_transfer()
#   [B] rarefaction_depth parameter threaded through to
#       .FEAST_internal()
# ============================================================

.validate_inputs <- function(abundance_table, metadata, stage_order, source_types) {
  if (!is.matrix(abundance_table) && !is.data.frame(abundance_table))
    stop("abundance_table must be a matrix or data.frame.")
  if (is.null(rownames(abundance_table)))
    stop("abundance_table must have row names (sample IDs).")
  if (!is.data.frame(metadata))
    stop("metadata must be a data.frame.")
  required_cols <- c("Env", "SourceSink", "id")
  missing <- setdiff(required_cols, colnames(metadata))
  if (length(missing) > 0)
    stop("metadata missing columns: ", paste(missing, collapse = ", "))
  bad_ss <- setdiff(unique(metadata$SourceSink), c("Source", "Sink"))
  if (length(bad_ss) > 0)
    stop("SourceSink contains invalid values: ", paste(bad_ss, collapse = ", "))
  missing_samples <- setdiff(rownames(abundance_table), rownames(metadata))
  if (length(missing_samples) > 0)
    stop("Some rows of abundance_table missing from metadata: ",
         paste(utils::head(missing_samples, 5), collapse = ", "))
  all_envs <- unique(metadata$Env)
  bad_stages <- setdiff(stage_order, all_envs)
  if (length(bad_stages) > 0)
    stop("stage_order values not in metadata$Env: ",
         paste(bad_stages, collapse = ", "))
  bad_sources <- setdiff(source_types, all_envs)
  if (length(bad_sources) > 0)
    stop("source_types values not in metadata$Env: ",
         paste(bad_sources, collapse = ", "))
  invisible(TRUE)
}

# ── Initial stage (unchanged logic) ──────────────────────────

.analyse_initial_stage <- function(abundance_table, metadata,
                                   first_stage, source_types,
                                   n_restarts, max_iter, tol,
                                   rarefaction) {
  sink_rows   <- rownames(metadata)[metadata$Env == first_stage &
                                      metadata$SourceSink == "Sink"]
  source_rows <- rownames(metadata)[metadata$Env %in% source_types &
                                      metadata$SourceSink == "Source"]
  if (length(sink_rows) == 0)
    stop("No Sink samples found for stage '", first_stage, "'.")
  if (length(source_rows) == 0)
    stop("No Source samples found for source_types: ",
         paste(source_types, collapse = ", "))

  per_sample <- list()
  for (s in sink_rows) {
    use_rows <- c(source_rows, s)
    sub_C    <- abundance_table[use_rows, , drop = FALSE]
    sub_meta <- metadata[use_rows, , drop = FALSE]
    sub_meta$SourceSink[rownames(sub_meta) %in% source_rows] <- "Source"
    sub_meta$SourceSink[rownames(sub_meta) == s]             <- "Sink"
    tryCatch({
      res <- .FEAST_internal(sub_C, sub_meta,
                             n_restarts  = n_restarts,
                             max_iter    = max_iter,
                             tol         = tol,
                             rarefaction = rarefaction)
      per_sample[[s]] <- res[1, ]
    }, error = function(e) {
      warning("Initial stage '", first_stage, "', sample '", s,
              "' failed: ", conditionMessage(e))
      per_sample[[s]] <<- NULL
    })
  }

  valid <- Filter(Negate(is.null), per_sample)
  if (length(valid) == 0) {
    warning("All samples failed in stage '", first_stage, "'.")
    dummy <- rep(NA, length(source_types) + 1)
    names(dummy) <- c(source_types, "Unknown")
    return(list(per_sample = per_sample, mean_profile = dummy))
  }
  profile_mat  <- do.call(rbind, valid)
  mean_profile <- colMeans(profile_mat, na.rm = TRUE)
  list(per_sample = per_sample, mean_profile = mean_profile)
}

# ── [A] Transfer stage with allow_continuous_sources ─────────

#' Analyse one transfer step in the Markov chain
#'
#' @param allow_continuous_sources Logical. If \code{TRUE}, the original
#'   \code{source_types} samples are included alongside the previous-stage
#'   samples as known sources. This models persistent environmental input
#'   (e.g. ongoing Daqu inoculation) rather than a single initial seeding.
#'   Default \code{FALSE} (original v1 behaviour).
#' @keywords internal
.analyse_transfer <- function(abundance_table, metadata,
                              previous_stage, current_stage,
                              source_types,
                              allow_continuous_sources,
                              n_restarts, max_iter, tol,
                              rarefaction) {

  prev_rows <- rownames(metadata)[metadata$Env == previous_stage]
  curr_rows <- rownames(metadata)[metadata$Env == current_stage &
                                    metadata$SourceSink == "Sink"]
  if (length(prev_rows) == 0)
    stop("No samples found for previous stage '", previous_stage, "'.")
  if (length(curr_rows) == 0)
    stop("No Sink samples found for stage '", current_stage, "'.")

  # [A] optional continuous source rows
  cont_rows <- character(0)
  if (allow_continuous_sources) {
    cont_rows <- rownames(metadata)[metadata$Env %in% source_types &
                                      metadata$SourceSink == "Source"]
    if (length(cont_rows) == 0)
      warning("allow_continuous_sources = TRUE but no Source rows found ",
              "for source_types. Falling back to standard transfer.")
  }

  inherited_label <- paste0(previous_stage, "_Inherited")
  per_sample <- list()

  for (s in curr_rows) {
    use_rows <- c(prev_rows, cont_rows, s)
    sub_C    <- abundance_table[use_rows, , drop = FALSE]
    sub_meta <- metadata[use_rows, , drop = FALSE]

    # Label previous stage as "Inherited" source
    sub_meta$SourceSink[rownames(sub_meta) %in% prev_rows] <- "Source"
    sub_meta$Env[rownames(sub_meta) %in% prev_rows]        <- inherited_label

    # [A] continuous sources keep their original Env labels + SourceSink=Source
    if (length(cont_rows) > 0) {
      sub_meta$SourceSink[rownames(sub_meta) %in% cont_rows] <- "Source"
      # Env labels already set in original metadata
    }

    sub_meta$SourceSink[rownames(sub_meta) == s] <- "Sink"

    tryCatch({
      res <- .FEAST_internal(sub_C, sub_meta,
                             n_restarts  = n_restarts,
                             max_iter    = max_iter,
                             tol         = tol,
                             rarefaction = rarefaction)
      vec <- res[1, ]
      per_sample[[s]] <- vec
    }, error = function(e) {
      warning("Transfer '", previous_stage, "'->'", current_stage,
              "', sample '", s, "' failed: ", conditionMessage(e))
      per_sample[[s]] <<- NULL
    })
  }

  valid <- Filter(Negate(is.null), per_sample)
  if (length(valid) == 0) {
    warning("All samples failed in transfer to stage '", current_stage, "'.")
    dummy_names <- c(inherited_label,
                     if (allow_continuous_sources) source_types else NULL,
                     "Unknown")
    dummy <- rep(0, length(dummy_names))
    dummy[length(dummy)] <- 1
    names(dummy) <- dummy_names
    return(list(per_sample = per_sample, mean_profile = dummy))
  }
  profile_mat  <- do.call(rbind, valid)
  mean_profile <- colMeans(profile_mat, na.rm = TRUE)
  list(per_sample = per_sample, mean_profile = mean_profile)
}

# ── Cumulative update ─────────────────────────────────────────

.update_cumulative <- function(prev_cumulative, transfer_profile,
                               current_stage, source_types,
                               allow_continuous_sources) {

  inherited_key <- grep("_Inherited$", names(transfer_profile), value = TRUE)
  p_inherited   <- if (length(inherited_key) > 0) transfer_profile[[inherited_key]] else 0
  p_unknown     <- transfer_profile[["Unknown"]]
  if (is.null(p_unknown) || is.na(p_unknown)) p_unknown <- 0

  # Multiply existing contributions by the inherited fraction
  new_cumulative <- prev_cumulative * p_inherited

  # New environmental contribution at this stage
  env_key <- paste0("Env_", current_stage)
  new_cumulative[[env_key]] <- p_unknown

  # [A] If continuous sources present, record their direct contributions
  if (allow_continuous_sources) {
    for (src in source_types) {
      if (src %in% names(transfer_profile)) {
        cont_key <- paste0("Cont_", src, "_", current_stage)
        new_cumulative[[cont_key]] <- transfer_profile[[src]]
      }
    }
  }

  new_cumulative
}

# ── Save helper ───────────────────────────────────────────────

.save_results <- function(per_stage_results, cumulative_contributions,
                          output_dir) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  for (stage in names(per_stage_results)) {
    valid <- Filter(Negate(is.null), per_stage_results[[stage]])
    if (length(valid) == 0) next
    df <- do.call(rbind, lapply(names(valid), function(nm) {
      data.frame(Sample = nm, t(valid[[nm]]),
                 check.names = FALSE, stringsAsFactors = FALSE)
    }))
    utils::write.csv(df,
                     file.path(output_dir, paste0("per_sample_", stage, ".csv")),
                     row.names = FALSE)
  }

  cumul_df <- do.call(rbind, lapply(names(cumulative_contributions), function(st) {
    data.frame(Stage = st, t(cumulative_contributions[[st]]),
               check.names = FALSE, stringsAsFactors = FALSE)
  }))
  utils::write.csv(cumul_df,
                   file.path(output_dir, "cumulative_contributions.csv"),
                   row.names = FALSE)
}

# ── Main exported function ────────────────────────────────────

#' Run dynamic time-series microbial source tracking (FEASTime)
#'
#' Extends the FEAST EM algorithm to time-series data using a
#' Markov chain-based iterative framework.
#'
#' @param abundance_table Numeric matrix (samples x OTUs) with sample
#'   IDs as row names.  Integer count data recommended.
#' @param metadata Data frame with columns \code{Env},
#'   \code{SourceSink} ("Source"/"Sink"), and \code{id}.
#'   Row names must match \code{abundance_table}.
#' @param stage_order Character vector. Ordered stage labels,
#'   e.g. \code{c("S0","S1","S2")}. First element = initial sink stage.
#' @param source_types Character vector. \code{Env} labels of the
#'   original external source groups, e.g. \code{c("Daqu1","Daqu2")}.
#' @param output_dir Character or \code{NULL}. Directory to save CSV
#'   results. Default \code{NULL}.
#' @param n_restarts Integer. EM random restarts per sample (default 10).
#' @param max_iter Integer. Max EM iterations per restart (default 1000).
#' @param tol Numeric. |delta loglik| convergence threshold (default 1e-6).
#' @param rarefaction_depth Integer or \code{NULL}. If provided, all
#'   samples are rarefied to this sequencing depth before analysis.
#'   \code{NULL} (default) skips rarefaction.
#'   Set to \code{min(rowSums(abundance_table))} as a safe default when
#'   samples have unequal sequencing depths.
#' @param allow_continuous_sources Logical. If \code{TRUE}, the original
#'   \code{source_types} samples are included as additional sources at
#'   every transfer step (models ongoing environmental inoculation).
#'   Default \code{FALSE}.
#' @param verbose Logical. Print progress messages. Default \code{TRUE}.
#'
#' @return A named list:
#' \describe{
#'   \item{per_stage_results}{Per-sample FEAST results at each stage.}
#'   \item{mean_profiles}{Mean source proportions per stage.}
#'   \item{cumulative_contributions}{Cumulative contributions traced to
#'     original sources at each stage.}
#'   \item{stage_order}{Input stage_order.}
#'   \item{source_types}{Input source_types.}
#'   \item{call_info}{Analysis parameters and timestamp.}
#' }
#' @export
#' @examples
#' \dontrun{
#' set.seed(42)
#' abu <- matrix(sample(0:500, 20 * 50, replace = TRUE), 20, 50)
#' rownames(abu) <- c(paste0("Src1_", 1:3), paste0("Src2_", 1:3),
#'                    paste0("S0_", 1:4), paste0("S1_", 1:4),
#'                    paste0("S2_", 1:4), paste0("S3_", 1:4))
#' colnames(abu) <- paste0("OTU_", 1:50)
#' meta <- data.frame(
#'   Env        = c(rep("Source1", 3), rep("Source2", 3),
#'                  rep("S0", 4), rep("S1", 4), rep("S2", 4), rep("S3", 4)),
#'   SourceSink = c(rep("Source", 6), rep("Sink", 16)),
#'   id         = rownames(abu),
#'   row.names  = rownames(abu), stringsAsFactors = FALSE
#' )
#' # Standard run
#' result <- run_feastime(abu, meta, c("S0","S1","S2","S3"),
#'                        c("Source1","Source2"), n_restarts = 3)
#'
#' # With rarefaction and continuous sources
#' result2 <- run_feastime(abu, meta, c("S0","S1","S2","S3"),
#'                         c("Source1","Source2"), n_restarts = 3,
#'                         rarefaction_depth = 300,
#'                         allow_continuous_sources = TRUE)
#' }
run_feastime <- function(abundance_table, metadata,
                         stage_order, source_types,
                         output_dir              = NULL,
                         n_restarts              = 10,
                         max_iter                = 1000,
                         tol                     = 1e-6,
                         rarefaction_depth       = NULL,
                         allow_continuous_sources = FALSE,
                         verbose                 = TRUE) {

  .validate_inputs(abundance_table, metadata, stage_order, source_types)
  if (is.data.frame(abundance_table))
    abundance_table <- as.matrix(abundance_table)

  .msg <- function(...) if (verbose) message(...)

  per_stage_results        <- list()
  mean_profiles            <- list()
  cumulative_contributions <- list()

  # ── Step 1: original sources -> first sink stage ──
  first_stage <- stage_order[1]
  .msg("[FEASTime] Step 1: ", paste(source_types, collapse = "+"),
       " -> ", first_stage)

  init_res <- .analyse_initial_stage(abundance_table, metadata,
                                     first_stage, source_types,
                                     n_restarts, max_iter, tol,
                                     rarefaction_depth)
  per_stage_results[[first_stage]]        <- init_res$per_sample
  mean_profiles[[first_stage]]            <- init_res$mean_profile
  cumulative_contributions[[first_stage]] <- init_res$mean_profile
  .msg("[FEASTime]   Done.")

  # ── Step 2: sequential Markov transfers ──
  if (length(stage_order) > 1) {
    .msg("[FEASTime] Step 2: sequential transfers",
         if (allow_continuous_sources) " [continuous sources ON]" else "")

    current_cumulative <- cumulative_contributions[[first_stage]]

    for (i in 2:length(stage_order)) {
      prev_stage <- stage_order[i - 1]
      curr_stage <- stage_order[i]
      .msg("[FEASTime]   ", prev_stage, " -> ", curr_stage)

      tr_res <- .analyse_transfer(abundance_table, metadata,
                                  prev_stage, curr_stage,
                                  source_types          = source_types,
                                  allow_continuous_sources = allow_continuous_sources,
                                  n_restarts            = n_restarts,
                                  max_iter              = max_iter,
                                  tol                   = tol,
                                  rarefaction           = rarefaction_depth)

      per_stage_results[[curr_stage]] <- tr_res$per_sample
      mean_profiles[[curr_stage]]     <- tr_res$mean_profile

      new_cumulative <- .update_cumulative(current_cumulative,
                                           tr_res$mean_profile,
                                           curr_stage,
                                           source_types,
                                           allow_continuous_sources)
      cumulative_contributions[[curr_stage]] <- new_cumulative
      current_cumulative <- new_cumulative

      inherited_key <- grep("_Inherited$", names(tr_res$mean_profile), value = TRUE)
      p_inh <- if (length(inherited_key) > 0) round(tr_res$mean_profile[[inherited_key]], 3) else NA
      .msg("[FEASTime]   P_inherited=", p_inh,
           "  P_unknown=", round(tr_res$mean_profile[["Unknown"]], 3))
    }
  }

  if (!is.null(output_dir)) {
    .save_results(per_stage_results, cumulative_contributions, output_dir)
    .msg("[FEASTime] Results saved to: ", output_dir)
  }
  .msg("[FEASTime] Complete.")

  list(
    per_stage_results        = per_stage_results,
    mean_profiles            = mean_profiles,
    cumulative_contributions = cumulative_contributions,
    stage_order              = stage_order,
    source_types             = source_types,
    call_info = list(
      n_restarts               = n_restarts,
      max_iter                 = max_iter,
      tol                      = tol,
      rarefaction_depth        = rarefaction_depth,
      allow_continuous_sources = allow_continuous_sources,
      timestamp                = Sys.time()
    )
  )
}
