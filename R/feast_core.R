# ============================================================
# feast_core.R  --  FEASTime v2
#
# Changes vs v1:
#   [1] .rarefy(): multinomial rarefaction to equalise depth
#   [2] .init_unknown_bioinformed(): replaces runif() Unknown init
#       with sink-source residual logic (FEAST unknown_initialize_1)
#   [3] .init_params(): now receives y_sink + x_source for [2]
#   [4] .M_step(): zero-weight guard — retain p_unknown_prev
#       instead of resetting to uniform (prevents oscillation)
#   [5] .run_EM(): convergence on |delta loglik|, not delta_rho
#   [6] .FEAST_internal(): rarefaction parameter added
# ============================================================

# ── [1] Rarefaction ──────────────────────────────────────────

#' Rarefy a count matrix to a target sequencing depth
#'
#' Each sample with rowSum > target_depth is randomly sub-sampled
#' (multinomial, weights = observed proportions) down to target_depth.
#' Samples already at or below target_depth are left unchanged.
#' Mirrors FEAST_rarefy() from Shenhav et al.
#'
#' @param C            Integer matrix, samples x OTUs.
#' @param target_depth Integer. Target sequencing depth.
#' @return Integer matrix of the same dimensions as C.
#' @keywords internal
.rarefy <- function(C, target_depth) {
  n_otus <- ncol(C)
  for (i in seq_len(nrow(C))) {
    row_total <- sum(C[i, ])
    if (row_total > target_depth) {
      s       <- sample(n_otus, size = target_depth,
                        prob = C[i, ] / row_total, replace = TRUE)
      C[i, ] <- tabulate(s, nbins = n_otus)
    }
  }
  C
}

# ── [2] Biologically-informed Unknown initialisation ─────────

#' Initialise the Unknown source profile from sink-source residuals
#'
#' Mirrors FEAST unknown_initialize_1 (Shenhav et al.):
#' \itemize{
#'   \item Per-OTU unknown = max(sink_j - sum_of_sources_j, 0)
#'   \item For "core" OTUs (present in >= 50\% of sources):
#'         unknown = min(source_abundance_across_sources) / 2
#' }
#' On highly sparse OTU tables (>80\% zeros), this is far more stable
#' than uniform random initialisation and accelerates EM convergence.
#'
#' @param y_sink   Numeric vector, length n_otus. Sink proportions.
#' @param x_source Matrix, n_otus x K. Source probability profiles.
#' @return Numeric vector, length n_otus, sums to 1.
#' @keywords internal
.init_unknown_bioinformed <- function(y_sink, x_source) {
  n_otus    <- length(y_sink)
  n_sources <- ncol(x_source)

  # Per-OTU excess in sink over all known sources combined
  sources_sum  <- rowSums(x_source)
  unknown_init <- pmax(y_sink - sources_sum, 0)

  # Soften estimates for core OTUs (present in >= 50% of sources)
  core_mask <- rowSums(x_source > 0) >= (0.5 * n_sources)
  if (any(core_mask) && n_sources > 1) {
    min_source              <- apply(x_source[core_mask, , drop = FALSE], 1, min)
    unknown_init[core_mask] <- min_source / 2
  }

  total <- sum(unknown_init)
  if (total <= 0) return(rep(1 / n_otus, n_otus))   # rare edge-case fallback
  unknown_init / total
}

# ── [3] Initialise EM parameters ─────────────────────────────

#' Initialise EM parameters
#'
#' Mixing proportions (rho) are random uniform.
#' Unknown source profile uses biologically-informed initialisation [2].
#'
#' @param n_sources Integer. Number of known sources.
#' @param n_otus    Integer. Number of OTUs.
#' @param y_sink    Numeric vector, length n_otus. Sink proportions.
#' @param x_source  Matrix, n_otus x K. Source profiles.
#' @return List: \code{rho} (length K+1), \code{p_unknown} (length n_otus).
#' @importFrom stats runif
#' @keywords internal
.init_params <- function(n_sources, n_otus, y_sink, x_source) {
  rho <- stats::runif(n_sources + 1)
  rho <- rho / sum(rho)

  p_unknown <- .init_unknown_bioinformed(y_sink, x_source)

  list(rho = rho, p_unknown = p_unknown)
}

# ── E-step ────────────────────────────────────────────────────

#' E-step: compute posterior latent source assignments
#'
#' @param y_sink    Numeric vector, length n_otus. Sink proportions.
#' @param x_source  Matrix, n_otus x K. Source profiles.
#' @param rho       Numeric vector, length K+1. Mixing proportions.
#' @param p_unknown Numeric vector, length n_otus. Unknown source profile.
#' @return Matrix, n_otus x (K+1). Row-normalised posterior weights.
#' @keywords internal
.E_step <- function(y_sink, x_source, rho, p_unknown) {
  full_source <- cbind(x_source, p_unknown)        # n_otus x (K+1)
  Z           <- sweep(full_source, 2, rho, `*`)   # unnormalised
  row_sums    <- rowSums(Z)
  row_sums[row_sums == 0] <- 1                      # guard against /0
  Z / row_sums
}

# ── [4] M-step with zero-weight guard ────────────────────────

#' M-step: update mixing proportions and Unknown source profile
#'
#' Zero-weight guard: when Unknown receives zero expected weight in an
#' iteration, \code{p_unknown_prev} is retained rather than resetting to
#' a uniform distribution.  Resetting to uniform causes severe oscillation
#' on the sparse OTU tables typical of microbiome data.
#'
#' @param y_sink         Numeric vector, length n_otus.
#' @param Z              Matrix, n_otus x (K+1). E-step posterior weights.
#' @param p_unknown_prev Numeric vector, length n_otus. Previous iteration.
#' @return List: updated \code{rho} and \code{p_unknown}.
#' @keywords internal
.M_step <- function(y_sink, Z, p_unknown_prev) {
  n_full          <- ncol(Z)
  expected_counts <- colSums(y_sink * Z)
  total_counts    <- sum(expected_counts)
  if (total_counts == 0) total_counts <- 1

  rho_new <- expected_counts / total_counts

  unknown_total <- sum(y_sink * Z[, n_full])
  p_unknown_new <- if (unknown_total > 0) {
    (y_sink * Z[, n_full]) / unknown_total
  } else {
    p_unknown_prev    # [4] retain previous; no uniform reset
  }

  list(rho = rho_new, p_unknown = p_unknown_new)
}

# ── [5] Log-likelihood helper ─────────────────────────────────

#' Compute mixture log-likelihood
#'
#' Used as the convergence criterion in .run_EM() [5].
#'
#' @param y_sink    Numeric vector.
#' @param x_source  Matrix, n_otus x K.
#' @param rho       Numeric vector, length K+1.
#' @param p_unknown Numeric vector, length n_otus.
#' @return Scalar log-likelihood value.
#' @keywords internal
.compute_loglik <- function(y_sink, x_source, rho, p_unknown) {
  full_source <- cbind(x_source, p_unknown)
  mix_prob    <- as.numeric(full_source %*% rho)
  mix_prob[mix_prob <= 0] <- .Machine$double.eps
  sum(y_sink * log(mix_prob))
}

# ── [5] Full EM with loglik convergence ───────────────────────

#' Run EM algorithm for one sink sample
#'
#' Multiple random restarts; the solution with the highest log-likelihood
#' is returned.
#'
#' Convergence criterion: \code{|loglik_new - loglik_old| < tol}.
#' This is statistically correct for EM — v1 used \code{delta_rho},
#' which can declare convergence while \code{p_unknown} is still changing.
#'
#' @param y_sink     Numeric vector, length n_otus. Sink proportions.
#' @param x_source   Matrix, n_otus x K. Source profiles.
#' @param n_restarts Integer. Number of random restarts (default 10).
#' @param max_iter   Integer. Max EM iterations per restart (default 1000).
#' @param tol        Numeric. |delta loglik| convergence threshold (default 1e-6).
#' @return List: \code{rho}, \code{p_unknown}, \code{iterations}, \code{loglik}.
#' @keywords internal
.run_EM <- function(y_sink, x_source, n_restarts = 10,
                    max_iter = 1000, tol = 1e-6) {

  n_otus    <- length(y_sink)
  n_sources <- ncol(x_source)

  best_loglik <- -Inf
  best_result <- NULL

  for (restart in seq_len(n_restarts)) {
    # [3] informed init
    params    <- .init_params(n_sources, n_otus, y_sink, x_source)
    rho       <- params$rho
    p_unknown <- params$p_unknown

    prev_loglik <- -Inf

    for (iter in seq_len(max_iter)) {
      Z          <- .E_step(y_sink, x_source, rho, p_unknown)
      new_params <- .M_step(y_sink, Z, p_unknown)   # [4] pass prev p_unknown
      rho        <- new_params$rho
      p_unknown  <- new_params$p_unknown

      # [5] loglik convergence
      cur_loglik <- .compute_loglik(y_sink, x_source, rho, p_unknown)
      if (abs(cur_loglik - prev_loglik) < tol) break
      prev_loglik <- cur_loglik
    }

    loglik <- .compute_loglik(y_sink, x_source, rho, p_unknown)
    if (loglik > best_loglik) {
      best_loglik <- loglik
      best_result <- list(rho        = rho,
                          p_unknown  = p_unknown,
                          iterations = iter,
                          loglik     = loglik)
    }
  }

  best_result
}

# ── Source profile preparation ────────────────────────────────

#' Prepare source OTU probability profiles
#'
#' Groups source samples by Env label, averages raw counts within each
#' group, then normalises to probability profiles (each column sums to 1).
#'
#' @param C_sources   Matrix, n_source_samples x n_otus. Raw counts.
#' @param source_envs Character vector. Group label per source sample.
#' @return Matrix, n_otus x K. Probability profiles, columns sum to 1.
#' @keywords internal
.prepare_source_profiles <- function(C_sources, source_envs) {
  unique_sources <- unique(source_envs)
  n_otus         <- ncol(C_sources)

  profile_matrix <- matrix(0, nrow = n_otus, ncol = length(unique_sources),
                           dimnames = list(colnames(C_sources), unique_sources))

  for (j in seq_along(unique_sources)) {
    src         <- unique_sources[j]
    idx         <- which(source_envs == src)
    mean_counts <- colMeans(C_sources[idx, , drop = FALSE])
    total       <- sum(mean_counts)
    profile_matrix[, j] <- if (total > 0) {
      mean_counts / total
    } else {
      rep(1 / n_otus, n_otus)
    }
  }
  profile_matrix
}

# ── [6] Single-call FEAST wrapper with rarefaction ───────────

#' Run FEAST source tracking for one or more sink samples
#'
#' Internal workhorse called by the dynamic analysis functions.
#' For each sink sample the EM algorithm is run against the provided
#' source profiles.
#'
#' @param C           Matrix, n_samples x n_otus. OTU count table.
#'   Row names must match \code{metadata} row names.
#' @param metadata    Data frame. Must contain columns \code{Env},
#'   \code{SourceSink} ("Source"/"Sink"), and \code{id}.
#' @param n_restarts  Integer. EM random restarts (default 10).
#' @param max_iter    Integer. Max EM iterations (default 1000).
#' @param tol         Numeric. Convergence tolerance (default 1e-6).
#' @param rarefaction Integer or \code{NULL}. If an integer, all samples
#'   are rarefied to this depth before analysis.
#'   \code{NULL} (default) skips rarefaction.
#'   Recommended: set to \code{min(rowSums(C))} when samples differ in
#'   sequencing depth.
#' @return Matrix, n_sinks x (K+1). Each row sums to 1;
#'   last column is "Unknown".
#' @keywords internal
.FEAST_internal <- function(C, metadata,
                            n_restarts  = 10,
                            max_iter    = 1000,
                            tol         = 1e-6,
                            rarefaction = NULL) {

  # ── input validation ──
  if (!all(rownames(C) == rownames(metadata)))
    stop("Row names of C and metadata must be identical and in the same order.")

  missing_cols <- setdiff(c("Env", "SourceSink", "id"), colnames(metadata))
  if (length(missing_cols) > 0)
    stop("metadata is missing columns: ", paste(missing_cols, collapse = ", "))

  sources_idx <- which(metadata$SourceSink == "Source")
  sinks_idx   <- which(metadata$SourceSink == "Sink")
  if (length(sources_idx) == 0) stop("No Source samples found in metadata.")
  if (length(sinks_idx)   == 0) stop("No Sink samples found in metadata.")

  # ── [6] rarefaction ──
  if (!is.null(rarefaction)) {
    # Require integer-like input
    if (any(C != floor(C))) {
      warning("Non-integer values detected. Converting via ceiling() before rarefaction.")
      C <- ceiling(C)
    }
    row_depths <- rowSums(C)
    if (any(row_depths < rarefaction))
      warning("Some samples have fewer reads than rarefaction_depth (", rarefaction,
              "). Those samples will not be rarefied further.")
    C <- .rarefy(C, rarefaction)
  } else {
    # Auto-detect relative abundance input and warn
    if (max(C, na.rm = TRUE) <= 1 && min(C, na.rm = TRUE) >= 0)
      message("[FEASTime] Input appears to be relative abundance (all values 0-1). ",
              "Provide integer count data and set rarefaction_depth for best results.")
  }

  C_sources   <- C[sources_idx, , drop = FALSE]
  C_sinks     <- C[sinks_idx,   , drop = FALSE]
  source_envs <- metadata$Env[sources_idx]

  source_profiles <- .prepare_source_profiles(C_sources, source_envs)
  unique_sources  <- colnames(source_profiles)
  n_sources <- length(unique_sources)
  n_sinks   <- nrow(C_sinks)

  results <- matrix(0, nrow = n_sinks, ncol = n_sources + 1,
                    dimnames = list(rownames(C_sinks),
                                    c(unique_sources, "Unknown")))

  for (i in seq_len(n_sinks)) {
    y_sink     <- as.numeric(C_sinks[i, ])
    total_sink <- sum(y_sink)
    if (total_sink == 0) {
      warning("Sink sample '", rownames(C_sinks)[i],
              "' has zero total counts. Filling with NA.")
      results[i, ] <- NA
      next
    }
    y_prop <- y_sink / total_sink

    tryCatch({
      em_out       <- .run_EM(y_prop, source_profiles,
                              n_restarts = n_restarts,
                              max_iter   = max_iter,
                              tol        = tol)
      results[i, ] <- em_out$rho
    }, error = function(e) {
      warning("EM failed for sink '", rownames(C_sinks)[i],
              "': ", conditionMessage(e))
      results[i, ] <- NA
    })
  }

  results
}
