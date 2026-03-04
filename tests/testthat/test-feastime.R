library(testthat)
library(FEASTime)

# ---- helper: small reproducible dataset --------------------
make_test_data <- function(seed = 1) {
  set.seed(seed)
  n_otus <- 30
  abu <- matrix(sample(0:200, 14 * n_otus, replace = TRUE),
                nrow = 14, ncol = n_otus)
  rownames(abu) <- c("Src1_1","Src1_2","Src2_1","Src2_2",
                     "S0_1","S0_2","S0_3",
                     "S1_1","S1_2","S1_3",
                     "S2_1","S2_2","S2_3","S2_4")
  colnames(abu) <- paste0("OTU_", seq_len(n_otus))

  meta <- data.frame(
    Env        = c(rep("Source1",2), rep("Source2",2),
                   rep("S0",3), rep("S1",3), rep("S2",4)),
    SourceSink = c(rep("Source",4), rep("Sink",10)),
    id         = rownames(abu),
    row.names  = rownames(abu),
    stringsAsFactors = FALSE
  )
  list(abu = abu, meta = meta)
}

# ---- tests -------------------------------------------------
test_that("run_feastime returns correct structure", {
  d <- make_test_data()
  res <- run_feastime(d$abu, d$meta,
                      stage_order  = c("S0","S1","S2"),
                      source_types = c("Source1","Source2"),
                      n_restarts = 2, verbose = FALSE)

  expect_type(res, "list")
  expect_true(all(c("per_stage_results","mean_profiles",
                    "cumulative_contributions",
                    "stage_order","source_types","call_info") %in% names(res)))
  expect_equal(length(res$cumulative_contributions), 3)
  expect_true(all(c("S0","S1","S2") %in% names(res$cumulative_contributions)))
})

test_that("cumulative contributions at S0 sum to ~1", {
  d <- make_test_data()
  res <- run_feastime(d$abu, d$meta,
                      stage_order  = c("S0","S1","S2"),
                      source_types = c("Source1","Source2"),
                      n_restarts = 2, verbose = FALSE)
  s0_sum <- sum(res$cumulative_contributions[["S0"]], na.rm = TRUE)
  expect_true(abs(s0_sum - 1) < 0.05)
})

test_that("input validation catches missing metadata column", {
  d <- make_test_data()
  bad_meta <- d$meta[, c("Env","SourceSink")]   # drop 'id'
  expect_error(
    run_feastime(d$abu, bad_meta, c("S0","S1","S2"),
                 c("Source1","Source2"), verbose = FALSE),
    "missing columns"
  )
})

test_that("input validation catches bad SourceSink value", {
  d <- make_test_data()
  d$meta$SourceSink[1] <- "Neither"
  expect_error(
    run_feastime(d$abu, d$meta, c("S0","S1","S2"),
                 c("Source1","Source2"), verbose = FALSE),
    "invalid values"
  )
})

test_that("summarise_feastime returns a data.frame", {
  d <- make_test_data()
  res <- run_feastime(d$abu, d$meta,
                      stage_order  = c("S0","S1","S2"),
                      source_types = c("Source1","Source2"),
                      n_restarts = 2, verbose = FALSE)
  s <- summarise_feastime(res)
  expect_s3_class(s, "data.frame")
  expect_true(all(c("Stage","Source","Mean","SD") %in% colnames(s)))
})
