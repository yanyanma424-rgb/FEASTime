# FEASTime

![R](https://img.shields.io/badge/R-%3E%3D3.6.0-blue)
![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

FEASTime extends the [FEAST](https://github.com/cozygene/FEAST) microbial source-tracking algorithm to **longitudinal time-series data** using a Markov chain-based iterative EM framework. Instead of a single static attribution, FEASTime tracks how microbial communities are inherited, replaced, and newly introduced across successive time points — making it well-suited for fermentation, gut development, or any process with ordered sampling stages.

If you use this package, please cite [Ma Y. (2025).](#citation) pending publication.

## Installation

The easiest way to install FEASTime is via `devtools`:
```r
devtools::install_github("MaYanyan/FEASTime")
```

## Quick Start

FEASTime requires two inputs: an OTU count table and a metadata table.
```r
library(FEASTime)

# Load example data
abu <- read.csv(system.file("extdata", "example_abundance.csv", package = "FEASTime"), row.names = 1)
meta <- read.csv(system.file("extdata", "example_metadata.csv", package = "FEASTime"), row.names = 1)

# Run FEASTime
result <- run_feastime(
  abundance_table   = abu,
  metadata          = meta,
  stage_order       = c("S0", "S1", "S2"),
  source_types      = c("Source1", "Source2"),
  rarefaction_depth = 100,
  n_restarts        = 10
)

# Visualise
plot_contributions(result)
plot_trends(result)

# Summary table
summarise_feastime(result)
```

## File Formats

### abundance_table

A samples × OTUs integer count matrix with sample IDs as row names. Values must be non-negative integers.
```
        OTU_1  OTU_2  OTU_3  OTU_4  OTU_5
Src1_1    300    120      0     45    200
Src1_2    280    110      0     50    210
Src2_1     10     20    400    300      5
S0_1      150     80    100    120    130
S1_1      100     50    150    130     90
S2_1       60     30    200    180     60
```

### metadata

A data frame with sample IDs as row names and three required columns:
```
        Env      SourceSink   id
Src1_1  Source1  Source       Src1_1
Src1_2  Source1  Source       Src1_2
Src2_1  Source2  Source       Src2_1
S0_1    S0       Sink         S0_1
S1_1    S1       Sink         S1_1
S2_1    S2       Sink         S2_1
```

- `Env`: environment/group label matching `stage_order` or `source_types`
- `SourceSink`: must be exactly `"Source"` or `"Sink"`
- `id`: sample ID, must match row names

Row names of `abundance_table` and `metadata` must be identical and in the same order.

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `stage_order` | required | Ordered vector of sink stage labels, e.g. `c("S0","S1","S2")` |
| `source_types` | required | Env labels of original external sources, e.g. `c("Source1","Source2")` |
| `rarefaction_depth` | `NULL` | Rarefy all samples to this depth before EM. Recommended: `min(rowSums(abundance_table))` |
| `allow_continuous_sources` | `FALSE` | If `TRUE`, include original sources at every transfer step to model ongoing environmental inoculation |
| `n_restarts` | `10` | Number of EM random restarts per sample |
| `max_iter` | `1000` | Maximum EM iterations per restart |
| `tol` | `1e-6` | \|Δ log-likelihood\| convergence threshold |
| `output_dir` | `NULL` | Directory to save per-stage CSV results (created if needed) |
| `verbose` | `TRUE` | Print progress messages |

## Output

`run_feastime()` returns a named list:

- **`per_stage_results`** — per-sample source proportions at each stage
- **`mean_profiles`** — mean source proportions per stage
- **`cumulative_contributions`** — cumulative contributions traced back to original sources at each stage
- **`stage_order`** — input stage order
- **`source_types`** — input source types
- **`call_info`** — analysis parameters and timestamp

If `output_dir` is set, two types of CSV files are saved:
- `per_sample_<stage>.csv` — per-sample results at each stage
- `cumulative_contributions.csv` — cumulative contributions across all stages

## Theory

FEASTime models microbial community succession as a Markov chain. At each time point T:

1. **Initial stage** — the EM algorithm estimates what proportion of the sink community at T=0 comes from each original source (e.g. Daqu types in baijiu fermentation).
2. **Transfer stages** — for each subsequent stage, the previous stage's community is treated as the known source. The EM estimates what fraction was inherited from T-1 versus newly introduced from the environment (Unknown).
3. **Cumulative tracing** — by multiplying inherited fractions across stages, FEASTime traces what proportion of the final community can be attributed to each original source.

This differs from standard FEAST, which performs independent static attribution at a single time point. FEASTime fills the gap for longitudinal study designs.

## Usage Examples

**Standard run:**
```r
result <- run_feastime(abu, meta,
                       stage_order  = c("S0", "S1", "S2"),
                       source_types = c("Source1", "Source2"),
                       n_restarts   = 10)
```

**With rarefaction:**
```r
result <- run_feastime(abu, meta,
                       stage_order       = c("S0", "S1", "S2"),
                       source_types      = c("Source1", "Source2"),
                       rarefaction_depth = min(rowSums(abu)),
                       n_restarts        = 10)
```

**With continuous source exposure** (models ongoing environmental inoculation at every stage):
```r
result <- run_feastime(abu, meta,
                       stage_order              = c("S0", "S1", "S2"),
                       source_types             = c("Source1", "Source2"),
                       allow_continuous_sources = TRUE,
                       n_restarts               = 10)
```

**Save results to disk:**
```r
result <- run_feastime(abu, meta,
                       stage_order  = c("S0", "S1", "S2"),
                       source_types = c("Source1", "Source2"),
                       output_dir   = "my_results/")
```

## Citation

If you use FEASTime, please cite:

> Ma Y. (2025). *FEASTime: Extending FEAST for microbial source tracking in time-series samples.* R package version 0.1.0. https://github.com/yanyanma424-rgb/FEASTime

Also cite the original FEAST paper:

> Shenhav L. et al. (2019). FEAST: fast expectation-maximization for microbial source tracking. *Nature Methods*, 16, 627–632.

## License

MIT © 2025 Ma Yanyan