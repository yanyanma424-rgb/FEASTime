\# FEASTime <img src="man/figures/logo.png" align="right" height="120" alt="" />



<!-- badges: start -->

!\[R](https://img.shields.io/badge/R-%3E%3D3.6.0-blue)

!\[License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)

<!-- badges: end -->



\*\*FEASTime\*\* extends the \[FEAST](https://github.com/cozygene/FEAST) microbial

source-tracking algorithm to \*\*longitudinal time-series data\*\* using a

Markov chain-based iterative EM framework.



Instead of a single static attribution, FEASTime tracks how microbial

communities are inherited, replaced, and newly introduced across successive

time points — making it well-suited for fermentation, gut development, or

any process with ordered sampling stages.



---



\## Key features



\- \*\*Dynamic Markov chain tracking\*\* — each stage is simultaneously a sink

&nbsp; (receiving microbes) and a source for the next stage

\- \*\*Cumulative contribution tracing\*\* — quantify how much of the final

&nbsp; community traces back to each original source

\- \*\*Rarefaction support\*\* — optional depth normalisation before EM

&nbsp; (`rarefaction\_depth` parameter)

\- \*\*Biologically-informed Unknown initialisation\*\* — mirrors the

&nbsp; `unknown\_initialize\_1` heuristic from FEAST for faster, more stable

&nbsp; convergence on sparse OTU tables

\- \*\*Continuous source exposure\*\* — optional `allow\_continuous\_sources`

&nbsp; flag models ongoing environmental inoculation at every transfer step

\- \*\*Built-in visualisation\*\* — stacked bar and trend-line plots



---



\## Installation

```r

\# Install from GitHub (requires devtools)

devtools::install\_github("MaYanyan/FEASTime")

```



---



\## Quick start

```r

library(FEASTime)



\# abundance\_table: samples x OTUs integer count matrix (row names = sample IDs)

\# metadata:        data.frame with columns Env, SourceSink ("Source"/"Sink"), id



result <- run\_feastime(

&nbsp; abundance\_table  = abu,

&nbsp; metadata         = meta,

&nbsp; stage\_order      = c("S0", "S1", "S2", "S3"),

&nbsp; source\_types     = c("Daqu1", "Daqu2"),

&nbsp; rarefaction\_depth = 5000,   # recommended when depths differ

&nbsp; n\_restarts       = 10

)



\# Visualise

plot\_contributions(result)

plot\_trends(result)



\# Summary table

summarise\_feastime(result)

```



---



\## Input format



\### `abundance\_table`



| | OTU\_1 | OTU\_2 | … |

|---|---|---|---|

| Daqu1\_rep1 | 120 | 0 | … |

| S0\_rep1 | 45 | 230 | … |



Rows = samples, columns = OTUs. Integer counts recommended.



\### `metadata`



| | Env | SourceSink | id |

|---|---|---|---|

| Daqu1\_rep1 | Daqu1 | Source | Daqu1\_rep1 |

| S0\_rep1 | S0 | Sink | S0\_rep1 |



Row names must match `abundance\_table` row names exactly.



---



\## Parameters



| Parameter | Default | Description |

|---|---|---|

| `stage\_order` | required | Ordered vector of sink stage labels |

| `source\_types` | required | Env labels of original external sources |

| `rarefaction\_depth` | `NULL` | Rarefy all samples to this depth before EM |

| `allow\_continuous\_sources` | `FALSE` | Include original sources at every transfer step |

| `n\_restarts` | `10` | EM random restarts per sample |

| `tol` | `1e-6` | \\|Δ log-likelihood\\| convergence threshold |

| `output\_dir` | `NULL` | Directory to save per-stage CSV results |



---



\## Citation



If you use FEASTime, please cite:



> Ma Y. et al. (2025). FEASTime: Extending FEAST for microbial source

> tracking in multi-source and multi-sink time-series samples. \*(submitted)\*



> Shenhav L. et al. (2019). FEAST: fast expectation-maximization for

> microbial source tracking. \*Nature Methods\*, 16, 627–632.



---



\## License



MIT © 2025 Ma Yanyan

