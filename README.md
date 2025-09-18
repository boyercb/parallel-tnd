# parallel-tnd

Replication files for:

Boyer, C., Li, K.Q., Shi, X., & Tchetgen Tchetgen, T. J. (2025). “Identification and estimation of vaccine effectiveness in the test-negative design under equi-confounding”. arXiv. https://doi.org/10.48550/arXiv.2504.20360

# Abstract
The test-negative design (TND) is widely used to evaluate vaccine effectiveness in real-world settings. In a TND study, individuals with similar symptoms who seek care are tested, and effectiveness is estimated by comparing vaccination histories of test-positive cases and test-negative controls. The TND is often justified on the grounds that it reduces confounding due to unmeasured health-seeking behavior, although this has not been formally described using potential outcomes. At the same time, concerns persist that conditioning on test receipt can introduce selection bias. We provide a formal justification of the TND under an assumption of odds ratio equi-confounding, where unmeasured confounders affect test-positive and test-negative individuals equivalently on the odds ratio scale. Health-seeking behavior is one plausible example. We also show that these results hold under the outcome-dependent sampling used in TNDs. We discuss the design implications of the equi-confounding assumption and provide alternative estimators for the marginal risk ratio among the vaccinated under equi-confounding, including outcome modeling and inverse probability weighting estimators as well as a semiparametric estimator that is doubly robust. When equi-confounding does not hold, we suggest a straightforward sensitivity analysis that parameterizes the magnitude of the deviation on the odds ratio scale. A simulation study evaluates the empirical performance of our proposed estimators under a wide range of scenarios. Finally, we discuss broader uses of test-negative outcomes to de-bias cohort studies in which testing is triggered by symptoms.

## Requirements

This code requires R version 4.0 or higher. The following R packages are required:

**Core simulation and analysis:**
- `data.table` - for data manipulation and simulation framework
- `progressr` - for progress tracking during simulations
- `readr` - for reading/writing simulation results
- `lmtest` - for linear model testing (coefci function)
- `sandwich` - for robust standard errors (vcovHC function)
- `survival` - for survival analysis functions
- `numDeriv` - for numerical derivatives in estimating equations

**Table and figure generation:**
- `ggplot2` - for plotting simulation results
- `patchwork` - for combining plots
- `tidyverse` - collection of tidy data packages (includes dplyr, tidyr, stringr)
- `dplyr` - for data manipulation in table generation
- `tidyr` - for data tidying
- `stringr` - for string manipulation
- `kableExtra` - for LaTeX table formatting

## Installation

1. Clone this repository:
```bash
git clone https://github.com/boyercb/parallel-tnd.git
cd parallel-tnd
```

2. Install required R packages:
```r
install.packages(c("data.table", "ggplot2", "progressr", "survival", 
                   "lmtest", "sandwich", "readr", "numDeriv",
                   "dplyr", "tidyr", "stringr", "kableExtra", 
                   "patchwork", "tidyverse"))
```

3. Create the data directory (if it doesn't exist):
```r
dir.create("data", showWarnings = FALSE)
```

## Repository Structure

- `code/` - Contains all R scripts for the simulation study
  - `run.R` - Main script that runs all simulation scenarios
  - `datagen.R` - Data generation functions
  - `estimators.R` - Implementation of all estimators
  - `sim.R` - Simulation wrapper functions
  - `table.R` - Code to generate results tables
  - `plot.R` - Code to generate figures
- `data/` - Data files (in this project only saved simulation results)
- `manuscript/` - LaTeX source files for the manuscript
- `results/` - Output files (tables and figures)

## How to run

### Running the full simulation study

To reproduce all simulation results from the paper, run the following command in R from the project root directory:

```r
source("code/run.R")
```

**Note:** The simulation study runs:
- **Scenarios 1-7:** 1,000 replications each with sample size N=15,000 
- **Scenario 8:** 2,000 replications each with sample size N=15,000 across 4 sub-scenarios

This will take several hours to complete (estimated 6-8 hours depending on your system).

### Generating tables and figures

After running the simulations, you can generate the tables and figures using:

```r
# Load simulation results
sims <- readr::read_rds("data/sims.rds")

# Generate tables (creates LaTeX files in results/ directory)
source("code/table.R")

# Generate figures (creates PDF files in results/ directory)
source("code/plot.R")
```

### Individual components

You can also run individual parts of the analysis:

```r
# Load required functions
source("code/datagen.R")
source("code/estimators.R") 
source("code/sim.R")

# Run a single scenario (much faster for testing)
# See run.R for all scenario definitions
```

## Simulation Scenarios

The simulation study evaluates 8 different scenarios:

1. **No unmeasured confounding** - Baseline scenario where TND assumptions hold
2. **Equi-confounding** - Main scenario where proposed methods should work
3. **Direct effect of vaccination on test-negative infection** - Vaccination affects test-negative outcomes (exclusion restriction violated)
4. **Equi-confounding violated** - Confounding differs between test-positive and test-negative
5. **Equi-selection violated** - Selection bias differs between test-positive and test-negative
6. **Equal effect of vaccination on testing** - Equal effects of vaccination on testing behavior
7. **Unequal effect of vaccination on testing** - Unequal effects of vaccination on testing behavior
8. **Effect heterogeneity** - Scenarios with treatment effect modification
   - 8a: Both models correctly specified
   - 8b: Propensity score model misspecified
   - 8c: Outcome model misspecified  
   - 8d: Both models misspecified

Each scenario compares multiple estimators:
- **TND estimators:**
  - Logistic regression (`logit_reg`) - traditional TND approach
  - Risk ratio among vaccinated - outcome modeling (`rrv_om`)
  - Risk ratio among vaccinated - inverse probability weighting (`rrv_ipw`) 
  - Risk ratio among vaccinated - doubly-robust (`rrv_dr`)
- **Cohort estimators:** 
  - Cohort with unmeasured confounders (`cohort_reg_U`) - oracle estimator
  - Cohort without unmeasured confounders (`cohort_reg_noU`) - naive estimator
  - Difference-in-differences (`did_reg`) - equivalent to TND under equi-confounding

## Output

The simulation generates:
- **Performance metrics:** Bias, coverage probability, and confidence interval length for each estimator
- **LaTeX tables:** Saved to `results/` directory:
  - `sims.tex` - Main simulation results
  - `sims_dr*.tex` - Additional robustness results
- **Figures:** PDF files showing estimator performance across scenarios (`sims1.pdf`, `sims2.pdf`)
- **Raw data:** Complete simulation results saved as `data/sims.rds` (R data format)
- **Sample sizes:** Actual TND sample sizes achieved for each simulation

## Citation

If you use this code, please cite:

Boyer, C., Li, K.Q., Shi, X., & Tchetgen Tchetgen, T. J. (2025). "Identification and estimation of vaccine effectiveness in the test-negative design under equi-confounding". arXiv. https://doi.org/10.48550/arXiv.2504.20360

