# parallel-tnd

Replication files for:

Boyer, C., Li, K.Q., Shi, X., & Tchetgen Tchetgen, T. J. (2025). “Identification and estimation of vaccine effectiveness in the test-negative design under equi-confounding”. arXiv. https://doi.org/10.48550/arXiv.2504.20360

# Abstract
The test-negative design (TND) is frequently used to evaluate vaccine effectiveness in real-world settings. In a TND study, individuals with similar symptoms who seek care are tested for the disease of interest, and vaccine effectiveness is estimated by comparing the vaccination history of test-positive cases and test-negative controls. The design has previously been justified on the grounds that it reduces confounding due to unmeasured health-seeking behavior, although this has not been formally described using potential outcomes. However, it is also widely acknowledged that, by conditioning participation on receipt of a test, the TND risks inducing selection bias. In this paper, we propose a formal justification for the TND based on the assumption of \textit{odds ratio equi-confounding}, where unmeasured confounders influence test-positive and test-negative individuals equivalently on the odds ratio scale, with health-seeking behavior being just one plausible example. We also show that these results hold under outcome-dependent sampling design of the TND. We discuss the implications of the equi-confounding assumption for TND design and provide alternative estimators for the marginal risk ratio among the vaccinated under equi-confounding, including estimators based on outcome modeling and inverse probability weighting as well as a semiparametric estimator that is doubly-robust.  When the equi-confounding assumption does not hold, we suggest a straightforward sensitivity analysis that parameterizes the magnitude of the deviation on the odds ratio scale. We conduct a simulation study to evaluate the empirical performance of our proposed estimators under a wide range of scenarios. Finally, we also discuss how test-negative outcomes may be used more broadly to de-bias estimates from cohort studies where testing is symptom-triggered.

## Requirements

This code requires R version 4.0 or higher. The following R packages are required:

- `data.table` - for data manipulation
- `ggplot2` - for plotting
- `progressr` - for progress tracking during simulations
- `survival` - for survival analysis functions
- `lmtest` - for linear model testing
- `sandwich` - for robust standard errors
- `readr` - for reading/writing data files
- `dplyr` - for data manipulation
- `tidyr` - for data tidying
- `stringr` - for string manipulation
- `kableExtra` - for table formatting
- `patchwork` - for combining plots
- `tidyverse` - collection of tidy data packages

## Installation

1. Clone this repository:
```bash
git clone https://github.com/boyercb/parallel-tnd.git
cd parallel-tnd
```

2. Install required R packages:
```r
install.packages(c("data.table", "ggplot2", "progressr", "survival", 
                   "lmtest", "sandwich", "readr", "dplyr", "tidyr", 
                   "stringr", "kableExtra", "patchwork", "tidyverse"))
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
- `manuscript/` - LaTeX source files for the manuscript
- `results/` - Output files (tables and figures)

## How to run

### Running the full simulation study

To reproduce all simulation results from the paper, run the following command in R from the project root directory:

```r
source("code/run.R")
```

**Note:** The simulation study runs 1,000 replications across 8 different scenarios with a sample size of 15,000 per simulation. This will take several hours to complete (estimated 4-6 hours depending on your system).

### Generating tables and figures

After running the simulations, you can generate the tables and figures using:

```r
# Load simulation results
sims <- readr::read_rds("data/sims.rds")

# Generate tables
source("code/table.R")

# Generate figures  
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

1. **No unmeasured confounding** - Baseline scenario where traditional TND assumptions hold
2. **Equi-confounding** - Main scenario where all novel assumptions hold
3. **Exclusion restriction violated** - Vaccination affects test-negative outcomes
4. **Equi-confounding is violated** - Confounding differs between test-positive cases and test-negative controls
5. **Equi-selection is violated** - Testing behavior differs for test-positive cases and test-negative controls
6. **Equal effects of vaccination on testing** - Vaccination affects testing behavior equally for test-positive cases and test-negative controls 
7. **Unequal effects of vaccination on testing** - Vaccination affects testing behavior differently for test-positive cases and test-negative controls wrong
8. **Effect heterogeneity** - Same as scenario 2 except vaccine has heterogeneous effects

Each scenario compares multiple estimators for risk ratio among vaccinated:
- TND standard logistic regression estimator
- TND outcome modeling (OM) estimator 
- TND inverse probability weighting (IPW) estimator  
- TND doubly-robust (DR) estimator
- Cohort study regression estimators (with and without unmeasured confounders)
- Cohort study difference-in-differences estimator

## Output

The simulation generates:
- Bias, coverage probability, and confidence interval length for each estimator
- Results tables in LaTeX format (saved to `results/` directory)
- Figures showing estimator performance across scenarios
- Raw simulation results saved as `data/sims.rds`

## Citation

If you use this code, please cite:

Boyer, C., Li, K.Q., Shi, X., & Tchetgen Tchetgen, T. J. (2025). "Identification and estimation of vaccine effectiveness in the test-negative design under equi-confounding". arXiv. https://doi.org/10.48550/arXiv.2504.20360

