# parallel-tnd

Replication files for:

Boyer, C., Li, K.Q., Shi, X., & Tchetgen Tchetgen, T. J. (2025). “Identification and estimation of vaccine effectiveness in the test-negative design under equi-confounding”. arXiv. https://doi.org/10.48550/arXiv.2504.20360

# Abstract
The test-negative design (TND) is frequently used to evaluate vaccine effectiveness in real-world settings. In a TND study, individuals with similar symptoms who seek care are tested for the disease of interest, and vaccine effectiveness is estimated by comparing the vaccination history of test-positive cases and test-negative controls. The design has previously been justified on the grounds that it reduces confounding due to unmeasured health-seeking behavior, although this has not been formally described using potential outcomes. However, it is also widely acknowledged that, by conditioning participation on receipt of a test, the TND risks inducing selection bias. In this paper, we propose a formal justification for the TND based on the assumption of \textit{odds ratio equi-confounding}, where unmeasured confounders influence test-positive and test-negative individuals equivalently on the odds ratio scale, with health-seeking behavior being just one plausible example. We also show that these results hold under outcome-dependent sampling design of the TND. We discuss the implications of the equi-confounding assumption for TND design and provide alternative estimators for the marginal risk ratio among the vaccinated under equi-confounding, including estimators based on outcome modeling and inverse probability weighting as well as a semiparametric estimator that is doubly-robust.  When the equi-confounding assumption does not hold, we suggest a straightforward sensitivity analysis that parameterizes the magnitude of the deviation on the odds ratio scale. We conduct a simulation study to evaluate the empirical performance of our proposed estimators under a wide range of scenarios. Finally, we also discuss how test-negative outcomes may be used more broadly to de-bias estimates from cohort studies where testing is symptom-triggered.

# How to run

