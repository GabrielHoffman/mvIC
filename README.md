
##### Evaluate BIC for multivariate regression to enable model selection with multiple response variables.

mvBIC extends the standard the standard Bayesian Information Criterion (BIC) to the case of multivariate regression.  The model fit across many response variables is evaluated and the criterion explicitly considers correlation between reponses.  mvBIC is appiciable to linear and linear mixed models.
```r
devtools::install_github("GabrielHoffman/mvBIC", repos=BiocManager::repositories())
```
This automatically installs dependencies from [Bioconductor](https://bioconductor.org)

See [manual](https://users.hpc.mssm.edu/~hoffmg01/software/mvBIC/mvBIC-manual.pdf) for examples and documentation.