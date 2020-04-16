
<div align="right">
<img src="https://users.hpc.mssm.edu/~hoffmg01/software/mvBIC/mvBIC_logo.png" alt="mvBIC logo" width="100px"><br>
</div>

# Evaluate BIC for multivariate model selection

mvBIC extends the standard the standard Bayesian Information Criterion (BIC) to the case of multivariate regression.  The model fit across many response variables is evaluated and the criterion explicitly considers correlation between reponses.  mvBIC is appiciable to linear and linear mixed models.

Forward stepwise regression with the mvBIC criterion enables automated variable selection for high dimensional datasets.


## Installation
```r
devtools::install_github("GabrielHoffman/mvBIC", repos=BiocManager::repositories())
```
This automatically installs dependencies from [Bioconductor](https://bioconductor.org)

## Examples

#### Multivariate linear model
```r
# Predict Sepal width and Length given Species
# Evaluate model fit
fit1 = lm( cbind(Sepal.Width, Sepal.Length) ~ Species, data=iris)
score = mvBIC( fit1 )
```
The variable `score` now contains the `mvBIC` value and method used. It also contains values for `n` samples, `p` response variables and `m` parameters.
```r
print(score)
[1] 979.6136
attr(,"nparamsMethod")
[1] "lm"
attr(,"params")
    n p m
1 150 2 4
```

```r
# add Petal width and length
# smaller mvBIC means better model
fit2 = lm( cbind(Sepal.Width, Sepal.Length) ~ Petal.Width + Petal.Length + Species, data=iris)
mvBIC( fit2 )
```

#### Forward stepwise regression
```r
# Combine respones on *rows*
Y = with(iris, rbind(Sepal.Width, Sepal.Length))

# variables to consider in the model
variables = c("Petal.Length", "Petal.Width", "Species")

# fit forward stepwise regression starting with model: ~1. 
bestModel = mvForwardStepwise( Y, ~ 1, data=iris, variables=variables)
```

Categorical variables can be modeled as random effects using the `(1|x)` syntax.
```r
# model Species as a random effect
variables = c("Petal.Length", "Petal.Width", "(1|Species)")

bestModel = mvForwardStepwise( Y, ~ 1, data=iris, variables=variables)
```
If a random effect is specified, a linear mixed model is fit and the number of parameter equal to the effective degrees of freedom of the model fit.  Note that using a linear mixed model is more computationally demanding, but 
- prevents overcorrection with variables with many categories
- regularizes the effect of estimate for each category
- gracefully handles colinearity between categorical variables


## Documentation
See [manual](https://users.hpc.mssm.edu/~hoffmg01/software/mvBIC/mvBIC-manual.pdf) for examples and documentation.