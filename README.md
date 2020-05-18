
<div align="right">
<img src="https://deepfigv.mssm.edu/img/software/mvIC/mvIC_logo.png" alt="mvIC logo" width="100px"><br>
</div>

![Image of Yaktocat](https://deepfigv.mssm.edu/img/software/mvIC/mvIC_logo.png)

# Evaluate information criteria for multivariate model selection

mvIC extends the standard the standard Akaike or Bayesian Information Criterion (AIC, BIC) to the case of multivariate regression.  The model fit across many response variables is evaluated and the criterion explicitly considers correlation between reponses.  mvIC is applicable to linear and linear mixed models.

Forward stepwise regression with the mvIC criterion enables automated variable selection for high dimensional datasets.


## Installation
```r
devtools::install_github("GabrielHoffman/mvIC", repos=BiocManager::repositories())
```
This automatically installs dependencies from [Bioconductor](https://bioconductor.org)

## Examples

#### Multivariate linear model
```r
# Predict Sepal width and Length given Species
# Evaluate model fit
fit1 = lm( cbind(Sepal.Width, Sepal.Length) ~ Species, data=iris)
score = mvIC( fit1 )
```

```r
# add Petal width and length
# smaller mvIC means better model
fit2 = lm( cbind(Sepal.Width, Sepal.Length) ~ Petal.Width + Petal.Length + Species, data=iris)
mvIC( fit2 )
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
See [manual](http://deepfigv.mssm.edu/img/software/mvIC/mvIC-manual.pdf) for examples and documentation.