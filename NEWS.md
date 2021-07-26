# mvIC 1.6.1
July 26, 2021

- fix bug in calculating effective degrees of freedom: 
 - report m = mean edf per response
 - for fixed effect, now report m instead of m + 1
- Add more unit tests
- pcTransform() varFrac is 1 but default
-  when criterion = "BIC", shrink = "EB" and n > p, fastApprox is exact with only fixed effects, but approximate with random effects
 - with p > n use original value of p in eclairs to calculate logML

# mvIC 1.6.0  
## Major changes
July 13, 2021

- add `fastApprox` argument that is exact for fixed effects and approximate if random effects are included
	 - include `pcTransform()` function
	 - this can be 100x faster for large datasets
- better estimate of lambda and logML 
	- use `lambda = 0.01` by default
- defaults are new `criterion = "BIC"`, `shrink.method = "EB"`, `nparamsMethod = "edf"`
- remove scaling of responses to allow use of `pcTransform()` on responses
- Using effect degrees of freedom is faster since it is built into `dream()`
	- fixed bugs handling edf for each regression
	- remove redundant calculations
-  Add more unit tests

# mvIC 1.5.1  
October 26, 2020
- enforce compatability with variancePartition >= 1.19.20

# mvIC 1.5.0
  May 5, 2020
- add back gdf

# mvIC 1.4.0
  May 1, 2020
- Add bestSubsetSearch()
- Add estimateMVN_EB for Empirical Bayes estimation from Multivariate Normal Inverse Wishart model (MNIWEB)


# mvIC 1.3.0
  May 1, 2020
- add rough code for empirical Bayes covariance estimation 


# mvIC 1.1.4
  April 27th, 2020
- Add other selection criteria


# mvIC 1.1.4
  April 26th, 2020
- Add fast shrinkage of eigen values for Touloumis_equal


# v1.1.3
  April 26th, 2020
- add generalized degrees of freedom for covariance estimate

# v1.1.1
  April 16th, 2020
- Improve documentation
- add naive BIC method
- add print and show methods


# v1.1   
April 15th, 2020
- define number of parameters for BIC as nparamsMethod = c("edf", "countLevels", "lme4")
- set default to effective degrees of freedom: "edf" this reduce the penalty for random effects compared to counting the levels
- set stopping criterion to delta of 10