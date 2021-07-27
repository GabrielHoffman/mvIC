# Gabriel Hoffman
# April 12, 2020

# The values of BIC and evalCriterion are different
# but the *difference* in values for two models are the same
test_mvBIC = function(){

	# library(RUnit)

	fit1 <- lm(dist ~ 1, cars)
	fit2 = lm(dist ~ speed, cars)

	# standard R functions
	a = BIC(fit2) - BIC(fit1)

	# mvIC: compare two models
	b = as.numeric(mvIC(fit2)) - as.numeric(mvIC(fit1))

	checkEqualsNumeric( a, b)
}


test_compare = function(){

	set.seed(1)
	df_iris = data.frame(iris, test = rnorm(nrow(iris)))

	fit1 = lm( cbind(Sepal.Width, Sepal.Length) ~ Species, data=df_iris)
	score1 = mvIC( fit1 )

	fit2 = lm( cbind(Sepal.Width, Sepal.Length) ~ Species + test, data=df_iris)
	score2 = mvIC( fit2 )

	d = as.numeric(score1) - as.numeric(score2)

	Y = with(iris, t(cbind(Sepal.Width, Sepal.Length)))
	res = mvForwardStepwise( Y, ~ Species, df_iris, variables= "test")

	checkEqualsNumeric( res$trace$delta[2], -d)
}


test_multiple_mbBIC_fixed = function(){

	# check for fixed effects models
	Y = with(iris, rbind(Sepal.Width, Sepal.Length))

	fit = lm( cbind(Sepal.Width, Sepal.Length) ~ Species, data=iris)

	scoreA = mvIC( fit )
	scoreB = mvIC_fit( Y, ~ Species, data=iris)

	checkEqualsNumeric( scoreA, scoreB)
}



test_multiple_mbBIC_random = function(){

	# y = rnorm(nrow(iris))
	# fit = lme4::lmer(y ~ (1|Species), data=iris )

	# check for mixed effects models
	Y = with(iris, rbind(Sepal.Width, Sepal.Length))

	fit1 = lme4::lmer(Sepal.Width ~ (1|Species), data=iris, REML=FALSE )
	fit2 = lme4::lmer(Sepal.Length ~ (1|Species), data=iris, REML=FALSE )

	scoreC = mvIC( list(fit1, fit2) )
	scoreD = mvIC_fit( Y, ~ (1|Species), data=iris, fastApprox=FALSE)

	# scoreC@params$dataTerm
	# scoreD@params$dataTerm

	# mvIC:::nparam(list(fit1, fit2))

	checkEqualsNumeric( scoreC, scoreD)
}

# July 13, 2021
test_new_settings = function(){ 

	Y = with(iris, cbind(Sepal.Width, Sepal.Length))
	rownames(Y) = rownames(iris)

	fit1 = lm( Y ~ Species, data=iris)
	a = mvIC( fit1, criterion="BIC", shrink.method="EB") 
	     
	fit = lm( pcTransform(Y) ~ Species, data=iris)
	b = mvIC( fit, criterion="BIC", shrink.method="EB")	

	fit = lm( pcTransform(Y) ~ Species, data=iris)
	c = mvIC( fit)	

	fit = lm( Y ~ Species, data=iris)
	d = mvIC( fit, fastApprox=TRUE)	

	checkEqualsNumeric( a, b) & 
	checkEqualsNumeric( a, c) & 
	checkEqualsNumeric( a, d)
}


test_new_settings_mvIC_fit = function(){ 

	Y = with(iris, cbind(Sepal.Width, Sepal.Length))
	rownames(Y) = rownames(iris)

	fit1 = lm( Y ~ Species, data=iris)
	a = mvIC( fit1) 

	b = mvIC_fit( t(Y), ~ Species, data=iris)

	c = mvIC_fit( t(Y), ~ Species, data=iris, fastApprox=TRUE)
	   
	checkEqualsNumeric( a, b) & 
	checkEqualsNumeric( a, c)
}


test_new_settings_mvIC_fit_random = function(){ 

	Y = with(iris, cbind(Sepal.Width, Sepal.Length))
	rownames(Y) = rownames(iris)

	# test with fixed effects
	variables = "Species"
	bm = mvForwardStepwise( t(Y), ~ 1, data=iris, variables=variables, verbose=FALSE)

	a = mvIC_fit( t(Y), ~ 1, data=iris)
	b = mvIC_fit( t(Y), ~ Species, data=iris)

	a_score = with(a@params, dataTerm + penalty)
	b_score = with(b@params, dataTerm + penalty)

	checkEqualsNumeric(bm$trace$score, c(a_score, b_score))

	# test with random effects
	variables = "(1|Species)"
	bm = mvForwardStepwise( t(Y), ~ 1, data=iris, variables=variables, verbose=FALSE)

	a = mvIC_fit( t(Y), ~ 1, data=iris)
	b = mvIC_fit( t(Y), ~ (1|Species), data=iris)

	a_score = with(a@params, dataTerm + penalty)
	b_score = with(b@params, dataTerm + penalty)

	checkEqualsNumeric(bm$trace$score, c(a_score, b_score))

	# Check formulas learned from exact and approximate methods
	# --------
	variables = c("Petal.Length", "Petal.Width", "(1|Species)")
	# exact
	bm2 = mvForwardStepwise( t(Y), ~ 1, data=iris, variables=variables, verbose=FALSE)

	# fastApprox is approximate for random effects
	bm3 = mvForwardStepwise( t(Y), ~ 1, data=iris, variables=variables, fastApprox=TRUE, verbose=FALSE)

	checkEquals(bm2$formula, bm3$formula)
}










	



