# Gabriel Hoffman
# April 12, 2020

# The values of BIC and evalCriterion are different
# but the *difference* in values for two models are the same
test_mvBIC = function(){

	# library(RUnit)

	n = 100
	p = 3

	X = matrix(rnorm(n*p),n,p)
	a = rnorm(n)
	y = X %*% rep(1,p) + rnorm(n)

	fit1 = lm( y ~ X)
	fit2 = lm( y ~ X + a)

	diff1 = BIC(fit2) - BIC(fit1)
	diff2 = evalCriterion( fit2) - evalCriterion( fit1)

	checkEqualsNumeric( diff1, diff2)
}



test_multiple_mbBIC = function(){

	# check for fixed effects models
	Y = with(iris, rbind(Sepal.Width, Sepal.Length))

	fit = lm( cbind(Sepal.Width, Sepal.Length) ~ Species, data=iris)

	scoreA = mvBIC( fit )
	scoreB = mvBIC_fit( Y, ~ Species, data=iris)

	checkEqualsNumeric( scoreA, scoreB)
}



test_multiple_mbBIC = function(){

	# check for mixed effects models
	Y = with(iris, rbind(Sepal.Width, Sepal.Length))

	fit1 = lme4::lmer(Sepal.Width ~ (1|Species), data=iris, REML=FALSE )
	fit2 = lme4::lmer(Sepal.Length ~ (1|Species), data=iris, REML=FALSE )

	scoreC = mvBIC( list(fit1, fit2) )
	scoreD = mvBIC_fit( Y, ~ (1|Species), data=iris)

	checkEqualsNumeric( scoreA, scoreB)
}


	



