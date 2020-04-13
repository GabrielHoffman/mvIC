# Gabriel Hoffman
# April 12, 2020

#' Multivariate forward stepwise regression 
#'
#' Multivariate forward stepwise regression evluated by multivariate BIC
#'
#' @param exprObj matrix of expression data (g genes x n samples), or ExpressionSet, or EList returned by voom() from the limma package
#' @param baseFormula specifies baseline variables for the linear (mixed) model.  Must only specify covariates, since the rows of exprObj are automatically used as a response. e.g.: \code{~ a + b + (1|c)} Formulas with only fixed effects also work, and \code{lmFit()} followed by \code{contrasts.fit()} are run.
#' @param data data.frame with columns corresponding to formula
#' @param variables array of variable names to be considered in the regression.  If variable should be considered as random effect, use '(1|A)'.
#' @param deltaCutoff stop interating of the model improvement is less than deltaCutoff 
#' @param countLevels default TRUE, count number of levels rather than number of variance compinents.  See description in \code{\link{nparam}}
#' @param verbose Default TRUE. Print messages
#' 
#' @return list with formula of final model, and trace of iterations during model selection
#' @examples
#' 
#' Y = with(iris, rbind(Sepal.Width, Sepal.Length))
#' 
#' # fit forward stepwise regression starting with model: ~1. Stop when model improvement is < 10
#' bestModel = mvForwardStepwise( Y, ~ 1, data=iris, variables=colnames(iris)[3:5], deltaCutoff=10)
#' 
#' bestModel
#' 
#' @import variancePartition
#' @importFrom stats as.formula
#' @export
mvForwardStepwise = function( exprObj, baseFormula, data, variables, deltaCutoff = 100, countLevels = TRUE, verbose=TRUE  ){

	# score base model
	suppressWarnings({
	baseScore = mvBIC_fit( exprObj, baseFormula, data, countLevels=countLevels, verbose=FALSE)
	})

	resultsList = list(	iteration 		= c(),
						variable 		= c(),
						delta 			= c(), 
						score 			= c(), 
						isBestVariable	= c(), 
						isAdded			= c())

	iteration = 1

	# run until break
	while(1){

		if( verbose ) message(paste0("Base model: ", paste0(baseFormula, collapse=' ')))

		# evaluate score of each potential model
		score = sapply( variables, function(feature){
			if( verbose ) message(paste("\tevaluating: +", feature))

			# formula of model to try
			form = as.formula(paste(paste0(baseFormula, collapse=' '), "+", feature))

			suppressWarnings({
			# evaluate multivariate BIC
			mvBIC_fit( exprObj, form, data, countLevels=countLevels, verbose=FALSE)
			})
			})

		# get index of minumum score
		i = which.min(score)

		# get difference between best and second best score
		delta = score[i] - baseScore
		if( verbose ) message("Best model delta: ", format(delta, digits=2))	

		# save results
		resultsList$iteration = c(resultsList$iteration, rep(iteration, length(score)))
		resultsList$variable = c(resultsList$variable, variables)
		resultsList$delta = c(resultsList$delta, score - baseScore)
		resultsList$score = c(resultsList$score, score)
		isBest = rep("", length(score))
		isBest[i] = "yes"
		resultsList$isBestVariable = c(resultsList$isBestVariable, isBest)
		isAdded = rep("", length(score))		

		iteration = iteration + 1

		# evaluate of best model is better than existing model
		if( delta < -deltaCutoff ){
			if( verbose ) message("Add variable to model: ", variables[i], '\n')

			# set new model, baseScore and possible variable list
			baseFormula = as.formula(paste(paste0(baseFormula, collapse=' '), "+", variables[i]))
			baseScore = score[i]
			variables = variables[-i]			
			isAdded[i] = "yes"
			resultsList$isAdded = c(resultsList$isAdded, isAdded)

			# if there are no additional variables to try
			if( length(variables) == 0){
				break
			}			

		}else{
			resultsList$isAdded = c(resultsList$isAdded, isAdded)

			if( verbose ){
				message(paste0("Final model: ", paste0(baseFormula, collapse=' ')))
				message("No additional variables selected")
			}
			break
		}
	}

	return( list(formula = baseFormula, 
				 trace = as.data.frame(resultsList) ) )
}








#' Evaluate multivariate BIC
#'
#' Evaluate multivariate BIC on matrix of response variables.  Smaller is better.
#'
#' @param exprObj matrix of expression data (g genes x n samples), or ExpressionSet, or EList returned by voom() from the limma package
#' @param formula specifies variables for the linear (mixed) model.  Must only specify covariates, since the rows of exprObj are automatically used as a response. e.g.: \code{~ a + b + (1|c)} Formulas with only fixed effects also work, and \code{lmFit()} followed by \code{contrasts.fit()} are run.
#' @param data data.frame with columns corresponding to formula
#' @param countLevels default TRUE, count number of levels rather than number of variance compinents.  See description in \code{\link{nparam}}
#' @param verbose Default TRUE. Print messages
#'
#' @description
#' Evaluate multivariate BIC while considering correlation between response variables. For n samples, p responses and m parameters for each model, evaluate the multivariate BIC as \deqn{n * logDet(\Sigma) + log(n) * (p*m + 0.5*p*(p+1))}
#' where \eqn{\Sigma} is the residual covariance matrix.  This formula extends the standard univariate BIC to the multivariate case.  
#' For one response the standard penalty is \eqn{log(n)*m}, this just adds a \eqn{log(n)} to that value, but only the differences between two models is important.  Estimating the \eqn{p x p} covariance matrix requires \eqn{0.5*p*(p+1)} parameters.
#' When \eqn{p > m} the residual covariance matrix Sigma is not full rank.  In this case the psudo-determinant is used instead.
#'
#' See References
#'
#' Pauler, DK. The Schwarz criterion and related methods for normal linear models. Biometrika (1998), 85, 1, pp. 13-27
#' 
#' Edward J. Bedrick and Chih-Ling Tsai. Model Selection for Multivariate Regression in Small Samples.  Biometrics,  50:1 1994 226-231
#'
#' TJ Wu, P Chen, Y Yan. The weighted average information criterion for multivariate regression model selection. Signal Processing 93.1 (2013): 49-55.
#'
#' @return multivariate BIC value
#' @examples
#' 
#' # create matrix of responses
#' Y = with(iris, rbind(Sepal.Width, Sepal.Length))
#' 
#' # Evaluate model 1
#' mvBIC_fit( Y, ~ Species, data=iris)
#' 
#' # Evaluate model 2
#' # smaller mvBIC means better model
#' mvBIC_fit( Y, ~ Petal.Width + Petal.Length + Species, data=iris)
#' 
#' @import variancePartition 
#'
#' @export
mvBIC_fit = function( exprObj, formula, data, countLevels = TRUE, verbose=FALSE ){

	# fit model and compute residuals
	suppressWarnings({
	modelFit = dream( exprObj, formula, data, REML=FALSE, computeResiduals=TRUE, quiet=!verbose)
	})

	# extract residuals
	residMatrix = residuals( modelFit, modelFit )

	# for a fixed effect model, cat fitVarPartModel() fails
	
	fit = tryCatch({
	    # get one full model fit
	    suppressWarnings({
		fitVarPartModel( exprObj[1,,drop=FALSE], formula, data, showWarnings=FALSE, quiet=!verbose)[[1]]
		})
	}, error = function(e) {
		return("isError")
	})

	if( is(fit, "character") ){
		# fixed effect
	   	m <- ncol(coef(modelFit)) + 1
	}else{
		# mixed model
		# extract number of parameters
		m <- nparam( fit, countLevels=countLevels)
	}
	
	mvBIC_from_residuals( residMatrix, m)
}	



#' residuals for MArrayLM
#'
#' residuals for MArrayLM
#'
#' @param object MArrayLM object from dream
#' @param ... other arguments, currently ignored
#'
#' @return results of residuals
#' @importFrom limma residuals.MArrayLM
#' @rdname residuals-method
#' @aliases residuals,MArrayLM-method
setMethod("residuals", "MArrayLM",
function( object, ...){
	if( is.null(object$residuals) ){
		# use residuals computed by limma
		res = residuals.MArrayLM( object, ...)
	}else{
		# use precomputed residuals
		res = object$residuals
	}
	res
})

#Matrix of residual from list of model fits
#
# Matrix of residual from list of model fits
#
# @param fitList list of model fits from \code{lm()} or \code{lmer()}
#
getResids = function(fitList){

	if(is(fitList, "list")){
		residMatrix = lapply(fitList, residuals)
		residMatrix = do.call(rbind, residMatrix)
	}else if(is(fitList, "mlm")){
		residMatrix = t(residuals( fitList))
	}else{
		stop("Fit must be list or mlm")
	}
	residMatrix
}

#' Number of parameters in model
#' 
#' Number of parameters in model from \code{lm()} or \code{lmer()}
#' 
#' @param object model fit by \code{lm()} or \code{lmer()}
#' @param countLevels default TRUE, count number of levels rather than number of variance compinents.  See description
#'
#' @description 
#' In the case of \code{lm()}, the result is the nubmer of coefficients + 1 for the variance term.  In the case of \code{lmer()}, there are two options.  For \code{countLevels = TRUE}, returns the number of fixed effects + number of levels in random effects + 1 for residual variance term.  This treats each level of a random effect as a parameter. For \code{countLevels = FALSE}, returns number of fixed effects + number of variance components.  Here a random effect with 10 levels is only counted as 1 parameter.  This tends to underpenalize. 
#' @return number of parameters
#' @importFrom stats coef
#' @importFrom methods is
nparam = function( object, countLevels = TRUE){

	# if object is not any of these
	if( !is(object, "lm") & !is(object, "mlm") & !is(object, 'merMod') ){
		# see if element of list is valid model fit
		object = object[[1]]
	}

	if( is(object, "mlm") ){
		# must be evaluated first because if object is 'mlm', it is also 'lm'
		# need 'mlm' to take presidence
		m = nrow(coef(object)) + 1
	}else if( is(object, "lm") ){
		m = length(coef(object)) + 1
	}else if( is(object, "merMod") ){

	    if( countLevels ){
	   	 	# fixed + number of random levels
	    	m = length(object@beta) + object@devcomp[["dims"]][['q']] + object@devcomp[["dims"]][["useSc"]]
	    }else{
	    	# lme4:::npar.merMod
			# counts each random effect as a single parameter
		    m = length(object@beta) + length(object@theta) + object@devcomp[["dims"]][["useSc"]]
	    }
	}else{
		stop("object is not a valid model fit from lm() or lmer()")
	}

    m
}

#' Evaluate multivariate BIC
#'
#' Evaluate multivariate BIC from a list of regression fits
#'
#' @param fitList list of model fits with \code{lm()} or \code{lmer()}.  All models must have same data, response and formula.
#' @param countLevels default TRUE, count number of levels rather than number of variance compinents.  See \code{\link{nparam}} for details
#'
#' @description
#' Evaluate multivariate BIC while considering correlation between response variables. For n samples, p responses and m parameters for each model, evaluate the multivariate BIC as \deqn{n * logDet(\Sigma) + log(n) * (p*m + 0.5*p*(p+1))}
#' where \eqn{\Sigma} is the residual covariance matrix.  This formula extends the standard univariate BIC to the multivariate case.  
#' For one response the standard penalty is \eqn{log(n)*m}, this just adds a \eqn{log(n)} to that value, but only the differences between two models is important.  Estimating the \eqn{p x p} covariance matrix requires \eqn{0.5*p*(p+1)} parameters.
#' When \eqn{p > m} the residual covariance matrix Sigma is not full rank.  In this case the psudo-determinant is used instead.
#'
#' See References
#'
#' Pauler, DK. The Schwarz criterion and related methods for normal linear models. Biometrika (1998), 85, 1, pp. 13-27
#' 
#' Edward J. Bedrick and Chih-Ling Tsai. Model Selection for Multivariate Regression in Small Samples.  Biometrics,  50:1 1994 226-231
#'
#' TJ Wu, P Chen, Y Yan. The weighted average information criterion for multivariate regression model selection. Signal Processing 93.1 (2013): 49-55.
#'
#' @return multivariate BIC value
#' @examples
#' # Predict Sepal width and Length given Species 
#' # Evaluate model fit
#' fit1 = lm( cbind(Sepal.Width, Sepal.Length) ~ Species, data=iris)
#' mvBIC( fit1 )
#' 
#' # add Petal width and length
#' # smaller mvBIC means better model
#' fit2 = lm( cbind(Sepal.Width, Sepal.Length) ~ Petal.Width + Petal.Length + Species, data=iris)
#' mvBIC( fit2 )
#' 
#' @importFrom methods is
#' @export
mvBIC = function( fitList, countLevels = TRUE ){

	# get residuals for 'mlm', 'lm', or list of 'lm' or 'lmer'
	residMatrix = getResids( fitList )

	# get number of parameters for multiple forms of fitList
	m = nparam( fitList, countLevels=countLevels )

	mvBIC_from_residuals( residMatrix, m )
}



##' Evaluate multivariate BIC from matrix of residuals
##'
##' Evaluate multivariate BIC from matrix of residuals
##' 
##' @param residMatrix matrix of residuals where rows are features
##' @param m number of parameters for each model
##' 
mvBIC_from_residuals = function( residMatrix, m){

	n = ncol(residMatrix) # number of samples
	p = nrow(residMatrix) # number of response variables

	# compute log determinant explicitly
	# slower and not defined for low rank matrices
	# dataTerm = n * determinant(crossprod(residMatrix), log=TRUE)$modulus[1]

	# Compute log det from singular values of residual matrix
	# Much faster for large data
	# When p > n, only considers the non-zero singular values	
	# this is the "pseudo-determinant" as corallary to the "pseudo-inverse"
	evalues = svd(residMatrix, nu=0, nv=0)$d^2
	evalues = evalues[evalues > 1e-10]
	# p = length(evalues)
	dataTerm = n * sum(log(evalues))

	# Penalty term for BIC
	# For one response the standard penalty is log(n)*m
	# 	this just adds a log(n) to that value
	#	but only the differences between two models is important
	# Estimating the p x p covariance matrix requires 0.5*p*(p+1) parameters
	penalty = log(n) * (p*m + 0.5*p*(p+1)) #- log(2*pi)*(p*m + 0.5*p*(p+1)) 

	# retrun data term plus penalty
	res = dataTerm + penalty
	attr(res, 'params') = data.frame(n=n, p=p, m=m)
	res
}






