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
#' @param nparamsMethod "edf": effective degrees of freedom. "countLevels" count number of levels in each random effect.  "lme4" number of variance compinents, as used by lme4.  See description in \code{\link{nparam}}
#' @param verbose Default TRUE. Print messages
#' @param ... additional arguements
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
mvForwardStepwise = function( exprObj, baseFormula, data, variables, deltaCutoff = 10, nparamsMethod = c("edf", "countLevels", "lme4"), verbose=TRUE,...  ){

	nparamsMethod = match.arg(nparamsMethod)

	if( ! is.data.frame(data) ){
		data = as.data.frame(data, stringsAsFactors=FALSE)
	}

	# scale features (i.e. rows)
	exprObj = scale_features( exprObj )

	# score base model
	suppressWarnings({
	baseScore = mvBIC_fit( exprObj, baseFormula, data, nparamsMethod=nparamsMethod, verbose=FALSE,...)
	})

	resultsList = list(	iter 		= c(),
						variable 	= c(),
						delta 		= c(), 
						score 		= c(),
						nparams 	= c(), 
						isBestVar	= c(), 
						isAdded		= c())

	iteration = 1

	# run until break
	while(1){

		if( verbose ) message(paste0("Base model: ", paste0(baseFormula, collapse=' ')))

		# evaluate score of each potential model
		score = lapply( variables, function(feature){
			if( verbose ) message(paste("\r\tevaluating: +", feature), appendLF=FALSE)

			# formula of model to try
			form = as.formula(paste(paste0(baseFormula, collapse=' '), "+", feature))

			suppressWarnings({
			# evaluate multivariate BIC
			mvBIC_fit( exprObj, form, data, nparamsMethod=nparamsMethod, verbose=FALSE, ...)
			})
			})

		# get index of minumum score
		i = which.min(unlist(score))

		# get difference between best and second best score
		delta = as.numeric(score[[i]] - baseScore)
		if( verbose ) message("\nBest model delta: ", format(delta, digits=2))	

		# save results
		resultsList$iter = c(resultsList$iter, rep(iteration, length(score)))
		resultsList$variable = c(resultsList$variable, variables)
		resultsList$delta = c(resultsList$delta, as.numeric(unlist(score) - baseScore) )
		resultsList$score = c(resultsList$score, unlist(score))
		resultsList$nparams = c(resultsList$nparams, unlist(lapply(score, function(x) round(attr(x, 'param')$m, 2))) )
		isBest = rep("", length(score))
		isBest[i] = "yes"
		resultsList$isBestVar = c(resultsList$isBestVar, isBest)
		isAdded = rep("", length(score))		

		iteration = iteration + 1

		# evaluate of best model is better than existing model
		if( delta < -deltaCutoff ){
			if( verbose ) message("Add variable to model: ", variables[i], '\n')

			# set new model, baseScore and possible variable list
			baseFormula = as.formula(paste(paste0(baseFormula, collapse=' '), "+", variables[i]))
			baseScore = score[[i]]
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
#' @param nparamsMethod "edf": effective degrees of freedom. "countLevels" count number of levels in each random effect.  "lme4" number of variance compinents, as used by lme4.  See description in \code{\link{nparam}}
#' @param verbose Default TRUE. Print messages
#' @param ... additional arguements
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
mvBIC_fit = function( exprObj, formula, data, nparamsMethod = c("edf", "countLevels", "lme4"), verbose=FALSE, ... ){

	nparamsMethod = match.arg(nparamsMethod)

	if( ! is.data.frame(data) ){
		data = as.data.frame(data, stringsAsFactors=FALSE)
	}

	# fit model and compute residuals
	suppressWarnings({
	modelFit = dream( exprObj, formula, data, REML=FALSE, computeResiduals=TRUE, quiet=!verbose)
	})

	# extract residuals
	residMatrix = residuals( modelFit, modelFit )

	# effect fixed 
	if( modelFit$method == "ls"){		
	   	m <- ncol(coef(modelFit)) + 1
	}else{

		if( nparamsMethod == "edf" ){
		    k = min(100, nrow(exprObj))
		}else{
			k = 1
		}

		fitList = fitVarPartModel( exprObj[1:k,,drop=FALSE], formula, data, showWarnings=FALSE, quiet=!verbose)

		# get number of parameters
		m <- nparam( fitList, nparamsMethod=nparamsMethod)
	}
	
	mvBIC_from_residuals( residMatrix, m, ... )
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
		residMatrix = t(t(residuals(fitList)))
	}
	residMatrix
}

#' Number of parameters in model
#' 
#' Number of parameters in model from \code{lm()} or \code{lmer()}
#' 
#' @param object model fit by \code{lm()} or \code{lmer()}
#' @param nparamsMethod "edf": effective degrees of freedom. "countLevels" count number of levels in each random effect.  "lme4" number of variance compinents, as used by lme4.  See description in \code{\link{nparam}}
#'
#' @description 
#' In the case of \code{lm()}, the result is the number of coefficients, + 1 for the variance term.  For a linear mixed model fit with \code{lmer()} there are 3 options.  "edf": effective degrees of freedom as computed by sum of diagonal values of the hat matrix return by \code{lmer()} . "countLevels", returns the number of fixed effects + number of levels in random effects + 1 for residual variance term.  This treats each level of a random effect as a parameter. "lme4", returns number of fixed effects + number of variance components.  Here a random effect with 10 levels is only counted as 1 parameter.  This tends to underpenalize. 
#' @return number of parameters
#' @importFrom stats coef
#' @importFrom methods is
#' @importFrom stats hatvalues
nparam = function( object, nparamsMethod = c("edf", "countLevels", "lme4")){

	nparamsMethod = match.arg(nparamsMethod)

	if( is(object, "list") & nparamsMethod == "edf" ){

		if( all(sapply(object, function(fit) is(fit, "merMod"))) ){
			# if object is a list of merMod fits, compute the mean effective degrees of freedom
			# this is more stable than code below using edf of just one model
			m = mean(sapply( object, function(fit) sum(hatvalues(fit)) )) + 1
			nparamsMethod = "already estimated"
		}
	}

	# if object is not any of these
	if( !is(object, "lm") & !is(object, "mlm") & !is(object, 'merMod') ){
		# see if element of list is valid model fit
		object = object[[1]]
	}

	if( is(object, "mlm") ){
		# must be evaluated first because if object is 'mlm', it is also 'lm'
		# need 'mlm' to take presidence
		m = nrow(coef(object)) + 1
		attr(m, "nparamsMethod") = "lm"
	}else if( is(object, "lm") ){
		m = length(coef(object)) + 1
		attr(m, "nparamsMethod") = "lm"
	}else if( is(object, "merMod") ){

		m = switch( nparamsMethod, 
			# effective degrees of freedom
			# Add term for residual variance
			"edf" = sum(hatvalues(object)) + 1,

	   	 	# fixed + number of random levels
			"countLevels" = length(object@beta) + object@devcomp[["dims"]][['q']] + object@devcomp[["dims"]][["useSc"]],

			# lme4:::npar.merMod
			# counts each random effect as a single parameter
			"lme4" = length(object@beta) + length(object@theta) + object@devcomp[["dims"]][["useSc"]],

			"already estimated" = m
			)

		attr(m, "nparamsMethod") = nparamsMethod

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
#' @param nparamsMethod "edf": effective degrees of freedom. "countLevels" count number of levels in each random effect.  "lme4" number of variance compinents, as used by lme4.  See description in \code{\link{nparam}}
#' @param ... additional arguements
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
mvBIC = function( fitList, nparamsMethod = c("edf", "countLevels", "lme4"), ...){

	nparamsMethod = match.arg(nparamsMethod)

	# get residuals for 'mlm', 'lm', or list of 'lm' or 'lmer'
	residMatrix = getResids( fitList )

	# get number of parameters for multiple forms of fitList
	m = nparam( fitList, nparamsMethod=nparamsMethod )

	mvBIC_from_residuals( residMatrix, m, ...)
}



#' Evaluate multivariate BIC from matrix of residuals
#'
#' Evaluate multivariate BIC from matrix of residuals
#' 
#' @param residMatrix matrix of residuals where rows are features
#' @param m number of parameters for each model
#' @param useMVBIC (development only, do not change)
#' @param logDetMethod method to compute logDet of correlation matrix for finite sample size
#' 
#' @importFrom methods new
mvBIC_from_residuals = function( residMatrix, m, useMVBIC = TRUE, logDetMethod = c("Touloumis_equal", "Touloumis_unequal", "Strimmer", "pseudodet")  ){

	logDetMethod = match.arg( logDetMethod )

	n = ncol(residMatrix) # number of samples
	p = nrow(residMatrix) # number of response variables

	if( useMVBIC ){
		# compute log determinant explicitly
		# slower and not defined for low rank matrices
		# dataTerm = n * determinant(crossprod(residMatrix), log=TRUE)$modulus[1]

		# Evaluate logDet based on logDetMethod
		logDet = rlogDet( residMatrix, logDetMethod )
		dataTerm = n * logDet
		
		# Penalty term for BIC
		# For one response the standard penalty is log(n)*m
		# 	this just adds a log(n) to that value
		#	but only the differences between two models is important
		# Estimating the p x p covariance matrix requires 0.5*p*(p+1) parameters
		penalty = log(n) * (p*m + 0.5*p*(p+1)) #- log(2*pi)*(p*m + 0.5*p*(p+1)) 

		# retrun data term plus penalty
		res = dataTerm + penalty
		attr(res, 'params') = data.frame(n=n, p=p, m=as.numeric(m), dataTerm=dataTerm, penalty=penalty, lambda = attr(logDet, "lambda"))

		if( ! is.null( attr(m, 'nparamsMethod') )){
			attr(res, 'nparamsMethod') = attr(m, 'nparamsMethod')
		}else{
			attr(res, 'nparamsMethod') = "lm"
		}
	}else{
		# Naive metric summing BIC from all models independently
		rss = apply(residMatrix, 1, function(x) sum(x^2))
		dataTerm = n*sum(log(rss/n))
		penalty = log(n) * m*p

		# retrun data term plus penalty
		res = dataTerm + penalty
		attr(res, 'params') = data.frame(n=n, p=p, m=as.numeric(m), dataTerm=dataTerm, penalty=penalty, lambda = NA)
		attr(res, 'nparamsMethod') = "naive"
	}

	new("mvBIC", as.numeric(res), 
				 nparamsMethod = attr(res, 'nparamsMethod'),
				 params = attr(res, 'params'))
}


#' Class mvBIC
#'
#' Class stores mvBIC score, method and parameter values
#'
#' @name mvBIC-class
#' @rdname mvBIC-class
#' @exportClass mvBIC
setClass("mvBIC", representation(nparamsMethod = "character", params="data.frame"), contains="numeric")



# Print mvBIC object
#
# Print mvBIC object
# 
# @param x mvBIC object
# @export
setMethod("print", "mvBIC", function( x ){

	cat("\t\tMultivariate BIC score\n\n")
	cat(paste("  Method:\t", x@nparamsMethod), "\n")
	cat(paste("  Samples:\t", x@params$n, "\n"))	
	cat(paste("  Responses:\t", x@params$p, "\n"))	
	cat(paste("  Parameters:\t", x@params$m, "\n"))	
	cat(paste("  lambda:\t", format(x@params$lambda, digits=3), "\n"))	
	cat("  Score:\t", as.numeric(x), digits=3, "\n\n")
})



# Show mvBIC object
#
# Show mvBIC object
# 
# @param object mvBIC object
# @export
setMethod("show", "mvBIC", function( object ){
	print( object )
})








