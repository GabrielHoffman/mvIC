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
#' @param criterion multivariate criterion ('AIC', 'BIC') or summing score assumign independence of reponses ('sum AIC', 'sum BIC')
#' @param shrink.method Shrink covariance estimates to be positive definite. Using "var_equal" assumes all variance on the diagonal are equal.  This method is the fastest because it is linear time.  Using "var_unequal" allows each response to have its own variance term, however this method is quadratic time.  Using "none" does not apply shrinkge, but is only valid when there are very few responses
#' @param nparamsMethod "edf": effective degrees of freedom. "countLevels" count number of levels in each random effect.  "lme4" number of variance compinents, as used by lme4.  See description in \code{\link{nparam}}
#' @param deltaCutoff stop interating of the model improvement is less than deltaCutoff. default is 5 
#' @param verbose Default TRUE. Print messages
#' @param ... additional arguements passed to logDet
#' 
#' @return list with formula of final model, and trace of iterations during model selection
#' @examples
#' 
#' Y = with(iris, rbind(Sepal.Width, Sepal.Length))
#' 
#' # fit forward stepwise regression starting with model: ~1. 
#' bestModel = mvForwardStepwise( Y, ~ 1, data=iris, variables=colnames(iris)[3:5])
#' 
#' bestModel
#' 
#' @import variancePartition
#' @importFrom stats as.formula
#' @importFrom dplyr as_tibble
#' @export
mvForwardStepwise = function( exprObj, baseFormula, data, variables, criterion = c("AIC", "BIC", "AICC", "CAIC", "sum AIC", "sum BIC"), shrink.method = c("var_equal", "var_unequal", "EB", "none"), nparamsMethod = c("edf", "countLevels", "lme4"), deltaCutoff = 5, verbose=TRUE,...  ){

	criterion = match.arg(criterion)
	shrink.method  = match.arg(shrink.method)
	nparamsMethod = match.arg(nparamsMethod)
	baseFormula = as.formula( baseFormula )

	if( ! is.data.frame(data) ){
		data = as.data.frame(data, stringsAsFactors=FALSE)
	}
	data = droplevels(data)

	# score base model
	suppressWarnings({
	baseScore = mvIC_fit( exprObj, baseFormula, data, criterion=criterion, shrink.method=shrink.method , nparamsMethod=nparamsMethod, verbose=FALSE,...)
	})

	resTrace = data.frame(	iter 		= 0,
							variable 	= as.character(baseFormula)[-1],
							delta 		= NA, 
							score 		= as.numeric(baseScore),								 
							isBest		= "yes", 
							isAdded		= "yes", 
							stringsAsFactors=FALSE)
	resTrace = cbind(resTrace, baseScore@params)

	iteration = 1

	# run until break
	while(1){

		if( verbose ) message(paste0("Base model: ", paste0(baseFormula, collapse=' ')))

		# evaluate score of each potential model
		score = lapply( variables, function(feature){
			if( verbose ) message(paste("\r\tevaluating: +", feature, '                  '), appendLF=FALSE)

			# formula of model to try
			form = as.formula(paste(paste0(baseFormula, collapse=' '), "+", feature))

			suppressWarnings({
			# evaluate multivariate score
			mvIC_fit( exprObj, form, data, criterion=criterion, shrink.method=shrink.method, nparamsMethod=nparamsMethod, verbose=FALSE, ...)
			})
			})

		# get index of minumum score
		i = which.min(unlist(score))

		# get difference between best and second best score
		delta = as.numeric(score[[i]] - baseScore)
		if( verbose ) message("\nBest model delta: ", format(delta, digits=2))	


		isBest = rep("", length(score))
		isAdded = rep("", length(score))	
		isBest[i] = "yes"
		if( delta < -deltaCutoff ){				
			isAdded[i] = "yes"
		}
		resNew = data.frame(iter 		= iteration,
							variable 	= variables,
							delta 		= as.numeric(unlist(score) - baseScore), 
							score 		= unlist(score),								 
							isBest		= isBest, 
							isAdded		= isAdded, 
							stringsAsFactors=FALSE)

		# get summary stats from model fits
		params = lapply(score, function(x) x@params)
		params = do.call(rbind, params)

		# combine results from this iteration
		resNew = cbind(resNew, params)

		# combine with results from previous interations
		resTrace = rbind(resTrace, resNew)

		iteration = iteration + 1

		# evaluate of best model is better than existing model
		if( delta < -deltaCutoff ){
			if( verbose ) message("Add variable to model: ", variables[i], '\n')

			# set new model, baseScore and possible variable list
			baseFormula = as.formula(paste(paste0(baseFormula, collapse=' '), "+", variables[i]))
			baseScore = score[[i]]

			variables = variables[-i]		

			# if there are no additional variables to try
			if( length(variables) == 0){
				break
			}			

		}else{
			if( verbose ){
				message(paste0("\nFinal model:\n  ", paste0(baseFormula, collapse=' ')))
			}
			break
		}
	}

	# remove some columns from resTrace that are constant
	idx = colnames(resTrace) %in% c("n", 'p', 'criterion', 'shrink.method')

	# return as mvIC_result object
	new("mvIC_result", 	list(formula = baseFormula, 
						settings= resTrace[1,idx],
						trace 	= as_tibble(resTrace[,!idx]) ))
}








#' Evaluate multivariate BIC
#'
#' Evaluate multivariate BIC on matrix of response variables.  Smaller is better.
#'
#' @param exprObj matrix of expression data (g genes x n samples), or ExpressionSet, or EList returned by voom() from the limma package
#' @param formula specifies variables for the linear (mixed) model.  Must only specify covariates, since the rows of exprObj are automatically used as a response. e.g.: \code{~ a + b + (1|c)} Formulas with only fixed effects also work, and \code{lmFit()} followed by \code{contrasts.fit()} are run.
#' @param data data.frame with columns corresponding to formula
#' @param criterion multivariate criterion ('AIC', 'BIC') or summing score assumign independence of reponses ('sum AIC', 'sum BIC')
#' @param shrink.method Shrink covariance estimates to be positive definite. Using "var_equal" assumes all variance on the diagonal are equal.  This method is the fastest because it is linear time.  Using "var_unequal" allows each response to have its own variance term, however this method is quadratic time.  Using "none" does not apply shrinkge, but is only valid when there are very few responses
#' @param nparamsMethod "edf": effective degrees of freedom. "countLevels" count number of levels in each random effect.  "lme4" number of variance compinents, as used by lme4.  See description in \code{\link{nparam}}
#' @param verbose Default TRUE. Print messages
#' @param ... additional arguements passed to logDet
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
#' #' Edward J. Bedrick and Chih-Ling Tsai. Model Selection for Multivariate Regression in Small Samples.  Biometrics,  50:1 1994 226-231
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
#' mvIC_fit( Y, ~ Species, data=iris)
#' 
#' # Evaluate model 2
#' # smaller mvIC means better model
#' mvIC_fit( Y, ~ Petal.Width + Petal.Length + Species, data=iris)
#' 
#' @import variancePartition 
#'
#' @export
mvIC_fit = function( exprObj, formula, data, criterion = c("AIC", "BIC", "AICC", "CAIC", "sum AIC", "sum BIC"), shrink.method = c("var_equal", "var_unequal", "EB", "none"), nparamsMethod = c("edf", "countLevels", "lme4"), verbose=FALSE, ... ){

	criterion = match.arg(criterion)
	shrink.method  = match.arg(shrink.method)
	nparamsMethod = match.arg(nparamsMethod)
	formula = as.formula(formula)

	if( ! is.data.frame(data) ){
		data = as.data.frame(data, stringsAsFactors=FALSE)
	}

	# scale features (i.e. rows)
	exprObj = scale_features( exprObj )

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
	
	mvIC_from_residuals( residMatrix, m, criterion=criterion, shrink.method=shrink.method,... )
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
#' @param criterion multivariate criterion ('AIC', 'BIC') or summing score assumign independence of reponses ('sum AIC', 'sum BIC')
#' @param shrink.method Shrink covariance estimates to be positive definite. Using "var_equal" assumes all variance on the diagonal are equal.  This method is the fastest because it is linear time.  Using "var_unequal" allows each response to have its own variance term, however this method is quadratic time.  Using "none" does not apply shrinkge, but is only valid when there are very few responses
#' @param nparamsMethod "edf": effective degrees of freedom. "countLevels" count number of levels in each random effect.  "lme4" number of variance compinents, as used by lme4.  See description in \code{\link{nparam}}
#' @param ... additional arguements passed to logDet
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
#' mvIC( fit1 )
#' 
#' # add Petal width and length
#' # smaller mvIC means better model
#' fit2 = lm( cbind(Sepal.Width, Sepal.Length) ~ Petal.Width + Petal.Length + Species, data=iris)
#' mvIC( fit2 )
#' 
#' @importFrom methods is
#' @export
mvIC = function( fitList, criterion =c("AIC", "BIC", "AICC", "CAIC", "sum AIC", "sum BIC"), shrink.method = c("var_equal", "var_unequal", "EB", "none"), nparamsMethod = c("edf", "countLevels", "lme4"), ...){

	criterion = match.arg(criterion)
	shrink.method  = match.arg(shrink.method)
	nparamsMethod = match.arg(nparamsMethod)

	# get residuals for 'mlm', 'lm', or list of 'lm' or 'lmer'
	residMatrix = getResids( fitList )

	# get number of parameters for multiple forms of fitList
	m = nparam( fitList, nparamsMethod=nparamsMethod )

	mvIC_from_residuals( residMatrix, m, criterion=criterion, shrink.method=shrink.method,...)
}



#' Evaluate multivariate BIC from matrix of residuals
#'
#' Evaluate multivariate BIC from matrix of residuals
#' 
#' @param residMatrix matrix of residuals where rows are features
#' @param m number of parameters for each model
#' @param criterion multivariate criterion ('AIC', 'BIC') or summing score assumign independence of reponses ('sum AIC', 'sum BIC')
#' @param shrink.method Shrink covariance estimates to be positive definite. Using "var_equal" assumes all variance on the diagonal are equal.  This method is the fastest because it is linear time.  Using "var_unequal" allows each response to have its own variance term, however this method is quadratic time.  Using "none" does not apply shrinkge, but is only valid when there are very few responses
#' @param ... other arguments passed to logDet
#' 
#' @importFrom methods new
mvIC_from_residuals = function( residMatrix, m, criterion =c("AIC", "BIC", "AICC", "CAIC", "sum AIC", "sum BIC"), shrink.method = c("var_equal", "var_unequal", "EB", "none"), ... ){

	criterion = match.arg(criterion)
	shrink.method  = match.arg(shrink.method)

	n = ncol(residMatrix) # number of samples
	p = nrow(residMatrix) # number of response variables

	if( criterion %in% c("AIC", "BIC", "AICC", "CAIC") ){

		if( criterion %in% c("AICC", "CAIC") & n < p){
			stop(paste("Criterion", criterion, "cannot be evaluated when n < p"))
		}
		
		# compute log determinant explicitly
		# slower and not defined for low rank matrices
		# dataTerm = n * determinant(crossprod(residMatrix), log=TRUE)$modulus[1]

		if( shrink.method == "EB"){

			# responses are *rows*
			# res = eb_cov_est( t(residMatrix) )
			# res = eb_cov_est2( t(residMatrix) )
			res = estimateMVN_EB( t(residMatrix) )
			lambda = res$alpha

            dataTerm = -2*res$logLik            
			gdf_cov = p + 1#(1-lambda)*p*(p-1)/2

		}else{
			# Evaluate logDet based on shrink.method
			logDet = rlogDet( residMatrix, shrink.method,... )
			dataTerm = n * logDet
			
			# get effective number of parameter used to estimate covariance by shrinkage
			gdf_cov = attr(logDet, "param")$gdf
			lambda = attr(logDet, "param")$lambda
		}

		# see Yanagihara, et al. 2015
		# doi:10.1214/15-EJS1022
		penalty = switch( criterion, 
						"AIC" 	= 2 * (p*(m-1) + gdf_cov),
						"BIC" 	= log(n) * (p*(m-1) + gdf_cov),
						"AICC"	= 2 * n*(p*(m-1) + gdf_cov) / (n-(m-1) - p - 1),
						"CAIC"	= (1+log(n)) * (p*(m-1) + gdf_cov))

		# retrun data term plus penalty
		res = dataTerm + penalty
		attr(res, 'params') = data.frame(	n 			= n, 
											p 			= p,
											m 			= as.numeric(m), 
											dataTerm 	= dataTerm, 
										 	penalty 	= penalty, 
										 	lambda 		= lambda,
										  	df_cov 		= gdf_cov, 
										  	criterion 	= criterion, 
										  	shrink.method = shrink.method,
										  	stringsAsFactors=FALSE)

		if( ! is.null( attr(m, 'nparamsMethod') )){
			attr(res, 'nparamsMethod') = attr(m, 'nparamsMethod')
		}else{
			attr(res, 'nparamsMethod') = "lm"
		}
	}else{
		# Naive metric summing BIC from all models independently
		rss = apply(residMatrix, 1, function(x) sum(x^2))
		dataTerm = n*sum(log(rss/n))

		penalty = switch(criterion, 
						"sum AIC" = 2 * m*p,
						"sum BIC" = log(n) * m*p)
		

		# retrun data term plus penalty
		res = dataTerm + penalty
		attr(res, 'params') = data.frame(	n 			= n, 
											p 			= p, 
											m 			= as.numeric(m),
											dataTerm 	= dataTerm, 
											penalty 	= penalty, 
											lambda 		= NA, 
											df_cov 		= NA,
											criterion 	= criterion, 
											shrink.method = "none",
											stringsAsFactors=FALSE)
		attr(res, 'nparamsMethod') = "naive"
	}

	new("mvIC", as.numeric(res), 
				 nparamsMethod = attr(res, 'nparamsMethod'),
				 params = attr(res, 'params'))
}


#' Class mvIC
#'
#' Class stores mvIC score, method and parameter values
#'
#' @name mvIC-class
#' @rdname mvIC-class
#' @exportClass mvIC
setClass("mvIC", representation(nparamsMethod = "character", params="data.frame"), contains="numeric")



# Print mvIC object
#
# Print mvIC object
# 
# @param x mvIC object
# @export
setMethod("print", "mvIC", function( x ){

	cat("\t\tMultivariate IC score\n\n")
	cat(paste("  Samples:\t", x@params$n, "\n"))	
	cat(paste("  Responses:\t", x@params$p, "\n"))	
	cat(paste("  Coef param:\t", round(x@params$m, digits=1), "\n"))
	cat(paste("  Cov param:\t", round(x@params$df_cov, digits=1), "\n"))	
	cat(paste("  Regression:\t", x@nparamsMethod), "\n")		
	cat("  Shrink method:", x@params$shrink.method, "\n")	
	cat(paste("  lambda:\t", format(x@params$lambda, digits=3), "\n"))	
	cat("  Criterion:\t", x@params$criterion, "\n")
	cat("  Score:\t", as.numeric(x), "\n\n")
})



# Show mvIC object
#
# Show mvIC object
# 
# @param object mvIC object
# @export
setMethod("show", "mvIC", function( object ){
	print( object )
})






#' Class mvIC_result
#'
#' Class stores result of \link{\code{mvForwardStepwise}}
#'
#' @name mvIC_result-class
#' @rdname mvIC_result-class
#' @exportClass mvIC_result
# setClass("mvIC_result", representation(formula = "formula", settings="data.frame", trace="data.frame"))
setClass("mvIC_result", contains="list")



# Print mvIC_result object
#
# Print mvIC_result object
# 
# @param x mvIC_result object
# @export
setMethod("print", "mvIC_result", function( x ){

	cat("\t\tMultivariate IC forward stepwise regression\n\n")
	cat("  Samples:\t", x$settings$n, '\n')
	cat("  Responses:\t", x$settings$p, '\n')
	cat("  Shrink method:", x$settings$shrink.method, '\n')
	cat("  Criterion:\t", x$settings$criterion, '\n')
	cat("  Iterations:\t", max(x$trace$iter), "\n\n")	
	cat('  Best model:', paste(as.character(x$formula), collapse=" "), '\n\n')
})



# Show mvIC_result object
#
# Show mvIC_result object
# 
# @param object mvIC_result object
# @export
setMethod("show", "mvIC_result", function( object ){
	print( object )
})







