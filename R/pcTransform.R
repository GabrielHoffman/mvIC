


#' PC transformation
#'
#' Transform data using PCA and retaining given fraction of variance
#'
#' @param X data matrix
#' @param varFrac fraction of variance to retain
#' @export
pcTransform = function(X, varFrac = 1){

	if( ncol(X) == 1){
		return(X)
	}

	# SVD
	dcmp = svd(X)

	# if the number of columns is larger then nrow/5, 
	# then truncate spectrum
	if( ncol(X) > nrow(X)/5 ){
		# number of PC's that explain at least varFrac of the variance
		# with(dcmp, cumsum(d^2/sum(d^2)))
		k = with(dcmp, sum(cumsum(d^2/sum(d^2)) <= varFrac))
	}else{
		k = ncol(X)
	}

	d = u = NULL # PASS R check

	# created rotated data
	X_est = with(dcmp, u[,seq(1,k),drop=FALSE] %*% diag(d[seq(1,k)]))
	rownames(X_est) = rownames(X)
	X_est
}    