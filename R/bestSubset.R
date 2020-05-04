# Gabriel Hoffman
# May 4, 2020

#' Evaluate information criteria for best subset selection
#'
#' Evaluate information criteria for best subset selection
#'
#' @param exprObj matrix of expression data (g genes x n samples), or ExpressionSet, or EList returned by voom() from the limma package
#' @param data data.frame with columns corresponding to formula
#' @param variables array of variable names to be considered in the regression.  If variable should be considered as random effect, use '(1|A)'.
#' @param maxk evaluate all combinations of variables up to size maxk
#'
#' @import data.table
#' @importFrom utils combn
#' @export
bestSubsetSearch = function(Y, data, variables, maxk=5){

  res = lapply( 1:maxk, function(k){
    idx = combn(length(variables), k)

    formArray = apply(idx, 2, function(i) paste('Y ~', paste(variables[i], collapse=' + ')))

    res = lapply(formArray, function(form){
      fit = lm(as.formula(form), data=data)

      res = lapply(c("EB", "none"), function(method){
        res = lapply(c("AIC", "BIC", "sum AIC", "sum BIC"), function(criterion){
          if( length(grep("sum", criterion)) > 0 & method != "none" ) return(NULL)

          score = mvIC(fit, shrink.method=method, criterion=criterion)
          data.frame(k=k,form=form, method=method, criterion=criterion, score = score[1], stringsAsFactors=FALSE)
        })
        do.call(rbind, res)
      })
      do.call(rbind, res)
    } )
    res = do.call(rbind, res)
  })
  data.table(do.call(rbind, res))
}