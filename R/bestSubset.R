# Gabriel Hoffman
# May 4, 2020

#' Evaluate information criteria for best subset selection
#'
#' Evaluate information criteria for best subset selection
#'
#' @param Y response matrix with responses ad columns and samples as rows
#' @param data data.frame with columns corresponding to formula
#' @param variables array of variable names to be considered in the regression.  Random effects are not allowed
#' @param maxk evaluate all combinations of variables up to size maxk
#'
#' @importFrom dplyr as_tibble
#' @importFrom utils combn
#' @importFrom stats lm
#' @export
bestSubsetSearch = function(Y, data, variables, maxk=5){

  res = lapply( 1:maxk, function(k){
    idx = combn(length(variables), k)

    formArray = apply(idx, 2, function(i) paste('Y ~', paste(variables[i], collapse=' + ')))

    res = lapply(formArray, function(form){
      fit = lm(as.formula(form), data=data)

      n = nrow(residuals(fit))
      p = ncol(residuals(fit))

      res = lapply(c("EB", "none"), function(method){
        res = lapply(c("AIC", "BIC", "AICC", "CAIC", "sum AIC", "sum BIC"), function(criterion){
         
          # if( length(grep("sum", criterion)) > 0 & method != "none" ) return(NULL)
          if( criterion %in% c("AICC", "CAIC") & (method != "none") ) return(NULL)
          if( (method == "none") & !(criterion %in% c("sum AIC", "sum BIC")) & (p > n) ) return(NULL)     
          if( (method == "EB") & (criterion %in% c("sum AIC", "sum BIC")) ) return(NULL)     

          score = mvIC(fit, shrink.method=method, criterion=criterion)

          data.frame(k=k,form=form, method=method, criterion=criterion, score = score[1], stringsAsFactors=FALSE)
        })
        do.call(rbind, res)
      })
      do.call(rbind, res)
    } )
    res = do.call(rbind, res)
  })
  as_tibble(do.call(rbind, res))
}


