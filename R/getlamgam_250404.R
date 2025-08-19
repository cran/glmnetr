###################################################################################################
###################################################################################################

#' Get numerical values for lam and gam 
#' 
#' @description 
#' This functin derives the numerical values for lam and gam (lambda and gamma) for usage in plot and 
#' predict() functions.  If the input variables lam and gam are unspecified then 
#' the cross validation informed lambda and gamma values ('lambda.min' and 'gamma.mim) 
#' which minimize the cross validation deviance, are returned.  One may also give
#' lam='lambda.1se' and gam='gamma.1se' to identify the corresponding numerical 
#' values.  Importantly one may also simply specify gam=1 (and lam=NULL) to get the best lasso fit 
#' from the unrelaxed lasso model (gamma=1), or gam=0 (and lam=NULL) to get the best lasso fit
#' from the fully relaxed lasso model (gamma=0).  
#' 
#' @param object glmnetr object as input
#' @param lam value for lam, may be NULL, typically NULL.
#' @param gam value for gam, may be MULL, typically NULL, 0 or 1. 
#' @param comment Default of TRUE to write to console information on lam and gam selected for output.
#' FALSE will suppress this write to console.    
#'
#' @return numerical values for lam and gam for use in plot(), predict() and summary() functions 
#'
#' @noRd
#' 
getlamgam = function(object, lam, gam, comment=0) {
  #  object, lam, gam, comment
  #  object = nested_elastic_fit1$cv_lasso_fit ; xs_new=NULL ; lam=NULL ; gam=NULL ; comment=TRUE ;
  #  object = nested_elastic_fit1$cv_lasso_fit ; xs_new=NULL ; lam=NULL ; gam=0 ; comment=TRUE ;  
  #  object = nested_elastic_fit1$cv_lasso_fit ; xs_new=NULL ; lam=NULL ; gam=1 ; comment=TRUE ;  
  lambda_index = NULL 
  lambdas = object$lambda 
  if ((is.null(lam)) & (is.null(gam))) {
    lam = object$relaxed$lambda.min 
    lambda_index = object$relaxed$index[1,1] 
    gam  = object$relaxed$gamma.min
    if (comment==TRUE) cat(paste0(" (lambda, gamma) pair minimizing CV average deviance is used \n"))
  } else if (inherits(lam,"character") & inherits(gam,"character")) {
    if ((lam=='lambda.1se') & (gam=="gamma.1se")) {
      lam = object$relaxed$lambda.1se 
      lambda_index = object$relaxed$index[2,1] 
      gam  = object$relaxed$gamma.1se
    }  else {
      lam = object$relaxed$lambda.min 
      lambda_index = object$relaxed$index[1,1] 
      gam  = object$relaxed$gamma.min
      if (comment==TRUE) cat(paste0(" (lambda, gamma) pair minimizing CV average deviance is used \n"))
    }
  } else if (!inherits(lam,"numeric") & inherits(gam,"numeric")) {
    if (gam==1) { 
      lam = object$lambda.min 
      lambda_index = object$index[1] 
      if (comment==TRUE) cat(paste0(" lambda minimizing CV average deviance among gamma=1 is used \n"))
    } else if (gam==0) { 
      lambda_index = object$relaxed$index.min.g0[1]  
      lam          = object$relaxed$lambda.min.g0  
      if (comment==TRUE) cat(paste0(" lambda minimizing CV average deviance among gamma=0 is used \n"))
    }
  } else if (inherits(lam,"numeric") & inherits(gam,"numeric")) {
    lambda_index = max(c(1:length(lambdas))[ (lam<=lambdas) ])
  }  else {cat(paste0( "class(lam)=", class(lam), " class(gam)=",class(gam), "\n"))}
  #  returnlist = list(lam=lam, gam=gam, lambda_index=lambda_index, lambdas=list(lambdas))            
  returnlist = list(lam=lam, gam=gam, lambda_index=lambda_index, lambdas=(lambdas))          
  return(returnlist)
}

###############################################################################################################
###############################################################################################################
