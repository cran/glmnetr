############################################################################################################
#' Get predicteds or coefficients using a glmnetr output object
#' 
#' @description 
#' Give predicteds based upon a glmnetr() output object. Because the glmnetr() function 
#' has no cross validation information, lambda and gamma must be specified.  To choose 
#' lambda and gamma based upon cross validation one may use the cv.glmnetr() or nested.glmnetr() 
#' and the corresponding predict() functions.  
#'
#' @param object A glmnetr output object
#' @param xs_new A desing matrix for predictions
#' @param lam The value for lambda for determining the lasso fit.  Required.  
#' @param gam The value for gamma for determining the lasso fit. Required.   
#' @param ... Additional arguments passed to the predict function.  
#'
#' @return Coefficients or predictions using a glmnetr output object.  When
#' outputting coefficients (beta), creates a list with the first element, beta_, 
#' including 0 and non-0 terms and the second element, beta, including only 
#' non 0 terms.
#' 
#' @seealso
#'   \code{\link{nested.glmnetr}} 
#' 
#' @noRd
#' 
predict.glmnetr_0_5_5 = function( object, xs_new=NULL, lam=NULL, gam=NULL, ...) {
  if (is.null(gam)) {gam=1} 
  if (inherits(lam,"numeric") & inherits(gam,"numeric")) {
    if (gam < 0) { gam = 0 
    } else if (gam > 1) { gam = 1 }
    
    lambda_index = sum(object$lambda>=lam)
    lambda_index = min(lambda_index, dim(object$betag0)[2])

    betag0 = object$betag0[,lambda_index]
    betag1 = object$betag1[,lambda_index]
    beta_   = (1-gam) * betag0 + gam * betag1
    beta_
    if (!is.null(xs_new)) { 
      score_new = xs_new %*% beta_ 
      outputobject = score_new 
    } else { 
      beta = beta_[(beta_!=0)]
      outputobject=list(beta_=beta_, beta=beta)
    }
    return( outputobject ) 
  } else { print( " Check specification of lamb(da) and gam(ma) ")}
}

###############################################################################################################
###############################################################################################################
