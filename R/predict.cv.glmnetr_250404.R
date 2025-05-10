###############################################################################################################
###############################################################################################################

#' Give predicteds for elastic net models form a nested.glmnetr() output object.
#' 
#' @description 
#' Give predicteds based upon a cv.glmnetr() output object. By default lambda and gamma 
#' are chosen as the minimizing values for the relaxed lasso model.  If gam=1 and lam=NULL 
#' then the best unrelaxed lasso model is chosen and if gam=0 and lam=NULL then
#' the best fully relaxed lasso model is selected.  
#' 
#' @param object  A cv.glmnetr (or nested.glmnetr) output object.
#' @param xs_new The predictor matrix.  If NULL, then betas are provided.  
#' @param alpha A specific value of alpha for plotting.  Used only when type is 
#' set to "elastic". Specifies which alpha is to be used for deviance plots.
#' Default is "alpha.min", else must be an element of the alpha vector used in 
#' running the elastic net model. This can be reviewed using summary(fit) 
#' where fit is a nested.glmnetr() output object. Note, alpha is 1 for the 
#' lasso model and alpha is 0 for the ridge model.
#' @param gamma The gamma value for choice of beta.  If NULL, then 
#' gamma.min is used from the cross validated tuned relaxed model.  We
#' use the term gam instead of gamma as gamma usually denotes a vector 
#' in the package.    
#' @param lambda The lambda value for choice of beta.  If NULL, then 
#' lambda.min is used from the cross validated tuned relaxed model.  We
#' use the term lam instead of lambda as lambda usually denotes a vector 
#' in the package.    
#' @param type type of model on which to base predictds. One of "lasso", "ridge" 
#' and "elastic" if elastic net model is fit. 
#' @param comment Default of TRUE to write to console information on lam and gam selected for output.
#' FALSE will suppress this write to console.   
#' @param ... Additional arguments passed to the predict function.  
#'
#' @return Either predicteds (xs_new*beta estimates based upon the predictor matrix xs_new)  
#' or model coefficients, based upon a cv.glmnetr() output object.  When 
#' outputting coefficients (beta), creates a list 
#' with the first element, beta_, including 0 and non-0 terms and the 
#' second element, beta, including only non 0 terms. 
#' 
#' @seealso
#'   \code{\link{summary.cv.glmnetr}} , \code{\link{nested.glmnetr}} 
#' 
#' @export
#'
predict.cv.glmnetr = function( object, xs_new=NULL, alpha=NULL, gamma=NULL, lambda=NULL, type="lasso", comment=TRUE, ...) {
#  cat( "  in predict.cv.glmnetr   class(object) = ", class(object), "\n")
  if (inherits(object,"nested.glmnetr")) { 
    if ( substr(object$version[2],1,21) == "glmnetr version 0.6-1" ) { 
      predict.cv.glmnetr_0_6_1( object, xs_new=xs_new, alpha=alpha, gamma=gamma, lambda=lambda, comment=comment, type=type, ...)
    } else {
      predict.cv.glmnetr_0_5_5( object, xs_new=xs_new, lam=lambda, gam=gamma, comment=comment, ...)
    }
    #  } else if (is.null(object$vals.elastic)) {
  } else if (inherits(object,"cv.glmnetr")) {
    type = "lasso"
    predict.cv.glmnetr_0_6_1( object, xs_new=xs_new, gamma=gamma, lambda=lambda, comment=comment, type=type, ...)
  } else if (inherits(object,"cv.glmnetr.el")) {
    type = "elastic"
    predict.cv.glmnetr_0_6_1( object, xs_new=xs_new, gamma=gamma, lambda=lambda, comment=comment, type=type, ...)
  } else if (inherits(object,"cv.glmnetr.list")) {
    type = "elastic"
    predict.cv.glmnetr_0_6_1( object, xs_new=xs_new, alpha=alpha, gamma=gamma, lambda=lambda, comment=comment, type=type, ...)
  } else {
    type = "lasso"
    predict.cv.glmnetr_0_6_1( object, xs_new=xs_new, gamma=gamma, lambda=lambda, comment=comment, type=type, ...)
  }
}

###############################################################################################################
###############################################################################################################
