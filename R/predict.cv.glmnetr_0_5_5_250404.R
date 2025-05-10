###############################################################################################################
###############################################################################################################

#' Give predicteds based upon a cv.glmnetr() output object.
#' 
#' @description 
#' Give predicteds based upon a cv.glmnetr() output object. By default lambda and gamma 
#' are chosen as the minimizing values for the relaxed lasso model.  If gam=1 and lam=NULL 
#' then the best unrelaxed lasso model is chosen and if gam=0 and lam=NULL then
#' the best fully relaxed lasso model is selected.  
#' 
#' @param object  A cv.glmnetr (or nested.glmnetr) output object.
#' @param xs_new The predictor matrix.  If NULL, then betas are provided.  
#' @param lam The lambda value for choice of beta.  If NULL, then 
#' lambda.min is used from the cross validated tuned relaxed model.  We
#' use the term lam instead of lambda as lambda usually denotes a vector 
#' in the package.    
#' @param gam The gamma value for choice of beta.  If NULL, then 
#' gamma.min is used from the cross validated tuned relaxed model.  We
#' use the term gam instead of gamma as gamma usually denotes a vector 
#' in the package.    
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
#' @noRd
#' 
predict.cv.glmnetr_0_5_5 = function( object, xs_new=NULL, lam=NULL, gam=NULL, comment=TRUE, ...) {
  if (inherits(object,"nested.glmnetr")) { object = object$cv.glmnet.fit }
  #  object=cv.cox.fit ; xs_new=NULL; lam=NULL; gam=NULL; comment=TRUE ;
  #  object=cv.cox.fit ; xs_new=NULL; lam=NULL; gam=0; comment=TRUE ;
  #  object = cv.glmnetr.fit.hp
  #  lam=NULL ; gam=NULL 
  #  lam="lambda.1se" ; gam="gamma.1se" 
  
  #  family = object$family 
  family = object$sample$family 
  lamgam = getlamgam(object, lam, gam, comment)
  lam    = lamgam$lam
  gam    = lamgam$gam
  lambda_index = lamgam$lambda_index
  lambdas      = lamgam$lambdas
  lamgam$lambdas[lamgam$lambda_index]
  
  if (inherits(lam,"numeric") & inherits(gam,"numeric")) {
    if (is.null(lambda_index)) { 
      #      lambda_index = c(1:length(lambdas))[lambdas==lam]  
      lambda_index = sum(lambdas>=lam)
    }
    #    cat(paste0(" lambda = " , lam, "   lambda_index = " , lambda_index, "   df = " , object$nzero[lambda_index] , "   gamma = " , gam, "\n")) 
    
    if (gam < 0) { gam = 0 
    } else if (gam > 1) { gam = 1 }
    
    if (family %in% c("binomial","gaussian")) {
      a0g1 = object$glmnet.fit$a0[lambda_index]
      a0g0 = object$glmnet.fit$relaxed$a0[lambda_index]
      a0   = (1-gam) * a0g0 + gam * a0g1
    } else { a0 = 0 }
    betag0 = object$glmnet.fit$relaxed$beta[,lambda_index]
    betag1 = object$glmnet.fit$beta[,lambda_index]
    beta_   = (1-gam) * betag0 + gam * betag1
    beta_
    if (!is.null(xs_new)) { 
      score_new = as.numeric( xs_new %*% beta_ ) 
      if (family %in% c("binomial","gaussian")) {
        score_new=score_new + a0 
      }
      outputobject = score_new
    } else { 
      if (family %in% c("binomial","gaussian")) { 
        beta_ = c(a0,beta_) 
        names(beta_)[1] = "Int" 
      }
      beta = beta_[(beta_!=0)]
      outputobject=list(beta_=beta_, beta=beta)
    }
    return( outputobject ) 
  } else { print( " Check specification of lamb(da) and gam(ma) ")}
}

###############################################################################################################
###############################################################################################################
