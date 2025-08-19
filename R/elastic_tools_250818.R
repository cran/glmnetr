##### elstic_tools_250522

# #' Title
# #'
# #' @param x 
# #' @param alpha 
# #'
# #' @return predictions 
# #'
# #' @noRd
# #' 
# pred_elastic_alp_fix = function(x, alpha, xs=NULL) {
#   alpha_i = c(1:length(x$alpha))[(alpha == round(x$alpha,digits=6))]    ## tests did not work without round() 
#   object  = x$cv_glmnet_fit_[[alpha_i]]  
  
#   family = object$family                       ## 250818
# ##  object2 = object[[alpha_i]]                 ## 250818 
  
#  #        alpha_i = object$index.elastic[1]
#  gam_ = object$relaxed$gamma.min
#  lambda_v = object$lambda
#  lambda.min = object$relaxed$lambda.min
#  beta0   = object$glmnet.fit$relaxed$beta # 
#  beta1   = object$glmnet.fit$beta # [,lambda_i]
#  beta = ( (1-gam_) * beta0 + gam_ * beta1 )
#  if (family == "cox") {                                                       
#    a00 = 0 
#    a01 = 0 
#  } else { 
#    a00 = object$glmnet.fit$relaxed$a0 
#    a01 = object2$glmnet.fit$a0
#  }
#  intercept = ( (1-gam_) * a00 + gam_ * a01 )           
#  if (!is.null(xs)) {
#    xb = xs %*% beta     
#    xb = xs %*% beta    
#    if (family != "cox") xb = xb + intercept      
#  } else {
#    xb = NULL 
#  }
#  nzero = object$nzero
#  return(list(intercept=intercept, beta=beta, nzero=nzero, gam=gam_,xb=xb)) 
# }

#' Title
#'
#' @param x 
#' @param gam 
#'
#' @return predictions 
#' 
#' @noRd
#' 
pred_elastic_gam_fix = function(x, gam, xs=NULL, track=0) {
  if (track >= 2) { cat("  In pred_elastic_gam_fix() \n") } 
  object = x$cv_glmnet_fit_
  family = object$family
#  class(object)
#  class(object[[1]])
#  names(object)
#  names(object[[1]]$relaxed)
  gam_ = gam 
  gamma = object[[1]]$relaxed$gamma 
  lgamma = length( gamma ) 
  gamma_i = c(1:lgamma)[(gam == round(gamma,digits=6))]
  alpha_v = object$alpha 
  lalpha = length(alpha_v)
  # lambda_v = object[[1]]$relaxed$statlist[[gamma_i]]$lambda
  for ( i_ in c(1:lalpha)) {
#    names( object[[i_]] )
#    names( object[[i_]]$glmnet.fit ) 
#    names( object[[i_]]$relaxed$statlist[[gamma_i]] ) 
    mincvm = min( object[[i_]]$relaxed$statlist[[gamma_i]]$cvm )
    lambda_it = which.min( object[[i_]]$relaxed$statlist[[gamma_i]]$cvm )
    lambda_t = object[[i_]]$relaxed$statlist[[gamma_i]]$lambda[lambda_it]
    if (i_ == 1) {
      measure.min = mincvm 
      alpha_i = i_ 
      lambda_i = lambda_it
      lambda = lambda_t 
    } else if (  mincvm < measure.min ) {
      measure.min = mincvm 
      alpha_i = i_             
      lambda_i = lambda_it
      lambda = lambda_t 
    }
  }
  alpha_ = alpha_v[alpha_i] 
  if (track >= 2) { print( c(alpha_i, alpha_, lambda_i, lambda) ) } 
  object2 = object[[alpha_i]]
  nzero = object2$nzero 
#  names(object2)
#  names(object2$glmnet.fit$relaxed)
  beta0 = object2$glmnet.fit$relaxed$beta[,lambda_i] # 
  if (family == "cox") { a00 = 0 } else { a00 = object2$glmnet.fit$relaxed$a0[lambda_i] } 
  beta1 = object2$glmnet.fit$beta[,lambda_i] 
  if (family == "cox") { a01 = 0 } else { a01 = object2$glmnet.fit$a0[lambda_i] }         
  beta = ( (1-gam) * beta0 + gam * beta1 )
  intercept = ( (1-gam) * a00 + gam * a01 )                                   
  if (!is.null(xs)) {
    xb = xs %*% beta                                                          
    if (family != "cox") xb = xb + intercept                                  
  } else {
    xb = NULL 
  }
  return(list(intercept=intercept, beta=beta, nzero=nzero, alp=alpha_, alpha_i=alpha_i, lambda=lambda, lambda_i=lambda_i, xb=xb)) 
}
