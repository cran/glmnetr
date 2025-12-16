###############################################################################################################
###############################################################################################################
#' Give predicteds for elastic net models form a nested..glmnetr() output object.
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
#' @param lambda The lambda value for choice of beta.  If NULL, then 
#' lambda.min is used from the cross validated tuned relaxed model.  We
#' use the term lam instead of lambda as lambda usually denotes a vector 
#' in the package.    
#' @param gamma The gamma value for choice of beta.  If NULL, then 
#' gamma.min is used from the cross validated tuned relaxed model.  We
#' use the term gam instead of gamma as gamma usually denotes a vector 
#' in the package.    
#' @param comment Default of TRUE to write to console information on lam and gam selected for output.
#' FALSE will suppress this write to console.   
#' @param type type of model on which to base predictds. One of "lasso", "ridge" 
#' and "elastic" if elastic net model is fit. 
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
#' @noRd
#' 
predict.cv.glmnetr_0_6_1 = function( object, xs_new=NULL, alpha=NULL, gamma=NULL, lambda=NULL, comment=TRUE, type="lasso", track=0, ...) {
#  object=nested_elastic_fit2 ; xs_new=NULL ; alpha= 0 ;  gamma= 1 ; lambda=0.04 ; comment=TRUE ; type="elastic" ;
  if (track >= 2) { cat( "  in predict.cv.glmnetr_0_6_1   class(object) = ", class(object), "\n") } 
  
  if (type == "elastic") {
    if ( (!is.null(alpha)) & (is.character(alpha)) ) {
      if (alpha != "alpha.min") { if (comment) { cat('  Only character value option for alpha is "alpha.min", alpha is set to default of NULL\n' ) } }
      alpha = NULL 
    }
  }
  
  if ( (!is.null(gamma)) & (is.character(gamma)) ) {
    if (gamma != "gamma.min") { if (comment) { cat('  Only character value option for gamma is "gamma.min", gamma is set to default of NULL\n' ) } }
    gamma = NULL 
  }
  
  if ( (!is.null(lambda)) & (is.character(lambda)) ) {
    if (lambda != "lambda.min") { if (comment) { cat('  Only character value option for lambda is "lambda.min", lambda is set to default of NULL\n' ) } } 
    lambda = NULL 
  }
   
  object_ = object 
  
  if (inherits(object,"nested.glmnetr")) { 
    if (!(type %in% c("lasso", "elastic", "ridge"))) { type = "lasso" }
    if (type == "lasso"  )         { 
      object = object$cv_lasso_fit 
    } else if ( (type == "elastic") & is.null(alpha) & is.null(gamma) ) { 
      object = object$cv_elastic_fit 
      
    } else if ( (type == "elastic") & (!is.null(alpha)) & (!is.null(gamma)) ) { 
      object = object$cv_glmnet_fit_
 
    } else if ( (type == "elastic") & (!is.null(alpha)) & is.null(gamma) ) { 
      if (!is.null(lambda)) {
        if (comment) { cat("  As gamma is NULL, lambda is also set to NULL\n") } 
        lambda = NULL 
      }
      object = object$cv_glmnet_fit_
      
    } else if ( (type == "elastic") & (is.null(alpha)) & (!is.null(gamma)) ) {
      if (!is.null(lambda)) {
        if (comment) { cat("  As alpha is NULL, lambda is also set to NULL\n") } 
        lambda = NULL 
      }
      object = object$cv_glmnet_fit_
      
    } else if (type == "elastic") { 
      object = object$cv_glmnet_fit_ 
      
    } else if (type == "ridge"  )  { 
      object = object$cv_ridge_fit
    }
  } else if (inherits(object,"cv.glmnetr"     )) { type = "lasso" 
  } else if (inherits(object,"cv.glmnetr.el"  )) { type = "elastic" ; alpha = NULL ; 
  } else if (inherits(object,"cv.glmnetr.list")) { type = "elastic"
  } else if (inherits(object,"cv.glmnet"      )) { type = "ridge"  }
  
  #=============================================================================
  
  if (type == "ridge") {
    if (is.null(lambda)) { 
      lam = object$lambda.min
      lambda_index = object$index[1]
    } else {
      lam = lambda 
      lambda_index = sum(object$lambda >= lam)
    }
    #----------
    if (!is.null(xs_new)) {  
      return( predict(object, newx=xs_new, s=lam) ) 
    } else {
      beta_ = as.matrix( predict(object, newx=NULL, s=lam, type = "coefficients" )  )  
      #      beta = beta_[ (beta_!=0)  , 1, drop=FALSE]
      beta_names = row.names(beta_)
      beta_ = as.numeric(beta_)
      names(beta_) = beta_names 
      beta = beta_[ beta_!=0 ] 
      outputobject = list(beta_=beta_, beta=beta)
      return( outputobject )
    }
    
  } else if (type == "lasso") {
    if (is.null(gamma)) {
      names(object$relaxed)
      gam = object$relaxed$gamma.min
      lam = object$relaxed$lambda.min
      gamma_index = object$relaxed$index[1,2]
      lambda_index = object$relaxed$index[1,1]
      cvm_ = object$relaxed$statlist[[ gamma_index ]]$cvm[ lambda_index ]
#      print("HERE 1")
    } else if (is.null(lambda)) {
      gamma_v  = object$relaxed$gamma 
      gamma_index  = c(1:length(gamma_v))[ (round(gamma,digits=6) == round(gamma_v,digits=6)) ]
      cvm_v = object$relaxed$statlist[[ gamma_index ]]$cvm
      lambda_v = object$relaxed$statlist[[ gamma_index ]]$lambda
      lambda_index = which.min( cvm_v )
      lambda.min   = lambda_v[ lambda_index ] 
      cvm_ = cvm_v[ lambda_index ] 
      gam = gamma
      lam = lambda.min
#      print("HERE 2")
    } else {                                    ## ((!isnull(gamma)) & (!is.null(lambda)) )
      gamma_v  = object$relaxed$gamma 
      gamma_index  = c(1:length(gamma_v))[ (round(gamma,digits=6) == round(gamma_v,digits=6)) ]
      lambda_v = object$relaxed$statlist[[ gamma_index ]]$lambda
      # lambda_index  = c(1:length(lambda_v))[ (round(lam,digits=6) == round(lambda_v,digits=6)) ]
      lambda_index = which.min( abs(lambda_v - lambda) ) 
      cvm_ = object$relaxed$statlist[[ gamma_index ]]$cvm[ lambda_index ]
      gam = gamma
      lam = lambda
#      print("HERE 3")
    }
#    cat(" Class(object) =", class(object),   "  gam = ", gam, "  lam=", lam, "\n")
#    cat("   gamma_index = ", gamma_index, " lambda_index= ", lambda_index, "\n")
    #----------
    if (!is.null(object$relaxed)) { class(object) = c("cv.relaxed", "cv.glmnet")
    } else { class(object) = c("cv.glmnet") }  
    #----------
    if (!is.null(xs_new)) {

      predxs = predict(object, newx=xs_new, s=lam, gamma=gam,...) 
      if (comment) {
        if (is.null(gamma)) {         cat("   gamma.min =", gam, "  lambda.min = ", lam,  "  deviance =", round(cvm_, digits=4) , "\n" )
        } else if (is.null(lambda)) { cat("   gamma     =", gam, "  lambda.min = ", lam,  "  deviance =", round(cvm_, digits=4) , "\n" )
        } else {                      cat("   gamma     =", gam, "  lambda = "    , lam,  "  deviance =", round(cvm_, digits=4) , "\n" )  }
      }
      return( predxs )
    } else if (inherits(gam,"numeric") &inherits(lam,"numeric") ) {
      beta_ = as.matrix( predict(object, newx=NULL, s=lam, gamma=gam, type = "coefficients" )  )  
#      beta_
#      class(beta_)
      beta_names = row.names(beta_)
      beta_ = as.numeric(beta_)
      names(beta_) = beta_names 
      beta = beta_[ beta_!=0 ] 
#      betanames = names( beta )
      if (comment) {
        if (is.null(gamma)) {         cat("   gamma.min =", gam, "  lambda.min = ", lam,  "  df =", length(beta), "  deviance =", round(cvm_, digits=4) , "\n" )
        } else if (is.null(lambda)) { cat("   gamma     =", gam, "  lambda.min = ", lam,  "  df =", length(beta), "  deviance =", round(cvm_, digits=4) , "\n" )
        } else {                      cat("   gamma     =", gam, "  lambda = "    , lam,  "  df =", length(beta), "  deviance =", round(cvm_, digits=4) , "\n" )  }
      }
      outputobject = list(beta_=beta_, beta=beta)
      return( outputobject )
    } else { if (comment) { cat( " Check specification of lambda and gamma \n") } } 
    
  } else if ( (type == "elastic") & is.null(alpha) & is.null(gamma) ) {
    alpha.min  = object$vals.elastic[1]
    gamma.min  = object$vals.elastic[2]
    lambda.min = object$vals.elastic[3]
    cvm_       = object$vals.elastic[4]
    #----------
    if (!is.null(object$relaxed)) { class(object) = c("cv.relaxed", "cv.glmnet")
    } else { class(object) = c("cv.glmnet") }  
    #----------
    if (!is.null(xs_new)) {
      predxs = predict(object, newx=xs_new, s=lambda.min, gamma=gamma.min,...) 
      if (comment) { cat("  aplha.min =", alpha.min, "  gamma.min =", gamma.min, "  lambda.min = ", lambda.min, "  deviance =", round(cvm_, digits=4) , "\n" ) } 
      return( predxs )
    } else {
      beta_ = as.matrix( predict(object, newx=NULL, s=lambda.min, gammma=gamma.min, type = "coefficients" )  )  
      beta_names = row.names(beta_)
      beta_ = as.numeric(beta_)
      names(beta_) = beta_names 
      beta = beta_[ beta_!=0 ] 
      if (comment) { cat("  aplha.min =", alpha.min, "  gamma.min =", gamma.min, "  lambda.min = ", lambda.min, "  df =", length(beta), "  deviance =", round(cvm_, digits=4) , "\n" ) } 
      outputobject = list(beta_=beta_, beta=beta)
      return( outputobject )
    }

  } else if ( (type == "elastic") & (!is.null(alpha)) & (!is.null(gamma)) ) {
    #      names( object[[1]]$relaxed ) 
    alpha_v = object$alpha
    lalpha  = length(alpha_v) 
    gamma_v = object[[1]]$relaxed$gamma      
    lgamma  = length(gamma_v) 
    alpha_i = c(1:lalpha)[alpha==round(alpha_v,digits=4)]
    object2 = object[[alpha_i]] 
    #----------
    if (is.null(lambda)) {
      gamma_i = c(1:lgamma)[gamma==round(gamma_v,digits=4)]
      #        str(object2)
      #        names(object2)
      #        names(object2$glmnet.fit)
      #        names(object2$glmnet.fit$relaxed$statlist)
      #        names(object2$relaxed$statlist)
      cvm_ = object2$relaxed$statlist[[gamma_i]]$cvm
      lambda.min = object2$relaxed$statlist[[gamma_i]]$lambda[ which.min(cvm_) ]
#      cat("  lambda.min = ", lambda.min)
      lam = lambda.min 
    } else {
      lam = lambda 
    }
    #----------
    if (!is.null(object2$relaxed)) { class(object2) = c("cv.relaxed", "cv.glmnet")
    } else { class(object2) = c("cv.glmnet") }  
    #----------
    if (!is.null(xs_new)) {
      predxs = predict(object2, newx=xs_new, s=lam, gamma=gamma,...) 
      if (comment) {
        if (is.null(lambda)) {
          cat("  aplha =", alpha, "  gamma =", gam, "  lambda.min = ", lam, "  deviance =", round(object2$measure.min, digits=4) , "\n" )
        } else {
          cat("  aplha =", alpha, "  gamma =", gam, "  lambda = ", lam, "  deviance =", round(object2$measure.min, digits=4) , "\n" ) 
        }
      }
      return( predxs )
    } else {
      beta_ = as.matrix( predict(object2, newx=NULL, s=lam, gammma=gamma, type = "coefficients" )  )  
      beta_names = row.names(beta_)
      beta_ = as.numeric(beta_)
      names(beta_) = beta_names 
      beta = beta_[ beta_!=0 ] 
      outputobject = list(beta_=beta_, beta=beta)
      if (comment) {
        if (is.null(lambda)) {
          cat("  aplha =", alpha, "  gamma =", gamma, "  lambda.min = ", lam, "  df =", length(beta), "  deviance =", round(object2$measure.min, digits=4) , "\n" )
        } else {
          cat("  aplha =", alpha, "  gamma =", gamma, "  lambda = ", lam, "  df =", length(beta), "  deviance =", round(object2$measure.min, digits=4) , "\n" )
        }
      }
      return( outputobject )
    }
    
  } else if ( (type == "elastic") & (!is.null(alpha)) & is.null(gamma) ) { 
    
    #      names( object[[1]]$relaxed ) 
    alpha_v = object$alpha
    lalpha  = length(alpha_v) 
    gamma_v = object[[1]]$relaxed$gamma      
    lgamma  = length(gamma_v) 
    alpha_i = c(1:lalpha)[alpha==round(alpha_v,digits=4)]
    
    object2 = object[[alpha_i]]
    
    gamma.min = object2$relaxed$gamma.min
    lambda.min = object2$relaxed$lambda.min
    
    if (!is.null(object2$relaxed)) { class(object2) = c("cv.relaxed", "cv.glmnet")
    } else { class(object2) = c("cv.glmnet") }  
    
    if (!is.null(xs_new)) {
      predxs = predict(object2, newx=xs_new, s=lambda.min, gamma=gamma.min,...) 
      if (comment) { cat("  aplha =", alpha, "  gamma.min =", gamma.min, "  lambda.min = ", lambda.min, "  deviance =", round(object2$measure.min, digits=4) , "\n" ) } 
      return( predxs )
    } else {
      beta_ = as.matrix( predict(object2, newx=NULL, s=lambda.min, gammma=gamma.min, type = "coefficients" )  )  
      beta_names = row.names(beta_)
      beta_ = as.numeric(beta_)
      names(beta_) = beta_names 
      beta = beta_[ beta_!=0 ] 
      if (comment) { cat("  aplha =", alpha, "  gamma.min =", gamma.min, "  lambda.min = ", lambda.min, "  df =", length(beta), "  deviance =", round(object2$measure.min, digits=4) , "\n" ) } 
      outputobject = list(beta_=beta_, beta=beta)
      return( outputobject )
    }
    
  } else if ( (type == "elastic") & (is.null(alpha)) & !is.null(gamma) ) { 
    names(object)
    names( object[[1]]$relaxed )
    names( object[[1]]$relaxed$statlist )
    names( object[[1]]$relaxed$statlist[[3]] )
    
    alpha_v = object$alpha 
    gamma_v = object[[1]]$relaxed$gamma
    lalpha  = length(alpha_v) 
    lgamma  = length(gamma_v) 
    gamma_i = c(1:lgamma)[gamma==round(gamma_v,digits=4)]
#    statlist = list()
    for ( i_ in c(1:length(object$alpha))) {
      lambda_ = object[[i_]]$relaxed$statlist[[gamma_i]]$lambda
      cvm_    = object[[i_]]$relaxed$statlist[[gamma_i]]$cvm
      cvsd_   = object[[i_]]$relaxed$statlist[[gamma_i]]$cvsd
      cvmup_  = object[[i_]]$relaxed$statlist[[gamma_i]]$cvup
      cvmlo_  = object[[i_]]$relaxed$statlist[[gamma_i]]$cvlo
      nzero_  = object[[i_]]$relaxed$statlist[[gamma_i]]$nzero
      
      alpha.min.cvm_ = min( cvm_ )  
      lambda.min_        = lambda_[ which.min(cvm_) ]
#      df.lambda.min_     = nzero_[ which.min(cvm_) ]
      
      if (i_ == 1) {
        alpha.min_i = 1
        alpha.min.cvm = alpha.min.cvm_ 
        lambda.min = lambda.min_
#        df.lambda.min = df.lambda.min_ 
      } else if (alpha.min.cvm_ <= alpha.min.cvm) {
        alpha.min_i = i_ 
        alpha.min.cvm = alpha.min.cvm_ 
        lambda.min = lambda.min_
#        df.lambda.min = df.lambda.min_ 
      }
#      statlist[[i_]] = object[[i_]]$relaxed$statlist[[gamma_i]]
    }
#    names(statlist) = paste0("a:",alpha_v)
#    names(statlist)
    object2 = object[[ alpha.min_i ]]
    alpha.min = alpha_v[ alpha.min_i ]
    
    if (!is.null(object2$relaxed)) { class(object2) = c("cv.relaxed", "cv.glmnet")
    } else { class(object2) = c("cv.glmnet") }  
    
    if (!is.null(xs_new)) {
      #      predxs = predict(object, newx=xs_new, s=lam, gamma=gam) 
      if (comment) { cat("  aplha.min =", alpha.min, "  gamma =", gamma, "  lambda.min = ", lambda.min, "  deviance =", round(alpha.min.cvm,digits=4) , "\n" ) } 
      predxs = predict(object2, newx=xs_new, s=lambda.min, gamma=gamma,...) 
      return( predxs )
    } else {
      beta_ = as.matrix( predict(object2, newx=NULL, s=lambda.min, gammma=gamma, type = "coefficients" )  )  
      beta_names = row.names(beta_)
      beta_ = as.numeric(beta_)
      names(beta_) = beta_names 
      beta = beta_[ beta_!=0 ] 
      df_ = length(beta)  
      if (comment) { cat("  aplha.min =", alpha.min, "  gamma =", gamma, "  lambda.min = ", lambda.min, "  df =", df_, "  deviance =", round(alpha.min.cvm,digits=4) , "\n" ) }
      outputobject = list(beta_=beta_, beta=beta)
      return( outputobject )
    }

  } else if (type == "elastic") {
    if (comment) { cat('  Not yet implemented for type "elastic" \n') } 
  }
}

###############################################################################################################
###############################################################################################################
