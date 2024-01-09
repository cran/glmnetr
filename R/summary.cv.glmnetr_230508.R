###############################################################################################################
###############################################################################################################
#' Output summary of a cv.glmnetr() output object.  
#' 
#' @description 
#' Summarize the cross-validation informed model fit.  The fully penalized
#' (gamma=1) beta estimate will not be given by default but can too be output  
#' using printg1=TRUE.  
#'
#' @param object a cv.glmnetr() output object.  
#' @param printg1 TRUE to also print out the fully penalized lasso beta, else FALSE to suppress.  
#' @param orderall By default (orderall=FALSE) the order terms enter into the lasso model 
#' is given for the number of terms that enter in lasso minimizing loss model.  If 
#' orderall=TRUE then all terms that are included in any lasso fit are described.    
#' @param ... Additional arguments passed to the summary function.  
#'
#' @return Coefficient estimates (beta) 
#' 
#' @seealso
#'   \code{\link{cv.glmnetr}} , \code{\link{nested.glmnetr}} 
#' 
#' @export
#'
#' @examples
#' # set seed for random numbers, optionally, to get reproducible results
#' set.seed(82545037)
#' sim.data=glmnetr.simdata(nrows=100, ncols=100, beta=NULL)
#' xs=sim.data$xs 
#' y_=sim.data$y_ 
#' event=sim.data$event
#' # for this example we use a small number for folds_n to shorten run time 
#' cv.glmnetr.fit = cv.glmnetr(xs, NULL, y_, NULL, family="gaussian", folds_n=3, limit=2) 
#' summary(cv.glmnetr.fit)
#' 
summary.cv.glmnetr = function(object, printg1="FALSE", orderall=FALSE, ...) {
  if (inherits(object,"nested.glmnetr")) { object = object$cv.glmnet.fit }
  
  gamma = object$relaxed$gamma 
  
  lambda.min.g1   = object$lambda.min  
  lambda_index.g1 = object$index[1] 
  df.min.g1       = object$nzero[lambda_index.g1]
  deviance.g1     = object$relaxed$statlist[[length(gamma)]][[2]][lambda_index.g1]
  beta.min.g1     = object$glmnet.fit$beta[ , lambda_index.g1 ]
  beta.min.g1     = beta.min.g1 [ beta.min.g1 != 0 ]
  #  beta.min.g1
  #  length( beta.min.g1 ) 
  
  lambda.min.r   = object$relaxed$lambda.min  
  gamma.min.r    = object$relaxed$gamma.min
  lambda_index.r = object$relaxed$index[1,1] 
  gamma_index.r  = object$relaxed$index[1,2]
  df.min.r       = object$nzero[ lambda_index.r ]
  deviance.r     = object$relaxed$statlist[[gamma_index.r]][[2]][lambda_index.r]
  betag0 = object$glmnet.fit$relaxed$beta[,lambda_index.r]
  betag1 = object$glmnet.fit$beta[,lambda_index.r]
  beta_  = (1 - gamma.min.r) * betag0 + gamma.min.r * betag1
  beta.min.r = beta_[(beta_!=0)]
  #  length( beta.min.r )
  
  lambda.min.g0 = object$relaxed$lambda.min.g0 
  index.min.g0  = object$relaxed$index.g0[1] 
  df.min.g0     = object$nzero[index.min.g0]
  deviance.g0   = object$relaxed$statlist[[1]][[2]][index.min.g0]
  beta.min.g0   = object$glmnet.fit$relaxed$beta[,index.min.g0] * 1 
  beta.min.g0   = beta.min.g0[ beta.min.g0!=0 ]
  df.min.g0     = length( beta.min.g0 )                       ## lasso may inlcude linear dependent sets which drop out with unrestricted fit 
  #  beta.min.g0
  
  #  beta.min.g0   = object$glmnet.fit$beta[,index.min.g0] * 1 
  #  beta.min.g0   = beta.min.g0[ beta.min.g0!=0 ]
  #  beta.min.g0
  
  #  colSums( object$glmnet.fit$relaxed$beta != 0 )
  #  colSums( object$glmnet.fit$beta != 0 )
  #  colSums( object$glmnet.fit$relaxed$beta != 0 ) - colSums( object$glmnet.fit$beta != 0 )
  
  # object$glmnet.fit$relaxed$beta[,9]
  cat(paste0( "\n  The relaxed minimum is obtained for lambda = ", round(lambda.min.r,digits=8) , " , index = " , lambda_index.r , " and gamma = " , gamma.min.r , "\n"  )) 
  cat(paste0( "  with df (number of non-zero terms) = ", df.min.r , ", average deviance = " , round(deviance.r, digits=6) , " and beta = " , "\n" ))  
  print( beta.min.r ) 
  
  cat(paste0( "\n  The fully relaxed (gamma=0) minimum is obtained for lambda = ", round(lambda.min.g0,digits=8) , " and index = " , index.min.g0 , "\n"  )) 
  cat(paste0( "  with df (number of non-zero terms) = ", df.min.g0 , ", average deviance = " , round(deviance.g0, digits=6) ,  " and beta = " , "\n" ))  
  print( beta.min.g0 )  
  
  cat(paste0( "\n  The UNrelaxed (gamma=1) minimum is obtained for lambda = ", round(lambda.min.g1,digits=8) , " and index = " , lambda_index.g1 , "\n"  )) 
  cat(paste0( "  with df (number of non-zero terms) = ", df.min.g1 , ", average deviance = " , round(deviance.g1, digits=6) , "\n\n" ))  
  #  cat(paste0( "  with df (number of non-zero terms) = ", df.min.g1 , " and beta = " , "\n" ))  
  if (printg1==TRUE) { print( beta.min.g1 ) }

  cat(paste0( "\n  Order coefficients entered into the lasso model (1st to last):" , "\n"))
  betag0all = object$glmnet.fit$relaxed$beta 
  nrows = dim(betag0all)[1]
  ncols = dim(betag0all)[2]
  ne0.index = c(rep(0,nrows))
  for (k_ in c(1:nrows)) {
    ne0.index[k_] = min( c(1:(ncols+1))[(c(betag0all[k_,],1)!=0)] )
  }
  #ne0.index
  #order(ne0.index)
  #ne0.index[order(ne0.index)]
  # sort(ne0.index)
  if (orderall==TRUE) { print( rownames(betag0all)[order(ne0.index)] )
  } else { print( rownames(betag0all)[order(ne0.index)][1:df.min.g1] ) } 
}

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
#'   \code{\link{predict.glmnetr}} , \code{\link{cv.glmnetr}} , \code{\link{nested.glmnetr}} 
#' 
#' @export
#'
#' @examples
#' # set seed for random numbers, optionally, to get reproducible results
#' set.seed(82545037)
#' sim.data=glmnetr.simdata(nrows=200, ncols=100, beta=NULL)
#' xs=sim.data$xs 
#' y_=sim.data$y_ 
#' event=sim.data$event
#' # for this example we use a small number for folds_n to shorten run time 
#' cv.glmnetr.fit = cv.glmnetr(xs, NULL, y_, NULL, family="gaussian", folds_n=3, limit=2) 
#' predict(cv.glmnetr.fit)
#' 
predict.cv.glmnetr = function( object, xs_new=NULL, lam=NULL, gam=NULL, comment=TRUE, ...) {
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
      }
      beta = beta_[(beta_!=0)]
      outputobject=list(beta_=beta_, beta=beta)
    }
    return( outputobject ) 
  } else { print( " Check specification of lamb(da) and gam(ma) ")}
}

###################################################################################################

#' get numerical values for lam and gam 
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
getlamgam = function(object, lam, gam, comment) {
  lambda_index = NULL
  lambdas = object$lambda 
  if ((is.null(lam)) & (is.null(gam))) {
    lam = object$relaxed$lambda.min 
    lambda_index = object$relaxed$index[1,1] 
    gam  = object$relaxed$gamma.min
    if (comment==TRUE) cat(paste0("\n (lambda, gamma) pair minimizing CV average deviance is used \n\n"))
  } else if (inherits(lam,"character") & inherits(gam,"character")) {
    if ((lam=='lambda.1se') & (gam=="gamma.1se")) {
      lam = object$relaxed$lambda.1se 
      lambda_index = object$relaxed$index[2,1] 
      gam  = object$relaxed$gamma.1se
    }  else {
      lam = object$relaxed$lambda.min 
      lambda_index = object$relaxed$index[1,1] 
      gam  = object$relaxed$gamma.min
    }
  } else if (!inherits(lam,"numeric") & inherits(gam,"numeric")) {
    if (gam==1) { 
      lam = object$lambda.min 
      lambda_index = object$index[1] 
      if (comment==TRUE) cat(paste0("\n lambda minimizing CV average deviance among gamma=1 is used \n\n"))
    } else if (gam==0) { 
      lambda_index = object$relaxed$index.g0[1]  
      lam          = lambdas[lambda_index] 
      if (comment==TRUE) cat(paste0("\n lambda minimizing CV average deviance among gamma=0 is used \n\n"))
    }
  } else if (inherits(lam,"numeric") & inherits(gam,"numeric")) {
    lambda_index = max(c(1:length(lambdas))[ (lam<=lambdas) ])
  }  else {cat(paste0( "class(lam)=", class(lam), " class(gam)=",class(gam), "\n"))}
#  returnlist = list(lam=lam, gam=gam, lambda_index=lambda_index, lambdas=list(lambdas))                     ## CHANGE 28MAR23 
  returnlist = list(lam=lam, gam=gam, lambda_index=lambda_index, lambdas=(lambdas))                          ## CHANGE 28MAR23
  return(returnlist)
}

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
#'   \code{\link{glmnetr}} , \code{\link{cv.glmnetr}} , \code{\link{nested.glmnetr}} 
#' 
#' @export
#'
#' @examples
#' \donttest{
#' set.seed(82545037)
#' sim.data=glmnetr.simdata(nrows=200, ncols=100, beta=NULL)
#' xs=sim.data$xs 
#' y_=sim.data$yt
#' event=sim.data$event
#' glmnetr.fit = glmnetr( xs, NULL, y_, event, family="cox")
#' betas = predict(glmnetr.fit,NULL,exp(-2),0.5 )
#' betas$beta
#' }
#' 
predict.glmnetr = function( object, xs_new=NULL, lam=NULL, gam=NULL, ...) {
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
