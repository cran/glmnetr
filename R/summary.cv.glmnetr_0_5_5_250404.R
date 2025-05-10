################################################################################
##### summary.cv.glmnetr_yymmdd.R ##############################################
################################################################################
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
#'   \code{\link{predict.cv.glmnetr}} , \code{\link{nested.glmnetr}} 
#' 
#' @noRd
#' 
summary.cv.glmnetr_0_5_5 = function(object, printg1="FALSE", orderall=FALSE, ...) {
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
  cat(paste0( "  The relaxed minimum is obtained for lambda = ", round(lambda.min.r,digits=8) , " and gamma = " , gamma.min.r , "\n"  )) 
  cat(paste0( "  with df (number of non-zero terms) = ", df.min.r , ", average deviance = " , round(deviance.r, digits=6) , " and beta = " , "\n" ))  
  print( beta.min.r ) 
  
  cat(paste0( "\n  The fully relaxed (gamma=0) minimum is obtained at lambda = ", round(lambda.min.g0,digits=8) , "\n"  )) 
  cat(paste0( "  with df (number of non-zero terms) = ", df.min.g0 , ", average deviance = " , round(deviance.g0, digits=6) ,  " and beta = " , "\n" ))  
  print( beta.min.g0 )  
  
  cat(paste0( "\n  The UNrelaxed (gamma=1) minimum is obtained at lambda = ", round(lambda.min.g1,digits=8) , "\n"  )) 
  cat(paste0( "  with df (number of non-zero terms) = ", df.min.g1 , ", average deviance = " , round(deviance.g1, digits=6) , "\n" ))  
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
