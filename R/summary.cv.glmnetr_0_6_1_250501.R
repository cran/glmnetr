################################################################################
#' Output summary for elastic net models fit within a nested.glmnetr() output object.  
#' 
#' @description 
#' Summarize the cross-validation informed model fit.  The fully penalized
#' (gamma=1) beta estimate will not be given by default but can too be output  
#' using printg1=TRUE.  
#'
#' @param object a nested.glmnetr() output object.  
#' @param type one of c("lasso", "elastic", "ridge") to select for summarizing,  
#' with default of "lasso".
#' @param printg1 TRUE to also print out the fully penalized lasso beta, else FALSE to suppress.  
#' @param orderall By default (orderall=FALSE) the order terms enter into the lasso model 
#' is given for the number of terms that enter in lasso minimizing loss model.  If 
#' orderall=TRUE then all terms that are included in any lasso fit are described.
#' @param betatol beta values less than betatol are set to 0 when they are 
#' approaching the rounding error of the machine architecture. Default is set 
#' to 1e-14. 
 
#' @param ... Additional arguments passed to the summary function.  
#'
#' @return Coefficient estimates (beta) 
#' 
#' @seealso
#'   \code{\link{predict.cv.glmnetr}} , \code{\link{nested.glmnetr}} 
#' 
#' @export
#'
#' @noRd
#' 
summary.cv.glmnetr_0_6_1 = function(object, type=NULL, printg1="FALSE", orderall=FALSE, betatol=1e-14, ...) {
  # object = nested_elastic_fit1 ; type = "lasso" ; type = "elastic" ; betatol=1e-14 ;
  if (is.null(type)) { type = "lasso" }
  if (inherits(object,"nested.glmnetr")) { 
    if        (type == "lasso"  ) { 
      object = object$cv_lasso_fit
    } else if (type == "elastic") { 
      if ( is.null(object$cv_elastic_fit ) ) {
        type = "lasso" 
        object = object$cv_lasso_fit
#        cat(paste(" There is no object for the elastic net fits. Presuambly elastic net not run. Program will stop.\n\n"))
        stop(.call="There is no object for the elastic net fits. Presuambly elastic net not run. Program will stop.")
      } else { 
        object = object$cv_elastic_fit 
        vals.elastic = object$vals.elastic
        alpha.min.r = vals.elastic[1]
      }
    } else if (type == "ridge" ) { 
      object = object$cv_ridge_fit  
#      print(names(object))
    }
  }

  lambda.min.g1   = object$lambda.min  
  lambda_index.g1 = object$index[1] 
  df.min.g1       = object$nzero[lambda_index.g1]
  deviance.g1     = object$cvm[lambda_index.g1]
  beta.g1         = object$glmnet.fit$beta[ , lambda_index.g1 ]
  beta.g1[ abs(beta.g1) < betatol ] = 0 
  beta.g1         = beta.g1 [ beta.g1 != 0 ]
  df.beta.g1      = length( beta.g1 )    ## glmnet nzero may include linear dependent sets which drop out during fitting

  
  if (type == "ridge") { 
    lambda.min = object$lambda.min
  } 
  
  if (type %in% c("lasso","elastic")) { 
    gamma          = object$relaxed$gamma   
    lambda.min.r   = object$relaxed$lambda.min  
    gamma.min.r    = object$relaxed$gamma.min
    lambda_index.r = object$relaxed$index[1,1] 
    gamma_index.r  = object$relaxed$index[1,2]
    deviance.r     = object$relaxed$statlist[[gamma_index.r]][[2]][lambda_index.r]
    df.min.r       = object$nzero[ lambda_index.r ]
    betag0         = object$glmnet.fit$relaxed$beta[,lambda_index.r]
    betag1         = object$glmnet.fit$beta[,lambda_index.r]
    beta_          = (1 - gamma.min.r) * betag0 + gamma.min.r * betag1
    beta_[ abs(beta_) < betatol ] = 0 
    beta.r         = beta_[ (beta_!=0) ]
    #  length( beta.r )
    
    lambda.min.g0 = object$relaxed$lambda.min.g0 
    index.min.g0  = object$relaxed$index.min.g0
    df.min.g0     = object$nzero[index.min.g0]
    deviance.g0   = object$relaxed$statlist[[1]]$cvm[index.min.g0]
    beta.g0       = object$glmnet.fit$relaxed$beta[,index.min.g0] * 1 
    beta.g0[ abs(beta.g0) < betatol ] = 0 
    beta.g0       = beta.g0[ beta.g0 != 0 ]
    df.min.g0     = length( beta.g0 )      
  } 
  
  if         (type == "lasso" ) { 
    cat( paste( "  For the lasso model fits:\n"))
    cat(paste0( "  The minimum deviance is obtained at gamma = " , gamma.min.r, " and lambda = ", round(lambda.min.r,digits=6) , "\n"  )) 
    cat(paste0( "  with df (number of non-zero terms) = ", df.min.r , ", average deviance = " , round(deviance.r, digits=6) , " and beta = " , "\n" ))  
    print( beta.r ) 
    cat(paste0( "\n  The fully relaxed (gamma=0) minimum deviance is obtained at lambda = ", round(lambda.min.g0,digits=6) , "\n"  )) 
    cat(paste0( "  with df (number of non-zero terms) = ", df.min.g0 , ", average deviance = " , round(deviance.g0, digits=6) ,  " and beta = " , "\n" ))  
    print( beta.g0 ) 
    if (printg1 == 1) {
      cat(paste0( "\n  The UNrelaxed (gamma=1) minimum devaince is obtained at lambda = ", round(lambda.min.g1,digits=6) , "\n"  )) 
      cat(paste0( "  with df (number of non-zero terms) = ", df.beta.g1 , ", average deviance = " , round(deviance.g1, digits=6) , "\n" ))  
      print( beta.g1 ) 
    }
  } else if (type == "elastic") { 
    cat( paste( "  For the elastic net model fits:\n"))
    cat(paste0( "  The minimum devaince is obtained at alpha = ", alpha.min.r, ", gamma = " , gamma.min.r, " and lambda = ", round(lambda.min.r,digits=6) , "\n"  )) 
    cat(paste0( "  with df (number of non-zero terms) = ", df.min.r , ", average deviance = " , round(deviance.r, digits=6) , " and beta = " , "\n" ))  
    print( beta.r ) 
  } else if (type == "ridge" ) { 
    cat( paste( "  For the ridge regression model fits:\n")) 
    cat(paste0( "  The minimum devaince is obtained at lambda = ", round(lambda.min,digits=6) , "\n"  )) 
    cat(paste0( "  with df (number of non-zero terms) = ", df.beta.g1 , ", average deviance = " , round(deviance.g1, digits=6) , " and beta = " , "\n" ))  
    print( beta.g1 )  
  }

  printorder = 1 
  if (type == "ridge" ) {
    printorder = 0 
  } else if (type == "elastic") {
    if (vals.elastic[1] == 0) { printorder = 0 }
  }
  
  if ( printorder ) {
    if (type == "lasso") {
      cat(paste0( "\n  Order coefficients entered into the lasso model (1st to last):" , "\n"))    
    } else if (type == "elastic") {
      cat(paste0( "\n  Order coefficients entered into the elastic net model (1st to last):" , "\n"))    
    } 
    betag1all = as.matrix( object$glmnet.fit$beta ) 
    #  print( class(betag1all))
    nrows = dim(betag1all)[1]
    ncols = dim(betag1all)[2]
    ne0.index = c(rep(0,nrows))
    for (k_ in c(1:nrows)) {
      ne0.index[k_] = min( c(1:(ncols+1))[(c(betag1all[k_,],1)!=0)] )           ## add an extra element so those not entered get a large order number 
    }
    if (orderall==TRUE) { print( rownames(betag1all)[order(ne0.index)] )        ## table(table(order(ne0.index)))
    } else { print( rownames(betag1all)[order(ne0.index)][1:df.min.g1] ) }     
  }
}

###############################################################################################################
###############################################################################################################
