################################################################################
##### plot.cv.glmnetr_yymmdd.R #################################################
################################################################################
#' Plot the relaxed lasso coefficients.  
#' 
#' @description 
#' Plot the relaxed lasso coefficients from a nested.glmnetr() output 
#' object. One may specify a single value for gamma. If gamma is unspecified 
#' (NULL), then the gamma which minimizes loss is used.
#'
#' @param x Either a glmnetr, cv.glmnetr or a nested.glmnetr output object.  
#' @param type one of c("lasso", "elastic", "ridge") to plot the deviance 
#' curves of the respective model fit. Default is "lasso" for tuned relaxed lasso. 
#' @param alpha A specific value of alpha for plotting.  Used only when type is 
#' set to "elastic". Specifies which alpha is to be used for deviance plots.
#' Default is "alpha.min", else must be an element of the alpha vector used in 
#' running the elastic net model. This can be reviewed using summary(fit) 
#' where fit is a nested.glmnetr() output object. Note, alpha is 1 for the 
#' lasso model and alpha is 0 for the ridge model.
#' @param gamma A specific level of gamma for plotting.  By default gamma.min 
#' from the deviance minimizing (lambda.min, gamma.min) pair is used.  
#' @param lambda.lo A lower limit of lambda for plotting.  
#' @param title A title for the plot  
#' @param comment Default of TRUE to write to console information on lam and gamma selected for output.
#' FALSE will suppress this write to console.  
#' @param ... Additional arguments passed to the plot function.  
#'
#' @return This program returns a plot to the graphics window, and may provide 
#' some numerical information to the R Console.  If the input object is from a
#' nested.glmnetr or cv.glmnetr object, and gamma is not specified, 
#' then the gamma.min from
#' the deviance minimizing (lambda.min, gamma.min) pair will be used, and the
#' minimizing lambda.min will be indicated by a vertical line.  Also, if one 
#' specifies gamma=0, the lambda which minimizes deviance for the restricted 
#' set of models where gamma=0 will indicated by a vertical line.  
#' 
#' @seealso 
#'   \code{\link{plot.cv.glmnetr}} , \code{\link{plot.nested.glmnetr}} , \code{\link{glmnetr}} 
#' 
#' @export
#' 
#' @importFrom graphics abline axis lines
#'
#' @noRd
#'  
plot.glmnetr_0_6_1 = function(x, type="lasso", alpha=NULL, gamma=NULL, lambda.lo=NULL, title=NULL, comment=TRUE, ...) {
#  cat(paste("\n plot. 1  type = ", type, ", gamma = ", gamma, ", alpha = ", alpha, ", lambda.lo = ", lambda.lo, ", title = ", title, ", comment = ", comment, "\n"))
  #  type="lasso" ; type="lasso" ; type="elastic" ;gamma=NULL ; alpha=NULL ; lambda.lo=NULL ; title=NULL ; comment=TRUE ;
  #  x = nested_elastic_fit1 ; alpha = NULL ; gamma = 0.6 ; 
  #  x = nested_elastic_fit1 ; type="lasso" ; gamma = 1 ; 
  #  alpha = NULL ; gamma = NULL ; title = NULL 
  #  alpha = 0.2  ; gamma = NULL ; title = NULL
  #  alpha = NULL ; gamma = 0.5  ; title = NULL
  #  alpha = 0.2  ; gamma = 0.5  ; title = NULL
  gam = gamma 
  object = x 
  
  if (!(inherits(object,"nested.glmnetr"))) { 
    print( "    plot.glmnetr is only intended for nested.glmnetr() output objects.\n",
           "    Program will termiante. \n")
    stop
  } 
  
  if ( (type != "elastic") & (!is.null(alpha)) ) {
    cat(paste('  Becuase type != "elastic", non NULL value for alpha will be ignored'))
    alpha = NULL 
  }
  
  if ( (type == "elastic") & (is.null(x$cv_glmnet_fit_)) ) {
    cat(paste('  Elastic models not fit. type = "lasso" will be used instead.'))
    type = "lasso"
  }
  
  ## make sure alpha value is proper numeric or NULL 
  ## alpha == NULL is functional as alpha == "alpha.min"
  if ((!is.null(alpha))) {
    if (is.character(alpha)) {
      if (alpha == "alpha.min") { 
        cat(paste( 'alpha should be either NULL, "alpha.min" or numeric\n' ))
        cat(paste( 'alpha is changed to "alpha.min" \n' ))
        alpha = NULL
      }
    } else if (is.null(object$alpha)) {
      alpha = NULL                                                              ## redundant to if ( (type != "elastic") & (!is.null(alpha)) ) alpha = NULL 
    } else if (!(alpha %in% round(object$alpha,digits=6))) { 
      cat(paste('  The numeric value for alpha was not used when fitting the elastec net model.\n'))
      cat(paste('    alpha is changed to "alpha.min" \n' ))
      alpha = NULL 
    }
  }
  
  if ( (type != "ridge") & (!is.null(gam)) ) {
    if (inherits(object,"nested.glmnetr")) {
      if (!(gam %in% round(object$cv_lasso_fit$relaxed$gamma,digits=6)) ) { 
        cat(paste('  The numeric value for gamma was not used when fitting the elastec net model.\n'))
        cat(paste('    gam is changed to "gamma.min" \n' ))
        gam = NULL 
      }
    }
  }

  dolasso = x$fits[1]
  if (dolasso==1) { 
    if (type == "ridge") {
        object = x$cv_ridge_fit
        lambda_v = object$lambda
        lambda.min = object$lambda.min
        beta0 = NULL 
        beta1  = object$glmnet.fit$beta
        beta  = beta1 
        nzero = object$nzero
#        dim(beta) ; length(nzero) ; length(lambda_v)
        title1 = paste0("Ridge fit at log(gamma.min) = ", round(log(lambda.min),digits=4) )
        cat(paste0( "  For ridge fit log(lambda.min) = ", round(log(lambda.min),digits=4) , "\n")) 
    } else if (type == "lasso"  ) { 
      if (is.null(gam)) {
        object = x$cv_lasso_fit
        gam_ = object$relaxed$gamma.min
        lambda_v = object$lambda
        lambda.min = object$relaxed$lambda.min
        beta0   = object$glmnet.fit$relaxed$beta # 
        beta1   = object$glmnet.fit$beta # [,lambda_i]
        beta = ( (1-gam_) * beta0 + gam_ * beta1 )
        nzero = object$nzero
        title1 = paste0("Lasso fit at gamma.min = ", gam_ )
        cat(paste0( "  For lasso fit at gamma.min = ", gam_, ",\n",
                    "  log(lambda.min) = ", round(log(lambda.min),digits=4) , "\n")) 
      } else { 
        object = x$cv_lasso_fit
        lgamma = length(object$relaxed$gamma )
        gam_i = c(1:lgamma)[(gam == object$relaxed$gamma)]
        gam_  = gam 
        lambda_v = object$lambda
        
        names(object$relaxed$statlist)
        statlist = object$relaxed$statlist[[gam_i]] 
        lambda_i = which.min(statlist$cvm)
        lambda.min = lambda_v[lambda_i]
        
        beta0   = object$glmnet.fit$relaxed$beta # 
        beta1   = object$glmnet.fit$beta # [,lambda_i]
        beta = ( (1-gam_) * beta0 + gam_ * beta1 )
        nzero = object$nzero
        title1 = paste0("Relaxed lasso fit at gamma = ", gam_ )
        cat(paste0( "  For relaxed lasso fit at gamma = ", gam_, ",\n",
                    "  log(lambda.min) = ", round(log(lambda.min),digits=4) , "\n"))
      }
    } else if (type == "elastic") { 
      if ((is.null(alpha)) & (is.null(gam))) {
        vals.elastic = x$cv_elastic_fit$vals.elastic
        object  = x$cv_elastic_fit 
        alpha_i = object$index.elastic[1]
        vals.elastic = x$cv_elastic_fit$vals.elastic
        gam_ = vals.elastic[2]
        lambda_v = object$lambda
        lambda.min = object$relaxed$lambda.min
        beta0   = object$glmnet.fit$relaxed$beta # 
        beta1   = object$glmnet.fit$beta # [,lambda_i]
        beta = ( (1-gam) * beta0 + gam * beta1 )
        nzero = object$nzero
        title1 = paste0("Elastic net fit where alpha.min=", vals.elastic[1], " and gamma.min = ", gam_ )
        cat(paste0( "  For Elastic net fit where alpha.min = ", vals.elastic[1], " and gamma.min = ", vals.elastic[2], ",\n",
                    "  log(lambda.min) = ", round(log(vals.elastic[3]),digits=4) , "\n"))
      } else if (is.null(gam)) {
        alpha_i = c(1:length(x$alpha))[(alpha == round(x$alpha,digits=6))]    ## tests did not work without round() 
        object  = x$cv_glmnet_fit_[[alpha_i]]  
#        alpha_i = object$index.elastic[1]
        gam_ = object$relaxed$gamma.min
        lambda_v = object$lambda
        lambda.min = object$relaxed$lambda.min
        beta0   = object$glmnet.fit$relaxed$beta # 
        beta1   = object$glmnet.fit$beta # [,lambda_i]
        beta = ( (1-gam_) * beta0 + gam_ * beta1 )
        nzero = object$nzero
        title1 = paste0("Elastic net fit at alpha = ", alpha, ", and gamma.min = ", gam_ )
        cat(paste0( "  For elastic net fit at alpha = ", alpha, "\n   gamma.min = ", gam_, 
                    " and log(lambda.min) = ", round(log(lambda.min),digits=4) , "\n"))
      } else if (is.null(alpha)) {
        object = x$cv_glmnet_fit_
        class(object[[1]])
        names(object[[1]])
        gam_ = gam 
        lgamma = length( x$gamma ) 
        gamma_i = c(1:lgamma)[(gam == round(x$gamma,digits=6))]
        lalpha = length(x$alpha)
        for ( i_ in c(1:lalpha)) {
          names( object[[i_]] )
          names( object[[i_]]$glmnet.fit ) 
          names( object[[i_]]$relaxed$statlist[[gamma_i]] ) 
          mincvm = min( object[[i_]]$relaxed$statlist[[gamma_i]]$cvm )
          #          print(mincvm)
          lambda_it = which.min( object[[i_]]$relaxed$statlist[[gamma_i]]$cvm )
          lambda = object[[i_]]$relaxed$statlist[[gamma_i]]$lambda[lambda_it]
          log(lambda)
          if (i_ == 1) {
            measure.min = mincvm 
            alpha_i = i_ 
            lambda_i = lambda_it
          } else if (  mincvm < measure.min ) {
            measure.min = mincvm 
            alpha_i = i_             
            lambda_i = lambda_it
          }
        }
        c(alpha_i, x$alpha[alpha_i], lambda_i)
        alpha_ = x$alpha[alpha_i] 
        object2 = object[[alpha_i]]
        lambda.min = object2$lambda[lambda_i]
        lambda_v = object2$lambda 
        lambda_n = length( lambda_v ) 
        nzero = object2$nzero 
        names(object2)
        names(object2$glmnet.fit$relaxed)
        beta0 = object2$glmnet.fit$relaxed$beta # 
        beta1 = object2$glmnet.fit$beta # [,lambda_i]
        beta = ( (1-gam) * beta0 + gam * beta1 )
        title1 = paste0("Elastic net fit at gamma=", gam, " where alpha.min=", alpha_)
        cat(paste0( "  For elastic net fit at gamma = ", gam, "\n   alpha.min=", alpha_,
                    " and log(lambda.min) = ", round(log(lambda.min),digits=4) , "\n"))
      } else {
        # gam = 0.5 ; alpha = 0.4 ; 
        gam_ = gam 
        gamma_v = round(x$gamma, digits=6)
        lgamma = length(gamma_v) 
        object = x$cv_glmnet_fit_
        alpha_i = c(1:length(x$alpha))[(alpha == round(x$alpha,digits=6))]
        gamma_i = c(1:lgamma)[(gam == gamma_v)]
        nzero = object[[alpha_i]]$nzero
        object2 = x$cv_glmnet_fit_[[ alpha_i ]]
        lambda_v = object2$glmnet.fit$lambda
#        which.min( nested_elastic_fit0$cv_glmnet_fit_[[1]]$cvm )
        lambda_i = which.min( object[[1]]$cvm )
        lambda.min = lambda_v[ lambda_i ]
        beta0   = object2$glmnet.fit$relaxed$beta # 
        beta1   = object2$glmnet.fit$beta # [,lambda_i]
        beta = ( (1-gam) * beta0 + gam * beta1 )
        title1 = paste0("Elastic net fit at alpha=", alpha, ", and gamma = ", gam)
        cat(paste0( "  For elastic net fit at alpha=", alpha, ", and gamma = ", gam, ",\n",
                    "  log(lambda.min) = ", round(log(lambda.min),digits=4) , "\n"))
      }
    }
 
    if (!is.null(beta0)) { 
      beta = gam_*beta1 + (1-gam_)*beta0
    } else { 
      beta = beta1
    } 

    if (!is.null(lambda.lo)) {
      beta = as.matrix(beta)[ , (lambda_v > lambda.lo) ]
      lambda_v = lambda_v[ lambda_v > lambda.lo ]
    }
    
    loglambda_v = log(lambda_v)
    
    maxyterm = ceiling(10*max(beta))/10
    minyterm = floor  (10*min(beta))/10
    
    maxxterm = ceiling(10*max(loglambda_v))/10
    minxterm = floor(10*min(loglambda_v))/10
    
    if (is.null(title)) { title = title1 }      
    plot(x=max(loglambda_v),0, xlim=c(minxterm, maxxterm), ylim=c(minyterm, maxyterm) , 
         main=title, xlab="Log Lambda", ylab="Coefficients" )
#     axis(at = log(lambda.min), tick=FALSE, line = -0.5, side = 3 , cex.axis=1) 
#     axis(2, labels=nzero) 
    if (!is.null(lambda.min)) { abline(a=NULL, b=NULL, v=log(lambda.min), h=NULL) }
    for (i_ in 1:dim(beta)[1]) {
      lines(loglambda_v, beta[i_,], pch=16, col=i_)
    }
  }
}

############################################################################################################
############################################################################################################
##### plot relaxed CV deviances ############################################################################
#' Plot cross-validation deviances, or model coefficients.   
#' 
#' @description 
#' By default, with coefs=FALSE, plots the average deviances as function of 
#' alpha, gamma and lambda and also indicates the alpha, gamma and lambda 
#' which minimize deviance based upon a cv.glmnetr() output object.
#' Optionally, with coefs=TRUE, plots the relaxed lasso coefficients.  
#'
#' @param x a cv.glmnetr()  output object.  
#' @param type one of c("lasso", "elastic", "ridge") to plot the deviance 
#' curves of the respective model fit. Default is "lasso" for tuned relaxed lasso. 
#' @param alpha A specific value of alpha for plotting.  Used only when type is 
#' set to "elastic". Specifies which alpha is to be used for deviance plots.
#' Default is "alpha.min", else must be an element of the alpha vector used in 
#' running the elastic net model. This can be reviewed using summary(fit) 
#' where fit is a nested.glmnetr() output object. Note, alpha is 1 for the 
#' lasso model and alpha is 0 for the ridge model.
#' @param gamma a specific level of gamma for plotting.  By default gamma.min will be used.  
#' @param lambda.lo a lower limit of lambda when plotting.  
#' @param plup   an indicator to plot the upper 95 percent two-sided confidence limits.  
#' @param title  a title for the plot.  
#' @param coefs  default of FALSE plots deviances, option of TRUE plots coefficients.  
#' @param comment default of TRUE to write to console information on lam and gamma selected for output.
#' FALSE will suppress this write to console. 
#' @param lty line type for the cross validated deviances. Default is 1.
#' @param track 2 to track progress by printing to console, 0 (default) to 
#' not track.
#' @param ... Additional arguments passed to the plot function.  
#'
#' @return This program returns a plot to the graphics window, and may provide 
#' some numerical information to the R Console.  If gamma is not specified, then 
#' then the gamma.min from the deviance minimizing (lambda.min, gamma.min) pair 
#' will be used, and the corresponding lambda.min will be indicated by a vertical
#' line, and the lambda minimizing deviance under the restricted set of models 
#' where gamma=0 will be indicated by a second vertical line.    
#' 
#' @seealso
#'   \code{\link{plot.glmnetr}} , \code{\link{plot.nested.glmnetr}} , \code{\link{cv.glmnetr}} 
#'   
#' @export 
#' 
#' @importFrom graphics abline axis lines
#'
#'@noRd
#'
plot.cv.glmnetr_0_6_1 = function(x, type="lasso", alpha=NULL, gamma=NULL, lambda.lo=NULL, plup=0, title=NULL, coefs=FALSE, comment=TRUE, lty=1, track=0, ...) {
  #  x = nested_elastic_fit1 ; type="lasso"   ; gamma=NULL ; alpha=NULL ; lambda.lo=NULL ; plup=0 ; title=NULL ; coefs=FALSE ; comment=TRUE ; lty=1 ;
  #  x = nested_elastic_fit1 ; type="lasso"   ; gamma=0.5  ; alpha=NULL ; lambda.lo=NULL ; plup=0 ; title=NULL ; coefs=FALSE ; comment=TRUE ; lty=1 ;
  #  x = nested_elastic_fit1 ; type="lasso"   ; gamma=0    ; alpha=NULL ; lambda.lo=NULL ; plup=0 ; title=NULL ; coefs=FALSE ; comment=TRUE ; lty=1 ;
  #  x = nested_elastic_fit1 ; type="elastic" ; gamma=1    ; alpha=NULL ; lambda.lo=NULL ; plup=0 ; title=NULL ; coefs=FALSE ; comment=TRUE ; lty=1 ;
  #  type = "ridge"
  #  cat(paste("\n plot.cv. 1 type = ", type, ", gamma = ", gamma, ", alpha = ", alpha, ", lambda.lo = ", lambda.lo, ", title = ", title, ", comment = ", comment,", coefs = ", coefs,  "\n"))
  #    cat(paste("\n  type = ", type, "; alpha = ", alpha, "; gamma = ", gamma, "; lambda.lo = ", lambda.lo, "; coefs = ", coefs, "; title = ", title, "; comment = ", comment, " ; track =", track,  ";\n"))
    if (track >= 2) { 
      cat( "  in plot.cv.glmnetr_0_6_1   class(x) = ", class(x), "\n")
    }
#  x = nested_elastic_fit1 ;  type =  "lasso" ; alpha = NULL  ; gamma =  0.25 ; lambda.lo = NULL  ; plup = 0 ;  coefs =  FALSE ; title = NULL  ; comment =  TRUE  ; track = 3 ; title = NULL ; lty = 1 ; 
  object = x 
#  gam = gamma 
  
  if (!(inherits(object,"nested.glmnetr"))) {
    cat("  plot.cv.glmnetr_0_6_1 expects a nested.glmnetr output object.\n. Function will stop.")
  }
  
  if (!(type %in% c("ridge", "lasso", "elastic"))) {
    cat( '  type not in c("ridge", "lasso", "elastic"), is set to "lasso" ')
    type = "lasso"
  }
  
  if (coefs == 1) {
    plot.glmnetr_0_6_1(x, type=type , gamma=gamma, alpha=alpha, lambda.lo=lambda.lo, title=title, comment=comment, ...)
  } else if ( x$fits[1] == 1 ) {
    
    if ( (type != "elastic") & (!is.null(alpha)) ) {
      cat(paste('  Becuase type != "elastic", non NULL value for alpha will be ignored'))
      alpha = NULL 
    }
    
    if ( (type == "elastic") & (is.null(x$cv_glmnet_fit_)) ) {
      cat(paste('  Elastic models not fit. type = "lasso" will be used instead.\n\n'))
      type = "lasso"
    }
    
    ## make sure alpha value is proper numeric or NULL 
    ## alpha == NULL is functionl as alpha == "alpha.min"
    if ((!is.null(alpha))) {
      if (is.character(alpha)) {
        if (alpha == "alpha.min") { 
          cat(paste( 'alpha should be either NULL, "alpha.min" or numeric\n' ))
          cat(paste( 'alpha is changed to "alpha.min" \n' ))
          alpha = NULL
        }
      } else if (!(alpha %in% round(object$alpha,digits=6))) {
        cat(paste('  The numeric value for alpha was not used when fitting the elastic net model.\n'))
        cat(paste('    alpha is changed to "alpha.min" \n' ))
        alpha = NULL 
      }
    }
    
    if ( (type != "ridge") & (!is.null(gamma)) ) {
      if (!(gamma %in% round(object$cv_lasso_fit$relaxed$gamma,digits=6)) ) { 
        cat(paste('  The numeric value for gamma was not used when fitting the elastic net model.\n'))
        cat(paste('    gamma is changed to "gamma.min" \n' ))
        gamma = NULL 
      }
    }
    
    ## for elastic 
    if (!is.null(alpha)) { alpha_  = alpha 
    } else { alpha_  = object$cv_elastic_fit$vals.elastic[1] } 
    vals.elastic = x$cv_elastic_fit$vals.elastic                                ## check ?? 
    
    if        (type == "lasso"  ) { 
      object = object$cv_lasso_fit
    } else if (type == "ridge"  ) { 
      object = object$cv_ridge_fit 
      gamma = 1 
    } else if (type == "elastic") { 
      if ( (is.null(alpha)) & (is.null(gamma)) ) { 
        alpha_i = object$cv_elastic_fit$index.elastic[1]
        alpha_  = object$cv_elastic_fit$vals.elastic[1]
        object  = object$cv_elastic_fit
        #          alpha_  = object$vals.elastic[1] 
      } else if ( (!is.null(alpha)) & (!is.null(gamma)) ) {
        alpha_i = c(1:length(x$alpha))[(alpha == round(x$alpha,digits=6))]    ## tests did not work without round() 
        object = object$cv_glmnet_fit_[[ alpha_i ]]
        gamma_v = object$relaxed$gamma 
        lgamma = length(gamma_v) 
        gamma_i = c(1:lgamma)[ gamma ==round(gamma_v,digits=4)]
        class(object) = "cv.glmnetr.el" 
        #          alpha_ = alpha 
      } else if (!is.null(alpha)) {
        alpha_i = c(1:length(x$alpha))[(alpha == round(x$alpha,digits=6))]    ## tests did not work without round() 
        object = object$cv_glmnet_fit_[[ alpha_i ]]
        class(object) = "cv.glmnetr.el" 
        #          alpha_ = alpha 
      } 
    }
    
    if ( (type=="ridge") | 
         ((type=="lasso"  ) & (!is.null(gamma)) ) | 
         ((type=="elastic") & (!is.null(alpha)) & (!is.null(gamma))) ) {
      #==== single deviance curve with upper and lower "95%" band ==============
      if (type=="ridge") { 
        statlist1 = list()
        statlist1$lambda = object$lambda
        statlist1$cvm    = object$cvm
        statlist1$cvsd   = object$cvsd 
        statlist1$cvup   = object$cvup
        statlist1$cvlo   = object$cvlo
        statlist1$nzero  = object$nzero
      } else if ( (type=="lasso"  ) & (!is.null(gamma)) ) {
        gamma_i = c(1:length(object$relaxed$gamma))[gamma==round(object$relaxed$gamma,digits=4)]
        statlist1 = object$relaxed$statlist[[gamma_i]]
      } else if ((type=="elastic") & (!is.null(alpha)) & (!is.null(gamma))) {
        statlist1 = object$relaxed$statlist[[gamma_i]] 
      }
      
      #        lines(log(statlist[[k_]][[1]]), statlist[[k_]][[j_]], pch=16, lty=lty, col=k_)
      
      deviance.min = min( statlist1$cvm )
      lambda.min = statlist1$lambda[ which.min( statlist1$cvm ) ]
      df.lambda.min = statlist1$nzero[ which.min( statlist1$cvm ) ]
      
      if (!is.null(lambda.lo)) {
        if (lambda.min < lambda.lo) { lambda.lo = lambda.min }
        index = ( statlist1$lambda >= lambda.lo )
        statlist1$lambda = statlist1$lambda[index]
        statlist1$cvm    = statlist1$cvm[index]
        statlist1$cvsd   = statlist1$cvsd[index]
        statlist1$cvup   = statlist1$cvup[index]
        statlist1$cvlo   = statlist1$cvlo[index]
        statlist1$nzero  = statlist1$nzero[index] 
      }
      
      if ( (type=="elastic") & (!is.null(alpha)) & (!is.null(gamma)) ) {
        title1 = paste0("Elastic net fit at alpha=",  alpha, " and gamma=", gamma )
      } else if ((type == "lasso") & (gamma == 0)) {
        title1 = paste0("Fully Relaxed lasso fit at gamma ", gamma ) 
      } else if ((type == "lasso") & (gamma == 1)) {
        title1 = paste0("Fully penalized lasso fit at gamma ", gamma ) 
      } else if ( (type == "lasso") & (!is.null(gamma)) ) {
        title1 = paste0("Relaxed lasso fit at assinged gamma ", gamma ) 
      } else if (type == "ridge") {
        title1 = paste0("Ridge fit")  
      }
      
      if (is.null(title)) {
        title = title1 
      }
      
      loglambda = log( statlist1$lambda )
      maxxterm  = ceiling(10*max(loglambda))/10
      minxterm  = floor(10*min(loglambda))/10
      
      lambda_n = length(statlist1$lambda)
      
      maxcvm = max ( statlist1$cvm )
      mincvm = min ( statlist1$cvm )
      maxyterm = ceiling(1000*maxcvm)/1000
      minyterm = floor(1000*mincvm)/1000
      
      if (x$sample[1]=="cox") { ylab="deviance/event" } else { ylab="deviance/record" }
      plot(log(statlist1$lambda), statlist1$cvm, type="l", lty=1, xlim=c(minxterm, maxxterm) , ylim=c(minyterm, maxyterm), 
           main=NULL, xlab="log(lambda)", ylab=ylab ) 
      axis(at = loglambda[1:length(statlist1$lambda)], labels = statlist1$nzero, tick=FALSE, line = -0.5, side = 3 , cex.axis=1)
      abline(a=NULL, b=NULL, v=log(statlist1$lambda[which.min(statlist1$cvm)]), h=NULL)
      
      title(main=title)
      lines(log(statlist1$lambda), statlist1$cvup, lty=2, col=4)
      lines(log(statlist1$lambda), statlist1$cvlo, lty=2, col=4)
      
      if (comment==TRUE) { 
        if ((type == "lasso") & (gamma == 0)) {
          pretext = paste0(" Fully relaxed lasso (gamma=0) at log(lambda) = "      , (ceiling(log(lambda.min)*1000)/1000), ", df = " , df.lambda.min, ", deviance = ", round(deviance.min, digits=4) , "\n")
          cat(pretext)
        } else if ( (type == "lasso") & (!is.null(gamma)) ) {
          pretext = paste0(" Relaxed lasso for assinged gamma = ", gamma, ", at log(lambda) = ", ceiling(log(lambda.min)*1000)/1000, ", df = " , df.lambda.min, ", deviance = ", round(deviance.min, digits=4) , "\n")
          cat(pretext)
        } else if (type == "elastic") {
          if ( (!is.null(alpha)) & (!is.null(gamma)) ) { 
            pretext = paste0( "Elastic Net at alpha = ", alpha_ , " and gamma = ", gamma,  " tuned for lambda minimizing CV average deviance (maximizing log likelihood)\n",
                              "   log(lambda) = ", ceiling(log(lambda.min)*1000)/1000, ", df = " , df.lambda.min, ", deviance = ", round(deviance.min, digits=4), "\n" )
            cat( pretext )
          } 
        } else if (type == "ridge") {
          pretext = paste0( " Ridge tuned for lambda", 
                            ", to minimize CV average deviance (maximizing log likelihood)\n",  "   lambda.min= ", ceiling(log(lambda.min)*1000)/1000, ", df = " , df.lambda.min, ", deviance = ", round(deviance.min, digits=4), "\n" ) 
          cat(pretext) 
        }
      }
      
    } else { 
      #===== MULTIPLE deviance curves ==========================================
      
      plotrv = function(statlist, ktop, j_, lty=1) {  ## j_=2 cvm, 4 cvup, 5 cvlo , k = gamma_n or alpha_n
        for (k_ in 1:ktop) { 
          lines(log(statlist[[k_]]$lambda), statlist[[k_]][[j_]], pch=16, lty=lty, col=k_)
        }
      }
      
      if ( (type == "elastic") & (is.null(alpha)) & (!is.null(gamma)) ) {
        #===== MULTIPLE deviance curves -- Elastic with non null gamma ===========
        
        alpha_v = object$alpha 
        gamma_v = object$gamma 
        lgamma  = length(gamma_v) 
        gamma_i = c(1:lgamma)[gamma==round(gamma_v,digits=4)]
        # gamma_i = c(1:length(object$relaxed$gamma))[gamma==round(object$relaxed$gamma,digits=4)] 
        statlist = list()
        for ( i_ in c(1:length(object$alpha))) {
          lambda_ = object$cv_glmnet_fit_[[i_]]$relaxed$statlist[[gamma_i]]$lambda
          cvm_    = object$cv_glmnet_fit_[[i_]]$relaxed$statlist[[gamma_i]]$cvm
          nzero_  = object$cv_glmnet_fit_[[i_]]$relaxed$statlist[[gamma_i]]$nzero
          alpha.measure.min_ = min( cvm_ )  
          lambda.min_        = lambda_[ which.min(cvm_) ]
          df.lambda.min_     = nzero_[ which.min(cvm_) ]
          
          if (i_ == 1) {
            alpha.measure.min_i = 1
            alpha.measure.min = alpha.measure.min_ 
            lambda.min = lambda.min_
            df.lambda.min = df.lambda.min_ 
#            statlist_t = object$cv_glmnet_fit_[[i_]]$relaxed$statlist[[gamma_i]]
          } else if (alpha.measure.min_ < alpha.measure.min) {
            alpha.measure.min_i = i_ 
            alpha.measure.min = alpha.measure.min_ 
            lambda.min = lambda.min_
            df.lambda.min = df.lambda.min_ 
          }
          statlist[[i_]] = object$cv_glmnet_fit_[[i_]]$relaxed$statlist[[gamma_i]]
        }
        
        for ( i_ in 1:length(object$alpha) ) {
          maxlambda_ = max( statlist[[i_]]$lambda ) 
          minlambda_ = min( statlist[[i_]]$lambda ) 
          maxcvm_ = max( statlist[[i_]]$cvm ) 
          mincvm_ = min( statlist[[i_]]$cvm ) 
          if (i_ == 1) {
            maxcvm = maxcvm_ 
            mincvm = mincvm_ 
            maxlambda = maxlambda_ 
            minlambda = minlambda_ 
          }
          if ( maxcvm_ > maxcvm ) { maxcvm = maxcvm_ } 
          if ( mincvm_ < mincvm ) { mincvm = mincvm_ } 
          if ( maxlambda_ > maxlambda ) { maxlambda = maxlambda_ } 
          if ( minlambda_ < minlambda ) { minlambda = minlambda_ } 
        }
        
        maxyterm = ceiling(1000*maxcvm)/1000
        minyterm = floor(1000*mincvm)/1000
        
        maxxterm = ceiling(10*log(maxlambda))/10
        minxterm = floor(  10*log(minlambda))/10
        
        if (x$sample[1]=="cox") { ylab="deviance/event" } else { ylab="deviance/record" }
        plot(log(statlist[[1]]$lambda[1]), statlist[[1]]$cvm[1], xlim=c(minxterm, maxxterm) , ylim=c(minyterm, maxyterm), 
             main=NULL, xlab="log(lambda)", ylab=ylab ) 
        abline(a=NULL, b=NULL, v=log(lambda.min), h=NULL)
        
        plotrv(statlist, length(alpha_v), 2)
        if (plup) { plotrv(statlist, length(alpha_v), 4, 3) } 
        
        #  for (i_ in 1:length(alpha_v)) {
        #    lines(log(statlist[[i_]]$lambda), statlist[[i_]]$cvm, pch=16, lty=lty, col=i_)
        #  }
        
        title1 = paste0("Elastic net fit at gamma = ", gamma, " where alpha.min = ", alpha_v[alpha.measure.min_i])
        title(main=title1)
        
        pretext = paste0( "Elastic Net at gamma = ", gamma, " tuned for alpha and lambda minimizing CV average deviance (maximizing log likelihood)\n",
                          "   alpha.min = ", alpha_v[alpha.measure.min_i], ", log(lambda) = ", ceiling(log(lambda.min)*1000)/1000, ", df = " , df.lambda.min, ", deviance = ", round(mincvm, digits=4), "\n" )
        cat( pretext )
      } else { 
        #===== MULTIPLE deviance curves other than elastic with non null gamma =======
        
        #        nzero              = object$nzero 
        gamma_v            = object$relaxed$gamma 
        
        lambda.min.g1      = object$lambda.min 
        df.lambda.min.g1   = object$nzero[ object$index[1] ] 
        df.lambda.min.g1   = object$nzero.min
        deviance.g1        = object$cvm[ object$index[1] ]
        
        index.r            = object$relaxed$index[1,1]
        lambda.min         = object$relaxed$lambda.min
        gamma.min          = object$relaxed$gamma.min
        df.lambda.min      = object$relaxed$nzero.min                           ## get beta and see that beta > betatol of 1e-14
        df.lambda.min      = object$nzero[ index.r ]                                   ## CHECK ?? 
        index.lam.r        = object$relaxed$index[1,1]   
        index.gam.r        = object$relaxed$index[1,2]   
        deviance.r         = object$relaxed$statlist[[ index.gam.r ]]$cvm[ index.lam.r ]
        
        index.g0           = object$relaxed$index.min.g0
        lambda.min.g0      = object$relaxed$lambda.min.g0 
        df.lambda.min.g0   = object$relaxed$nzero.min.g0                        ## CHECK ?? 
        deviance.g0        = object$relaxed$statlist[[1]]$cvm[ index.g0 ]
        
        if        ( (type == "elastic") & (is.null(alpha)) & (is.null(gamma)) ) {
          title1 = paste0("Elastic net fit where alpha.min=", object$vals.elastic[1], ", gamma.min=", object$vals.elastic[2]) 
#          , ", log(l)=", round(log(lambda.min), digits=3)
        } else if ( (type == "elastic") & (is.null(gamma))   ) {
          indx = object$relaxed$index
          title1 = paste0("Elastic net fit at alpha=", alpha, " where gamma.min=", object$relaxed$gamma.min ) 
        } else if ( (type == "lasso") & (is.null(gamma)) ) {
          title1 = paste0("Relaxed lasso fit at gamma.min = ", object$relaxed$gamma.min ) 
        } 
        
        if (is.null(title)) {
          title = title1 
        }
        
        if (comment==TRUE) { 
          if (type == "lasso") {
            pretext = paste0( " Relaxed lasso minimizing CV average deviance (maximizing log likelihood)\n", 
                              "     at gamma.min = ", gamma.min, ", log(lambda) = ", ceiling(log(lambda.min   )*1000)/1000, ", df = " , df.lambda.min, ", deviance = ", round(deviance.r, digits=4), "\n",  
                              "   fully relaxed (gamma=0) at log(lambda) = ", ceiling(log(lambda.min.g0)*1000)/1000, ", df = " , df.lambda.min.g0, ", deviance = ", round(deviance.g0, digits=4) , "\n",
                              "   fully penalized (gamma=1) at log(lambda) = ", ceiling(log(lambda.min.g1)*1000)/1000, ", df = " , df.lambda.min.g1, ", deviance = ", round(deviance.g1, digits=4) , "\n" )
            cat(pretext)
          } else if (type == "elastic") {
            if ( (is.null(alpha)) & (is.null(gamma)) ) { 
              pretext = paste0( "Elastic Net tuned for alpha, gamma and lambda minimizing CV average deviance (maximizing log likelihood)\n",
                                "   alpha.min= ", alpha_ , ", gamma.min = ", gamma.min, ", log(lambda) = ", ceiling(log(lambda.min   )*1000)/1000, ", df = " , df.lambda.min, ", deviance = ", round(deviance.r, digits=4), "\n" )
              cat( pretext )
            } else if (is.null(gamma)) { 
              pretext = paste0( "Elastic Net at alpha = ", alpha, " tuned for gamma and lambda minimizing CV average deviance (maximizing log likelihood)\n",
                                "   gamma.min = ", gamma.min, ", log(lambda) = ", ceiling(log(lambda.min)*1000)/1000, ", df = " , df.lambda.min, ", deviance = ", round(deviance.r, digits=4), "\n" )
              cat( pretext )
            } 
          } else if (type == "ridge") {
            pretext = paste0( " Ridge tuned for lambda", 
                              ", to minimize CV average deviance (maximizing log likelihood)\n",  "   lambda.min= ", ceiling(log(lambda.min.g1)*1000)/1000, ", df = " , df.lambda.min.g1, ", deviance = ", round(deviance.g1, digits=4), "\n" ) 
            cat(pretext) 
          }
        }
        
        statlist=object$relaxed$statlist 
        #  statlist[[5]]$cvm[1:8]
        #  index
        
        gamma_n = length(statlist)
        
        if (!is.null(lambda.lo)) {                                              ## check how done elsewhere CHECK 
          index = statlist[[1]]$lambda >= lambda.lo 
          nzero     = object$nzero[index] 
          for (k_ in 1:gamma_n) {
            statlist[[k_]]$lambda = statlist[[k_]]$lambda[index]                ## lambda 
            statlist[[k_]]$cvm = statlist[[k_]]$cvm[index]                      ## cvm 
            statlist[[k_]]$cvup = statlist[[k_]]$cvup[index]                    ## cvup 
            statlist[[k_]]$cvlo = statlist[[k_]]$cvlo[index]                    ## nzero 
          }
        }
        
        loglambda = log(statlist[[1]]$lambda)
        
        maxxterm = ceiling(10*max(loglambda))/10
        minxterm = floor(10*min(loglambda))/10
        
        lambda_n = length(statlist[[1]]$lambda)
        cvmg0 = statlist[[1]]$cvm ; cvmg0 ; 
        cvmg0 = cvmg0[1:(lambda_n-1)] ; cvmg0 ;                                 ## CHECK this
        cvmg1 = statlist[[gamma_n]]$cvm
        cvmg1 = cvmg1[1:(lambda_n-1)]
        cvms = c()
        for (k_ in 1:gamma_n) {
          cvms = cbind(cvms,statlist[[k_]]$cvm)
        }
        maxcvm = max ( cvms )
        mincvm = min ( cvms )
        
        maxyterm = ceiling(1000*maxcvm)/1000
        minyterm = floor(1000*mincvm)/1000
        
        # llmin = log(lambda.min)
        if (x$sample[1]=="cox") { ylab="deviance/event" } else { ylab="deviance/record" }
        plot(log(statlist[[1]]$lambda[1]), statlist[[1]]$cvm[1], xlim=c(minxterm, maxxterm) , ylim=c(minyterm, maxyterm), 
             main=NULL, xlab="log(lambda)", ylab=ylab ) 
        #plot(log(statlist[[1]]$lambda[1]), statlist[[1]]$cvm[1], xlim=c(minxterm, maxxterm) , ylim=c(minyterm, maxyterm), 
        #     xlab="", type="o", xaxt = "n", ylab="log likelihood / subject" ) 
        # axis(1, at = 1:lambda_n, labels = loglambda, side = 1)
        #axis(1, at = 1:lambda_n, labels = object$nzero, line = 2.5, side = 3 )
        axis(at = loglambda[1:lambda_n], labels = object$nzero, tick=FALSE, line = -0.5, side = 3 , cex.axis=1)
        abline(a=NULL, b=NULL, v=log(lambda.min), h=NULL)
        abline(a=NULL, b=NULL, v=log(lambda.min.g0), h=NULL, lty=2 )
        abline(a=NULL, b=NULL, v=log(lambda.min.g1), h=NULL, lty=2 )
        #title(main="Relaxed lasso CV", sub="sub-title", xlab="log(lambda)", ylab="log likelihood / subject")
        #     xlab="log(lambda)", ylab="log likelihood / subject" )
        object$val
        
        title(main=title)
        
        plotr = function(statlist, j_, lty=1) {  ## j_=2 cvm, 4 cvup, 5 cvlo 
          for (k_ in 1:gamma_n) {
            lines(log(statlist[[k_]]$lambda), statlist[[k_]][[j_]], pch=16, lty=lty, col=k_)
          }
        }

#        plotr(statlist,2)
#        if (plup) { plotr(statlist,4,3) } 
        
        plotrv(statlist, gamma_n, 2)
        if (plup) { plotrv(statlist, gamma_n, 4, 3) } 
      }
    }
  } 
}

####################################################################################################
####################################################################################################