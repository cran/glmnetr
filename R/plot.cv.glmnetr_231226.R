############################################################################################################
############################################################################################################
#' Plot the relaxed lasso coefficients.  
#' 
#' @description 
#' Plot the relaxed lasso coefficients from either a glmnetr(), cv.glmnetr() or nested.glmnetr() output object.
#' One may specify gam, single value for gamma.  If gam is unspecified (NULL), then cv.glmnetr and 
#' nested.glmnetr() will use the gam which minimizes loss, and glmentr() will use gam=1.
#'
#' @param x Either a glmnetr, cv.glmnetr or a nested.glmnetr output object.  
#' @param gam A specific level of gamma for plotting.  By default gamma.min from
#' the deviance minimizing (lambda.min, gamma.min) pair will be used.  
#' @param lambda.lo A lower limit of lambda for plotting.  
#' @param title A title for the plot  
#' @param comment Default of TRUE to write to console information on lam and gam selected for output.
#' FALSE will suppress this write to console.  
#' @param ... Additional arguments passed to the plot function.  
#'
#' @return This program returns a plot to the graphics window, and may provide 
#' some numerical information to the R Console.  If the input object is from a
#' nested.glmnetr or cv.glmnetr object, and gamma is not specified, 
#' then the gamma.min from
#' the deviance minimizing (lambda.min, gamma.min) pair will be used, and the
#' minimizing lambda.min will be indicated by a vertical line.  Also, if one 
#' specifies gam=0, the lambda which minimizes deviance for the restricted 
#' set of models where gamma=0 will indicated by a vertical line.  
#' 
#' @seealso 
#'   \code{\link{plot.cv.glmnetr}} , \code{\link{plot.nested.glmnetr}} , \code{\link{glmnetr}} 
#' 
#' @export
#' 
#' @importFrom graphics abline axis lines
#'
#' @examples
#' \donttest{
#' set.seed(82545037)
#' sim.data=glmnetr.simdata(nrows=200, ncols=100, beta=NULL)
#' xs=sim.data$xs 
#' y_=sim.data$yt
#' event=sim.data$event
#' glmnetr.fit = glmnetr( xs, NULL, y_, event, family="cox")
#' plot(glmnetr.fit)
#' }
#' 
plot.glmnetr = function(x, gam=NULL, lambda.lo=NULL, title=NULL,comment=TRUE, ...) {
  object = x 
  if (inherits(object,"nested.glmnetr")) {
    dolasso = object$fits[1]
    if (dolasso==1) { 
      object = object$cv_glmnet_fit 
    } else { 
      cat(paste0(" No relaxed lasso for plotting" , "\n")) 
      stop 
    }
  } 
  
  if        ( inherits(object,"cv.glmnetr") ) { 
    glmnet.fit = object$glmnet.fit ; 
    beta0  = glmnet.fit$relaxed$beta  ; dim(beta0)
    beta1  = glmnet.fit$beta          ; dim(beta1)
    nzero  = object$nzero
    cv_ = 1 ; 
  } else if ( inherits(object,"glmnetr") ) { 
    glmnet.fit = object            ; 
    beta0  = glmnet.fit$betag0  ; dim(beta0)
    beta1  = glmnet.fit$betag1  ; dim(beta1)
    nzero  = glmnet.fit$cv.glmnet.fit$nzero
    cv_ = 0 ; 
  } else { stop } 

  lambda = glmnet.fit$lambda
  lambda_n = length(lambda)
  
  if (!is.null(lambda.lo)) {
    beta0  = beta0[,c(lambda>=lambda.lo)]  
    beta1  = beta1[,c(lambda>=lambda.lo)]  
    lambda = lambda[c(lambda>=lambda.lo)]  
  }

  if (cv_) {
    lambda.min = object$relaxed$lambda.min
    gamma.min  = object$relaxed$gamma.min
    df.lambda.min = nzero[object$relaxed$index[1,1]]
    lambda.min.g0 = object$relaxed$lambda.min.g0
    df.lambda.min.g0 = nzero[object$relaxed$index.g0[1]]
    lambda.min.g1 = object$lambda.min 
    df.lambda.min.g1 = nzero[object$index[1]]
#  cat(paste0(" Number of model terms at lambda.min = " , df.lambda.min, "\n"))
  }
  
#  if (is.null(gam)) { 
#    if (cv_) { 
#      gam=gamma.min 
#      if (is.null(title)) {title = paste0("Relaxed lasso fit at gamma.min = ", gam , "\n") } 
#    } else { 
#      if (is.null(gam)) { gam=1 } 
#      if (is.null(title)) { title = paste0("lasso fit, gamma = ", gam , "\n") } 
#    }
#  } else { if (is.null(title)) { title = paste0("Relaxed lasso fit, gamma = ", gam , "\n") } }
  
  lambda.vl = NULL 
  if (cv_) { 
    if (is.null(gam)) { 
      gam=gamma.min 
      if (is.null(title)) {title = paste0("Relaxed lasso fit at minimizing gamma ", gam , "\n") } 
      lambda.vl = log(lambda.min) 
      if (comment==TRUE) { cat(paste0(" min CV average deviance (max log likelihood) ", "\n", 
                 "  at log(lambda.min) = ", round(log(lambda.min),digits=3), ", gamma.min = ", gamma.min, ", df = " , df.lambda.min, "\n")) }
    } else if (gam==0) { 
      if (is.null(title)) {title = paste0("Fully relaxed lasso fit for gamma = 0", "\n") }
      lambda.vl = log(lambda.min.g0) 
      if (comment==TRUE) { cat(paste0(" Fully relaxed min CV average deviance (max log likelihood) ", "\n", 
                 "  at log(lambda.min) = ", round(log(lambda.min.g0),digits=3), ", df = " , df.lambda.min.g0, "\n")) }
    } else if (gam==1) { 
      if (is.null(title)) {title = paste0("Fully penalized lasso fit for gamma = 1", "\n") }
      lambda.vl = log(lambda.min.g1) 
      if (comment==TRUE) { cat(paste0(" Fully penalized min CV average deviance (max log likelihood) ", "\n", 
                 "  at log(lambda) = ", round(log(lambda.min.g1),digits=3), ", df = " , df.lambda.min.g1, "\n")) }
    } else { 
      if (is.null(title)) {title = paste0("Relaxed lasso fit for gamma = ", gam, "\n") }
    }
  } else {
    if (is.null(gam)) { gam=1 }
    if (is.null(title)) { title = paste0("Relaxed lasso fit for gamma = ", gam , "\n") }
  }

  loglambda = log(lambda) 

  maxxterm = ceiling(10*max(loglambda))/10
  minxterm = floor(10*min(loglambda))/10

#  dim(beta)
#  length(lambda)
#  length(beta[1,])
#  length(loglambda)

  maxyterm = ceiling(10*max(cbind(beta0,beta1)))/10
  minyterm = floor(10*min(cbind(beta0,beta1)))/10

  beta = gam*beta1 + (1-gam)*beta0 

  maxyterm = ceiling(10*max(beta))/10
  minyterm = floor  (10*min(beta))/10

  plot(x=max(loglambda),0, xlim=c(minxterm, maxxterm), ylim=c(minyterm, maxyterm) , 
       main=title, xlab="Log Lambda", ylab="Coefficients" )
  axis(at = loglambda[1:lambda_n], labels = nzero, tick=FALSE, line = -0.5, side = 3 , cex.axis=1)
#  axis(2, labels=nzero)
  if (!is.null(lambda.vl)) { abline(a=NULL, b=NULL, v=lambda.vl, h=NULL) }
#   plot(loglambda, beta[1,], ylim=c(minyterm, maxyterm) )
#   plot(loglambda, beta[1,], ylim=c(minyterm, maxyterm) )
  for (i_ in 1:dim(beta)[1]) {
    lines(loglambda, beta[i_,], pch=16, col=i_)
  }
}

############################################################################################################
############################################################################################################
##### plot relaxed CV deviances ############################################################################
#' Plot cross-validation deviances, or model coefficients.   
#' 
#' @description 
#' By default, with coefs=FALSE, plots the average deviances as function of lam (lambda) and gam (gamma), and also 
#' indicates the gam and lam which minimize deviance based upon a cv.glmnetr() output object.
#' Optionally, with coefs=TRUE, plots the relaxed lasso coefficients.  
#'
#' @param x a cv.glmnetr()  output object.  
#' @param gam a specific level of gamma for plotting.  By default gamma.min will be used.  
#' @param lambda.lo a lower limit of lambda when plotting.  
#' @param plup   an indicator to plot the upper 95 percent two-sided confidence limits.  
#' @param title  a title for the plot.  
#' @param coefs  default of FALSE plots deviances, option of TRUE plots coefficients.  
#' @param comment default of TRUE to write to console information on lam and gam selected for output.
#' FALSE will suppress this write to console. 
#' @param ... Additional arguments passed to the plot function.  
#'
#' @return This program returns a plot to the graphics window, and may provide 
#' some numerical information to the R Console.  If gam is not specified, then 
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
#' @examples 
#' # set seed for random numbers, optionally, to get reproducible results
#' set.seed(82545037)
#' sim.data=glmnetr.simdata(nrows=100, ncols=100, beta=NULL)
#' xs=sim.data$xs 
#' y_=sim.data$y_ 
#' event=sim.data$event
#' # for this example we use a small number for folds_n to shorten run time 
#' cv_glmnetr_fit = cv.glmnetr(xs, NULL, y_, NULL, family="gaussian", folds_n=3, limit=2) 
#' plot(cv_glmnetr_fit)
#' plot(cv_glmnetr_fit, coefs=1)
#' 
plot.cv.glmnetr = function(x, gam=NULL, lambda.lo=NULL, plup=0, title=NULL, coefs=FALSE, comment=TRUE, ...) {
  object = x 
  method=0 
  if (coefs==TRUE) {
    plot.glmnetr(object, gam=gam, lambda.lo=lambda.lo, title=NULL) 
  } else {
  gamma            = object$relaxed$gamma
  nzero            = object$nzero
  lambda.min       = object$relaxed$lambda.min
  gamma.min        = object$relaxed$gamma.min
  df.lambda.min    = nzero[object$relaxed$index[1,1]]
  lambda.min.g0    = object$relaxed$lambda.min.g0
  df.lambda.min.g0 = nzero[object$relaxed$index.g0[1]]
  lambda.min.g1    = object$lambda.min 
  df.lambda.min.g1 = nzero[object$index[1]]
  if (comment==TRUE) { 
    cat(paste0(" min CV average deviance (max log likelihood) for ", "\n", 
               "   relaxed at log(lambda)  = ",        ceiling(log(lambda.min   )*1000)/1000, ", gamma.min = ", gamma.min, ", df = " , df.lambda.min, "\n",
               "   fully relaxed at log(lambda)   = ", ceiling(log(lambda.min.g0)*1000)/1000, ", df = " , df.lambda.min.g0, "\n",
               "   fully penalized at log(lambda) = ", ceiling(log(lambda.min.g1)*1000)/1000, ", df = " , df.lambda.min.g1, "\n")) 
  }
  
  statlist=object$relaxed$statlist 
  statlist[[5]][[2]][1:8]
#  index

  gamma_n = length(statlist)

  if (!is.null(lambda.lo)) {
    index = statlist[[1]][[1]] >= lambda.lo 
    nzero     = nzero[index] 
    for (k_ in 1:gamma_n) {
      statlist[[k_]][[1]] = statlist[[k_]][[1]][index]
      statlist[[k_]][[2]] = statlist[[k_]][[2]][index]
      statlist[[k_]][[4]] = statlist[[k_]][[4]][index]
      statlist[[k_]][[5]] = statlist[[k_]][[5]][index]
    }
  }
  
  loglambda = log(statlist[[1]][[1]])
  
  maxxterm = ceiling(10*max(loglambda))/10
  minxterm = floor(10*min(loglambda))/10
  
  lambda_n = length(statlist[[1]][[2]])
  cvmg0 = statlist[[1]][[2]] ; cvmg0 ; 
  cvmg0 = cvmg0[1:(lambda_n-1)] ; cvmg0 ; 
  cvmg1 = statlist[[gamma_n]][[2]]
  cvmg1 = cvmg1[1:(lambda_n-1)]
  cvms = c()
  for (k_ in 1:gamma_n) {
    cvms = cbind(cvms,statlist[[k_]][[2]])
  }
  maxcvm = max ( cvms )
  mincvm = min ( cvms )

  maxyterm = ceiling(1000*maxcvm)/1000
  minyterm = floor(1000*mincvm)/1000

  # llmin = log(lambda.min)
  family = object$sample$family 
  if (family=="cox") { ylab="deviance/event" } else { ylab="deviance/record" }
  plot(log(statlist[[1]][[1]][1]), statlist[[1]][[2]][1], xlim=c(minxterm, maxxterm) , ylim=c(minyterm, maxyterm), 
       main=NULL, xlab="log(lambda)", ylab=ylab ) 
  #plot(log(statlist[[1]][[1]][1]), statlist[[1]][[2]][1], xlim=c(minxterm, maxxterm) , ylim=c(minyterm, maxyterm), 
  #     xlab="", type="o", xaxt = "n", ylab="log likelihood / subject" ) 
  # axis(1, at = 1:lambda_n, labels = loglambda, side = 1)
  #axis(1, at = 1:lambda_n, labels = nzero, line = 2.5, side = 3 )
  axis(at = loglambda[1:lambda_n], labels = nzero, tick=FALSE, line = -0.5, side = 3 , cex.axis=1)
  abline(a=NULL, b=NULL, v=log(lambda.min), h=NULL)
  abline(a=NULL, b=NULL, v=log(lambda.min.g0), h=NULL)
  #title(main="Relaxed lasso CV", sub="sub-title", xlab="log(lambda)", ylab="log likelihood / subject")
  #     xlab="log(lambda)", ylab="log likelihood / subject" )

  main = "Relaxed lasso CV Deviance" 
  
  if (is.null(gam)) {main = paste0(main, ", gamma.min=", gamma.min) 
  } else { main = paste0(main, ", gamma=", gam)  }
  
  title(main=main)

  plotr = function(statlist, j_,lty=1) { 
    for (k_ in 1:gamma_n) {
      lines(log(statlist[[k_]][[1]]), statlist[[k_]][[j_]], pch=16, lty=lty, col=k_)
    }
  }

#  dosingle = as.numeric(c(0)) ; 
  dosingle = c(0) ; 
  if (!is.null(gam)) {
    if (gam %in% gamma) { dosingle = 1 } 
  }
#  print(gam) 
  
  if (dosingle) {
    col_ = c(1:gamma_n)[(gam==gamma)]
    lines(log(statlist[[col_]][[1]]), statlist[[col_]][[2]], pch=16, lty=1, col=1)
    lines(log(statlist[[col_]][[1]]), statlist[[col_]][[4]], pch=16, lty=3, col=1)
    lines(log(statlist[[col_]][[1]]), statlist[[col_]][[5]], pch=16, lty=3, col=1)
  } else { 
    plotr(statlist,2)
    if (plup) { plotr(statlist,4,3) } 
  }
  }
}

####################################################################################################
####################################################################################################

#' Plot the cross validated relaxed lasso deviances or coefficients from a nested.glmnetr call.  See 
#' plot.cv.glmnetr().    
#'
#' @param x A nested.glmnetr output object 
#' @param gam A specific level of gamma for plotting.  By default gamma.min will be used.  
#' @param lambda.lo A lower limit of lambda when plotting.  
#' @param title A title
#' @param plup Plot upper 95 percent two-sided confidence intervals for the deviance plots.  
#' @param coefs Default is FALSE to plot deviances.  Option of TRUE to plot coefficients.  
#' @param comment Default of TRUE to write to console information on lam and gam selected for output.
#' FALSE will suppress this write to console.  
#' @param ... Additional arguments passed to the plot function.  
#'
#' @return This program returns a plot to the graphics window, and may provide 
#' some numerical information to the R Console.  
#' 
#' @seealso
#'   \code{\link{plot.glmnetr}} , \code{\link{plot.cv.glmnetr}} , \code{\link{nested.glmnetr}} 
#'
#' @author Walter Kremers (kremers.walter@mayo.edu)
#' 
#' @export
#' 
#' @examples
#' \donttest{
#' sim.data=glmnetr.simdata(nrows=1000, ncols=100, beta=NULL)
#' xs=sim.data$xs 
#' y_=sim.data$yt
#' event=sim.data$event
#' # for this example we use a small number for folds_n to shorten run time 
#' fit3 = nested.glmnetr(xs, NULL, y_, event, family="cox", folds_n=3) 
#' plot(fit3)
#' plot(fit3, coefs=TRUE)
#' }
#' 
plot.nested.glmnetr = function(x, gam=NULL, lambda.lo=NULL, title=NULL, plup=0, coefs=FALSE, comment=TRUE, ... ) {
  object = x 
  dolasso = object$fits[1]
  if (dolasso==1) {
    cv_glmnetr_fit = object$cv_glmnet_fit
    plot(cv_glmnetr_fit, gam=gam, lambda.lo=lambda.lo, plup=plup, title=title, coefs=coefs,comment=comment, ... ) 
    #    else { plot.glmnetr(cv_glmnetr_fit, gam=gam, lambda.lo=lambda.lo, title=title) }
  } else { 
    cat(paste0(" No relaxed lasso for plotting" , "\n")) 
  }
}

####################################################################################################
####################################################################################################


