################################################################################
##### plot.cv.glmnetr_yymmdd.R #################################################
################################################################################
#' Plot the relaxed lasso coefficients.  
#' 
#' @description 
#' Plot the relaxed lasso, elastic net or ridge model coefficients from a nested.glmnetr() output object.
#' One may specify a value for gamma.  If gamma is unspecified (NULL), then 
#' the plot will be for the gamma which minimizes loss. 
#'
#' @param x A nested.glmnetr output object.  
#' @param type one of c("lasso", "elastic", "ridge") to plot the deviance 
#' curves of the respective model fit. Default is "lasso". 
#' @param alpha A specific level of alpha for plotting. By default alpha.min
#' will be used such that the triplet (alpha.min, gamma.min, lambda.min) 
#' minimizes the model deviance.
#' @param gamma A specific level of gamma for plotting.  By default gamma.min will 
#' be used such that the pair (gamma.min, lambda.min) minimzes the model 
#' deviance.
#' @param lambda.lo A lower limit of lambda for plotting.  
#' @param title A title for the plot  
#' @param comment Default of TRUE to write to console information on gamma and 
#' lambda selected for output. FALSE will suppress this write to console.  
#' @param ... Additional arguments passed to the plot function.  
#'
#' @return This program returns a plot to the graphics window, and may provide 
#' some numerical information to the R Console.  If gamma is not specified, 
#' then the gamma.min from
#' the deviance minimizing (lambda.min, gamma.min) pair will be used, and the
#' minimizing lambda.min will be indicated by a vertical line.  Also, if one 
#' specifies gamma=0, the lambda which minimizes deviance for the restricted 
#' set of models where gamma=0 will indicated by a vertical line.  
#' 
#' @seealso 
#'   \code{\link{plot.cv.glmnetr}} , \code{\link{plot.nested.glmnetr}}  
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
#' yt=sim.data$yt
#' yg=sim.data$yt
#' event=sim.data$event
#' glmnetr.fit = nested.glmnetr(xs, start=NULL, yg, event=event, family="gaussian",
#'   resample=0, folds_n=4)
#' plot(glmnetr.fit, type="lasso")
#' }
#' 
plot.glmnetr = function(x, type="lasso", alpha=NULL, gamma=NULL, lambda.lo=NULL, title=NULL,comment=TRUE, ...) {
  if        ( substr(x$version[2],1,21) == "glmnetr version 0.6-3" ) {
    plot.glmnetr_0_6_1(x, type=type, alpha=alpha, gamma=gamma, lambda.lo=lambda.lo, title=title, comment=comment ,  ...) 
  } else if ( substr(x$version[2],1,21) == "glmnetr version 0.6-2" ) {
    plot.glmnetr_0_6_1(x, type=type, alpha=alpha, gamma=gamma, lambda.lo=lambda.lo, title=title, comment=comment ,  ...) 
  } else if ( substr(x$version[2],1,21) == "glmnetr version 0.6-1" ) {
    plot.glmnetr_0_6_1(x, type=type, alpha=alpha, gamma=gamma, lambda.lo=lambda.lo, title=title, comment=comment ,  ...) 
  } else {
    plot.glmnetr_0_5_5(x, gam=gamma, lambda.lo=lambda.lo, title=title, comment=comment , ...)  
  }
}

############################################################################################################
############################################################################################################
##### plot relaxed CV deviances ############################################################################
#' Plot cross-validation deviances, or model coefficients.   
#' 
#' @description 
#' By default, with coefs=FALSE, plots the average deviances as function of lambda and gamma, and also 
#' indicates the gamma and lambda which minimize deviance for the lasso, elastic net or ridge model. 
#' Optionally, with coefs=TRUE, plots the relaxed lasso coefficients.  
#'
#' @param x a nested.glmnetr()  output object. 
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
#' @param comment default of TRUE to write to console information on lambda and gamma selected for output.
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
#'   \code{\link{plot.glmnetr}} , \code{\link{plot.nested.glmnetr}}  
#'   
#' @export 
#' 
#' @importFrom graphics abline axis lines 
#'
plot.cv.glmnetr = function(x, type="lasso", alpha=NULL, gamma=NULL, lambda.lo=NULL, plup=0, title=NULL, coefs=FALSE, comment=TRUE, lty=1, track=0, ...) {
  if (track >= 2) { cat( "  in plot.cv.glmnetr   class(x) = ", class(x), "\n") }
  if        ( substr(x$version[2],1,21) == "glmnetr version 0.6-3" ) {
    plot.cv.glmnetr_0_6_1(x, type=type, alpha=alpha, gamma=gamma, lambda.lo=lambda.lo, plup=plup, title=title, coefs=coefs, comment=comment, lty=lty, track=track, ...) 
  } else if ( substr(x$version[2],1,21) == "glmnetr version 0.6-2" ) {
    plot.cv.glmnetr_0_6_1(x, type=type, alpha=alpha, gamma=gamma, lambda.lo=lambda.lo, plup=plup, title=title, coefs=coefs, comment=comment, lty=lty, track=track, ...) 
  } else if ( substr(x$version[2],1,21) == "glmnetr version 0.6-1" ) {
    plot.cv.glmnetr_0_6_1(x, type=type, alpha=alpha, gamma=gamma, lambda.lo=lambda.lo, plup=plup, title=title, coefs=coefs, comment=comment, lty=lty, track=track, ...) 
  }else {
    plot.cv.glmnetr_0_5_5(x,                         gam  =gamma, lambda.lo=lambda.lo, plup=plup, title=title, coefs=coefs, comment=comment, ...)
  }
}

####################################################################################################
####################################################################################################
