####################################################################################################
####################################################################################################

#' Plot results from a nested.glmnetr() output
#' 
#' @description Plot the nested cross validation performance numbers, cross 
#' validated relaxed lasso deviances or coefficients from a nested.glmnetr() call.  
#'
#' @param x A nested.glmnetr output object
#' @param type type of plot to be produced from the (nested) cross validation 
#' model fits and evaluations. One of c("devrat", "devian", "agree", 
#' "intcal", "lincal") to plot estimates of one these performance measures. 
#' One of c("lasso", "elastic", "ridge") 
#' to plot model coefficients or deviances as function of lambda 
#' and to some degree gamma and alpha. Default is "devrat", 
#' the fractional reduction in deviance relative to the null model deviance. 
#' Use, "devrat" to plot deviance ratios, "devain" to plot devainces, 
#' "agree" to plot agreement e.g. R-square or concordance), "lincal" to plot 
#' the linear calibration slope coefficients, "intcal" to plot the 
#' linear calibration intercept coefficients or "devian" to plot the 
#' deviances from the nested cross validation. 
#' For each performance measure estimates from the individual (outer) 
#' cross validation fold are depicted by thin lines of different colors and styles, 
#' while the composite value from all folds is depicted by a thicker black line, 
#' and the performance measures naively calculated on the all data using the model 
#' derived from all data is depicted in a thicker red line.  
#' @param alpha A specific value of alpha for plotting.  Used only when type is 
#' set to "elastic". Specifies which alpha is to be used for deviance plots.
#' Default is "alpha.min", else must be an element of the alpha vector used in 
#' running the elastic net model. This can be reviewed using summary(fit) 
#' where fit is a nested.glmnetr() output object. Note, alpha is 1 for the 
#' lasso model and alpha is 0 for the ridge model.  
#' @param gamma A specific level of gamma for plotting.  By default gamma.min will 
#' be used.  Applies only for types in c("lasso", "elastic"). 
#' @param lambda.lo A lower limit of lambda when plotting.  Applies only for type = "lasso". 
#' @param title A title
#' @param plup Plot upper 95 percent two-sided confidence intervals for the deviance 
#' plots.  Applies only for type = "lasso". 
#' @param coefs 1 (TRUE) to plot coefficients, else 0 (FALSE) to plot deviances 
#' as function of tuning paramters. Only applies for type in 
#' c("devrat", "devian", "agree", "intcal", "lincal").  See option 'type'.  
#' @param comment Default of TRUE to write to console information on lam and gam selected for output.
#' FALSE will suppress this write to console.  Applies only for type = "lasso". 
#' @param pow Power to which agreement is to be raised when the "gaussian" model 
#' is fit, i.e. 2 for R-square, 1 for correlation.  Does not apply to type = "lasso". 
#' @param ylim y axis limits for model performance plots, i.e. does not apply to 
#' type = "lasso".  The ridge model may calibrate very poorly obscuring plots for 
#' type of "lincal" or "intcal", so one may specify the ylim value.  If ylim is 
#' set to 1, then the program will derive a reasonable range for ylim.  If ylim is 
#' set to 0, then the entire range for all models will be displayed.  Does not 
#' apply to type = "lasso". 
#' @param plot By default 1 to produce a plot, 0 to return the data used in the 
#' plot in the form of a list. 
#' @param fold By default 1 to display model performance estimates form 
#' individual folds (or replicaitons for boostrap evaluations) when type of 
#' "agree", "intcal", "lincal", "devrat" for "devian". If 0 then the individual 
#' fold calculations are not displayed. When there are many replications as
#' sometimes the case when using bootstrap, one may specify the number of 
#' randomly selected lines for plotting.   
#' @param xgbsimple 1 (default) to include results for the untuned XGB model, 0 to not include.  
#' @param track 2 to track progress by printing to console, 0 (default) to 
#' not track.
#' @param ... Additional arguments passed to the plot function.  
#'
#' @return This program returns a plot to the graphics window, and may provide 
#' some numerical information to the R Console.  
#' 
#' @seealso
#'   \code{\link{plot_perf_glmnetr}} , \code{\link{calplot}} , \code{\link{plot.cv.glmnetr}} , \code{\link{nested.glmnetr}}
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
#' yg=sim.data$y_
#' event=sim.data$event
#' # for this example we use a small number for folds_n to shorten run time 
#' fit3 = nested.glmnetr(xs, NULL, yg, event, family="gaussian", folds_n=3, resample=0) 
#' plot(fit3)
#' plot(fit3, coefs=1)
#' }
#' 
plot.nested.glmnetr = function(x, type="devrat", alpha=NULL, gamma=NULL, lambda.lo=NULL, title=NULL, plup=0, coefs=0, comment=TRUE,
                               pow=2, ylim=1, plot=1, fold=1, xgbsimple=0, track=0, ... ) {
  if (track >= 2) { cat( "  in plot.nested.glmnetr   class(x) = ", class(x), "\n") }
  if (x$resample== 0){
    if (!(type %in% c('lasso', 'elastic', 'ridge'))) { type = "lasso" }
    if (!(coefs %in% (c(0,1)))) { coefs = 0 } 
  }
  
  if (is.null(type)) { type = "devrat" } 
  
  if (substr(type,1,4) == "coef") { 
    type = "lasso"
    coefs = TRUE  
  }
  
  if ((coefs == 1) & (type == "devrat")) { cat(paste('  Type devrat takes priority over coefs of 1' )) }
  
  object = x 
  if ( type[1] %in% c('lasso', 'elastic', 'ridge') ) {
    dolasso = object$fits[1]
    if (dolasso==1) {
      plot.cv.glmnetr(object, type=type, alpha=alpha, gamma=gamma, lambda.lo=lambda.lo, plup=plup, title=title, coefs=coefs, comment=comment, ... )
    } else { 
      cat(paste0(" Specifid analysis not found. \n\n")) 
    }
  } else if (type %in% c( "devian", "devrat", "agree", "lincal", "intcal")) {
    #    type="agree", pow=2, ylim=1, plot=1
    if (plot == 0) {
      rlist = plot_perf_glmnetr( object, type=type, pow=pow, ylim=ylim, plot=0, track=track) 
      return(rlist)
    } else {
      plot_perf_glmnetr( object, type=type, pow=pow, ylim=ylim, plot=plot, fold=fold, xgbsimple=xgbsimple, track=track) 
    }
  } else {
    cat(paste("  Invalid value for type, ", type, "\n") )
  }
}

####################################################################################################
####################################################################################################
