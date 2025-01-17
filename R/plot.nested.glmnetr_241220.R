####################################################################################################
####################################################################################################

#' Plot results from a nested.glmnetr() output
#' 
#' @description Plot the nested cross validation performance numbers, cross 
#' validated relaxed lasso deviances or coefficients from a nested.glmnetr() call.  
#'
#' @param x A nested.glmnetr output object
#' @param type type of plot to be produced form the (nested) cross validation 
#' performance measures, and the lasso model tuning or lasso model 
#' coefficients. For the lasso model the options include "lasso" to plot deviances
#' informing hyperparmeter choice or "coef" to plot lasso parameter estimates.  
#' Else nested cross validation performance measures are 
#' plotted.  
#' To show cross validation performance measures the options include 
#' "devrat" to plot deviance ratios, i.e. the 
#' fractional reduction in deviance relative to the null model deviance, 
#' "agree" to plot agreement, "lincal" to plot the linear
#' calibration slope coefficients, "intcal" to plot the linear calibration intercept 
#' coefficients or "devian" to plot the deviances from the nested cross 
#' validation. 
#' For each performance measure estimates from the individual (outer) 
#' cross validation fold are depicted by thin lines of different colors and styles, 
#' while the composite value from all fol=ds is depicted by a thicker black line, 
#' and the performance measures naively calculated on the all data using the model 
#' derived from all data is depicted in a thicker red line.  
#' @param gam A specific level of gamma for plotting.  By default gamma.min will 
#' be used.  Applies only for type = "lasso".   
#' @param lambda.lo A lower limit of lambda when plotting.  Applies only for type = "lasso". 
#' @param title A title
#' @param plup Plot upper 95 percent two-sided confidence intervals for the deviance 
#' plots.  Applies only for type = "lasso". 
#' @param coefs Depricated.  See option 'type'.  To plot coefficients specify type = "coef". 
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
#' event=sim.data$event
#' # for this example we use a small number for folds_n to shorten run time 
#' fit3 = nested.glmnetr(xs, NULL, y_, event, family="cox", folds_n=3) 
#' plot(fit3)
#' plot(fit3, type="coef")
#' }
#' 
plot.nested.glmnetr = function(x, type="devrat", gam=NULL, lambda.lo=NULL, title=NULL, plup=0, coefs=FALSE, comment=TRUE,
                               pow=2, ylim=1, plot=1, fold=1, xgbsimple=0, ... ) {
  
  if (is.null(type)) { type = "devrat" } 
  
  if (substr(type,1,4) == "coef") { 
    type = "coef"
    coefs = 1 
  }
  
  if ((coefs == 1) & (type == "devrat")) {
    type = "coef" 
    warning('  The coefs option is depricated.  Use instead type="coef"')
  }
  
  object = x 
  if ( type %in% c('lasso', 'coef') ) {
    dolasso = object$fits[1]
    if (dolasso==1) {
      cv_glmnetr_fit = object$cv_glmnet_fit
      plot(cv_glmnetr_fit, gam=gam, lambda.lo=lambda.lo, plup=plup, title=title, coefs=coefs, comment=comment, ... ) 
      #    else { plot.glmnetr(cv_glmnetr_fit, gam=gam, lambda.lo=lambda.lo, title=title) }
    } else { 
      cat(paste0(" No relaxed lasso for plotting" , "\n")) 
    }
  } else if (type %in% c( "devian", "agree", "lincal", "intcal", "devrat")) {
    #    type="agree", pow=2, ylim=1, plot=1
    if (plot == 0) {
      rlist = plot_perf_glmnetr( object, type=type, pow=pow, ylim=ylim, plot=0) 
      return(rlist)
    } else {
      plot_perf_glmnetr( object, type=type, pow=pow, ylim=ylim, plot=plot, fold=fold, xgbsimple=xgbsimple) 
    }
  } else {
    cat(paste("  Invalid value for type, ", type, "\n") )
  }
}

####################################################################################################
####################################################################################################
