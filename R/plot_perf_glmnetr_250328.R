################################################################################
##### plot.perf.glmnetr_yymmdd.R ###############################################
################################################################################
#' Plot nested cross validation performance summaries
#' 
#' @description 
#' This function plots summary information from a nested.glmnetr() output object, that 
#' is from a nested cross validation performance.  Alternamvely one can output the 
#' numbers otherwise displayed to a list for extraction or customized plotting.  Performance
#' measures for plotting include "devrat" the deviance ratio, i.e. the fractional 
#' reduction in deviance relative to the null model deviance, "agree" a measure 
#' of agreement, "lincal" the slope from a linear calibration and "intcal" the 
#' intercept from a linear calibration.  Performance measure estimates 
#' from the individual (outer) cross validation fold are depicted by thin lines 
#' of different colors and styles, while the composite value from all folds is 
#' depicted by a thicker black line, and the performance measures naively 
#' calculated on the all data using the model derived from all data is 
#' depicted by a thicker red line. 
#' 
#' @param x A nested.glmnetr output object
#' @param type determines what type of nested cross validation performance measures are 
#' plotted.  Possible values are 
#' "devrat" to plot the deviance ratio, i.e. the fractional reduction 
#' in deviance relative to the null model deviance,
#' "devian" to plot deviance, 
#' "agree" to plot agreement in terms of concordance, correlation or R-square, 
#' "lincal" to plot the linear calibration slope coefficients, 
#' "intcal" to plot the linear calibration intercept coefficients, 
#' from the (nested) cross validation. 
#' @param pow Power to which agreement is to be raised when the "gaussian" model 
#' is fit, i.e. 2 for R-square, 1 for correlation.  Does not apply to type = "lasso". 
#' @param ylim y axis limits for model perforamnce plots, i.e. does not apply to 
#' type = "lasso".  The ridge model may calibrate very poorly obscuring plots for 
#' type of "lincal" or "intcal", so one may specify the ylim value.  If ylim is 
#' set to 1, then the program will derive a reasonable range for ylim.  If ylim is 
#' set to 0, then the entire range for all models will be displayed.  Does not 
#' apply to type = "lasso". 
#' @param fold By default 1 to display using a spaghetti the performance as 
#' calculated from the individual folds, 0 to display using dots only the composite 
#' values calculated using all folds.  
#' @param xgbsimple 1 to include results for the untuned XGB model, 0 (default) to not include.  
#' @param plot By default 1 to produce a plot, 0 to return the data used in the 
#' plot in the form of a list. 
#' @param track 2 to track progress by printing to console, 0 (default) to 
#' not track.
#' 
#' @return This program returns a plot to the graphics window by default, and returns
#' a list with data used in teh plots if the plot=1 is specified.
#' 
#' @seealso
#'   \code{\link{plot.nested.glmnetr}} , \code{\link{nested.glmnetr}}
#'
#' @author Walter Kremers (kremers.walter@mayo.edu)
#' 
#' @export
#' 
plot_perf_glmnetr = function( x, type="devrat", pow=2, ylim=1, fold=1, xgbsimple=0, plot=1, track=0) {
  if (track >=3) { cat( "  in plot_perf_glmnetr   class(x) = ", class(x), "\n") }
  object = x 
  if (is.null(x$version[2])) {
    cat("     Output object 'glmnetr' version not identified. \n",
        "    Analysis will not be performed.") 
  } else { 
    if        (substr(x$version[2],1,21) == "glmnetr version 0.6-1") { 
      plot_perf_glmnetr_0_6_1(object, type=type, pow=pow, ylim=ylim, fold=fold, xgbsimple=xgbsimple, plot=plot, track=track) 
    } else { 
      plot_perf_glmnetr_0_5_5(object, type=type, pow=pow, ylim=ylim, fold=fold, xgbsimple=xgbsimple, plot=plot)
    }
  }
} 
