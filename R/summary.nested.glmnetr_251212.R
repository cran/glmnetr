################################################################################
##### summary.nested.glmnetr_yymmdd.R ##########################################
################################################################################
#' A redirect to the summary() function for nested.glmnetr() output objects 
#'
#' @param x a nested.glmnetr() output object.  
#' @param ... additional pass through inputs for the print function.
#' 
#' @return - a nested cross validation fit summary, or a cross validation model summary.  
#' 
#' @seealso
#'    \code{\link{summary.nested.glmnetr}} , \code{\link{nested.glmnetr}}
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
#' print(fit3)
#' }
#' 
print.nested.glmnetr = function(x, ...) {
  summary(x) 
}
#
################################################################################
#' round elements of a summary.glmnetr() output
#'
#' @param summdf a summary data frame from summary.nested.glmnetr() obtained using 
#' the option table=0 
#' @param digits the minimum number of decimals to display the elements of the data
#' frame 
#' @param resample 1 (default) if the summdf object is a summary for an analysis including
#' nested cross validation, 0 if only the full data models were fit.
#'
#' @return a data frame with same form as the input but with rounding for easier 
#' display
#' 
#' @seealso
#'   \code{\link{summary.nested.glmnetr}} , \code{\link{nested.glmnetr}} 
#' 
#' @export 
#'
roundperf = function(summdf, digits=3, resample=1) {
  if (resample == 1) {
    set1 = c(1,2,5,6,8) 
    set2 = c(3) 
    set3 = c(4:7)
  } else {
    set1 = c(1,3) 
    set2 = NULL 
    set3 = c(2)
  }
  for ( j_ in set1 ) {
    digitst = digits 
    for (i_ in c(0:3)) {
      if (sum(!is.na(summdf[j_])) == 0 )  {
        conditin = (min(abs(summdf[j_]), rm.na=TRUE) < 2*10^(-digits-i_)) 
        if (is.na(conditin)) { conditin = 0 } 
        if ( conditin ) { digitst = digits + i_ + 1 }
      }
    }
    #        print(c(j_,digitst))
    summdf[j_] = round(summdf[j_], digits=digitst)
  }
  for ( j_ in set2 ) {
    digitst = digits 
    for (i_ in c(0:3)) {
      #          print( summdf[j_] )
      #          print( summdf[j_] - 1)
      conditin = (min(abs(summdf[j_][!is.na(summdf[j_])]-1)) < 2*10^(-digits-i_)) 
      if (is.na(conditin)) { conditin = 0 }
      if (conditin) { digitst = digits + i_ + 1 }
      #          print(c(i_,digitst))
    }
    
    summdf[j_] = round(summdf[j_], digits=digitst)
  }
  for ( j_ in set3 ) {
    digitst = digits 
    for (i_ in c(0:3)) {
      conditin = (max(abs(summdf[j_])) > (1 -2*10^(-digits-i_)))
      if (is.na(conditin)) { conditin = 0 }
      if (conditin) { digitst = digits + i_ + 1 }
    }
    summdf[j_] = round(summdf[j_], digits=digitst)
  }
  return( summdf )
}

################################################################################  
#' Summarize a nested.glmnetr() output object
#'
#' @description 
#' Summarize the model fit from a nested.glmnetr() output object, i.e. the fit of 
#' a cross-validation informed relaxed lasso model fit, inferred by nested cross 
#' validation.  Else summarize the cross-validated model fit.    
#'
#' @param object a nested.glmnetr() output object.  
#' @param cvfit  default of FALSE to summarize fit of a cross validation informed 
#' relaxed lasso model fit, inferred by nested cross validation.  Option of TRUE 
#' will describe the cross validation informed relaxed lasso model itself. 
#' @param type When cvfit is TRUE, one of c("lasso", "elastic", "ridge") to 
#' select for summarizing, with default of "lasso". 
#' @param pow the power to which the average of correlations is to be raised.  Only 
#' applies to the "gaussian" model.  Default is 2 to yield R-square but can be on to 
#' show correlations.  Pow is ignored for the family of "cox" and "binomial".   
#' @param printg1 TRUE to also print out the fully penalized lasso beta, else to suppress.  
#' Only applies to cvfit=TRUE.
#' @param digits digits for printing of deviances, linear calibration coefficients 
#' and agreement (concordances and R-squares).
#' @param call 1 to print call used in generation of the object, 0 or NULL to not print 
#' @param onese 0 (default) to not include summary for 1se lasso fits in tables, 1 to include 
#' @param table 1 to print table to console, 0 to output the tabled information to a data frame
#' @param tuning 1 to print tuning parameters, 0 (default) to not print
#' @param width character width of the text body preceding the performance 
#' measures which can be adjusted between 60 and 120.
#' @param cal 1 print performance statistics for lasso 
#' models calibrated on training data, 2 to print performance statistics for 
#' lasso and random forest models calibrated on training data, 0 (default) to 
#' not print.  Note, despite any intuitive appeal these training data 
#' calibrated models may sometimes do rather poorly.
#' @param ... Additional arguments passed to the summary function.  
#' 
#' @return - a nested cross validation fit summary, or a cross validation model summary.  
#' 
#' @seealso
#'   \code{\link{nested.compare}} , \code{\link{nested.cis}} , \code{\link{summary.cv.glmnetr}} , \code{\link{roundperf}} , 
#'   \code{\link{plot.nested.glmnetr}} , \code{\link{calplot}} , \code{\link{nested.glmnetr}} 
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
#' summary(fit3)
#' }
#' 

summary.nested.glmnetr = function(object, cvfit=FALSE, type="lasso", pow=2, printg1=FALSE, 
                                  digits=4, call=NULL, onese=0, table=1, tuning=0, width=84, cal=0, ...) {
# cvfit=FALSE ; pow=2 ; printg1=FALSE ; digits=4 ; call=NULL ; onese=0 ; table=1 ; tuning=0 ; width=108  ; cal = 1 ; 
  if (is.null(object$version[2])) {
    cat("     Output object 'glmnetr' version not identified. \n",
        "    Analysis will not be performed.") 
  } else { 
    if        (substr(object$version[2],1,21) == "glmnetr version 0.6-3") {summary.nested.glmnetr_0_6_2(object, cvfit=cvfit, type=type, pow=pow, printg1=printg1, digits=digits, call=call, onese=onese, table=table, tuning=tuning, width=width, cal=cal, ...) 
    } else if (substr(object$version[2],1,21) == "glmnetr version 0.6-2") {summary.nested.glmnetr_0_6_2(object, cvfit=cvfit, type=type, pow=pow, printg1=printg1, digits=digits, call=call, onese=onese, table=table, tuning=tuning, width=width, cal=cal, ...) 
    } else if (substr(object$version[2],1,21) == "glmnetr version 0.6-1") {summary.nested.glmnetr_0_6_1(object, cvfit=cvfit, type=type, pow=pow, printg1=printg1, digits=digits, call=call, onese=onese, table=table, tuning=tuning, width=width, cal=cal, ...) 
    } else if (object$version[2] == "glmnetr version 0.5-5 (2024-12-28)") {summary.nested.glmnetr_0_5_3(object, cvfit=cvfit, pow=pow, printg1=printg1, digits=digits, call=call, onese=onese, table=table, tuning=tuning, width=width, cal=cal, ...) 
    } else if (object$version[2] == "glmnetr version 0.5-4 (2024-10-24)") {summary.nested.glmnetr_0_5_3(object, cvfit=cvfit, pow=pow, printg1=printg1, digits=digits, call=call, onese=onese, table=table, tuning=tuning, width=width, cal=cal, ...) 
    } else if (object$version[2] == "glmnetr version 0.5-3 (2024-08-28)") {summary.nested.glmnetr_0_5_3(object, cvfit=cvfit, pow=pow, printg1=printg1, digits=digits, call=call, onese=onese, table=table, tuning=tuning, width=width, cal=cal, ...) 
    } else if (object$version[2] == "glmnetr version 0.5-2 (2024-07-10)") {summary.nested.glmnetr_0_5_2(object, cvfit=cvfit, pow=pow, printg1=printg1, digits=digits, Call=call, onese=onese, table=table, tuning=tuning, width=width, cal=cal, ...) 
    } else if (object$version[2] == "glmnetr version 0.5-1 (2024-05-10)") {summary.nested.glmnetr_0_5_1(object, cvfit=cvfit, pow=pow, printg1=printg1, digits=digits, Call=call, onese=onese, table=table, tuning=tuning, width=width, cal=cal, ...) 
    } else if (object$version[2] == "glmnetr version 0.4-6 (2024-04-21)") {summary.nested.glmnetr_0_4_5(object, cvfit=cvfit, pow=pow, printg1=printg1, digits=digits, Call=call, onese=onese, table=table, tuning=tuning, width=width, cal=cal, ...) 
    } else if (object$version[2] == "glmnetr version 0.4-5 (2024-04-20)") {summary.nested.glmnetr_0_4_5(object, cvfit=cvfit, pow=pow, printg1=printg1, digits=digits, Call=call, onese=onese, table=table, tuning=tuning, width=width, cal=cal, ...) 
    } else if (object$version[2] == "0.4-5 dev 240410") {summary.nested.glmnetr_0_4_5(object, cvfit=cvfit, pow=pow, printg1=printg1, digits=digits, Call=call, onese=onese, table=table, tuning=tuning, width=width, cal=cal, ...) 
    } else if (object$version[2] == "0.4-4 dev 240322") {summary.nested.glmnetr_0_4_4(object, cvfit=cvfit, pow=pow, printg1=printg1, digits=digits, Call=call, onese=onese, table=table, tuning=tuning, width=width, cal=cal, ...) 
    } else if (object$version[2] == "0.4-3") {summary.nested.glmnetr_0_4_3(object, cvfit=cvfit, pow=pow, printg1=printg1, digits=digits, Call=call, onese=onese, table=table, tuning=tuning, width=width, cal=cal, ...) 
    } else if (object$version[2] == "0.4-2") {summary.nested.glmnetr_0_4_2(object, cvfit=cvfit, pow=pow, printg1=printg1, digits=digits, Call=call, onese=onese, table=table, tuning=tuning, width=width, cal=cal, ...) 
    }
  }
} 

####################################################################################################################################
####################################################################################################################################

