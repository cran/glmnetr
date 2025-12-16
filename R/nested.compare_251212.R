###############################################################################################################
###############################################################################################################
#' Compare cross validation fit performances from a nested.glmnetr output.
#'
#' @description 
#' Compare cross-validation model fits in terms of average performances from the 
#' nested cross validation fits. In general the standard deviations for the
#' performance measures evaluated on the leave-out samples may be biased. While 
#' the standard deviations of the paired within fold differences of 
#' performances intuitively might be less biased this has not been shown. See 
#' the package vignettes for more discussion.    
#'
#' @param object A nested.glmnetr output object.
#' @param type determines what type of nested cross validation performance measures are 
#' compared.  Possible values are "devrat" to compare the deviance ratios, i.e. the 
#' fractional reduction in deviance relative to the null model deviance, 
#' "agree" to compare agreement, "lincal" to compare the linear calibration 
#' slope coefficients, "intcal" to compare the linear calibration intercept 
#' coefficients, from the nested cross validation. 
#' @param digits digits for printing of z-scores, p-values, etc. with default of 4 
#' @param pow the power to which the average of correlations is to be raised.  Only 
#' applies to the "gaussian" model.  Default is 2 to yield R-square but can be on to 
#' show correlations.  pow is ignored for the family of "cox" and "binomial".  
#' @param table 1 to print table to console, 0 to output the tabled information to a data frame
#'  
#' @return A printout to the R console. 
#' 
#' @seealso
#'   \code{\link{nested.cis}} , \code{\link{summary.nested.glmnetr}} , \code{\link{nested.glmnetr}} 
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
#' nested.compare(fit3)
#' }
#' 
nested.compare = function( object, type="devrat", digits=4, pow=1, table=1 ) {
  if ( (type == "intcal") & (object$sample[1] == 'cox') ) { 
    cat("     There is no intercept for the Cox model to compare! \n",
        "    Analysis will not be performed.") 
  } else if (is.null(object$version[2])) {
    cat("     Output object 'glmnetr' version not identified. \n",
        "    Analysis will not be performed.") 
  } else {
    if        (substr(object$version[2],1,21) == "glmnetr version 0.6-3") { nested.compare_0_6_2(object, type=type, digits=digits, pow=pow, table=table ) 
    } else if (substr(object$version[2],1,21) == "glmnetr version 0.6-2") { nested.compare_0_6_2(object, type=type, digits=digits, pow=pow, table=table ) 
    } else if (substr(object$version[2],1,21) == "glmnetr version 0.6-1") { nested.compare_0_6_1(object, type=type, digits=digits, pow=pow ) 
    } else if (object$version[2] == "glmnetr version 0.5-5 (2024-12-28)") { nested.compare_0_5_3(object, type=type, digits=digits, pow=pow ) 
    } else if (object$version[2] == "glmnetr version 0.5-4 (2024-10-24)") { nested.compare_0_5_3(object, type=type, digits=digits, pow=pow ) 
    } else if (object$version[2] == "glmnetr version 0.5-3 (2024-08-28)") { nested.compare_0_5_3(object, type=type, digits=digits, pow=pow ) 
    } else if (object$version[2] == "glmnetr version 0.5-2 (2024-07-10)") { nested.compare_0_5_2(object, type=type, digits=digits, pow=pow ) 
    } else if (object$version[2] == "glmnetr version 0.5-1 (2024-05-10)") { nested.compare_0_5_1(object, type=type, digits=digits, pow=pow ) 
    } else if (object$version[2] == "glmnetr version 0.4-6 (2024-04-21)") { nested.compare_0_5_1(object, type=type, digits=digits, pow=pow ) 
    } else if (object$version[2] == "glmnetr version 0.4-5 (2024-04-20)") { nested.compare_0_5_1(object, type=type, digits=digits, pow=pow ) 
    } else if (object$version[2] == "0.4-5 dev 240410") { nested.compare_0_5_1(object, type=type, digits=digits, pow=pow ) 
    } else if (object$version[2] == "0.4-4 dev 240322") { nested.compare_0_5_1(object, type=type, digits=digits, pow=pow ) 
    } else if (object$version[2] == "0.4-3") { nested.compare_0_5_1(object, type=type, digits=digits, pow=pow ) 
    } else if (object$version[2] == "0.4-2") { 
      if (type == "intcal") { cat("\n intcal is not stored for most models")
      } else if (type == "devrat") { cat("\n     Needed info for devrat is not available\n",
                                         "    Analysis will not be performed.")
      } else { nested.compare_0_5_1(object, type=type, digits=digits, pow=pow ) 
      }
    } 
  }
}

###############################################################################################################
###############################################################################################################

#' A redirect to nested.compare
#'
#' @description 
#' See nested.compare(), as glmnetr.compcv() is depricated  
#'
#' @param object A nested.glmnetr output object.
#' @param digits digits for printing of z-scores, p-values, etc. with default of 4 
#' @param pow the power to which the average of correlations is to be raised.   
#' @param type determines what type of nested cross validation performance measures are 
#' compared.  Possible values are "devrat" to compare the deviance ratios, i.e. the 
#' fractional reduction in deviance relative to the null model deviance, 
#' "agree" to compare agreement, "lincal" to compare the linear calibration 
#' slope coefficients, "intcal" to compare the linear calibration intercept 
#' coefficients, from the nested cross validation. 
#'  
#' @return A printout to the R console. 
#' 
#' @seealso
#'   \code{\link{nested.compare}}
#' 
#' @export
#' 

glmnetr.compcv = function(object, digits=4, type="devrat", pow=1) {
  nested.compare(object, type, digits, pow) 
}
  
###############################################################################################################
###############################################################################################################
