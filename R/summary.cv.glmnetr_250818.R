################################################################################
##### summary.cv.glmnetr_yymmdd.R ##############################################
################################################################################
#' Output summary for elastic net models fit within a nested.glmnetr() output object.  
#' 
#' @description 
#' Summarize the cross-validation informed model fit.  The fully penalized
#' (gamma=1) beta estimate will not be given by default but can too be output  
#' using printg1=TRUE.  
#'
#' @param object a nested.glmnetr() output object.  
#' @param printg1 TRUE to also print out the fully penalized lasso beta, else FALSE to suppress.  
#' @param orderall By default (orderall=FALSE) the order terms enter into the lasso model 
#' is given for the number of terms that enter in lasso minimizing loss model.  If 
#' orderall=TRUE then all terms that are included in any lasso fit are described.    
#' @param type once of c("lasso", "elastic", "ridge") to select for summarizing,  
#' with default of "lasso". 
#' @param ... Additional arguments passed to the summary function.  
#'
#' @return Coefficient estimates (beta) 
#' 
#' @seealso
#'   \code{\link{predict.cv.glmnetr}} , \code{\link{nested.glmnetr}} 
#' 
#' @export
#'
summary.cv.glmnetr = function(object, printg1="FALSE", orderall=FALSE, type="lasso", ...) {
  if        ( substr(object$version[2],1,21) == "glmnetr version 0.6-2" ) { 
    summary.cv.glmnetr_0_6_2(object, printg1=printg1, orderall=orderall, type=type, ...)
  } else if ( substr(object$version[2],1,21) == "glmnetr version 0.6-1" ) { 
    summary.cv.glmnetr_0_6_1(object, printg1=printg1, orderall=orderall, type=type, ...)
  } else {
    summary.cv.glmnetr_0_5_5(object, printg1=printg1, orderall=orderall, ...) 
  }
}

###############################################################################################################
###############################################################################################################
