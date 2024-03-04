################################################################################

#' Fit a Random Forest model on data provided in matrix and vector formats.
#' 
#' @description Fit an Random Forest model using the rfsrc() function of the 
#' randomForestSRC package.  
#'
#' @param xs     predictor input - an n by p matrix, where n (rows) is sample size, and p (columns) 
#' the number of predictors.  Must be in matrix form for complete data, no NA's, no Inf's, etc.,
#' and not a data frame. 
#' @param start  an optional vector of start times in case of a Cox model. Class numeric of length same as number of patients (n)
#' @param y_     dependent variable as a vector: time, or stop time for Cox model, Y_ 0 or 1 for binomial (logistic), numeric for gaussian. 
#' Must be a vector of length same as number of sample size. 
#' @param event  event indicator, 1 for event, 0 for census, Cox model only.
#' Must be a numeric vector of length same as sample size.
#' @param family model family, "cox", "binomial" or "gaussian" (default) 
#' @param mtryc  a vector (numeric) of values to search over for optimization of the 
#' Random Forest fit.  This if for the mtry input variable of the rfsrc() program 
#' specifying the number of terms to consider in each step of teh Random Forest fit.
#' @param ntreec a vector (numeric) of 2 values, the first for the number of forests
#' (ntree from rfsrc()) to use when searhcing for a better bit and the second to use
#' when fitting the final model.  More trees should give a better fit but 
#' require more computations and storage for the final. 
#' model.   
#' @param seed a seed for set.seed() so one can reproduce the model fit.  If 
#' NULL the program will generate a random seed.  Whether specified or NULL, the 
#' seed is stored in the output object for future reference.  Note,
#' for the default this randomly generated seed depends on the seed in memory at that 
#' time so will depend on any calls of set.seed prior to the call of this function.  
#' @param track 1 to output a brief summary of the final selected model, 2 to 
#' output a brief summary on each model fit in search of a better model or 0 
#' (default) to not output this information.
#'
#' @return a Random Forest model fit 
#' 
#' @author Walter Kremers (kremers.walter@mayo.edu)
#' 
#' @importFrom randomForestSRC rfsrc 
#'
#' @export
#' 
rf_tune = function(xs, start=NULL, y_, event=NULL, family=NULL, mtryc=NULL, ntreec=NULL, seed=NULL, track=0) {
  
  if (is.null(seed)) { seed = round(runif(1)*1e9) } 
  set.seed( seed ) 
  
  if (is.null(family) == 1) {
    if (!is.null(event) == 1) {
      family = "cox"
    } else if ( (length(table(y_)) == 2) & (names(table(y_))[1] == 0) & (names(table(y_))[2] == 1) ) { family = "binomial" 
    } else { family = "gaussian"} 
    print(family) 
  }
  
  if (track >= 1) { 
    time_start = diff_time()
    time_split = time_start
  } 
  
  if (family == "cox") {
    df = as.data.frame(cbind(xs, y_=y_, event = event)) 
  } else {
    df = as.data.frame(cbind(xs, y_)) 
    names(df)
  }
#  df = as.data.frame(cbind(xs, y_=y_)) 
#  dim(df) 
  
#  print(mtryc)
#  print (family) 
  
  if (is.null(mtryc)) {
    if (family %in% c("cox","binomial")) {
      mtryc = round( sqrt(dim(xs)[2]) * c(0.67 , 1, 1.5, 2.25, 3.375, 5.0625) ) 
      mtryc = mtryc [mtryc < dim(xs)[2]] 
      #mtryc = mtryc [mtryc < 1] 
      if (length(mtryc) == 0) { mtryc = 1 }
    } else {
      mtryc = round( dim(xs)[2]*c(0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6) )  
#      mytryc = round( dim(xs)[2]*c(0.33, 0.5, 0.67) )  
    }
  }
  
#  print(mtryc)
  
  if (is.null(ntreec)) {
    ntreec = c(10,  50)
    if (family == "gaussian") { ntreec = c(50, 500) 
    } else { ntreec = c(25, 250) }
  }
  
  ##------------------------------------------------------------------------------
  
  ##mtryc = mtryc[1:3]
  ## k_ = 1 
  
  #  rfsrc does not analyze (start,stop) time survival data 
  if ((family == "cox") & (!is.null(start)==1)) { skip = 1 } else { skip = 0 } 
  
  if (skip == 0) {
    
    for ( k_ in c(1:length(mtryc))) {
      if (family == "cox") {
        rf = rfsrc( Surv(y_, event) ~ . , data=df, mtry=mtryc[k_], nsplit=8, 
                    ntree=ntreec[1], membership = TRUE, importance=TRUE) 
      } else {
        rf = rfsrc( y_ ~ . , data=df, mtry=mtryc[k_], nsplit=8, 
                    ntree=ntreec[1], membership = TRUE, importance=TRUE) 
      }
      
      if (k_ ==1) {
        k_best = 1 
        mtry_best = mtryc[k_] 
        err.rate = rf$err.rate[length(rf$err.rate)]                               ## should be OOB error rate 
        err.rate.best = err.rate
        err.ratev = err.rate
        rf_best = rf 
      }  else if (rf$err.rate[length(rf$err.rate)] < err.rate.best) {
        k_best = k_ 
        mtry_best = mtryc[k_] 
        err.rate.best = rf$err.rate[length(rf$err.rate)]
        rf_best = rf 
      }
      
      if (k_ > 1) { err.ratev = c(err.ratev, rf$err.rate[length(rf$err.rate)]) } 
      
      if (track >= 2) { 
        print( c(k_, mtryc[k_], mtry_best, rf$err.rate[length(rf$err.rate)] ) )
        time_split = diff_time(time_start, time_split)
      }
    }
    
    if (ntreec[2] > ntreec[1]) {
      if (family == "cox") {
        rf_tuned = rfsrc( Surv(y_, event) ~ . , data=df, mtry=mtryc[k_best], nsplit=8, 
                    ntree=ntreec[2], membership = TRUE, importance=TRUE) 
      } else {
        rf_tuned = rfsrc( y_ ~ . , data=df, mtry=mtryc[k_best], nsplit=8, 
                    ntree=ntreec[2], membership = TRUE, importance=TRUE) 
      }
    } else {
      rf_tuned = rf_best
    }
    
    if (track >= 1) {
      print( c(mtry_best, rf_tuned$err.rate[length(rf_tuned$err.rate)] ) )
      time_split = diff_time(time_start, time_split)
    }
    
    rffit = list(rf_tuned=rf_tuned, err.ratev=err.ratev, rf_tuned$err.rate[length(rf_tuned$err.rate)], 
                 mtryc=mtryc, ntreec=ntreec, seed=seed )

  } else { rffit = list(rffit="NONE") }
  
  class(rffit) <- c("rf_tune")
  
  return(rffit)
  
}

################################################################################
#' 
#' Summarize output from rf_tune() function
#'
#' @param object output from an rf_tune() function
#' @param ... optional pass through parameters to pass to summary.rfsrc() 
#'
#' @return summary to console
#' 
#' @export
#'
summary.rf_tune = function(object, ...) {
  cat(paste0("\n    search set for mtry and ntree : \n\n"))
  print(object$mtryc)
  print(object$ntreec)
  cat(paste0("\n    selected vale for for mtry from search: ", object$rf_tuned$mtry, "\n\n"))
  rf_tuned = object$rf_tuned 
  summary(rf_tuned)
}

################################################################################
#'
#' Print output from rf_tune() function
#' 
#' @param x output from an rf_tune() function
#' @param ... optional pass through parameters to pass to print.rfsrc() 
#'
#' @return summary to console
#' 
#' @export
#'
#'
print.rf_tune = function(x, ...) {
  rf_tuned = x$rf_tuned 
  print( rf_tuned ) 
}

################################################################################

