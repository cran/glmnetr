################################################################################
##### ObliqueRandomForest_yymmdd.R #############################################
################################################################################
#' Fit a Random Forest model on data provided in matrix and vector formats.
#' 
#' @description Fit an Random Forest model using the orsf() function of the 
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
#' Random Forest fit.  This if for the mtry input variable of the orsf() program 
#' specifying the number of terms to consider in each step of teh Random Forest fit.
#' @param ntreec a vector (numeric) of 2 values, the first for the number of forests
#' (ntree from orsf()) to use when searhcing for a better bit and the second to use
#' when fitting the final model.  More trees should give a better fit but 
#' require more computations and storage for the final. 
#' model.   
#' @param nsplitc This nsplit of orsf(), a non-negative integer for the  number of 
#' random splits for a predictor.
#' @param seed a seed for set.seed() so one can reproduce the model fit.  If 
#' NULL the program will generate a random seed.  Whether specified or NULL, the 
#' seed is stored in the output object for future reference.  Note,
#' for the default this randomly generated seed depends on the seed in memory at that 
#' time so will depend on any calls of set.seed prior to the call of this function. 
#' @param tol a small number, a lower bound to avoid division by 0  
#' @param track 1 to output a brief summary of the final selected model, 2 to 
#' output a brief summary on each model fit in search of a better model or 0 
#' (default) to not output this information.
#'
#' @return a Random Forest model fit 
#' 
#' @seealso 
#'    \code{\link{summary.orf_tune}} , \code{\link{rederive_orf}} , \code{\link{nested.glmnetr}}
#    , \code{\link{predict_nested_orf}}    
#' 
#' @author Walter Kremers (kremers.walter@mayo.edu)
#' 
#' @importFrom aorsf orsf 
#'
#' @export
#' 
# xs, start=NULL, y_, event=NULL, family=NULL, mtryc=NULL, ntreec=NULL, nsplitc=8, seed=NULL, track=0 
# xs, start=NULL, y_, event=NULL, family=NULL, mtryc=NULL, ntreec=NULL, nsplitc=8, seed=NULL, track=0 
# xs=hp.xs ; start=NULL ; y_=hp.tm2deep_infection ; event=hp.deep_infection ; family="cox" ;
# mtryc=NULL ; ntreec=NULL ; nsplitc=8 ; seed=NULL ; track=1 ;

orf_tune = function(xs, start=NULL, y_, event=NULL, family=NULL, mtryc=NULL, ntreec=NULL, nsplitc=8, seed=NULL, tol=1e-5, track=0) {
  
  if (is.matrix(y_)) { y_ = as.numeric(y_) }
  
  if (is.null(seed)) { seed = round(runif(1)*1e9) } 
  set.seed( seed ) 
  
  if (is.null(family) == 1) {
    if (!is.null(event) == 1) {
      family = "cox"
    } else if ( (length(table(y_)) == 2) & (names(table(y_))[1] == 0) & (names(table(y_))[2] == 1) ) { 
      family = "binomial" 
    } else { 
      family = "gaussian"
    } 
    print(family) 
  }
  
  if        (family == "binomial") { 
#    y_ = as.factor(y_)
    class(y_) 
    table(y_)
    split_rule="gini"  
#    split_rule="variance"                                                       ## this should not be !!
  } else if (family == "gaussian") { split_rule="variance"
  } else if (family == "cox"     ) { split_rule="logrank" } 
  
  if (track >= 2) { 
    time_start = diff_time()
    time_split = time_start
  } 
  
  xs = xs[,(diag(cov(xs)) != 0)]
  
  if (family == "cox") {
    df = as.data.frame(cbind(xs, y_=y_, event = event)) 
  } else if (family=="binomial") {
#    df = as.data.frame(cbind(xs, y_)) 
    y_ = factor(y_)
    df = as.data.frame(cbind(as.data.frame(xs), y_)) 
#    print(class(df$y_)) 
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
    ntreec = c(500,  500)
    if (family == "gaussian") { ntreec = c(100, 500) 
    } else { ntreec = c(100, 500) }
  }
  
  if (is.null(nsplitc)) {
    nsplitc = c(8)
  } else if ((is.na(nsplitc[1])) | (nsplitc[1] < 2)) {
    nsplitc[1] = 8 
  }
  
  ##------------------------------------------------------------------------------
  
  ##mtryc = mtryc[1:3]
  ## k_ = 1 
  
  #  rfsrc does not analyze (start,stop) time survival data == How about rosf ??
  if ((family == "cox") & (!is.null(start)==1)) { skip = 1 } else { skip = 0 } 
  
  if (skip == 0) {
    
    for ( k_ in c(1:length(mtryc))) {
      set.seed( seed ) 
      if (family == "cox") {
#        rf = rfsrc( Surv(y_, event) ~ . , data=df, mtry=mtryc[k_], nsplit=nsplitc[1], ntree=ntreec[1], 
#                    splitrule=splitrule, membership = TRUE, importance=TRUE) 
        orf = orsf( Surv(y_, event) ~ . , data=df, mtry=mtryc[k_], n_split=nsplitc[1], n_tree=ntreec[1], split_rule=split_rule) 
      } else {
        orf = orsf( y_ ~ . , data=df, mtry=mtryc[k_], n_split=nsplitc[1], n_tree=ntreec[1], split_rule=split_rule) 
      }
      
      if (family == "cox") { y__ = Surv(y_,event) } else { y__ = y_ } 
      orf_perf1  = orf_perform(orf_model=orf, dframe=as.data.frame(xs), ofst=NULL, y__=y__, family=family, tol=tol)

      if (k_ ==1) {
        k_best = 1 
        mtry_best = mtryc[k_] 
        oobag_fun = orf$eval_oobag$stat_values 
        oobag_fun.best = oobag_fun
        oobag_funv = oobag_fun
        orf_best = orf 
      }  else if (orf$eval_oobag$stat_values > oobag_fun.best) {                # C, AUC-ROC or R-square
        k_best = k_ 
        mtry_best = mtryc[k_] 
        oobag_fun.best = orf$eval_oobag$stat_values
        orf_best = orf 
      }
      
      if (k_ > 1) { oobag_funv[k_] = orf$eval_oobag$stat_values } 
      
      if (track >= 2) { 
        print( c(k_, mtryc[k_], mtry_best, orf$eval_oobag$stat_values ) )
        time_split = diff_time(time_start, time_split)
      }
    }
    
    if (ntreec[2] > ntreec[1]) {
      set.seed( seed ) 
      if (family == "cox") {
#        rf_tuned = rfsrc( Surv(y_, event) ~ . , data=df, mtry=mtryc[k_best], nsplit=nsplitc[1], ntree=ntreec[2], 
#                   splitrule=splitrule, membership = TRUE, importance=TRUE) 
        orf_tuned = orsf( Surv(y_, event) ~ . , data=df, mtry=mtryc[k_best], n_split=nsplitc[1], n_tree=ntreec[2], split_rule=split_rule) 
      } else {
        orf_tuned = orsf( y_ ~ . , data=df, mtry=mtryc[k_best], n_split=nsplitc[1], n_tree=ntreec[2], split_rule=split_rule) 
      }
    } else {
      orf_tuned = orf_best
    }
    
    if (track >= 2) {
      print( c(mtry_best, orf_tuned$eval_oobag$stat_values ) )
      time_split = diff_time(time_start, time_split)
    }
    
#    orffit = list(orf_tuned=orf_tuned, err.ratev=err.ratev, err.rate=orf_tuned$err.rate[length(orf_tuned$err.rate)], 
#                 mtryc=mtryc, ntreec=ntreec, seed=seed )

#    orf_tuned$err.ratev=err.ratev    
    orffit = list(orf_tuned=orf_tuned, oobag_funv=oobag_funv, mtryc=mtryc, ntreec=ntreec, seed=seed )

  } else { orf_tuned = list(orffit="NONE") }
  
  class(orffit) <- c("orf_tune")
  
  return(orffit)
  
}

################################################################################
#' 
#' Summarize output from rf_tune() function
#'
#' @param object output from an rf_tune() function
#' @param ... optional pass through parameters to pass to summary.orsf() 
#'
#' @return summary to console
#' 
#' @seealso 
#    \code{\link{predict_nested_rf}} ,
#'    \code{\link{rf_tune}} , \code{\link{nested.glmnetr}}
#' 
#' @export
#'
summary.orf_tune = function(object, ...) {
  cat(paste0("\n    search set for mtry : "))
  cat(object$mtryc)
  cat(paste0("\n    search set for ntree : "))
  cat(object$ntreec)
  cat(paste0("\n    selected value for mtry from search: ", object$orf_tuned$mtry, "\n\n"))
  orf_tuned = object$orf_tuned 
  print(orf_tuned) 
}

################################################################################
#'
#' Print output from orf_tune() function
#' 
#' @param x output from an orf_tune() function
#' @param ... optional pass through parameters to pass to print.orf() 
#'
#' @return summary to console
#' 
#' @seealso 
#    \code{\link{orf_xbhat}} , \code{\link{orf_perform}} ,   
#'    \code{\link{summary.orf_tune}} , \code{\link{orf_tune}} , \code{\link{nested.glmnetr}}
#' 
#' @export
#'
#'
print.orf_tune = function(x, ...) {
  orf_tuned = x$orf_tuned 
  print( orf_tuned ) 
}

################################################################################

