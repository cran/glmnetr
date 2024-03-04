################################################################################
#' Plot nested cross validation fit summaries
#' 
#' @description 
#' This function plots summary information from a nested.glmnetr() output object, that 
#' is from a nested cross validation performance.  Alternamvely one can output the 
#' numbers otherwise displayed to a list for extration or customized plotting. Performance
#' mesures for plotting include "agree" a measure of agreement, "lincal" the slope 
#' from a linear calibration, "intcal" the intercept from a lienar calibration, 
#' and "devrat" the deviance ratio, i.e. the fractional reduction in deviances 
#' relative compared to the null model deviances, Performance measure estimates 
#' from the individual (outer) cross validation fold are depicted by thin lines 
#' of different colors and styles, while the composite value from all folds is 
#' depicted by a thicker black line, and the performance measures naively 
#' calculated on the all data using the model derived from all data is 
#' depicted by a thicker red line. 
#' 
#' @param x A nested.glmnetr output object
#' @param type determines what type of nested cross validation performance measures are 
#' plotted.  Possible values are "agree" to plot agreement, "lincal" to plot the linear
#' calibration slope coefficients, "intcal" to plot the linear calibration intercept 
#' coefficients or "devrat" to plot the model reduction in deviances relative to 
#' the null model deviances, from the nested cross 
#' validation. 
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
#' @param pch A number indicating which symbol is to be used in place of the lines 
#' when fold=0. 
#' @param plot By default 1 to produce a plot, 0 to return the data used in the 
#' plot in the form of a list. 
#' 
#' @return This program returns a plot to the graphics window by default, and returns
#' a list with data used in teh plots if the plot=1 is specified.
#' 
#' @seealso
#'   \code{\link{plot.glmnetr}} , \code{\link{plot.cv.glmnetr}} , \code{\link{plot.nested.glmnetr}} 
#'
#' @author Walter Kremers (kremers.walter@mayo.edu)
#' 
#' @importFrom stats interaction.plot 
#' 
#' @export
#' 
plot_perf_glmnetr = function( x, type="agree", pow=2, ylim=1, fold=1, pch=20, plot=1 ) {

  object = x 
  
  if (is.null(type)) { type == "agree" }
  if (is.null(fold)) { fold = 1 
  } else if (fold != 0) { fold = 1 }

  family = object$sample[1] 
  
  if (is.null(pow)) { 
    if (family != "gaussian") { pow = 1
    } else {pow = 2} 
  } 

  if (type != "agree") { pow = 1} 
  
  if (type == "agree") { 
    if (family != "gaussian") { 
      pow = 1 
      ylab = "NCV Concordance"
    } else if (pow == 1) {
      ylab = "NCV Correlation"  
    } else if (pow==2) {
      ylab = "NCV R-square"  
    } else {
      ylab = paste0("Correlation^",pow)  
    }
  } else if (type == "intcal") { ylab = "NCV Calibration Intercept"
  } else if (type == "lincal") { ylab = "NCV Calibration Slope"
  } else if (type == "devian") { ylab = "NCV Deviance" 
  } else if (type == "devrat") { ylab = "NCV Deviance Reduction Fraction" 
  } else { ylab = "None" }
  
  ensemble = object$ensemble 
  ensemble
  
  m2.ll.null = object$null.m2LogLik.cv
  m2.ll.sat  = object$sat.m2LogLik.cv
  n__ = object$n.cv
  
  null.dev = as.numeric(object$sample[6])
  sat.m2LogLik = as.numeric(object$sample[8])
  null.m2LogLik = null.dev - sat.m2LogLik

  toplot1= c(NA,NA,NA)
  nms = NA
  frst = 1 
  i_ = 0 
  
  if (object$fits[1] == 1) {
    if         (type == "agree") {  object1 = object$lasso.agree.cv 
    } else if (type == "lincal") {  object1 = object$lasso.lincal.cv 
    } else if (type == "intcal") {  object1 = object$lasso.intcal.cv 
    } else if (type == "devian") {  object1 = object$lasso.devian.cv 
    } else if (type == "devrat") {
      object1 = object$lasso.devian.cv
      devratl = list()
      for (j_ in c(1:7)) {
        devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) ; object1[,j_] = devratl[[j_]][[1]]
      }
    }
    nfold = dim(object1)[1]
    toplot1= cbind( rep(1,nfold) , c(1:nfold), object1[,2]) 
    toplot1= rbind( toplot1, cbind( rep(2,nfold) , c(1:nfold), object1[,4]) )
    toplot1= rbind( toplot1, cbind( rep(3,nfold) , c(1:nfold), object1[,6]) )                
    toplot1= rbind( toplot1, cbind( rep(4,nfold) , c(1:nfold), object1[,7]) )                                

    if        (type ==  "agree") {  object2 = object$lasso.agree.naive 
    } else if (type == "lincal") {  object2 = object$lasso.lincal.naive 
    } else if (type == "intcal") {  object2 = object$lasso.intcal.naive 
    } else if (type == "devian") {  object2 = object$lasso.devian.naive
    } else if (type == "devrat") {
      object2 = object$lasso.devian.naive
      for (j_ in c(1:7)) {
        object2[j_] = devrat_(object2[j_], null.m2LogLik, sat.m2LogLik, 1 )[[2]]   
      }
    }
    toplot2 = cbind( 1 , 1, object2[2]) 
    toplot2 = rbind( toplot2, cbind( 2 , 1, object2[4]) )
    toplot2 = rbind( toplot2, cbind( 3 , 1, object2[6]) )                
    toplot2 = rbind( toplot2, cbind( 4 , 1, object2[7]) )                                

    object3 = colMeans(object1,na.rm=T)      
    if (type == "devrat") {
      for (j_ in c(1:7)) {
        object3[j_] = devratl[[j_]][[2]]
      }
    }
    toplot3 = cbind( 1 , 1, object3[2]) 
    toplot3 = rbind( toplot3, cbind( 2 , 1, object3[4]))
    toplot3 = rbind( toplot3, cbind( 3 , 1, object3[6]))                
    toplot3 = rbind( toplot3, cbind( 4 , 1, object3[7]))                                
    
    nms = c("Lasso", "Relaxed", "G0 Relax", "Ridge")
    i_ = 4 
  }
  
  toplot2
  
  ##### XGB ####################################################################
  if (object$fits[2] == 1) {
    
    if         (type == "agree") {  object1 = object$xgb.agree.cv 
    } else if (type == "lincal") {  object1 = object$xgb.lincal.cv 
    } else if (type == "intcal") {  object1 = object$xgb.intcal.cv 
    } else if (type == "devian") {  object1 = object$xgb.devian.cv 
    } else if (type == "devrat") {
      object1 = object$xgb.devian.cv
      devratl = list()
      for (j_ in c(1:6)) {
        devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) ; object1[,j_] = devratl[[j_]][[1]]
      }
    }
    nfold = dim(object1)[1]    
    
    if        (type ==  "agree") {  object2 = object$xgb.agree.naive 
    } else if (type == "lincal") {  object2 = object$xgb.lincal.naive 
    } else if (type == "intcal") {  object2 = object$xgb.intcal.naive 
    } else if (type == "devian") {  object2 = object$xgb.devian.naive 
    } else if (type == "devrat") {
      object2 = object$xgb.devian.naive
      for (j_ in c(1:6)) {
        object2[j_] = devrat_(object2[j_], null.m2LogLik, sat.m2LogLik, 1 )[[2]]   
      }
    }
    
    object3 = colMeans(object1,na.rm=T)      
    if (type == "devrat") {
      for (j_ in c(1:6)) {
        object3[j_] = devratl[[j_]][[2]]
      }
    }
  
    if ( sum(ensemble[c(1,5)]) >=1 ) { 
      i_= i_ + 1  
      toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,4]) )
      toplot2 = rbind( toplot2, cbind( i_ , 1, object2[4]) )
      toplot3 = rbind( toplot3, cbind( i_ , 1, object3[4]) )
      nms[i_] = "XGB"
    } 
    
    if ( sum(ensemble[c(2,6)]) >=1 ) { 
      i_ = i_ + 1 
      toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,5]) )
      toplot2 = rbind( toplot2, cbind( i_ , 1, object2[5]) )
      toplot3 = rbind( toplot3, cbind( i_, 1, object3[5]) )
      nms[i_] = "XGB feat."
    }
    
    if ( sum(ensemble[c(3,4,7,8)]) >=1 ) {   
      i_ = i_ + 1
      toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,6]) )
      toplot2 = rbind( toplot2, cbind( i_ , 1, object2[6]) )
      toplot3 = rbind( toplot3, cbind( i_, 1, object3[6]) )
      nms[i_] = "XGB offs."
    }
    frst = 0 
  }
  
  toplot2
  
  ##### Random Forest ##########################################################
  if (object$fits[3] ==1) {
    
    if         (type == "agree") {  object1 = object$rf.agree.cv 
    } else if (type == "lincal") {  object1 = object$rf.lincal.cv 
    } else if (type == "intcal") {  object1 = object$rf.intcal.cv 
    } else if (type == "devian") {  object1 = object$rf.devian.cv
    } else if (type == "devrat") {
      object1 = object$rf.devian.cv
      devratl = list()
      for (j_ in c(1:3)) {
        devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) ; object1[,j_] = devratl[[j_]][[1]]
      }
    }
    nfold = dim(object1)[1]    
    
    if        (type ==  "agree") {  object2 = object$rf.agree.naive 
    } else if (type == "lincal") {  object2 = object$rf.lincal.naive 
    } else if (type == "intcal") {  object2 = object$rf.intcal.naive 
    } else if (type == "devian") {  object2 = object$rf.devian.naive 
    } else if (type == "devrat") {
      object2 = object$rf.devian.naive
      for (j_ in c(1:3)) {
        object2[j_] = devrat_(object2[j_], null.m2LogLik, sat.m2LogLik, 1 )[[2]]   
      }
    }
    
    object3 = colMeans(object1,na.rm=T)      
    if (type == "devrat") {
      for (j_ in c(1:3)) {
        object3[j_] = devratl[[j_]][[2]]
      }
    }
    
    if ( sum(ensemble[c(1,5)]) >=1 ) { 
      i_ = i_ + 1
      toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,1]) )
      toplot2 = rbind( toplot2, cbind(  i_ , 1, object2[1]) )
      toplot3 = rbind( toplot3, cbind(  i_ , 1, object3[1]) )
      nms[i_] = "RF"
    }
    if ( sum(ensemble[c(2,6)]) >=1 ) { 
      i_ = i_ + 1
      toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,2]) )
      toplot2 = rbind( toplot2, cbind(  i_ , 1, object2[2]) )
      toplot3 = rbind( toplot3, cbind(  i_ , 1, object3[2]) )  
      nms[i_] = "RF feat."
    }
    if ((sum(ensemble[c(3,4,7,8)]) >=1) & (family == "gaussian")) { 
      i_ = i_ + 1
      toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,3]) )
      toplot2 = rbind( toplot2, cbind( i_ , 1, object2[3]) )
      toplot3 = rbind( toplot3, cbind( i_ , 1, object3[3]) )
      nms[i_] = "RF offs."
    }
    frst = 0 
  }
  
  toplot2
  
  ##### Neural Network #########################################################
  if (object$fits[5] == 1) {

    if        (type == "agree" ) {  object1 = object$ann.agree.cv 
    } else if (type == "lincal") {  object1 = object$ann.lincal.cv 
    } else if (type == "intcal") {  object1 = object$ann.intcal.cv 
    } else if (type == "devian") {  object1 = object$ann.devian.cv 
    } else if (type == "devrat") {
      object1 = object$ann.devian.cv
      devratl = list()
      for (j_ in c(1:8)) {
        devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) ; object1[,j_] = devratl[[j_]][[1]]
      }
    }
    nfold = dim(object1)[1]    
    
    if        (type ==  "agree") {  object2 = object$ann.agree.naive 
    } else if (type == "lincal") {  object2 = object$ann.lincal.naive 
    } else if (type == "intcal") {  object2 = object$ann.intcal.naive 
    } else if (type == "devian") {  object2 = object$ann.devian.naive 
    } else if (type == "devrat") {
      object2 = object$ann.devian.naive
      for (j_ in c(1:8)) {
        object2[j_] = devrat_(object2[j_], null.m2LogLik, sat.m2LogLik, 1 )[[2]]   
      }
    }
    
    object3 = colMeans(object1,na.rm=T)      
    if (type == "devrat") {
      for (j_ in c(1:8)) {
        object3[j_] = devratl[[j_]][[2]]
      }
    }
    
    for (j_ in c(1:8)) {
      if ( ensemble[j_]==1 ) { 
        i_ = i_ + 1
        toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,j_]) )
        toplot2 = rbind( toplot2, cbind( i_ , 1, object2[j_]) )    
        toplot3 = rbind( toplot3, cbind( i_ , 1, object3[j_]) )
        nms[i_] = paste0("ANN ", j_) 
      }
    }
  }
  
  toplot2
  
  ##### RPART ##################################################################
  if (object$fits[4] == 1) {

    if        (type == "agree" ) {  object1 = object$rpart.agree.cv 
    } else if (type == "lincal") {  object1 = object$rpart.lincal.cv 
    } else if (type == "intcal") {  object1 = object$rpart.intcal.cv 
    } else if (type == "devian") {  object1 = object$rpart.devian.cv 
    } else if (type == "devrat") {
      object1 = object$rpart.devian.cv
      devratl = list()
      for (j_ in c(1:9)) {
        devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) ; object1[,j_] = devratl[[j_]][[1]]
      }
    }
    nfold = dim(object1)[1]    
    
    if        (type ==  "agree") {  object2 = object$rpart.agree.naive 
    } else if (type == "lincal") {  object2 = object$rpart.lincal.naive 
    } else if (type == "intcal") {  object2 = object$rpart.intcal.naive 
    } else if (type == "devian") {  object2 = object$rpart.devian.naive 
    } else if (type == "devrat") {
      object2 = object$rpart.devian.naive
      for (j_ in c(1:9)) {
        object2[j_] = devrat_(object2[j_], null.m2LogLik, sat.m2LogLik, 1 )[[2]]   
      }
    }
    
    object3 = colMeans(object1,na.rm=T)      
    if (type == "devrat") {
      for (j_ in c(1:9)) {
        object3[j_] = devratl[[j_]][[2]]
      }
    }
      
    if ( sum(ensemble[c(1,5)]) >=1 ) { 
      i_ = i_ + 1 
      toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,1]) )
      toplot2 = rbind( toplot2, cbind( i_ , 1, object2[1]) )
      toplot3 = rbind( toplot3, cbind( i_ , 1, object3[1]) )
      nms[i_] = "RPART .00"
      i_ = i_ + 1     
      toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,2]) )
      toplot2 = rbind( toplot2, cbind( i_ , 1, object2[2]) )
      toplot3 = rbind( toplot3, cbind( i_ , 1, object3[2]) )
      nms[i_] = "RPART .01"
      i_ = i_ + 1     
      toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,3]) )
      toplot2 = rbind( toplot2, cbind( i_ , 1, object2[3]) )
      toplot3 = rbind( toplot3, cbind( i_ , 1, object3[3]) )
      nms[i_] = "RPART .02"
    }
    if ( sum(ensemble[c(2,6)]) >=1 ) { 
      i_ = i_ + 1 
      toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,4]) )
      toplot2 = rbind( toplot2, cbind( i_ , 1, object2[4]) )
      toplot3 = rbind( toplot3, cbind( i_ , 1, object3[4]) )
      nms[i_] = "RP Fe .00"
      i_ = i_ + 1     
      toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,5]) )
      toplot2 = rbind( toplot2, cbind( i_ , 1, object2[5]) )
      toplot3 = rbind( toplot3, cbind( i_ , 1, object3[5]) )
      nms[i_] = "RP Fe .01"
      i_ = i_ + 1     
      toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,6]) )
      toplot2 = rbind( toplot2, cbind( i_ , 1, object2[6]) )
      toplot3 = rbind( toplot3, cbind( i_ , 1, object3[6]) )
      nms[i_] = "RP Fe .02"
    }
    
    if ((sum(ensemble[c(3,4,7,8)]) >= 1) & (!(family %in% c("binomial")))) { 
      i_ = i_ + 1 
      toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,7]) )
      toplot2 = rbind( toplot2, cbind( i_ , 1, object2[7]) )
      toplot3 = rbind( toplot3, cbind( i_ , 1, object3[7]) )
      nms[i_] = "RP Of .00"
      i_ = i_ + 1     
      toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,8]) )
      toplot2 = rbind( toplot2, cbind( i_ , 1, object2[8]) )
      toplot3 = rbind( toplot3, cbind( i_ , 1, object3[8]) )
      nms[i_] = "RP Fe .01"
      i_ = i_ + 1     
      toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,8]) )
      toplot2 = rbind( toplot2, cbind( i_ , 1, object2[9]) )
      toplot3 = rbind( toplot3, cbind( i_ , 1, object3[9]) )
      nms[i_] = "RP Of .02"
    }
  }
  
  toplot2
  
  ##### Step Wise ##############################################################
  if ( (object$fits[6] == 1) | (object$fits[7] == 1) ) {
    object1 = object$step.agree.cv
    object2 = object$step.agree.naive
    object3 = object$step.agree.cv
    object3 = colMeans(object3,na.rm=T)
    
    if        (type == "agree" ) {  object1 = object$step.agree.cv 
    } else if (type == "lincal") {  object1 = object$step.lincal.cv 
    } else if (type == "intcal") {  object1 = object$step.intcal.cv 
    } else if (type == "devian") {  object1 = object$step.devian.cv
    } else if (type == "devrat") {
      object1 = object$step.devian.cv
      devratl = list()
      for (j_ in c(1:3)) {
        devratl[[j_]] = devrat_(object1[,j_], m2.ll.null, m2.ll.sat, n__ ) ; object1[,j_] = devratl[[j_]][[1]]
      }
    }
    nfold = dim(object1)[1]    
    
    if        (type ==  "agree") {  object2 = object$step.agree.naive 
    } else if (type == "lincal") {  object2 = c(1,1,1)
    } else if (type == "intcal") {  object2 = c(0,0,0)
    } else if (type == "devian") {  object2 = object$step.devian.naive 
    } else if (type == "devrat") {
      object2 = object$step.devian.naive
      for (j_ in c(1:3)) {
        object2[j_] = devrat_(object2[j_], null.m2LogLik, sat.m2LogLik, 1 )[[2]]   
      }
    }
    
    object3 = colMeans(object1,na.rm=T)      
    if (type == "devrat") {
      for (j_ in c(1:3)) {
        object3[j_] = devratl[[j_]][[2]]
      }
    }
  }
    
#  if ( (object$fits[6] == 1) & (type %in% c("intcal","lincal"))) {
  if (object$fits[6] == 1) {
    i_ = i_ + 1 
    toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,1]) )
    toplot2 = rbind( toplot2, cbind( i_ , 1, object2[1]) )
    toplot3 = rbind( toplot3, cbind( i_ , 1, object3[1]) )
    nms[i_] = "Step df"
    
    i_ = i_ + 1 
    toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,2]) )
    toplot2 = rbind( toplot2, cbind( i_ , 1, object2[2]) )
    toplot3 = rbind( toplot3, cbind( i_ , 1, object3[2]) )
    nms[i_] = "Step p"
  }
  
  toplot2
  
#  if ( (object$fits[7] == 1) & (type %in% c("intcal","lincal"))) {
  if  (object$fits[7] == 1) {
    i_ = i_ + 1 
    toplot1= rbind( toplot1, cbind( rep(i_,nfold) , c(1:nfold), object1[,3]) )
    toplot2 = rbind( toplot2, cbind( i_ , 1, object2[3]) )
    toplot3 = rbind( toplot3, cbind( i_ , 1, object3[3]) )
    toplot3
    nms[i_] = "AIC"
    aic_which = i_ 
  }
  
  toplot2
  
  ################################################################################
  ################################################################################
#  plot(toplot1[,1], toplot1[,3], col=toplot1[,2], pch=20, ylab="Correlation", xlab="Model", ylim=c(0,1))
#  plot(toplot1[,1], (toplot1[,3])^2, col=toplot1[,2], pch=20, ylab="R-square", xlab="Model", ylim=c(0,1))
  
#  interaction.plot(toplot1[,1], toplot1[,2], (toplot1[,3])^pow, pch=20, col=toplot1[,2], 
#                   axes=F, ylab=ylab, xlab="Model",legend=F, ylim=c(0,1))

  if (is.null(ylim)) {  
    ylim = 1 
  } 
  if (length(ylim) == 1) {
    if (ylim == 0) { 
      limridge = 0
    } else { 
      ylim = 1 
      limridge = 1  
    }
  }
  
  if (length(ylim) != 2) {
    if (type == "agree") { 
      if (family == "gaussian") {
        ylim = c(0,1) 
      } else {
        ylim = c(0.5,1) 
      }
    } else if (type == "intcal") { 
      ylim = range(rbind(toplot1[,3],toplot2[,3],toplot3[,3]),na.rm=T)
      if (limridge == 1) {
        ylim2 = range(rbind(toplot1[(toplot1[,1]!=4),3],toplot2[ (toplot2[,1] != 4) ,3], toplot3[ (toplot3[,1]!=4) ,3]),na.rm=T) 
        ylim[2] = min( ylim[2], 2 * ylim2[2])
      }
    } else if (type == "lincal") { 
      ylim = range(rbind(toplot1[,3],toplot2[,3],toplot3[,3]), na.rm=T) 
      if (limridge == 1) {
        ylim2 = range(rbind(toplot1[(toplot1[,1]!=4),3],toplot2[ (toplot2[,1] != 4) ,3], toplot3[ (toplot3[,1]!=4) ,3]),na.rm=T) 
        ylim[2] = min( ylim[2], 2 * ylim2[2])
      }
      ylim[2] = max(ylim,1)
    } else if (type == "devian") { 
      if (family == "cox") {
        ylim = range(rbind(toplot1[,3],toplot3[,3]),na.rm=1) 
      } else {
        ylim = range(rbind(toplot1[,3],toplot2[,3],toplot3[,3]),na.rm=1) 
      }
      if ((limridge == 1) & (object$fits[7] == 1)) {
        ylim2 = range(rbind(toplot1[(toplot1[,1] != aic_which),3],toplot2[ (toplot2[,1] != aic_which) ,3], toplot3[ (toplot3[,1] != aic_which) ,3]),na.rm=1) 
        ylim[1] = max( ylim[1], 0.8 * ylim2[1])
        
        if (family == "cox") {
          ylim[2] = min( ylim[2], 1   * ylim2[2])
        } else {
          ylim[2] = min( ylim[2], 1.2 * ylim2[2])
        }
      }
#      ylim[1] = 0 ; 
    } else if (type == "devrat") {
      ylim = range(rbind(toplot1[,3],toplot2[,3],toplot3[,3]),na.rm=1) 
      ylim[1] = 0 
    } 
  }
  
#  print(ylim)

  if ((family == "cox") & (type == "intcal")) {
    cat(paste(" For the Cox model there is no intercept term and all calibration intercpet\n",
              " coefficients are by convention 0\n"))
  } else {
    if (plot != 0) {
      if (fold == 1) {
      interaction.plot(toplot1[,1], toplot1[,2], (toplot1[,3])^pow, pch=pch, col=toplot1[,2], 
                       ylab=ylab, xlab="",legend=F, ylim=ylim,
                       xaxt = "n")
      axis(1, at = c(1:length(nms)), labels = nms, las=2)
      if (!((type %in% c("devian","devrat")) & (family == "cox"))) { 
        lines(toplot2[,1], toplot2[,3]^pow,col=2,lwd=2) # pch=6,
      }
      lines(toplot3[,1], toplot3[,3]^pow,col=1,lwd=2) # pch=6,
      if (type == "lincal") { abline(h=1, lty=3, col=1, lwd=2) }
      if ((type == "intcal") | (type == "devrat")) { abline(h=0, lty=3, col=1, lwd=2) }
      } else {
        plot( toplot3[,1], (toplot3[,3])^pow, pch=pch, col=toplot3[,2], 
                         ylab=ylab, xlab="",ylim=ylim,
                         xaxt = "n")
        axis(1, at = c(1:length(nms)), labels = nms, las=2)
      }
    }
    if (plot %in% c(0,2)) { return( list(toplot1=toplot1, toplot2=toplot2, toplot3=toplot3, nms=nms) ) }
  }
  
}


