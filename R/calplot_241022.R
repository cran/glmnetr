################################################################################
 
#  library(glmnetr)
#  library(splines)
#  library(survival)
 
################################################################################
################################################################################

#' Construct a single calibration plot for a nested.glmnetr output object 
#'
#' @param object1 A nested.glmnetr() output object for calibration
#' @param wbeta  Which Beta should be plotted, an integer.  This will depend on
#' which machine learning models were run when creating the output object.   If 
#' unsure the user can run the function without specifying wbeta and a legend 
#' will be directed to the console.  
#' @param df The degrees of freedom for the spline function, default=5. 
#' @param resample 1 to base the splines on the leave out X*Beta's ($xbetas.cv), or 0 to 
#' use the naive X*Beta's ($xbetas).
#' @param oob 1 (default) to construct calibration plots using the out of bag 
#' data points, 0 to use inbag (including resampled data points) data 
#' points. Only applies when bootstrap is used instead of k-fold cross 
#' validation. 
#' @param bootci Calculate bootstrap confidence ientervals for calibration 
#' curves. This is for exploration only. The applicability of bootstrap 
#' confidence intervals for these calibration curves is questionable. 
#' @param trim the percent of top and bottom of the data to be trimmed away when 
#' producing plots.
#' @param plot 1 by default to plot the calibration
#' @param lwd line width for plots, 1 by default
#' @param title title for inclusion on plots, default NULL 
#' @param xlim xlim for the plots.  This does not effect the curves within 
#' the plotted region.  Caution, for the "cox" framework the xlim are specified
#' in terms of the X*beta and not the HR, even when HR is described on the axes.
#' @param ylim ylim for the plots, which will usually only be specified in a 
#' second run of for the same data.  This does not effect the curves within 
#' the plotted region.  Caution, for the "cox" framework the ylim are specified
#' in terms of the X*beta and not the HR, even when HR is described on the axes.
#' @param xlab a user specified label for the x axis
#' @param ylab a user specified label for the y axis
#' @param col.term a number for the line depicting the overall calibration estimates
#' @param col.se a number for the line depicting the +/- 2 * standard error 
#' lines for the overall calibration estimates 
#' @param rug 1 to plot a rug for the model x*betas, 0 (default) to not.
#' @param plothr  a power > 1 determining the spacing of the values 
#' on the axes, e.g. 2, exp(1), sqrt(10) or 10.  The default of 0 plots the 
#' X*Beta.  This only applies for "cox" survival data models.   
#' @param repn The replication or repeat number for (nested) CV or bootstrap.
#' @param knottype 1 (default) to use XBeta used for the spline fit to 
#' choose ns() knots for gaussian and binomial families, 2 to use the XBeta
#' from all re-samples to determine the knots.  
#' @param ... allowance to pass terms to the invoked plot function
#'
#' @return Plots are returned by default, and optionally data for plots are 
#' output to a list.
#' 
#' @importFrom stats sd quantile 
#' @importFrom utils str 
#' @importFrom graphics box rug title
#' @importFrom survival coxph pspline
#' @importFrom splines ns
#' 
#' @author Walter Kremers (kremers.walter@mayo.edu)
#' 
#' @noRd

calplot0 = function(object1, wbeta=NULL, df=5, resample=1, oob=1, bootci=0, repn=NULL, knottype=1,
                    trim=0, plot=1, lwd=1, title=NULL, xlim=NULL, ylim=NULL, xlab=NULL, ylab=NULL, 
                    col.term=1, col.se=2, rug=1, plothr=0, ... ) {
#  print(c(100, oob ))  
  ## resample replaces cv 
  ## repn replaces fold 
  #  object1 = object ; wbeta = 5 ; df = 3 ; resample = 1 ; trim = 0 ; repn = 1 ; plot = 1 ; lwd = 1 ; title=NULL ; xlim=NULL ; ylim=NULL ; xlab=NULL ; ylab=NULL ; col.term=1 ; col.se=2 ; rug=1 ; plothr=0 ; knottype=0 ; 
  #  plot=plotfold_ ; title=paste("repn",i_) ; 
  
#  object1 = object ; repn=NULL ; knottype=1 ; lwd=2 ;                           ## naive calibration plot test 
  
  family = object1$sample[1] 
  bootstrap =  as.numeric( object1$tuning[8] ) 
  folds_n = as.numeric( object1$tuning[1] )
  
  if (!(knottype %in% c(-1,0,1,2,3))) { knottype = 1 } 
#  if ( bootstrap == 0 ) { bootci = 0 }
#  if ( bootstrap == 1) {
#    if ( (bootci == 1) & (oob == 1) ) { cat("  For (bootci == 1) & (oob == 1), oob is set to 0") ; oob = 0 ; }
#  }
  
#  print(c(resample,bootstrap,oob,folds_n)) 
  
  ## check that resampling was done in original analysis 
  if (resample == 1) {
    if (bootstrap == 0) {
      if ( is.null(object1$xbetas.cv) ) {
        resample = 0 
        cat(paste("  xbeta.cv is NULL so resample set to 0\n"))
      }
    }
    if (bootstrap >= 1) {
      if (is.null(object1$xbetas.boot.oob) ) {
        resample = 0 
        cat(paste("  xbeta.boot is NULL so resample set to 0\n"))
      }
    }
  }
  
  stop_ = 0 
  
  ##### get range of full model data ###########################################
  minxb = min( object1$xbetas[,wbeta ] )                                        # object1$xbetas[1:10 , wbeta, drop=FALSE ] # xb.hat.p[1:10]
  maxxb = max( object1$xbetas[,wbeta ] )
  fully_ = object1$y_                                                           # length(fully_) # fully_[1:20] # yt[1:20] # (fully_[1:20]<0)
  
#    print( data.frame(cbind( varxb, num=1:length(varxb))) )
  
##---------------------
  
#  print(c(resample,bootstrap,oob)) 
  if (resample == 1) { 
    if (bootstrap == 0) {
      origxb = object1$xbetas.cv[,wbeta ]  
      xbname = colnames(object1$xbetas.cv)[wbeta] 
    } else if (bootstrap >= 1) {
      if (oob == 1) { 
        origxb = object1$xbetas.boot.oob[,wbeta+2 ]  
        xbname = colnames(object1$xbetas.boot.oob)[wbeta+2] 
      } else {
        origxb = object1$xbetas.boot.inb[,wbeta+2 ]  
        xbname = colnames(object1$xbetas.boot.inb)[wbeta+2] 
#        print("Here 241010 1")
      }
    }
  } else { ## if (resample == 0)
    origxb = object1$xbetas[,wbeta]                                             ## cbind(origxb , xb.hat.p )[1:20,]
    xbname = colnames(object1$xbetas)[wbeta]  
  }
  
  ##---------------------

  ## select subset for single fold or bootstrap sample, if applicable ----------
  if ( (resample == 1) & (!is.null(repn)) ) {
    if (bootstrap == 0) {
      foldid = object1$foldid 
      table( foldid ) 
      xb_ = origxb[ (foldid == repn) ] 
      y_  = fully_[ (foldid == repn) ] 
    } else if (bootstrap >= 1) {
      if (oob == 1) { 
        xb_ = origxb[ object1$xbetas.boot.oob[,1]==repn ]                       ## check this 
        oob_tmp = object1$xbetas.boot.oob[object1$xbetas.boot.oob[,1]==repn,2]  ## out of bag sequence IDs 
        y_  = fully_[ oob_tmp ]                                                 ## get correpsoing y_'s 
      } else { 
        xb_ = origxb[ object1$xbetas.boot.inb[,1]==repn ]                       ## check this 
        inb_tmp = object1$xbetas.boot.inb[object1$xbetas.boot.inb[,1]==repn,2]  ## out of bag sequence IDs 
        y_  = fully_[ inb_tmp ]                                                 ## get corresponding y_'s 
      } 
    }
  }
  ##---------------------
  
  ## data for all rows -- really only makes sense for cross validation ---------
  ## and only for inspection, not intended for inference -----------------------
  if ( (resample == 1) & (is.null(repn)) ) {
    if (bootstrap == 0) {
      y_  = fully_
      xb_ = origxb
    } 
    if (bootstrap >= 1) {                              
      xb_ = origxb 
      for ( i_ in c(1:object1$tuning[8])) {
        oob_tmp = object1$xbetas.boot.oob[object1$xbetas.boot.oob[,1]==i_,2]
        y_tmp  = fully_[oob_tmp]
        xb_tmp = origxb[oob_tmp]
        if (i_ == 1) { 
          resampno = rep(1,length(y_tmp))
          y_ = y_tmp
          xb_ = xb_tmp
        } else if (i_ != 1) { 
          resampno = rep(i_,length(y_tmp))
          y_ = c(y_,y_tmp)
          xb_ = c(xb_,xb_tmp)
        } 
      }
    }
  }
  ##---------------------
 
  ##----- calibrate naively without resampling, i.e. using model based upon all data ---------------------
  ##----- uses folds but more for inspection than inference ---------
  ##                                                                            ## check, avoid for bootstrap >=1, resample = 0 and !is.null(repn)
  if ( (resample == 0) & (!is.null(repn)) ) {
    foldid = object1$foldid 
    y_  = fully_[ (foldid == repn) ]  
    xb_ = origxb[ (foldid == repn) ] 
  }
  ##---------------------
  
  ##---------------------
  if ( (resample == 0) & (is.null(repn)) ) {
    y_  = fully_                                                                # length(y_) # y_[1:20] # yt[1:20]
    xb_ = origxb                                                                ## cbind(xb_ , xb.hat.p )[1:20,]
#    print(xb_[1:10])
  }
  ##---------------------
  
  c( length( y_ ) , length( xb_ ) ) 
  
  ## get range of actual data for plotting individual spline fits 
  minxb_ = min(xb_) 
  maxxb_ = max(xb_) 

#  print(c(bootci,resample,repn))

  ## reorder data in terms of xb_ 
  ordr = order(xb_)                                                             ## table(table(ordr)) ; summary(xb_) ; summary(y_) ; 
  xb_ = xb_[ordr] 
  y_  = y_ [ordr]                                                               ## summary(xb_) ; summary(y_) ; 
#  cat(paste("3 var xb_ = ", var(xb_), "\n"))
  
  if (family == "cox") {
 #   length(table(table(xb_))) 
    lxb_ = length(unique(xb_))
    if ( lxb_ < (df+1))  { df = lxb_ - 1 }
    if (df > 1) {
      object2 = coxph(y_ ~ pspline(xb_,df=df), data=as.data.frame(cbind(y_, xb_)) )                         # termplot(object2,1) ; abline(a=0,b=1) ; title("calplot")
#      object3 = coxph(Surv(yt,event) ~ pspline(xb_,df=df), data=as.data.frame(cbind(yt, event, xb.hat.p)) ) # termplot(object3,1) ; abline(a=0,b=1) ; title("predict()")
    } else {
      object2 = coxph(y_ ~ xb_, data=as.data.frame(cbind(y_, xb_)) ) 
    }
  } else {
    if (var(xb_)>0) { 
      #        cat( paste( "knottype = " , knottype, "\n" ) )   
      if        (knottype == 1) { myspline = ns(xb_, df=df, intercept=FALSE)  
      } else if (knottype == 2) { myspline = ns(origxb, df=df, intercept=FALSE)
      } else if (knottype == 3) { myspline = ns(object1$xbetas[,wbeta ], df=df, intercept=FALSE) 
      } else if (knottype <= 0) { myspline = ns(xb_, df=df, intercept=FALSE) } 
      if (knottype<=0) { 
        print("Outside spline") 
        print(attr(myspline,c("degree")))
        print(attr(myspline,c("knots")))
        print(attr(myspline,c("Boundary.knots")))
        print(attr(myspline,c("intercept")))
      }
      if (knottype %in% c(1,2,3)) { 
        Knots = attr(myspline,"knots")
        Boundary.Knots = attr(myspline,"Boundary.knots")
        object2 = glm(y_ ~ ns(xb_,knots=Knots, Boundary.knots=Boundary.Knots) , family=family ) # data=as.data.frame(cbind(y_, xb))  
      } 
      if (knottype <= 0) { 
        object2 = glm(y_ ~ ns(xb_,df=df) , family=family ) # data=as.data.frame(cbind(y_, xb))  
        if (knottype==-1) { str(object2$model) }
      } 
    } else {
      object2 = glm(y_ ~ 1 , family=family ) # data=as.data.frame(cbind(y_, xb))  
    }
    summary(object2)
  }
  
  plotxb = minxb + rep(0:100)*(maxxb-minxb)/100
  
  ## mean(predict(object2)) ## 9e-16
  meancalxb = mean(predict(object2,data.frame(xb_=origxb)))
  #      cat( paste(title, meancalxb, "\n" ) )  ## 9e-16
  meanxb_ = mean(xb_)
  dfr = data.frame(xb_=plotxb)                   ## title="derived range" 
#  dfr = data.frame(xb_=xb_)                     ## title="original data"
  est = predict(object2, newdata=dfr)   
  predicteds = predict(object2, newdata=data.frame(xb_=plotxb), se=TRUE)  
  est = predicteds$fit
  se = predicteds$se.fit
  estci = cbind(plotxb, est, se, est - 2 * se, est + 2 * se)
#  estci = cbind(plotxb, est, se, est - qnorm(0.975) * se, est + qnorm(0.975) * se)
  
  lty1 = 1 ; lty2 = 2 ;   
  if (is.null(ylab)) {
    if (family %in% c("cox") ) {
      if (plothr > 1) {
        ylab = "Calibrated HR"
      } else {
        ylab = "Calibrated log(HR)"
      }
    } else if (family %in% c("binomial","gaussian") ) {
      ylab = "Calibrated X*Beta"
    } 
  }
  if (is.null(xlab)) { xlab = myxlab(family, plothr, resample, oob, xbname, bootstrap, bootci) }
  xlims = quantile(xb_, (c(trim,(100-trim))/100))
  #    print(xlims)
  summary(estci[,1])
  #    estci_ = estci[((estci[,1]>=xlims[1]) & (estci[,1]<= xlims[2])),]
  #    summary(estci_[,1])
  if (is.null(xlim)) { xlim = c(minxb, maxxb) }
  if (is.null(ylim)) { ylim = c(min(estci[,2]), max(estci[,2])) }
  if (plot %in% c(1,2)) {
    if ((minxb < minxb_) | (maxxb > maxxb_)) { 
      indxlo = c(1:length(estci[,1])) [ estci[,1] < minxb_ ] 
      indxhi = c(1:length(estci[,1])) [ estci[,1] > maxxb_ ] 
      if ( length( c(indxlo,indxhi) ) > 0) {
        indxfold = c(1:length(estci[,1])) [-c(indxlo,indxhi) ] 
      } else {
        indxfold = c(1:length(estci[,1]))
      }
      plot(  estci[indxlo,1], estci[indxlo,2] , type="l", 
             lty=lty2, col=col.term, lwd=lwd, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, axes=0 , ... ) 
      lines( estci[indxlo,1], estci[indxlo,4], lty=lty2, lwd=lwd, col=col.se ) 
      lines( estci[indxlo,1], estci[indxlo,5], lty=lty2, lwd=lwd, col=col.se ) 
      lines( estci[indxfold,1], estci[indxfold,2], lty=lty1, lwd=lwd, col=col.term ) 
      lines( estci[indxfold,1], estci[indxfold,4], lty=lty1, lwd=lwd, col=col.se ) 
      lines( estci[indxfold,1], estci[indxfold,5], lty=lty1, lwd=lwd, col=col.se ) 
      lines( estci[indxhi,1], estci[indxhi,2], lty=lty2, lwd=lwd, col=col.term ) 
      lines( estci[indxhi,1], estci[indxhi,4], lty=lty2, lwd=lwd, col=col.se ) 
      lines( estci[indxhi,1], estci[indxhi,5], lty=lty2, lwd=lwd, col=col.se ) 
    } else {
      plot(  estci[,1], estci[,2] , type="l", col=col.term, lwd=lwd, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, axes=0 , ... ) 
      lines( estci[,1], estci[,4] , lty=lty1, lwd=lwd, col=col.se ) 
      lines( estci[,1], estci[,5] , lty=lty1, lwd=lwd, col=col.se )
    }
    
    if (family == "cox") {
      abline(a=-meanxb_, b=1, lty=3)
      abline(h=0, lty=1)
      myaxis(ylim, plothr, est , side=2)
      myaxis(xlim, plothr, plotxb, side=1)
    } else { 
      abline(a=0, b=1, lty=3)
      axis(side=2)
      axis(side=1)
    }
    box() 
    
    keep_  = ((xb_ >= xlim[1]) & (xb_ <= xlim[2]))
    if (sum(keep_==0) > 0) {
      cat(paste("  Due to user specified xlim ", sum(keep_==0) , "tick marks are not displayed repn (replicaiton number) ", repn, ")\n"))
    }
#    print(length(xb_))
    xb_ = xb_[ keep_ ]
    y_  = y_ [ keep_ ]
#    print(length(xb_))
    if (rug == 1) { myrug(xb_, y_, family) } 
    if (!is.null(title)) { title ( title ) }
    #        termplot(object2, terms=1,se=TRUE, xlab=xlab, ylab="Calibrated log(HR)", lwd.term=lwd, lwd.se=lwd )
  }
  if (plot %in% c(0,2)) {
    return( list( estci=estci, meancalxb=meancalxb, meanxb_=meanxb_, minmax=c(minxb_,maxxb_), xb_=xb_, y_=y_ ) )
  }
}

#
################################################################################
#
#' Construct calibration plots for a nested.glmnetr output object 
#' 
#' Using k-fold cross validation this function constructs calibration plots for 
#' a nested.glmnetr output object. Each hold out subset of the k-fold cross 
#' validation is regressed on the x*beta predicteds based upon the 
#' model fit using the non-hold out data using splines.  This yields k spline 
#' functions for evaluating model performance.  These k spline functions are 
#' averaged to provide an overall model calibration.  Standard deviations of 
#' the k spline fits are also calculated as a function of the predicted X*beta, 
#' and these are used to derive and plot approximate 95% confidence intervals 
#' (mean +/- 2 * SD/sqrt(k)).  Because regression equations can be unreliable 
#' when extrapolating beyond the data range used in model derivation, we 
#' display this overall calibration fit and CIs with solid lines only for the 
#' region which lies within the ranges of the predicted x*betas for  
#' all the k leave out sets.  The spline fits are made using the same framework 
#' as in the original machine learning model fits, i.e. one of "cox", "binomial" 
#' or "gaussian"family. For the "cox" famework the pspline() funciton is used, 
#' and for the "binomial" and "gaussian" frameworks the ns() function is 
#' used.  Predicted X*betas beyond the range of any of the hold 
#' out sets are displayed by dashed lines to reflect the lessor certainty when 
#' extrapolating even for a single hold out set.   
#' 
#' Optionally, for comparison,  
#' the program can fit a spline based upon the predicted x*betas ignoring the 
#' cross validation structure, or one can fit a spline using the x*betas 
#' calculated using the model based upon all data. 
#'
#' @param object A nested.glmnetr() output object for calibration
#' @param wbeta  Which Beta should be plotted, an integer.  This will depend on
#' which machine learning models were run when creating the output object.   If 
#' unsure the user can run the function without specifying wbeta and a legend 
#' will be directed to the console.  
#' @param df The degrees of freedom for the spline function 
#' @param resample 1 to base the splines on the leave out X*Beta's ($xbetas.cv 
#' or $xbetas.boot.oob), or 0 to use the naive X*Beta's ($xbetas).  This can be done 
#' to see biases associated with the naive approach. 
#' @param oob 1 (default) to construct calibration plots using the out-of-bag 
#' data points, 0 to use in bag (including resampled data points) data 
#' points. This option only applies when bootstrap is used instead of k-fold 
#' cross validation, and when resample is set to 1. For cross validation 
#' evaluations out-of-bag samples (folds) are always used for evaluation. The 
#' purpose of oob = 0 is to allow evaluation of the variability of 
#' bootstrap calibrations ignoring bias like done in Riley et al., 2023,
#' doi: 10.1186/s12916-023-03212-y and Austin and Steyerberg 2013, 
#' doi: 10.1002/sim.5941
#' @param bootci 1 to calculate bootstrap confidence intervals for calibration 
#' curves adjusting for bias, 0 (default) to simply plot the calibration curves 
#' based upon the inbag data. This is for exploration only, and only when 
#' bootstrap samples were used for model performance evaluation. The 
#' applicability of bootstrap confidence intervals for these calibration curves 
#' is questionable. If bootci is set to 1 then oob is set to 0.   
#' @param plot 1 by default to produce plots, 0 to output data for plots only, 
#' 2 to plot and output data. 
#' @param plotfold 0 by default to not plot the individual fold calibrations, 1 
#' to overlay the k leave out spline calibration fits in a single figure and 2 
#' to produce separate plots for each of the k hold out calibration curves. 
#' @param plothr  a power > 1 determining the spacing of the values 
#' on the axes, e.g. 2, exp(1), sqrt(10) or 10.  The default of 0 plots the 
#' X*Beta.  This only applies fore "cox" survival data models. 
#' @param knottype 1 (default) to use XBeta used for the spline fit to 
#' choose knots in ns() for gaussian and binomial families, 2 to use the XBeta
#' from all re-samples to determine the knots. 
#' @param trim the percent of top and bottom of the data to be trimmed away when 
#' producing plots.  The original data are still used used calcualting the curves 
#' for plotting. 
#' @param vref Similar to trim but instead of trimming the spline lines, plots 
#' vertical refence lines aht the top vref and bottom vref percent of the model
#' X*Betas's
#' @param xlim xlim for the plots.  This does not effect the curves within 
#' the plotted region.  Caution, for the "cox" framework the xlim are specified
#' in terms of the X*beta and not the HR, even when HR is described on the axes.
#' @param ylim ylim for the plots, which will usually only be specified in a 
#' second run of for the same data.  This does not effect the curves within 
#' the plotted region.  Caution, for the "cox" framework the ylim are specified
#' in terms of the X*beta and not the HR, even when HR is described on the axes.
#' @param xlab a user specified label for the x axis
#' @param ylab a user specified label for the y axis
#' @param col.term a number for the line depicting the overall calibration estimates
#' @param col.se a number for the line depicting the +/- 2 * standard error 
#' lines for the overall calibration estimates 
#' @param rug 1 to plot a rug for the model x*betas, 0 (default) to not.
#' @param seed an integer seed used to random select the multiple of X*Betas
#' to be used in the rug when using bootstraping for model evaluation as sample 
#' elements may be included multiple times as test (Out Of Bag) data.
#' @param cv Deprecated. Use resample option instead.
#' @param fold Deprecated. This term is now ignored.
#' @param ... allowance to pass terms to the invoked plot function  
#'
#' @return Calibration plots are returned by default, and optionally data for plots 
#' are output to a list.
#' 
#' @seealso
#'   \code{\link{plot.nested.glmnetr}} , \code{\link{summary.nested.glmnetr}} , \code{\link{nested.glmnetr}} 
#'   
#' @author Walter Kremers (kremers.walter@mayo.edu)
#'
#' @export
#'
calplot = function(object, wbeta=NULL, df=3, resample=NULL, oob=1, bootci=0, 
                   plot=1, plotfold=0, plothr=0, knottype=1, trim=0, vref=0, xlim=NULL, ylim=NULL, 
                   xlab=NULL, ylab=NULL, col.term=1, col.se=2, rug=1, seed=NULL, cv=NULL, fold=NULL, ... ) {
usefold = 1   
# print(c(101, oob ))  
#  cat( paste( " knottype 1 =", knottype, "\n")) 
# wbeta=5 ; df=3 ; resample=1 ; oob=1 ; bootci=0 ; trim=0 ; vref=0 ; plot=1 ; plotfold=0 ; xlim=NULL ; ylim=NULL ; xlab=NULL ; ylab=NULL ; col.term=1 ; col.se=2 ; rug=1 ; seed=NULL ; plothr=0 ; knottype=0 ; resample=1 ; fold=NULL ; cv=NULL ; 
# usefold=1 ; plot = 2 ; 
# resample=0 ; 

# object=nested.bin.boot.fit ; wbeta=16 ; bootci=0 ; oob=0 ; plotfold=1 ;

# object = nested.bin.fit3B 
# object = nested.bin.fit10Boot 
# object=nps_csf_glo_ratcap_fit
# object=nested_glmnetr_fit_hp_las
# object = nested.cox.fit 
# object1 = object ;
  
#  object1 = object ;  lwd=2 ;  trim = 0 ; usefold=NULL ;  
#  wbeta=5 ; resample=0  ; usefold=0 ; plotfold=0 ; ylim=c(-0.2,0.1) ; 
  
#  cat(paste("HERE1",ylim[1],ylim[2],"\n"))
  
  ## for random selection of subset of resamples for plotting 
  if (!is.null(seed)) { 
    seedr_ = seed
  } else { 
    x_ = object$seed$seedr  
    seedr_ = x_[length(x_)] 
  }
  
  set.seed( seedr_ ) 
  
##  if ( (!is.null(fold)) & (is.null(usefold)) ) { 
##    usefold = fold
##    cat(paste("  fold is deprecated and will be dropped from future versions. Use usefold instead\n", 
##              " as its meaning is clearer, though it also applies to bootstrap resamples." ))
##  } else if ( is.null(usefold) ) { usefold = 1 } 
  
  if ((!is.null(cv)) & (is.null(resample))) {
    resample = cv 
    cat(paste("  cv is deprecated and will be dropped from future versions. Use resample instead\n", 
                " as its meaning is clearer, as it also applies to bootstrap resamples." ))
  }
  if (is.null(resample)) { 
    resample = 1 
  } 
  
  #----------------------------------------
  
  if ( plothr == 1 ) { plothr = exp(1) }
  
  family = object$sample[1] 
  tuning = object$tuning 
  bootstrap = as.numeric( tuning[8] )
  unique    = as.numeric( tuning[9] )
  nobs = object$sample[2] 
  
  if ((bootstrap == 0) & (is.null(object$xbetas.cv))) {
    resample = 0 
    cat(paste("  xbeta.cv is NULL so resample set to 0\n"))
  } 
  if ((bootstrap >= 1) & (is.null(object$xbetas.boot.oob))) {
    resample = 0 
    cat(paste("  xbeta.boot.oob is NULL so resample set to 0\n"))
  } 
  
  if (!(knottype %in% c(-1,0,1,2,3))) { knottype = 1 } 
  if ( ( bootstrap == 0 ) & ( bootci == 1) ) { bootci = 0 ; cat("  bootci does not apply for cross validation evaluation so bootci is ignored") ; }
  if ((bootstrap >= 1) & (bootci == 1)) {
    if (resample != 1) { cat("  For (bootci == 1) & (resample == 0), resample is set to 1\n") ; resample = 1 ; }
    if (oob      != 0) { cat("  For (bootci == 1) & (oob == 1), oob is set to 0\n") ; oob = 0 ; }
    if (usefold  != 1) { usefold = 1 ; cat("  For (bootci == 1) & (usefold == 0), usefold is ignored\n") }
    if (knottype != 1) { knottype= 1 ; # cat("  For (bootci == 1) & (knottype == 0), usefold is ignored\n")
    }
  }
  
#  print(c(resample, bootstrap, usefold))
  
  if ( (resample == 0) & (usefold == 1) ) {
    usefold = 0 
#    cat(paste( "For (resample==0), usefold is set to 0\n"))
  }

  stop_ = 0 
  
  if (is.null(wbeta)) {
    cat(paste(" specify num for wbeta = \n"))
    if (resample == 1) {  
      if (bootstrap == 0) {
        colnames(object$xbetas.cv)[1] = "null"
        Var = diag(cov((object$xbetas.cv))) 
      } else if (bootstrap >= 1) {
#        colnames(object$xbetas.boot.oob)[c(1,2)] = c("rep","oobid") 
#        object$xbetas.boot.oob[,-1] 
#        cov(object$xbetas.boot.oob[,-1])
        Var = diag(cov(object$xbetas.boot.oob[,-c(1,2)])) 
      }
    } else { 
      colnames(object$xbetas)[1] = "null"
      Var = diag(cov((object$xbetas)))
    } 
    if (min(Var) > 0.01) { Var = round(Var,digits=2) 
    } else if (min(Var) > 0.001) { Var = round(Var,digits=3) 
    } else if (min(Var) > 0.0001) { Var = round(Var,digits=4) 
    } else if (min(Var) > 0.00001) { Var = round(Var,digits=5) 
    }
    print( data.frame(cbind( Var, num=1:length(Var))) )
    stop_ = 1 
  } 
  
  if (!(resample %in% c(-1,0,1))) { stop_ = 1 }                                 ## -1 ?? 
  
  if (resample == 0) { oob = 0 }
  
  if ( (stop_ == 0) & (usefold == 0) ) {
    estci = calplot0(object, wbeta=wbeta, df=df, resample=resample, oob=oob, repn=NULL, knottype=1, 
                     trim=trim, plot=plot, lwd=2, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                     col.term=col.term, col.se=col.se, rug=rug, plothr=plothr , ... )$estci 
    stop_ = 2 
#    print( estci )
  }
  
  if ( (stop_ == 0) & (bootci == 1) ) {
    estci_orig = calplot0(object, wbeta=wbeta, df=df, resample=0, bootci=1, oob=oob, repn=NULL, knottype=knottype, 
                     trim=trim, plot=0, lwd=2, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                     col.term=col.term, col.se=col.se, rug=rug, plothr=plothr , ... )$estci  
#    print( estci_orig[1:5,] )
  }
  
  if ( (stop_ == 0) & (usefold == 1) ) { 
    if (resample == 1) { 
      if (bootstrap == 0) {
        xbeta = object$xbetas.cv[,wbeta ]  
        xbname = colnames(object$xbetas.cv)[wbeta] 
      } else if ((bootstrap >= 1) & (oob == 1)) {
        xbeta = object$xbetas.boot.oob[,wbeta+2 ]  
        xbname = colnames(object$xbetas.boot.oob)[wbeta+2] 
      } else if (bootstrap >= 1) {
        xbeta = object$xbetas.boot.inb[,wbeta+2 ]  
        xbname = colnames(object$xbetas.boot.inb)[wbeta+2] 
      }
    } else {
      xbeta = object$xbetas[,wbeta]  
      xbname = colnames(object$xbetas)[wbeta] 
    } 
#    print(xbeta[1:10])
    y_ = object$y_ ; 
    
    if (var(xbeta) == 0) {
      warning("  There may be no varaibalbity in the model X*Beta's\n",
              "    The calibration plot may be nonsensical")
    }

    minxb = min(object$xbetas[,wbeta ]) 
    maxxb = max(object$xbetas[,wbeta ]) 
    if (bootci == 1) {
#      minxb = min( object$xbetas[,wbeta+2 ], object$xbetas.boot.inb[,wbeta+2 ] )
#      maxxb = max( object$xbetas[,wbeta+2 ], object$xbetas.boot.inb[,wbeta+2 ] )
    }
    if (is.null(xlim)) {  xlim = c( minxb, maxxb )  }
    
    if (is.null(ylab)) {
      if (family == "cox") {
        if (plothr > 1) { ylab = "Calibrated HR"
        } else { ylab = "Calibrated log(HR)" } 
      } else {
        ylab = "Calibrated X*Beta"
      }
    }
    
    if (is.null(xlab)) { xlab = myxlab(family, plothr, resample, oob, xbname, bootstrap, bootci) }
    
    folds_n = as.numeric( object$tuning[1] )
    if (!(plotfold %in% c(0,1,2))) { plotfold = 0 }                             ## save estimates
    plotfold_ = plotfold 
    if (plotfold == 2) { plotfold_ = 2 }                                        ## make individual plots
    if (plotfold == 1) { plotfold_ = 0 }                                        ## make single plot with folds overlaid 
    
    if (bootstrap == 0) { reps = folds_n
    } else if (bootstrap >= 1) { reps = bootstrap }
    
    if (usefold == 1) {
      est.resample = matrix(0, nrow=reps, ncol=101)
      se.resample = est.resample
      lower.resample = est.resample
      upper.resample = est.resample
      meancalxb  = rep(0,folds_n)
      meanxb_x = rep(0,folds_n)
      minplot = min(minxb)
      maxplot = max(maxxb)
      
      for (i_ in c(1:reps)) { 
        calplotout = calplot0(object, wbeta=wbeta, df=df, resample=resample, bootci=bootci, oob=oob, repn=i_, knottype=knottype, 
                              trim=trim, plot=plotfold_, lwd=2, title=paste("fold",i_), xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, 
                              col.term=col.term, col.se=col.se, rug=rug, plothr=plothr , ... )
        
        ##      object1=object ; repn=i_ ; plot=1 ; lwd=2 ; title=paste("fold",i_) ;
        ## get back the xbeta from the whole data set -- no predicteds used
        if (i_ == 1) { plotxb = calplotout$estci[,1] } 
        est.resample[i_,]  = calplotout$estci[,2]
        se.resample[i_,]   = calplotout$estci[,3]
        lower.resample[i_,] = calplotout$estci[,4]
        upper.resample[i_,] = calplotout$estci[,5]
        
        meancalxb[i_] = calplotout$meancalxb 
        meanxb_x[i_] = calplotout$meanxb_ 
        minplot = max( minplot, calplotout$minmax[1] ) 
        maxplot = min( maxplot, calplotout$minmax[2] ) 
        if (i_ == 1) {
          xb_x = calplotout$xb_ ;
          y_x = calplotout$y_
        } else {
          xb_x = c(xb_x, calplotout$xb_) ;
          y_x  = c(y_x , calplotout$y_ )
        }
      }
      
      #      cat( paste( "  Range of X*Beta for calibration: c( ", floor(minxb), " , ", ceiling(maxxb), " ) \n"))
      cat(paste( "  Range of X*Beta for calibration: \n"))
      cat(c(minxb,maxxb),"\n")
      cat(paste( "  Range of calibrated confidence intervals: \n")) 
      if (bootci == 1) {
        #        cat( c( min(estci_orig[,4]), max(estci_orig[,5]) ) ,"\n") 
      } else {
        cat( c( min(est.resample), max(est.resample) ) ,"\n") 
      }  
      #      sample(length(xb_x), nobs )
      #      plotxb.resample = calplotout$estci[,1]
      c(minplot, maxplot)
      ## get average and SE from the folds_n folds 
      
      if (bootci == 1) { 
        #        print( estci_orig[1:20,1] )
        #        print( plotxb )
        est.resample = t( 2 * estci_orig[,2] - t(est.resample) )
      }    ## check row or column subtraction 
      
      est = colMeans(est.resample)
      
      if (bootstrap == 0) {
        se = sqrt(apply(est.resample,2,var)/folds_n)    
      } else if (bootstrap >= 1) {
        se = sqrt(apply(est.resample,2,var))    
      }
      upper = est + 2 * se 
      lower = est - 2 * se 
      
      if (plot %in% c(1,2)) {
        lty1 = 1 ; lty2 = 2 ; 
        
        xlims = quantile(xbeta, (c(trim,(100-trim))/100))
        keep_ = ( ((plotxb >= xlims[1]) & (plotxb <= xlims[2])) & ((plotxb >= xlim[1]) & (plotxb <= xlim[2])) ) 
        if ( trim > 0) { cat("  Spline curves are trimed at x's of", xlims[1], "and", xlims[2], ", or", trim, "% hi and lo\n") } 
        
        plotxb_ = plotxb[keep_]
        est_ = est[keep_]
        lower_ = lower[keep_]
        upper_ = upper[keep_]
        
        if (vref > 0) { vrefs = quantile(xbeta, (c(vref,(100-vref))/100)) }
        
        if (plotfold == 1) { lwd = 3 } else { lwd = 2 }
        
        #          plot(plotxb_, est_, type="l", lty=lty1, col=col.term, lwd=lwd, xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, axes=0, ... )
        #          lines( plotxb_, lower_, lty=lty1, col=col.se, lwd=lwd) 
        #          lines( plotxb_, upper_, lty=lty1, col=col.se, lwd=lwd) 
        
        indxlo = c(1:length(plotxb_)) [ plotxb_ < minplot ] 
        indxhi = c(1:length(plotxb_)) [ plotxb_ > maxplot ] 
        if ( length( c(indxlo,indxhi) ) > 0) {
          indxfold = c(1:length(plotxb_)) [-c(indxlo,indxhi) ] 
        } else {
          indxfold = c(1:length(plotxb_))
        }
        
        xlim0 = c(min(plotxb_), max(plotxb_)) 
        if (!is.null(xlim)) { xlim0 = xlim }
        if (is.null(ylim)) {
          ylim = c(min(est_), max(est_))
          if ( length(indxfold) > 0) { 
            if ( min(lower_[indxfold] ) < ylim[1] ) { ylim[1] = min(lower_[indxfold] ) }
            if ( max(upper_[indxfold] ) > ylim[2] ) { ylim[2] = max(upper_[indxfold] ) }
          }
        }
        
        plot(  plotxb_[indxlo], est_[indxlo] , type="l", 
               lty=lty2, col=col.term, lwd=lwd, xlab=xlab, ylab=ylab, xlim=xlim0, ylim=ylim, axes=0 , ... )
        lines( plotxb_[indxlo]  , lower_[indxlo  ], lty=lty2, lwd=lwd, col=col.se ) 
        lines( plotxb_[indxlo]  , upper_[indxlo  ], lty=lty2, lwd=lwd, col=col.se ) 
        lines( plotxb_[indxfold], est_  [indxfold], lty=lty1, lwd=lwd, col=col.term ) 
        lines( plotxb_[indxfold], lower_[indxfold], lty=lty1, lwd=lwd, col=col.se ) 
        lines( plotxb_[indxfold], upper_[indxfold], lty=lty1, lwd=lwd, col=col.se ) 
        lines( plotxb_[indxhi]  , est_  [indxhi  ], lty=lty2, lwd=lwd, col=col.term ) 
        lines( plotxb_[indxhi]  , lower_[indxhi  ], lty=lty2, lwd=lwd, col=col.se ) 
        lines( plotxb_[indxhi]  , upper_[indxhi  ], lty=lty2, lwd=lwd, col=col.se ) 
        
        if (family == "cox") { 
          abline(a=-mean(xbeta), b=1, lty=3)
          abline(h=0, lty=1)
          myaxis(ylim, plothr, est)
          myaxis(xlim, plothr, NULL, side=1)
        } else { 
          abline(a=0, b=1, lty=3)
          axis(side=1)
          axis(side=2)
        }
        box()
        
        if( bootstrap >= 1) {
          xbeta_rug = object$xbetas[,wbeta]  
        } else {
          xbeta_rug = xbeta  
        }
        
        keep_rug  = ((xbeta_rug >= xlim[1]) & (xbeta_rug <= xlim[2]))
        if (sum(keep_rug==0) > 0) {
          cat(paste("  Due to user specified xlim ", sum(keep_rug==0) , "tick marks are not displayed in rug\n"))
          cat(paste("  min:max xb = ", min(xbeta), max(xbeta) , "\n"))
        }
        
        xbeta_rug = xbeta_rug[ keep_rug ]
        y_rug     = y_[ keep_rug ]
        
        #        c( min(xbeta), min(xbeta_) ) 
        if (rug == 1) { myrug(xbeta_rug, y_rug, family) } 
        if (vref > 0) { abline(v=vrefs, lty=3, col=2) }
      }
      
      if (plotfold == 1) {
        for (i_ in c(1:reps)) {  
          lines( plotxb, est.resample[i_,], lty=i_, col=i_+1,  lwd=1) 
        }
      }
      
    }
  }

  if (stop_ %in% c(0,2)) {
    if (plot %in% c(0,2)) {
      if (usefold == 1) {  
        return_list = list(estimates=cbind(plotxb, est, se, lower, upper), 
                           est.resample=est.resample, se.resample=se.resample, lower.resample=lower.resample, upper.resample=upper.resample)  
      } else {
        return_list = list( estci = estci )  
      }
      return( return_list )
    }
  }
} 

#===============================================================================
# x_ = rnorm(100)^2
# y_ = rbinnom(100,1,(1/(1+exp(x_)))) 

#' A customized rug 
#'
#' @param x_ the x variable to be included in the rug
#' @param y_ an event variable for placing the x_ values top for events and bottom
#' for non-events 
#' @param family one of "cox", "binomial" or "gaussian" 
#'
#' @return A rug to plot 
#'
#' @noRd

myrug = function(x_, y_=NULL, family="gaussian") {
  if (family == "cox") {
    if (dim(y_)[2]== 3) { e_ = y_[,3] 
    } else { e_ = y_[,2] }
#    ticksize = (0.1*y_[,1]/max(y_[,1])) 
#    ticksize[ticksize < 0.001] = 0.001
#    print(c(length(ticksize), length(y_)))
#    print(summary(ticksize))
#    ticksize0 = ticksize[e_==0]
#    ticksize1 = ticksize[e_==1]
    ticksize0 = 0.03 
    ticksize1 = 0.03 
    rug(x_[e_==0], side=1, col=1, ticksize=ticksize0)
    rug(x_[e_==1], side=3, col=1, ticksize=ticksize1)
  } else if (family == "binomial") {
    e_ = y_ 
    rug(x_[e_==0], side=1, col=1,ticksize=0.03)
    rug(x_[e_==1], side=3, col=1,ticksize=0.03)
  } else {
    rug(x_, side=1, col=1, ticksize=0.03)
  }
}

#===============================================================================
#' Un-log the log(HR)'s for plotting
#'
#' @param zlim the y limits on the log(HR) scale A nested.glmnetr() output object for calibration
#' @param plothr  a power of determining the placement of y axis 
#' values, e.g. 2, sqrt(10) or 10.  The default of 0 plots the X*Beta.
#' @param est values of X*Beta estimates to be plotted which will determine the 
#' range to plot in case zlim is NULL. 
#' @param side 1 for horizontal axis, 2 for vertical 
#' 
#' @return An axis statement 
#'  
#' @noRd

myaxis = function(zlim,plothr,est=NULL,side=2, digits=2) {
  if (is.null(plothr)) { plothr = 0 }
  if ((is.null(zlim)) & (!is.null(est))) {
    zlim = c(min(est),max(est)) 
  } else if ((is.null(zlim)) & (is.null(est))) {
    plothr = 0
  }
  if (plothr > 1) {
    xbs = zlim
    hrs = exp(zlim) 
    logXhrs = log(hrs) / log(plothr)
    logXhrs_ = c(ceiling((logXhrs[1])): floor(logXhrs[2]))
    hrs_ = plothr^(logXhrs_)
    loghrs_ = log(hrs_)
    labls = round(hrs_,digits=digits) 
    myaxis = axis(side, at=loghrs_, labels = labls, las=1 )  
  } else {
    myaxis = axis(side)  
  }
  return(myaxis)
}

#===============================================================================
#' Set up default lable for xaxis 
#'
#' @param family family
#' @param plothr  a power of determining the placement of y axis 
#' values, e.g. 2, sqrt(10) or 10.  The default of 0 plots the X*Beta.
#' @param resample resample 
#' @param xbname xbname
#' @param bootstrap bootstrap
#' @param bootci bootci
#' 
#' @return An axis label
#'  
#' @noRd

myxlab = function(family, plothr, resample, oob, xbname, bootstrap,bootci=0) {
  if (family %in% c("cox") ) {
    if (plothr > 1) {
      if (resample==1) { 
        if (bootstrap == 0) { 
          xlab = paste(xbname,"HR (Nested CV)") 
        } else if (bootstrap >= 1) {
          if (bootci == 1) { xlab = paste(xbname,"HR (Bootstrap In Bag, bias adjusted)")
          } else if (oob == 1) { xlab = paste(xbname,"HR (Bootstrap Out Of Bag)") 
          } else {  xlab = paste(xbname,"HR (Bootstrap In Bag)") }
        }
      } else { 
        xlab = paste(xbname,"HR (naive)")
      }
    } else {
      if (resample==1) { 
        if (bootstrap == 0) { 
          xlab = paste(xbname,"log(HR) (Nested CV)") 
        } else if (bootstrap >= 1) {
          if (bootci == 1) { xlab = paste(xbname,"HR (Bootstrap In Bag, bias adjusted)")
          } else if (oob == 1) { xlab = paste(xbname,"log(HR) (Bootstrap Out Of Bag)") 
          } else {  xlab = paste(xbname,"log(HR) (Bootstrap In Bag)") }
        }
      } else { 
        xlab = paste(xbname,"log(HR) (naive)")
      }
    }
  } else if (family %in% c("binomial","gaussian") ) {
    if (resample==1) { 
      if (bootstrap == 0) { 
        xlab = paste(xbname,"X*Beta (Nested CV)") 
      } else if (bootstrap >= 1) { 
        if (bootci == 1) { xlab = paste(xbname,"X*Beta (Bootstrap In Bag, bias adjusted)")
        } else if (oob == 1) { xlab = paste(xbname,"X*Beta (Bootstrap Out Of Bag)") 
        } else { xlab = paste(xbname,"X*Beta (Bootstrap In Bag)") }
      }
    } else { xlab = paste(xbname,"X*Beta (naive)") }
  }
  return(xlab)
}

################################################################################
################################################################################