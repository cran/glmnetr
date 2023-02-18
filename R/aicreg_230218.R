#' Identify model based upon AIC criteria from a stepreg() putput 
#'
#' @param xs predictor input - an n by p matrix, where n (rows) is sample size, and p (columns) 
#' the number of predictors.  Must be in matrix form for complete data, no NA's, no Inf's, etc.,
#' and not a data frame. 
#' @param start start time, Cox model only - class numeric of length same as number of patients (n)
#' @param y_ output vector: time, or stop time for Cox model, y_ 0 or 1 for binomial (logistic), numeric for gaussian. 
#' Must be a vector of length same as number of sample size. 
#' @param event event indicator, 1 for event, 0 for census, Cox model only.
#' Must be a numeric vector of length same as sample size.
#' @param steps_n maximum number of steps done in stepwise regression fitting 
#' @param family model family, "cox", "binomial" or "gaussian" 
#' @param object A stepreg() output.  If NULL it will be derived.  
#' @param track Indicate whether or not to update progress in the console.  Default of
#' 0 suppresses these updates.  The option of 1 provides these updates.  In fitting 
#' clinical data with non full rank design matrix we have found some R-packages to
#' take a very long time or possibly get caught in infinite loops.  Therefore we allow
#' the user to track the package and judge whether things are moving forward or 
#' if the process should be stopped.  
#'
#' @return The identified model in form of a glm() or coxph() output object, with an 
#' entry of the stepreg() output object.
#' 
#' @export
#' 
#' @examples
#' set.seed(18306296)
#' sim.data=glmnetr.simdata(nrows=100, ncols=100, beta=c(0,1,1))
#' # this gives a more intersting case but takes longer to run
#' xs=sim.data$xs           
#' # this will work numerically
#' xs=sim.data$xs[,c(2,3,50:55)] 
#' y_=sim.data$yt  
#' event=sim.data$event
#' cox.aic.fit = aicreg(xs, NULL, y_, event, family="cox", steps_n=40) 
#' summary(cox.aic.fit)
#' 
#' y_=sim.data$yt  
#' norm.aic.fit = aicreg(xs, NULL, y_, NULL, family="gaussian", steps_n=40) 
#' summary(norm.aic.fit)
#' 
# track=track added to stepreg call 221220
aicreg = function(xs, start, y_, event, steps_n=steps_n, family=family, object=NULL, track=0) {
  if (!is.null(object)) {
    #  stepreg.fit.all.best = as.matrix( cv.stepreg.fit.all$stepreg.fit.all.best ) 
    if (inherits(object,"cv.stepreg")) { 
      stepreg.fit.all.best = object$stepreg.fit.all.best 
    } else if (inherits(object,"nested.glmnetr")) { 
      stepreg.fit.all.best = object$stepreg.fit.all.best 
    }
  } else {
    if (track >= 1) { cat(paste0("\n", " ########## Derive stepwise model on all data ################################################" , "\n")) } 
    stepreg.fit.all  = stepreg(xs, start, y_, event, steps_n=steps_n, family=family, track=track)  
    class(stepreg.fit.all) = "data.frame" 
    stepreg.fit.all.best = stepreg.fit.all[(stepreg.fit.all$best==1) , ] 
  }
  aic = 2*stepreg.fit.all.best[,1] - 2*stepreg.fit.all.best[,5]                ;# plot(c(1:steps_n),aic)
  if (family == "binomial") { aic = aic + 2 
  } else if (family == "gaussian") { aic = aic + 4 }
  best.aic.all.index = which.min( aic )
  aic.fit.all.d = stepreg.fit.all.best[best.aic.all.index,] 
  temp_ = colnames(aic.fit.all.d)
  aic.fit.all.n = as.numeric(aic.fit.all.d)
  names(aic.fit.all.n) = temp_ 
  if ((length(temp_)-7)%%2 == 1) { 
    intercept = TRUE 
    nvars = (length(temp_)-8)/2
  } else { 
    intercept = FALSE 
    nvars = (length(temp_)-7)/2 
  } 
#  temp_[8:(7+nvars)]
  if (!intercept) { aic.fit.all.vars = temp_[(nvars+8):(length(temp_))][aic.fit.all.n[8:(nvars+7)]==1] 
  } else { aic.fit.all.vars = ( temp_[(nvars+9):(length(temp_))][aic.fit.all.n[8:(nvars+7)]==1] ) }
  #    if (family=="cox") {
  #      beta = stepreg.fit.all.best [ best.aic.all.index , (nvars+8):(2*nvars+7) ]        ## get beta, the regression estimate
  #      beta = matrix(beta, nrow=nvars, ncol=1) 
  #    } else {
  #      beta = stepreg.fit.all.best [ best.aic.all.index , (nvars+8):(2*nvars+8) ]        ## get beta, the regression estimate
  #      beta = matrix(beta, nrow=(nvars+1), ncol=1) 
  #    }
  if (family=="cox") {
    if (is.null(start)) { 
      dataf = as.data.frame( cbind( y_, event, xs  ) )
      form1 = formula( paste("Surv(" , colnames(dataf)[1] ,  ", " , colnames(dataf)[2] , ") ~ ", paste(aic.fit.all.vars, collapse = " + " ) ) )  
    } else { 
      dataf = as.data.frame( cbind( start, y_, event, xs  ) ) 
      form1 = formula( paste("Surv(" , colnames(dataf)[1] ,  ", " , colnames(dataf)[2],  ", " , colnames(dataf)[3] , ") ~ ", paste(aic.fit.all.vars, collapse = " + " ) ) )  
    }
    aic.fit.all = coxph(form1, dataf) 
  } else {
    dataf = as.data.frame( cbind( y_, xs  ) )
    form1 = formula( paste("y_ ~ " , paste(aic.fit.all.vars, collapse = " + " ) ) )  
    aic.fit.all = glm(form1, family=family, data=dataf)
  }
  aic.fit.all$stepaic
  aic.fit.all$stepaic = aic 
  aic.fit.all$aic.fit.n = aic.fit.all.n
  aic.fit.all$stepreg.fit.aic = stepreg.fit.all.best 
  summary(aic.fit.all)
  return(aic.fit.all)
}