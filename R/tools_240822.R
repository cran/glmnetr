################################################################################
##### tools_yymmdd.R ###########################################################
################################################################################
#' Get seeds to store, facilitating replicable results 
#'
#' @param seed The intput seed as a start, NULL, a vector of lenght 1 or 2, or 
#' a list with vectors of lenght 1 or the number of folds, $seedr for most models 
#' and $seedt for the ANN fits 
#' @param folds_n The number of folds in general
#' @param folds_ann_n The number of folds for the ANN fits 
#'
#' @return seed(s) in a list format for input to subsequent runs
#' 
#' @seealso
#'   \code{\link{nested.glmnetr}} 
#' 
#' @export
#'
glmnetr_seed = function(seed, folds_n=10, folds_ann_n=NULL) {
  if (is.null(folds_ann_n)) { folds_ann_n = folds_n } 
  
  if (is.null(seed)) { 
    seedr = round(runif(1)*1e9) 
    if (seedr == 0) { seedr = 1 } 
    set.seed( seedr ) 
    seedr = c( seedr, round(runif(folds_n)*1e9) )
    seedr = ifelse((seedr==0), 1, seedr)
    
    seedt = round(runif(1)*1e9) 
    if (seedt == 0) { seedt = 1 } 
    set.seed( seedt )
    seedt = c( seedt, round( runif(folds_ann_n)*1e9) )
    seedt = ifelse((seedt==0), 1, seedt) 
    
    seed = list(seedr=seedr, seedt=seedt)
  } else if (is.list(seed)) {           ############################ seed a list   
    
    if (is.null(seed$seedr)) { 
      seed$seedr = round(runif(1)*1e9)
    } else if(is.na(seed$seedr[1])) { 
      seed$seedr[1] = round(runif(1)*1e9)  
    } else if (seed$seedr[1]==0) { 
      seed$seedr[1] = round(runif(1)*1e9)  
    }
    seed$seedr = ifelse((seed$seedr==0), 1, seed$seedr)
    set.seed(seed$seedr[1])
    seedr_t = round(runif(folds_n)*1e9) 
    
    if (is.null(seed$seedt)) { 
      seed$seedt = round(runif(1)*1e9)
    } else if(is.na(seed$seedt[1])) { 
      seed$seedt[1] = round(runif(1)*1e9)  
    } else if (seed$seedt[1]==0) { 
      seed$seedt[1] = round(runif(1)*1e9)  
    } 
    seed$seedt = ifelse((seed$seedt==0), 1, seed$seedt)
    set.seed(seed$seedt[1])
    seedt_t = round(runif(folds_ann_n)*1e9) 
    
    for (i_ in c(1:folds_n)) {
      seedr_ = seed$seedr[i_+1] 
      if (is.null(seedr_)) { 
        seed$seedr[i_+1] = seedr_t[i_] 
      } else if (is.na(seedr_)) { 
        seed$seedr[i_+1] = seedr_t[i_] 
      } else if (seedr_==0) { 
        seed$seedr[i_+1] = seedr_t[i_] 
      } 
    }
    seed$seedr = seed$seedr %% 2.147e9 
    seed$seedr = ifelse((seed$seedr==0), 1, seed$seedr)

    for (i_ in c(1:folds_ann_n)) {
      seedt_ = seed$seedt[i_+1] 
      if (is.null(seedt_)) { 
        seed$seedt[i_+1] = seedt_t[i_] 
      } else if (is.na(seedt_)) { 
        seed$seedt[i_+1] = seedt_t[i_]  
      } else if (seedt_==0) { 
        seed$seedt[i_+1] = seedt_t[i_] 
      } 
    }
    seed$seedt = seed$seedt %% 2.147e9 
    seed$seedt = ifelse((seed$seedt==0), 1, seed$seedt)
    
  } else {                            ########################### seed a numeric   
    seedr = seed[1]
    seedt = seed[2]
    
    if (is.null(seedr))      { seedr = round(runif(1)*1e9) 
    } else if (is.na(seedr)) { seedr = round(runif(1)*1e9) 
    } else if (seedr==0)     { seedr = round(runif(1)*1e9) }
    seedr = ifelse((seedr==0), 1, seedr)
    set.seed(seedr)  
    seedr = c(seedr, round(runif(folds_n)*1e9) )  
    seedr = ifelse((seedr==0), 1, seedr)
    
    if (is.null(seedt))      { seedt = round(runif(1)*1e9) 
    } else if (is.na(seedt)) { seedt = round(runif(1)*1e9) 
    } else if (seedt==0)     { seedt = round(runif(1)*1e9) }  
    seedt = ifelse((seedt==0), 1, seedt)
    set.seed(seedt)  
    seedt = c(seedt, round(runif(folds_ann_n)*1e9) ) 
    seedt = ifelse((seedt==0), 1, seedt)
    
    seed = list(seedr=seedr, seedt=seedt)
  }
  
  return ( seed ) 
  
}

#temp1 = seedx(NULL)
#temp1
#seed = c( temp1$seedr[1], temp1$seedt[1] ) 
#seedl = list( seedr=temp1$seedr[1], seedt=temp1$seedt[1] ) 
#seed
#seedl
#seedx( seed ) 
#seedx( seed ) 

################################################################################
################################################################################

#' Generate foldid's by factor levels
#'
#' @param event the outcome variable in a vector identifying the different potential 
#' levels of the outcome
#' @param fold_n the numbe of folds to be constructed
#'
#' @return foldid's in a vector the same length as event
#' 
#' @seealso
#'   \code{\link{get.foldid}} , \code{\link{nested.glmnetr}} 
#'  
#' @export
#'
factor.foldid = function(event, fold_n=10) { 
  nobs    = length(event)
  nlevels = length(table(event))
  if (nlevels == 1) {
    if (is.factor(levels)) { levels = as.numeric(levels) } 
    #    foldid = sample( rep( 1:fold_n , ceiling(nobs/fold_n) )[ 1:nobs ] , nobs ) 
  } else if (nlevels == 2) {
    if (length(names(table(event))) == 2) { 
      event = as.factor(event) 
    }
  }
  if (is.factor(event)) {
    event_idm = matrix(rep(0,nobs*nlevels),nrow=nlevels,ncol=nobs)
    for (i_ in c(1:nlevels) ) { 
      event_idm[i_,] = 1*(event==levels(event)[i_]) 
    }
    event_idm
    nobslvl = rowSums(event_idm)
    nobslvl 
    for (i_ in c(1:nlevels) ) { 
      if (nobslvl[i_] < fold_n) {
        cat(paste("nlevels = ", nlevels, "i_ = ", i_, " nobslvl[i_] = " , nobslvl[i_] ))
        nobslvl[i_] = 1 
      }
      rep(1:fold_n, ceiling(nobslvl[i_]/fold_n)) 
      foldidi = sample( rep( 1:fold_n , ceiling(nobslvl[i_]/fold_n) )[ 1:nobslvl[i_] ] , nobslvl[i_] ) 
      length(foldidi)
      sum((event_idm[i_,]==1))
      event_idm[i_,][(event_idm[i_,]>=1)] = foldidi 
      event_idm
    }
    foldid = colSums(event_idm)
  } else {
    foldid = sample( rep( 1:fold_n , ceiling(nobs/fold_n) )[ 1:nobs ] , nobs ) 
  }
  return(foldid)
}

#' Generate foldid's by 0/1 factor for bootstrap like samples where unique option between 0 and 1
#'
#' @param event the outcome variable in a vector identifying the different potential 
#' levels of the outcome
#' @param fraction the fraction of the whole sample included in the bootstratp sample
#'
#' @return foldid's in a vector the same length as event
#' 
#' @seealso
#'   \code{\link{get.foldid}} , \code{\link{nested.glmnetr}} 
#'  
#' @export
#'

boot.factor.foldid = function(event, fraction) { 
    nobs = length(event)
    nobs0 = sum(event==0)
    nobs1 = sum(event==1)
    noob  = (fraction*nobs )%/%1 + rbinom(1,1,(fraction*nobs )%%1)  
    noob0 = (fraction*nobs0)%/%1 + rbinom(1,1,(fraction*nobs0)%%1)
    noob1 = noob - noob0 
    
    ids0 = c(1:nobs)[event==0]
    length(ids0)
    xids0inb = sort( sample( c(1:nobs0), noob0, replace=0 ) ) 
    xids0oob = c(1:nobs0)[ -xids0inb ]
    ids0inb = ids0[ xids0inb]
    ids0oob = ids0[-xids0inb]
    table(table(ids0inb))
    table(table(ids0oob))
    table(table(c( ids0inb , ids0oob)))
    
    ids1 = c(1:nobs)[event==1]
    length(ids1)
    xids1inb = sort( sample( c(1:nobs1), noob1, replace=0 ) ) 
    xids1oob = c(1:nobs1)[ -xids1inb ]
    ids1inb = ids1[ xids1inb]
    ids1oob = ids1[-xids1inb]
    table(table(ids1inb))
    table(table(ids1oob))
    table(table(c( ids1inb , ids1oob)))
    
    table(table( c(ids0inb, ids0oob, ids1inb, ids1oob) ))
    
    returnlist = list(ids0inb=ids0inb, ids0oob=ids0oob, ids1inb=ids1inb, ids1oob=ids1oob) 
    
    idsinb = sort( c(returnlist$ids0inb , returnlist$ids1inb) )
  #  table(table(idsinb))
  return(returnlist)
}

################################################################################
################################################################################

#' Get foldid's with branching for cox, binomial and gaussian models 
#'
#' @param y_ see help for cv.glmnetr() or nested.glmnetr() 
#' @param event see help for cv.glmnetr() or nested.glmnetr() 
#' @param family see help for cv.glmnetr() or nested.glmnetr() 
#' @param folds_n see help for cv.glmnetr() or nested.glmnetr() 
#' @param stratified see help for cv.glmnetr() or nested.glmnetr() 
#' 
#' @return A numeric vector with foldid's for use in a cross validation 
#' 
#' @seealso
#'   \code{\link{factor.foldid}} , \code{\link{nested.glmnetr}} 
#'  
#' @export
#' 
get.foldid = function( y_, event, family, folds_n , stratified = 1 ) {
  nobs = length( y_ )
  if (stratified >= 1) {
    if   (family == "binomial") { 
      foldid = factor.foldid(y_   , folds_n)    
    } else if (family == "cox") { 
      foldid = factor.foldid(event, folds_n)
    } else {
      foldid = sample( rep( 1:folds_n , ceiling(nobs/folds_n) )[ 1:nobs ] , nobs )   
    }
  } else {
    foldid = sample( rep( 1:folds_n , ceiling(nobs/folds_n) )[ 1:nobs ] , nobs )  
  }
  return( foldid )
}
################################################################################
################################################################################
#' Get foldid's when id variable is used to identify groups of dependent 
#' sampling units. With branching for cox, binomial and gaussian models 
#'
#' @param y_ see help for cv.glmnetr() or nested.glmnetr() 
#' @param event see help for cv.glmnetr() or nested.glmnetr() 
#' @param id see help for nested.glmnetr() 
#' @param family see help for cv.glmnetr() or nested.glmnetr() 
#' @param folds_n see help for cv.glmnetr() or nested.glmnetr() 
#' @param stratified see help for cv.glmnetr() or nested.glmnetr() 
#' 
#' @return A numeric vector with foldid's for use in a cross validation 
#' 
#' @importFrom stats aggregate 
#' 
#' @seealso
#'   \code{\link{factor.foldid}} , \code{\link{nested.glmnetr}} 
#'  
#' @export
#' 
get.id.foldid = function(y_, event, id, family, folds_n, stratified) {
  rder = c(1:length(y_))
  dftmp = as.data.frame( cbind(rder, id, y_) ) 
  if (family %in% c("cox")) { 
    dftmp$event = event 
    sumsby = aggregate(dftmp$event, list(dftmp$id), FUN=sum)
    colnames(sumsby) = c("id","event")
    sumsby$y_ = sumsby$event 
    sumsby$event = 1*(sumsby$event > 0)
  } else if (family %in% c("binomial")) { 
    sumsby = aggregate(dftmp$y_, list(dftmp$id), FUN=sum)
    colnames(sumsby) = c("id","y_")
    sumsby$y_ = 1*(sumsby$y_>0)
    if (length(unique(sumsby$event))==1) { stratified = 0 ; print(" stratified = 0 ")}
    sumsby$event = NULL
  } else {
    sumsby = aggregate(dftmp$y_, list(dftmp$id), FUN=mean)
    colnames(sumsby) = c("id","y_")
    sumsby$event = NULL
  }
#  print( sumsby )
#  print( sumsby$id )
#  print( sumsby$event )
  foldidtmp = get.foldid(sumsby$y_, sumsby$event, family, folds_n, stratified=stratified)    ##<<------------------
  foldidkey = as.data.frame( cbind(sumsby$id, foldidtmp) )
  names(foldidkey) = c("id","foldid")
  dftmp[1:10,]
  foldidkey[1:10,]
  dfmerged = merge(dftmp, foldidkey)
  dfmerged[1:10,]
  if (family != "cox") {
    dfmerged = dfmerged[order(dfmerged$rder) ,] [,c(2,1,3,4)]
  } else {
    dfmerged = dfmerged[order(dfmerged$rder) ,] [,c(2,1,4,5)]
  }
  foldid = dfmerged[,4] 
  return(list(foldid=foldid,foldidkey=foldidkey))
}

################################################################################
################################################################################

#' Calculate the CoxPH saturated log-likelihood
#'
#' @description Calculate the saturated log-likelihood for the Cox 
#' model using both the Efron and Breslow approximations for the case where all ties
#' at a common event time have the same weights (exp(X*B)).  For 
#' the simple case without ties 
#' the saturated log-likelihood is 0 as the contribution to the log-likelihood at 
#' each event time point can be made arbitrarily close to 1 by assigning a much larger 
#' weight to the record with an 
#' event.  Similarly, in the case of ties 
#' one can assign a much larger weight to be associated with one of the event times 
#' such that the associated record contributes a 1 to the likelihood.  Next one 
#' can assign a very large weight to a second tie, but smaller than the first tie 
#' considered, and this too will contribute a 1 to the 
#' likelihood.  Continuing in this way for this and all time points with ties, 
#' the partial log-likelihood is 0, just like for the no-ties case.  Note, this is 
#' the same argument with which we derive the log-likelihood of 0 for the no ties
#' case. Still, to be consistent with others we derive the saturated log-likelihood
#' with ties under the constraint that all ties at each event time carry the same weights.
#'
#' @param y_ Time variable for a survival analysis, whether or not there is a start time
#' @param e_ Event indicator with 1 for event 0 otherwise. 
#'
#' @return Saturated log likelihood for the Efron and Breslow approximations.
#' 
#' @seealso
#'   \code{\link{nested.glmnetr}} 
#' 
#' @export
#'
cox.sat.dev =  function(y_, e_) {
  times = y_[(e_==1)]
  ties = table(table(times))
  d_ = as.numeric(names(ties))
  loglik_sat_ef = -sum(ties*log(factorial(d_)))
  loglik_sat_br = -sum(ties*(d_*log(d_)))
  sat.dev = c(loglik_sat_ef , loglik_sat_br )
  names(sat.dev) = c("loglik_sat_efron", "loglik_sat_breslow")
  return( sat.dev )
}

###############################################################################################################
###############################################################################################################

#' Get elapsed time in c(hour, minute, secs) 
#'
#' @param time1 start time
#' @param time2 stop time
#'
#' @return Returns a vector of elapsed time in (hour, minute, secs) 
#' 
#' @seealso
#'   \code{\link{diff_time}} 
#'
#' @export
#' 
diff_time1 = function(time1, time2) {
  hour = difftime(time1, time2, units="hour") ;   
  class(hour) = NULL ; 
  hour   = floor( hour) 
  minute = difftime(time1, time2, units="min") ;   
  class(minute) = NULL ; 
  minute = floor( minute %% 60) 
  secs   = difftime(time1, time2, units="secs") ;   
  class(secs) = NULL ; 
  secs = floor( secs %% 60) 
  hr_mn_sc = c(hour, minute, secs)
#  print(paste( "hr_mn_sc = ", length( hr_mn_sc),"\n"))
  #names(hr_mn_sc) = c("Hour", "Min", "Sec") 
  return(hr_mn_sc)
}

###############################################################################################################
###############################################################################################################

#' Output to console the elapsed and split times  
#'
#' @param time_start beginning time for printing elapsed time 
#' @param time_last last time for calculating split time
#'
#' @return Time of program invocation 
#' 
#' @seealso
#'   \code{\link{diff_time}} , \code{\link{nested.glmnetr}} 
#' 
#' @export 
#'
#' @examples
#' time_start = diff_time()
#' time_last = diff_time(time_start)
#' time_last = diff_time(time_start,time_last)
#' time_last = diff_time(time_start,time_last)
#'
diff_time = function(time_start=NULL, time_last=NULL) {
  time_sys = Sys.time()
  time_elapse = diff_time1(time_sys, time_start)
  if (is.null(time_start)) { 
    cat(paste0("Start at Sys.time = ", time_sys,  " \n"))
  } else if (is.null(time_last)) { 
    cat(paste0("Sys.time = ", time_sys, 
               ",  elapsed time = " , time_elapse[1], ":", time_elapse[2], ":", time_elapse[3] , " h:m:s",  " \n"))
  } else { 
    time_split  = diff_time1(time_sys, time_last ) 
    cat(paste0("Sys.time = ", time_sys, 
               ",  elapsed time = " , time_elapse[1], ":", time_elapse[2], ":", time_elapse[3] , " h:m:s",
               ",  time split = "  , time_split[1] , ":", time_split[2] , ":", time_split[3] , " h:m:s", " \n"))
  } 
  time_last = time_sys ; 
  return(time_last)
}

###############################################################################################################
###############################################################################################################
#' Generate example data 
#' 
#' @description
#' Generate an example data set with specified number of observations, 
#' and predictors.  The first column in the design matrix is identically equal to 1
#' for an intercept.  Columns 2 to 5 are for the 4 levels of a character variable, 
#' 6 to 11 for the 6 levels of another character variable.  Columns 12 to 17 are for 3
#' binomial predictors, again over parameterized.  Such over parameterization
#' can cause difficulties with the glmnet() of the 'glmnet' package.  
#'
#' @param nrows Sample size (>=100) for simulated data, default=1000.
#' @param ncols Number of columns (>=17) in design matrix, i.e. predictors, default=100. 
#' @param beta Vector of length <= ncols for "left most" coefficients.  If beta has 
#' length < ncols, then the values at length(beta)+1 to ncols are set to 0.
#' Default=NULL, where a beta of length 25 is assigned standard normal values.
#' @param intr either NULL for no interactions or a vector of length 3 to impose a 
#' product effect as decribed by 
#' intr[1]*xs[,3]*xs[,8] + intr[2]*xs[,4]*xs[,16] + intr[3]*xs[,18]*xs[,19] + intr[4]*xs[,21]*xs[,22]
#' @param nid number of id levels where each level is associated with a random 
#' effect, of variance 1 for normal data.  
#'
#' @return A list with elements xs for desing matrix, y_ for a quantitative outcome, 
#' yt for a survival time, event for an indicator of event (1) or censoring (0), 
#' in the Cox proportional hazards survival model setting, yb for 
#' yes/no (binomial) outcome data, and beta the beta used in random number generation. 
#' 
#' @seealso
#'   \code{\link{nested.glmnetr}} 
#'
#' @importFrom stats rnorm runif rexp rbinom rbeta 
#' @importFrom Matrix rankMatrix
#' 
#' @export
#' 
#' @examples
#' sim.data=glmnetr.simdata(nrows=1000, ncols=100, beta=NULL)
#' # for Cox PH survial model data 
#' xs=sim.data$xs 
#' y_=sim.data$yt
#' event=sim.data$event
#' # for linear regression model data 
#' xs=sim.data$xs 
#' y_=sim.data$y_
#' # for logistic regression model data 
#' xs=sim.data$xs 
#' y_=sim.data$yb
#' 
glmnetr.simdata = function(nrows=1000, ncols=100, beta=NULL, intr=NULL, nid=NULL) {
  if (nrows<100) {nrows=100}
  if (ncols<17) {ncols=17}
  if (is.null(beta)) {
    beta = c(rep(0,ncols)) 
    beta[1:25]=rnorm(25)
  } else {
    beta_ = c(rep(0,ncols)) 
    beta_[1:length(beta)]=beta
    print(beta)
    print(beta_)
    beta=beta_
  }
  
  xs = matrix(rep(0,nrows*ncols), nrow=nrows, ncol=ncols)
  for (i_ in c(1:nrows)) { 
    runi = runif(5)
    #    print(i_)
    #    print(runi[1])
    xs[i_,1] = 1 
    if        (runi[1]< 0.3) {xs[i_,2] = 1
    } else if (runi[1]< 0.5) {xs[i_,3] = 1 
    } else if (runi[1]< 0.8) {xs[i_,4] = 1 
    } else                   {xs[i_,5] = 1 }
    if        (runi[2]< 0.1) {xs[i_,6] = 1
    } else if (runi[2]< 0.3) {xs[i_,7] = 1 
    } else if (runi[2]< 0.6) {xs[i_,8] = 1 
    } else if (runi[2]< 0.8) {xs[i_,9] = 1 
    } else if (runi[2]< 0.9) {xs[i_,10] = 1 
    } else                   {xs[i_,11] = 1 }
    if        (runi[3]< 0.1) {xs[i_,12] = 1
    } else                   {xs[i_,13] = 1 }
    if        (runi[2]< 0.2) {xs[i_,14] = 1
    } else                   {xs[i_,15] = 1 }
    if        (runi[2]< 0.3) {xs[i_,16] = 1
    } else                   {xs[i_,17] = 1 }
    xs[i_, c(18:ncols)] = rnorm(ncols-17)
  }
  
  #  xs[1:10,1:20]
  rankMatrix(xs)[1]
  #  summary(xs)
  score = xs %*% beta 
  summary(score)
  var(score)
  if (!is.null(intr)) { 
    #    score = score + intr[1]*xs[,3]*xs[,8]                           + intr[2]*xs[,18]*xs[,19] + intr[3]*xs[,21]*xs[,22]
    score = score + intr[1]*xs[,3]*xs[,8] + intr[2]*xs[,4]*xs[,16] + intr[3]*xs[,18]*xs[,19] + intr[4]*xs[,21]*xs[,22]
  }
  var(score)
  # normal 

  if (!is.null(nid)) {
    ids = rep(c(1:nid),ceiling(nrows/nid) ) [1:nrows]
    dfids = as.data.frame(cbind(1:nrows,ids))
    names(dfids) = c('order','id')
#    dfids

    dftmp2 = as.data.frame(cbind(c(1:nid),rnorm(nid)))
    names(dftmp2) = c('id','rnrm')
#    dftmp2
    
    dfran = merge(dfids,dftmp2)
#    dfran[1:10,]
    
    dfran = dfran[order(dfran$order) ,c(2,1,3)] 

    id = dfran$id
    score = score + dfran$rnrm 
  } else { id = NULL }
  
  y_ = as.vector( score ) + rnorm( nrows ) 
  
  # survival exponential 
  ct = 2*as.vector( -log(runif(nrows)) )
  rates = exp(score)
  yt = rexp(length(rates),rates) 
  event = as.vector( (yt < ct)*1 ) 
  yt    = as.vector( event*yt + (1-event)*ct ) 
  # start times
  runis = (runif(nrows) < 0.3)*1 
  start = as.vector( rbeta(nrows,2,4) * yt * runis )
  start[(start > yt-0.0001)] = yt[(start > yt-0.0001)]-0.0001
  sum((start > yt-0.0001))
  # logistic binomial 
  p_ = 1/(1+exp(-score)) 
  yb = rbinom(length(p_),1,p_)
  
  nameXS = c("X1")  
  for (i_ in c(2:dim(xs)[2])) {
    nameXS = c(nameXS, paste0("X", i_))
  } 
  colnames(xs) = nameXS 
  #  dim(xs)
  #  length(y_)
  #  length(event)
  #  length(yt)
  #  sum(event)
  return(list(xs=xs,start=start,y_=y_, yt=yt, id=id, event=event,yb=yb, beta=beta))
}

###############################################################################################################
###############################################################################################################
#' Calculate deviance ratios for CV based 
#'
#' @description
#' Calculate deviance ratios for individual folds and collectively.  Calculations
#' are based upon the average -2 Log Likelihoods calculated on each leave out test fold
#' data for the models trained on the other (K-1) folds.
#'
#' @param m2.ll.mod -2 Log Likelihoods calculated on the test data 
#' @param m2.ll.null -2 Log Likelihoods for the null models  
#' @param m2.ll.sat  -2 Log Likelihoods for teh saturated models
#' @param n__ sample zize for the indivual foles, or number of events for the Cox model 
#'
#' @return a list with devrat.cv for the deviance ratios for the indivual folds, 
#' and devrat, a single collective deviance ratio
#' 
#' @seealso
#'   \code{\link{nested.glmnetr}} 
#' 
#' @export
#'
devrat_ = function( m2.ll.mod, m2.ll.null, m2.ll.sat, n__ ) {
  devrat.cv = ((m2.ll.null - m2.ll.mod)/(m2.ll.null - m2.ll.sat )) 
  devrat = (sum(n__*(m2.ll.null - m2.ll.mod)) / sum(n__*(m2.ll.null - m2.ll.sat )) ) 
  return(list(devrat.cv=devrat.cv, devrat=devrat, avedevrat=mean(devrat.cv), 
              sedevrat = sqrt(var(devrat.cv)/length(devrat.cv)) ) )
}

###############################################################################################################
###############################################################################################################

get.DevRat = function( devian.cv, null.m2LogLik.cv, sat.m2LogLik.cv, n.cv ) {
  ncoldevian = dim(devian.cv)[2] 
  AllDevRat = rep(1,ncoldevian)
  for (j_ in c(1:ncoldevian)) {
    AllDevRat[j_] = devrat_(devian.cv[,j_], null.m2LogLik.cv, sat.m2LogLik.cv, n.cv )[[2]]
  }
  return( AllDevRat )
}

