###############################################################################################################
###############################################################################################################
## do cross validation to choose tuning paramter, i.e. lambda and gamma #######################################
#' Get a cross validated tuned relaxed lasso model fit.
#'
#' @description
#' This function derives a relaxed lasso model and derives hyperparmaters 
#' using cross validaiton.  It is analogous to the glmnet() function of the 
#' _glmnet_ package, but handles caes where glmnet() may run slowly when using the
#' relaxed-=TRUE option.   
#' 
#' @param xs Predictor matrix 
#' @param start Vector of start times or the Cox model.  May be NULL.
#' @param y_ outcome vector
#' @param event event vector in case of the Cox model.  May be NULL for other models.  
#' @param family Model family, one of "cox", "gaussian" or "binomial".  
#' @param lambda The lambda vector.  May be NULL.  
#' @param gamma The gamma vector.  Default is c(0,0.25,0.50,0.75,1). 
#' @param folds_n Number of folds for cross validation.  Default and recommended is 10.
#' @param limit limit the small values for lambda after the initial fit.  
#' This will calcualtions that have minimal impact on the cross validation.  
#' Default is 2 for moderate limitation, 1 for less limitation, 0 for none.     
#' @param fine Use a finer step in determining lambda.  Of little value unless one 
#' repeats the cross valiaiton many times to more finely tune the hyper parameters.  
#' See the _glmnet_ documentation.  
#' @param seed A seed for set.seed to assure one can get the same results twice.  If NULL 
#' the program will generate a random seed.  Whether specified or NULL, the seed is stored in the output
#' object for future reference.  
#' @param foldid A vector of integers to associate each record to a fold.  Should be integers between 1 and folds_n.
#' @param time Indicate whether or not to update progress in the console.  Default of
#' 0 suppresses these updates.  The option of 1 provides these updates.  In fitting 
#' clinical data with non full rank design matrix we have found some R-packages to
#' take a vary long time or seemingly be caught in infinite loops.  Therefore we allow
#' the user to track the package and judge whether things are moving forward or 
#' if the process should be stopped. 
#' 
#' @details
#' #' This is main program for model derivation.  As currently implemented the 
#' package requires the data to be input as #' vectors and matrices with no missing values 
#' (NA).  All data vectors and matrices must be numerical.  For categorical variables one
#' should first construct corresponding numerical variables to represent these 
#' categories.  To take advantage of the lasso model, one can use one hot coding
#' assigning an indicator for each level of each categorical varaible, or creating 
#' as well other contrasts variables suggested by the subject matter.
#'  
#' @return A cross validated relaxed lasso model fit.
#' 
#' @author Walter Kremers (kremers.walter@mayo.edu)
#' 
#' @seealso
#'   \code{\link{glmnetr}} , \code{\link{nested.glmnetr}}  , \code{\link{glmnetr.simdata}} 
#'   
#' @export
#' 
#' @importFrom stats runif 
#' @importFrom survival Surv 
#' @importFrom glmnet cv.glmnet 
#' @importFrom Matrix rankMatrix 
#'
#' @examples
#' # set seed for random numbers, optionally, to get reproducible results
#' set.seed(82545037)
#' sim.data=glmnetr.simdata(nrows=100, ncols=100, beta=NULL)
#' xs=sim.data$xs 
#' y_=sim.data$y_ 
#' event=sim.data$event
#' # for this example we use a small number for folds_n to shorten run time 
#' cv.glmnetr.fit = cv.glmnetr(xs, NULL, y_, NULL, family="gaussian", folds_n=3, limit=2) 
#' plot(cv.glmnetr.fit)
#' plot(cv.glmnetr.fit, coefs=1)
#' summary(cv.glmnetr.fit)
#' 
cv.glmnetr = function( xs, start, y_, event, family="cox", lambda=NULL, gamma=c(0,0.25,0.50,0.75,1), folds_n=10, limit=2, fine=0, time=0, seed=NULL, foldid=NULL) {

  if (is.null(gamma)) { gamma=c(0,0.25,0.50,0.75,1) }
  if (!(0 %in% gamma)) { gamma=c(0, gamma)}
  if (!(1 %in% gamma)) { gamma=c(gamma, 1)}
  gamma = gamma[(gamma >= 0)]
  gamma = gamma[(gamma <= 1)]
  gamma[order(gamma)]
  gamma_n = length(gamma)
  #  cat(paste0(" gamma_n = ", gamma_n, "\n")) 
  #  family = "cox" 
  
  if (time==1) { time_start = difftime2() }
  
  # coxcontrol = survival::coxph.control()                                        ## if switch to coxph.fit 
  # glmcontrol = glm.control()                                                    ## if switch to glm.fit 
  
  xs_ncol = dim(xs)[2]                                                          ## number of candidate predictor variables 
  nobs    = dim(xs)[1] 
  one     = c( rep(1, nobs) )
  
  if (is.null(foldid)) { 
    if (is.null(seed)) { seed = round(runif(1)*1000000000) }
    set.seed(seed) 
    foldid = sample( rep( 1:folds_n , ceiling(nobs/folds_n) )[ 1:nobs ] , nobs ) 
    if (family=="cox") {
      length(event)
      event_0_id  = c(1:length(event))[event==0]
      event_1_id  = c(1:length(event))[event==1]
      nobs_0 = length(event_0_id)
      nobs_1 = length(event_1_id)
    
      foldid_0 = sample( rep( 1:folds_n , ceiling(nobs_0/folds_n) )[ 1:nobs_0 ] , nobs_0 ) ;
      foldid_1 = sample( rep( 1:folds_n , ceiling(nobs_1/folds_n) )[ 1:nobs_1 ] , nobs_1 ) ;
      #    length(foldid_0)
      #    length(foldid_1)
      #    table(foldid_0)
      #    table(table(foldid_0))
      #    table(table(foldid_1))
    
      foldid_0id = cbind(event_0_id, foldid_0)
      foldid_1id = cbind(event_1_id, foldid_1)
      foldid_id = rbind(foldid_0id, foldid_1id)
      colnames(foldid_id) = c("id","fold")
      #    table(foldid_id[,1])
      #    table( table(foldid_id[,1]) )
    
      #    table(foldid_id[,2])
      #    table( table(foldid_id[,2]) )
    
      ord = order(foldid_id[,1])
      #    table(table(ord))
      #    ord[1:20]
      
      #    foldid_id[1:20,]
      foldid_id2 = foldid_id[ord,]
      table(table(foldid_id2[,1]))
      table(table(foldid_id2[,2]))
      #    foldid_2[1:20,]
      foldid = foldid_id2[,2]
    }
  }
  
  ##============================================================================  initial panalized fit 
  if (time==1) {
    cat(paste0("  Fitting cv.glmnet() with relax=FALSE to get range for lambda", "\n"))
  }
  if (family=="cox") {                                                          ## TODO1 - move to inside of glmnetr 
    if ( is.null(start) ) { 
      cv.glmnet.fit.f = cv.glmnet( xs, Surv(       y_, event), family="cox", lambda=lambda, relax=FALSE ) 
    } else {         
      cv.glmnet.fit.f = cv.glmnet( xs, Surv(start, y_, event), family="cox", lambda=lambda, relax=FALSE )   
    }
  } else {
    cv.glmnet.fit.f = cv.glmnet( xs, y_ , family=family , relax=FALSE) 
  }
  
  ## limit==1 is generally enough, probably ==0 nice plots for confirmation ## 
  ## limit==3 is just for testing purpose with some "difficult" datasets    ##
  if (fine==1) { limit = max(limit, 2) }
  if        (limit>=3) { lambda.lo = cv.glmnet.fit.f$lambda.min 
  } else if (limit>=2) { lambda.lo = cv.glmnet.fit.f$lambda.min^2 / cv.glmnet.fit.f$lambda.1se^1 
  } else if (limit>=1) { lambda.lo = cv.glmnet.fit.f$lambda.min^3 / cv.glmnet.fit.f$lambda.1se^2 
  } else               { lambda.lo = min(cv.glmnet.fit.f$lambda) }
  lambda    = cv.glmnet.fit.f$lambda[cv.glmnet.fit.f$lambda >= lambda.lo] 
  lambda_n = length(lambda)
  if (lambda_n < 10) {  
    lambda_n = min(10, length(cv.glmnet.fit.f$lambda))
    lambda   = cv.glmnet.fit.f$lambda[c(1:lambda_n)]
  }    
  #  print(lambda) 
  
  #  df.l.rank = sum(cv.glmnet.fit.f$df < rankMatrix(xs)[1])
  #  lamda = lambda[1:min(lambda_n, df.l.rank)]
  
  if (time==1) { time_last = difftime2(time_start) }
  
  if (fine==1) {
    ratio = sqrt(cv.glmnet.fit.f$lambda[2]/cv.glmnet.fit.f$lambda[1])
    lambda = sort( c(lambda, ratio*lambda) , decreasing=TRUE) 
    lambda_n = length(lambda)
    
    if (family=="cox") {                                                          ## TODO1 - move to inside of glmnetr 
      if ( is.null(start) ) { 
        cv.glmnet.fit.f = cv.glmnet( xs, Surv(       y_, event), family="cox", lambda=lambda, relax=FALSE ) 
      } else {         
        cv.glmnet.fit.f = cv.glmnet( xs, Surv(start, y_, event), family="cox", lambda=lambda, relax=FALSE ) 
      }
    } else {
      cv.glmnet.fit.f   = cv.glmnet( xs, y_ , family=family , relax=FALSE) 
    }
  } 
  
  glmnetr.fit = glmnetr( xs, start, y_, event, lambda=lambda, gamma=gamma, object=cv.glmnet.fit.f, time=0, family=family) 
  
#  if (family!="cox") { length(glmnetr.fit$a0g0) }
  
  if (time==1) { cat(paste0("Split for Beta g:0  ")) ; time_last = difftime2(time_start, time_last) }
  
  devratiolist = glmnetr_devratio( cv.glmnet.fit.f, glmnetr.fit, xs_new=xs, start_new=start, y_new=y_, event_new=event, family=family) 
  devratio     = devratiolist$devratio 
  
  if (time==1) { cat(paste0("Split for dev.ratio ")) ; time_last = difftime2(time_start, time_last) }
  
  neventsg0 = 0 
  
  cv_sum_ll_r   = matrix(rep(0,gamma_n*lambda_n), nrow=gamma_n, ncol=lambda_n)
  cv_sum_ll2_r  = matrix(rep(0,gamma_n*lambda_n), nrow=gamma_n, ncol=lambda_n)
  cv_sum_dev_r  = matrix(rep(0,gamma_n*lambda_n), nrow=gamma_n, ncol=lambda_n)
  cv_sum_dev2_r = matrix(rep(0,gamma_n*lambda_n), nrow=gamma_n, ncol=lambda_n)

  lambda_n 
  
  ##============================================================================  begin folds  
  for (i_ in 1:folds_n) {
    ##=== begin folds ==========================================================
    if (time==1) {
      cat(paste0("\n", " ########## Entering Cross Validation fold  ", i_, "  of  " , folds_n , "  ############################" , "\n"))
    }
    ##### set up train and test data sets in matrix form for glmnet & stepreg #####
    trainxs = xs[(foldid!=i_),]   
    testxs  = xs[(foldid==i_),]  
    #    test1 = c(rep(1,dim(testxs)[1])) 
    #    testxs1 = cbind(test1, testxs)                                              ## for models with Intercept 
    
    dim(trainxs)[1] + dim(testxs)[1] 
    dim(trainxs) 
    dim(testxs) 
    
    if ( is.null(start) ) {
      trainstart = NULL
      teststart  = NULL
    } else {
      trainstart = start[(foldid!=i_)]    
      teststart  = start[(foldid==i_)]   
    }
    
    trainy_ = y_[(foldid!=i_)]    
    testy_  = y_[(foldid==i_)]   
    
    if (family == "cox") {
      trainevent = event[(foldid!=i_)] ;  
      testevent  = event[(foldid==i_)] ; 
    } else {
      trainevent = NULL
      testevent  = NULL 
    }    
    
    if ((family=="cox") & is.null(start) ) {
      train_data = data.frame(trainy_, trainevent, trainxs)
      test_data  = data.frame(testy_ , testevent , testxs )
    } else if ((family=="cox") & (!is.null(start)) ) {
      train_data = data.frame(trainstart, trainy_, trainevent, trainxs)
      test_data  = data.frame(teststart , testy_ , testevent , testxs )
    } else {
      train_data = data.frame(trainy_, trainxs)
      test_data  = data.frame(testy_ , testxs )
    }
    
    glmnetr.fit.train = glmnetr( trainxs, trainstart, trainy_, trainevent, lambda, gamma, object=NULL, time=0, family=family) 
    #    names(glmnetr.fit.train) 
    #    colSums(glmnetr.fit.train$betag1!=0)
    #    sum( abs( glmnetr.fit.train$betag0[,9] - glmnetr.fit.train$betag0[,10] ) ) 
    #    sum( abs( glmnetr.fit.train$betag1[,9] - glmnetr.fit.train$betag1[,10] ) ) 
    # object=glmnetr.fit.train ; xs_new=testxs ; start_new=teststart ; y_new=testy_ ; event_new=testevent ; 
    ll_1fold.fit = glmnetrll_1fold( object=glmnetr.fit.train, 
                                    xs_new=testxs, start_new=teststart, y_new=testy_, event_new=testevent, 
                                    family=family, lambda_n=lambda_n, gamma=gamma) 
    
    #    print ( "here 01" )
    #    print ( ll_1fold.fit$cvloglik_r )
    #    print ( ll_1fold.fit$nobs ) 
    if (family %in% c("binomial","gaussian")) { 
      deviance = ll_1fold.fit$cvdevian_r / (ll_1fold.fit$nobs - 1)
    } else { 
      deviance = ll_1fold.fit$cvdevian_r / (ll_1fold.fit$nevent - 1)
    } 
    maxdev = 0 
    for (i_ in c(1:dim(deviance)[1])) {
      for (j_ in c(1:dim(deviance)[2])) {
        if (!is.na(deviance[i_,j_])) { maxdev = max(maxdev,deviance[i_,j_])}
      }
    }
    deviance[is.na(deviance)] = maxdev 
    cv_sum_dev_r  = cv_sum_dev_r  + deviance 
    cv_sum_dev2_r = cv_sum_dev2_r + deviance^2 
    cv_sum_ll_r   = cv_sum_ll_r   +  ll_1fold.fit[[1]] 
    cv_sum_ll2_r  = cv_sum_ll2_r  + (ll_1fold.fit[[1]])^2 
    
    if (family=="cox") {neventsg0 = neventsg0 + (ll_1fold.fit$nevent > 0) }
    
    if ( sum(deviance < 0) > 0 ) { print( sum(testevent)) }
    
    if (time==1) { time_last = difftime2(time_start, time_last) }
    ##=== end folds ============================================================
  } ## end of fold ########################################################################
  
#  print(class(glmnetr.fit$betag1))
#  print(glmnetr.fit$betag1)
#  dimbetag1 = dim( glmnetr.fit$betag1 ) 
  # matrix ( data=rep(0,(6*folds_n)), nrow = folds_n, ncol = 6 ) 

#  the followoing line works in program, but not examples using devtools::check 
#  nzero = colSums( (glmnetr.fit$betag1!=0)*1 )[1:lambda_n]                      ## MOD221129
## therefore the code - Start MOD221129 
  betag1 = as.matrix(glmnetr.fit$betag1)
#  betag1 = matrix(data=rep(0,(dimbetag1[1]*dimbetag1[2]), nrow=dimbetag1[1], ncol=dimbetag1[2])
#  dimbetag1 = dim( glmnetr.fit$betag1 ) 
#  betag1 = matrix(data=rep(0,(dimbetag1[1]*dimbetag1[2])), nrow=dimbetag1[1], ncol=dimbetag1[2])
#  for (im in c(1:dimbetag1[1])) {
#    for (jm in c(1:dimbetag1[2])) {
#      betag1[im,jm] = glmnetr.fit$betag1[im,jm]
#    }
#  }
  nzero = colSums( betag1!=0 )[1:lambda_n]  
# stop mod221129    
  names(nzero) = NULL
  nzero
  
  #=============================================================================
  
  cv.mean_dev_r  = cv_sum_dev_r   / folds_n                        ## / neventg0 
  cv.mean_dev2_r = cv_sum_dev2_r  / folds_n                        ## / neventg0 
  cv.var_dev_r   = cv.mean_dev2_r - cv.mean_dev_r^2 
  
  cvm.dev  = cv.mean_dev_r  
  cvsd.dev = sqrt(cv.var_dev_r/folds_n) 
  cvup.dev = cvm.dev + cvsd.dev 
  cvlo.dev = cvm.dev - cvsd.dev 
  
  statlist <- list()
  for(k_ in 1:gamma_n) {
    statlist_ = list(lambda=lambda, cvm =cvm.dev[k_,] , cvsd=cvsd.dev[k_,], 
                     cvup=cvup.dev[k_,], cvlo=cvlo.dev[k_,],  nzero=nzero)
    statlist[[paste0("g:", gamma[k_])]] <- statlist_
  }
  #  statlist
  
  #=============================================================================
  
  cv.mean_ll_r  = cv_sum_ll_r   / folds_n   
  cv.mean_ll2_r = cv_sum_ll2_r  / folds_n
  cv.var_ll_r   = cv.mean_ll2_r - cv.mean_ll_r^2 
  
  cvm.ll  = cv.mean_ll_r  
  cvsd.ll = sqrt(cv.var_ll_r/folds_n) 
  cvup.ll = cvm.ll + cvsd.ll 
  cvlo.ll = cvm.ll - cvsd.ll 
  
  statlist_ll <- list()
  for(k_ in 1:gamma_n) {
    statlist_ = list(lambda=lambda, cvm=cvm.ll[k_,]  , cvsd=cvsd.ll[k_,], 
                     cvup=cvup.ll[k_,], cvlo=cvlo.ll[k_,], nzero=nzero)
    statlist_ll[[paste0("g:", gamma[k_])]] <- statlist_
  }
  
  ##############################################################################
  ##############################################################################
  
  #  cvm=cvm.ll  ;  cvup = cvup.ll  ;  cvlo = cvlo.ll  ; 
  cvm=cvm.dev ;  cvsd = cvsd.dev ; cvup = cvup.dev ;  cvlo = cvlo.dev ;
  
  min.dev.index = which.min(cvm) 
  max.dev = max(cvm)
  min.dev = min(cvm)
  lambda.min.r.index = ceiling(min.dev.index/gamma_n)
  gamma.min.r.index  = min.dev.index %% gamma_n                                  ## a %% b - remainder of a/b 
  if (gamma.min.r.index == 0) { gamma.min.r.index = gamma_n }
  lambda.min.r = lambda[lambda.min.r.index]
  gamma.min.r  = gamma[gamma.min.r.index]
  nzero.min.r = nzero[lambda.min.r.index]
  
  ##========== get lambda.1se relaxed ==========================================
  
  upper.dev = cvup[ min.dev.index ] 
  
  cvm_t = cvm 
  
  candidatei = (cvm <= upper.dev) 
  candidatei
  
  candidatei[,min(lambda_n, (lambda.min.r.index+1)):lambda_n] = FALSE 
  candidatei
  if (gamma_n > 1) {
    if (gamma.min.r.index > 1) { candidatei[1:(gamma.min.r.index-1),] = FALSE } 
  }
  candidatei
  
  cvm_t[candidatei==FALSE] = min.dev
  cvm_t
  
  onese.index = which.max(cvm_t) 
  onese.index
  
  lambda.1se.r.index = ceiling(onese.index/gamma_n) 
  gamma.1se.r.index  = onese.index %% gamma_n 
  if (gamma.1se.r.index == 0) { gamma.1se.r.index = gamma_n }
  
  lambda.1se.r.index
  gamma.1se.r.index
  
  lambda.1se.r = lambda[lambda.1se.r.index]
  gamma.1se.r  = gamma[gamma.1se.r.index]
  
  nzero.1se.r = nzero[lambda.1se.r.index]
  
  index.r = matrix(rep(0,4), nrow=2, ncol=2)
  index.r[1,1] = lambda.min.r.index
  index.r[1,2] = gamma.min.r.index
  index.r[2,1] = lambda.1se.r.index
  index.r[2,2] = gamma.1se.r.index
  rownames(index.r) = c("min","1se")
  colnames(index.r) = c("Lambda","Gamma")
  index.r
  
  ########### unrelaxed (fully penalized) lambda.min and lambda.1se ############
  
  cvm.g1  = cvm [ gamma_n,]
  #  print(gamma_n) 
  #  print(cvm.g1) 
  lambda.min.g1.index = which.min(cvm.g1) 
  lambda.min.g1 = lambda[lambda.min.g1.index] 
  upper.g1.dev = cvup[gamma_n,lambda.min.g1.index] 
  upper.g1.dev
  
  candidatei = (cvm.g1 <= upper.g1.dev) 
  candidatei
  candidatei[min(lambda_n, (lambda.min.g1.index+1)):lambda_n] = FALSE 
  candidatei
  
  cvm.g1_t = cvm.g1 
  cvm.g1_t[candidatei==FALSE] = min.dev
  cvm.g1_t
  
  lambda.1se.g1.index = which.max(cvm.g1_t)
  lambda.1se.g1.index 
  lambda.1se.g1 = lambda[lambda.1se.g1.index]
  lambda.1se.g1 
  
  index.g1 = matrix(rep(0,2), nrow=2, ncol=1)
  index.g1[1,1] = lambda.min.g1.index
  index.g1[2,1] = lambda.1se.g1.index
  rownames(index.g1) = c("min","1se")
  colnames(index.g1) = c("Lambda")
  index.g1
  
  #  nzero.min.g1 = nzero[lambda.min.g1.index]
  #  nzero.1se.g1 = nzero[lambda.1se.g1.index]
  
  #  c(lambda.min.g1, lambda.min.g1.index, nzero.min.g1 )
  #  c(lambda.1se.g1, lambda.1se.g1.index, nzero.1se.g1 )
  
  ########### COMPLETELY relaxed lambda.min and lambda.1se #####################
  
  cvm.g0 = cvm [1,] 
  lambda.min.g0.index = which.min(cvm.g0) 
  lambda.min.g0 = lambda[lambda.min.g0.index] 
  upper.g0.dev = cvup[1,lambda.min.g0.index] 
  upper.g0.dev
  
  candidatei = (cvm.g0 <= upper.g0.dev) 
  candidatei
  candidatei[min(lambda_n, (lambda.min.g0.index+1)):lambda_n] = FALSE 
  candidatei
  
  cvm.g0_t = cvm.g0 
  cvm.g0_t[candidatei==FALSE] = min.dev
  cvm.g0_t
  
  lambda.1se.g0.index = which.max(cvm.g0_t)
  lambda.1se.g0.index 
  lambda.1se.g0 = lambda[lambda.1se.g0.index]
  lambda.1se.g0 
  
  index.g0 = matrix(rep(0,2), nrow=2, ncol=1)
  index.g0[1,1] = lambda.min.g0.index
  index.g0[2,1] = lambda.1se.g0.index
  rownames(index.g0) = c("min","1se")
  colnames(index.g0) = c("Lambda.g0")
  index.g0
  
  #  nzero.min.g0 = nzero[lambda.min.g0.index]
  #  nzero.1se.g0 = nzero[lambda.1se.g0.index]
  
  #  c(lambda.min.g0, lambda.min.g0.index, nzero.min.g0 )
  #  c(lambda.1se.g0, lambda.1se.g0.index, nzero.1se.g0 )
  
  #  table = cbind(object$nzero[1:lambda_n],devratio[,2],devratio[,1], lambda)
  #  rownames(table) = c(1:lambda_n)
  #  colnames(table) = c("Df", "% Dev", "% Dev R", "Lambda")
  #  table
  
  call <- match.call()
  
  glmnet.fit = list(a0       = glmnetr.fit$a0g1[1:lambda_n] ,
                    beta     = glmnetr.fit$betag1[,1:lambda_n] , 
                    df       = nzero , 
                    dim      = c(xs_ncol, lambda_n) , 
                    lambda   = lambda , 
                    dev.ratio = devratio[,2][1:lambda_n] , 
                    nulldev  = cv.glmnet.fit.f$glmnet.fit$nulldev , 
                    call     = call , 
                    npasses  = NULL ,
                    jerr     = NULL ,
                    offset   = NULL ,
                    nobs     = nobs )
  #  glmnet.fit$beta     = glmnetr.fit$betag1[,1:lambda_n]        ##  cv.glmnet.fit.f$glmnet.fit$beta[,1:lambda_n]    
  #  glmnet.fit$df       = nzero                                  ##  cv.glmnet.fit.f$glmnet.fit$df[1:lambda_n]  
  #  glmnet.fit$lambda   = lambda                                 ##  cv.glmnet.fit.f$glmnet.fit$lambda[1:lambda_n]  
  #  glmnet.fit$dev.ratio = devratio[,2][1:lambda_n]              ##  cv.glmnet.fit.f$glmnet.fit$dev.ratio[1:lambda_n]    
  
  ## the relaxed as func. of gamma and CV results  
  relaxed_base = list( statlist=statlist   , gamma=gamma, 
                       lambda.min = lambda.min.r, lambda.1se = lambda.1se.r, 
                       gamma.min  = gamma.min.r , gamma.1se  = gamma.1se.r , 
                       nzero.min  = nzero.min.r , nzero.1se  = nzero.1se.r , 
                       index=index.r, 
                       lambda.min.g0 = lambda.min.g0 , lambda.1se.g0 = lambda.1se.g0 , index.g0=index.g0, statlist_ll=statlist_ll )
  
  ## the fully relaxed fit    ##beta      = glmnetr.fit$betag0
  glmnet.fit$relaxed =  list( a0        = glmnetr.fit$a0g0[1:lambda_n] ,
                              beta      = glmnetr.fit$betag0[,1:lambda_n] ,
                              df        = nzero ,
                              dim       = dim(glmnetr.fit$betag0) ,
                              lambda    = lambda ,
                              dev.ratio = devratio[,1][1:lambda_n] ,
                              nulldev   = cv.glmnet.fit.f$glmnet.fit$nulldev ,
                              npasses   = cv.glmnet.fit.f$glmnet.fit$npasses ,
                              jerr      = cv.glmnet.fit.f$glmnet.fit$jett ,
                              offset    = cv.glmnet.fit.f$glmnet.fit$offset ,
                              call      = call , 
                              nobs      = nobs )
  
  #  print( lambda )
  #  print( glmnet.fit$lambda )
  #  print( cv.glmnet.fit.f$lambda )
  
  class(glmnet.fit) = "glmnetr" 
  class(cv.glmnet.fit.f) = "list" 
  
  name = "Partial Likelihood Deviance" ; 
  names(name) = "deviance"
  
  returnlist = list(lambda=lambda, cvm=cvm.dev[gamma_n,], cvsd=cvsd.dev[gamma_n,], cvup=cvup.dev[gamma_n,], cvlo=cvlo.dev[gamma_n,], nzero=nzero, 
                    call=call, name=name, glmnet.fit=glmnet.fit, 
                    lambda.min=lambda.min.g1, lambda.1se=lambda.1se.g1, index=index.g1, relaxed=relaxed_base,
                    cv.glmnet.fit.f=cv.glmnet.fit.f, xscolumns= dim(xs)[2], xsdf=rankMatrix(xs)[1], seed=seed, foldid=foldid, family=family)
  class(returnlist) = c("cv.glmnetr") 
  return(returnlist) 
}

###############################################################################################################
###############################################################################################################
## get relaxed glmnet model fit for (whole data) inlcuding for both gamma=0 and gamma=1 #######################
## object = cv.glmnet.fit.f 
## xs = trainxs ; start = trainstart ; y_ = trainy_ ; lambda=NULL ; gamma=c(0,0.25,0.50,0.75,1) ;
## xs_tmp = xs ; start_tmp = start ; y_tmp = y_ ; lambda=NULL ; gamma=c(0,0.25,0.50,0.75,1) ; object=NULL ; time=0 ; family="gaussian" 

#' Fit relaxed part of lasso model
#' 
#' @description
#' This function derives the relaxed lasso fits and optionally calls glmnet() to 
#' derive the fully penalized lasso fit.
#'
#' @param xs_tmp predictor (X) matrix
#' @param start_tmp start time in case Cox model and (Start, Stop) time for use in model  
#' @param y_tmp outcome (Y) variable, in case of Cox model (stop) time 
#' @param event_tmp event variablee in case of Cox model  
#' @param family one of "cox", "gaussian" or "binomial" 
#' @param lambda lambda vector, as in _glmnet_, default is NULL
#' @param gamma gamma vector, as with _glmnet_, default c(0,0.25,0.50,0.75,1) 
#' @param object an output object from _glmnet_ using relax=FALSE with the model fits for
#'    the fully penalized lass fits, i.e. gamma=1.  Defualt is NULL in which case these are derived 
#'    within the function. 
#' @param time Indicate whether or not to update progress in the console.  Default of
#' 0 suppresses these updates.  The option of 1 provides these updates.  In fitting 
#' clinical data with non full rank design matrix we have found some R-packages to
#' take a vary long time or seemingly be caught in infinite loops.  Therefore we allow
#' the user to track the package and judge whether things are moving forward or 
#' if the process should be stopped. 
#'
#' @return A list with two matrices, one for the model coefficients with
#'     gamma=1 and the other with gamma=0.  
#'     
#' @export 
#' 
#' @seealso
#'   \code{\link{cv.glmnetr}} , \code{\link{nested.glmnetr}}  , \code{\link{glmnetr.simdata}} 
#'   
#' @importFrom stats formula glm 
#' @importFrom survival Surv coxph 
#' @importFrom glmnet cv.glmnet 
#'
#' @examples
#' set.seed(82545037)
#' sim.data=glmnetr.simdata(nrows=200, ncols=100, beta=NULL)
#' xs=sim.data$xs 
#' y_=sim.data$yt
#' event=sim.data$event
#' glmnetr.fit = glmnetr( xs, NULL, y_, event, family="cox")
#' plot(glmnetr.fit)
#'
glmnetr = function(xs_tmp, start_tmp, y_tmp, event_tmp, family="cox", lambda=NULL, gamma=c(0,0.25,0.50,0.75,1), object=NULL, time=0) {

#  cv.glmnet.fit.f$glmnet.fit$a0
#  cv.glmnet.fit.f$glmnet.fit$relaxed$a0
  
  if (time==1) { time_start = difftime2() }
  
  xs_tmp_ncol = dim(xs_tmp)[2]
  
  if (!is.null(object)) { 
    cv.glmnet.fit=object
  } else {
    if (family=="cox") {
      if ( is.null(start_tmp) ) { 
        cv.glmnet.fit = cv.glmnet( xs_tmp, Surv(y_tmp, event_tmp), family="cox", lambda=lambda, gamma=gamma, relax=FALSE) 
      } else {
        cv.glmnet.fit = cv.glmnet( xs_tmp, Surv(start_tmp, y_tmp, event_tmp), family="cox", lambda=lambda, gamma=gamma, relax=FALSE ) 
      }
    } else {
      cv.glmnet.fit = cv.glmnet( xs_tmp, y_tmp , family=family, lambda=lambda, relax=FALSE) 
    }
  }

  lambda = cv.glmnet.fit$lambda 
  lambda_n = length( lambda )    
    
  if (time==1) { time_last = difftime2(time_start) }
  
  a0g1 = NULL 
  
  dim(cv.glmnet.fit$glmnet.fit$beta)
  betag1 = cv.glmnet.fit$glmnet.fit$beta[,1:lambda_n]                         ## beta for gamma=1, un-relaxed model 
  if (family != "cox") { a0g1 = cv.glmnet.fit$glmnet.fit$a0[1:lambda_n] }
      
  bnames=list()
  for (col_ in 1:lambda_n) { 
    bname = rownames(betag1)[ (betag1[,col_]!=0) ]
    if (col_==1) { bnames = list(bname) }                                       ## print(bnames)
    else       { bnames = c(bnames, list(bname)) }      
  }
    
  a0g0 = NULL
  
  #---------------------------------------------------------------------------
  for (col_ in 1:lambda_n) { 
    bname = bnames[[col_]]  
    bname
    if (length(bname) == 0) { 
      bnamenone = TRUE 
      beta = rep(0,dim(betag1)[1])
      if (family != "cox") { a0 = a0g1[col_] } 
    } else { bnamenone = FALSE }                                              # print("No predictors") 
    if (col_ == 1) { same = FALSE 
    } else {
      if ( length(bname) != length(bnamelast) ) {
        same = FALSE
      } else if (sum((bname != bnamelast) != 0)) {
        same = FALSE 
      }
    }
      
    if ((same==FALSE ) & (bnamenone==FALSE)) { 
      if (family=="cox") {
        if ( is.null(start_tmp) ) { 
          if (bnamenone) { 
            dataf = as.data.frame( cbind( y_tmp, event_tmp, one=c(1) ) )
            form2 = formula( paste("Surv(" , colnames(dataf)[1] ,  ", " , colnames(dataf)[2] , ") ~ one ") )  
          } else { 
            dataf = as.data.frame( cbind( y_tmp, event_tmp, xs_tmp ) ) 
            form2 = formula( paste("Surv(" , colnames(dataf)[1] ,  ", " , colnames(dataf)[2] , ") ~ ", paste(bname, collapse = " + " ) ) )  
          }
        } else { 
          if (bnamenone) { 
            dataf = as.data.frame( cbind( start_tmp, y_tmp, event_tmp, one=c(1) ) )
            form2 = formula( paste("Surv(" , colnames(dataf)[1] ,  ", " , colnames(dataf)[2],  ", " , colnames(dataf)[3]  , ") ~ one") ) 
          } else {
            dataf = as.data.frame( cbind( start_tmp, y_tmp, event_tmp, xs_tmp ) )
            form2 = formula( paste("Surv(" , colnames(dataf)[1] ,  ", " , colnames(dataf)[2],  ", " , colnames(dataf)[3]  , ") ~ ", paste(bname, collapse = " + " ) ) ) 
            }
          }
          fit2 = coxph( form2 , dataf ) 
      } else {
        if (bnamenone) { 
          dataf = as.data.frame( cbind( y_tmp, one=c(1) ) )
          form2 = formula( paste( colnames(dataf)[1] , " ~ one ") )  
        } else {
          dataf = as.data.frame( cbind( y_tmp,  xs_tmp ) ) 
          form2 = formula( paste( colnames(dataf)[1] ,  " ~ ", paste(bname, collapse = " + " ) ) )  
        }
        fit2 = glm( form2 , dataf, family=family) 
      }
      
       xs_tmpnames = colnames(xs_tmp)
       xs_tmpnames
#       xs_tmpk = length(xs_tmpnames)
       beta = rep(0,xs_tmp_ncol)  
       if (bnamenone == FALSE) {
         bnamep = as.numeric( (xs_tmpnames %in% bname ) )                       ## present(indicators)
         bnamep                                                                 ## same as riskvarsnewp
         bnamei = c(1:xs_tmp_ncol)[ (bnamep==1) ]                               ## numerical index
         bnamei
         if (family == "cox") { 
           beta[ bnamei ] = fit2$coefficients
         } else { 
           a0 = fit2$coefficients[ 1 ]
           beta[ bnamei ] = fit2$coefficients[-1]
         } 
         beta 
       }
    }  
    
    if  (col_ == 1) { 
      betag0 = beta  
      if (family != "cox") { a0g0 = a0 }
    } else { 
      betag0 = cbind(betag0, beta) 
      if (family != "cox") { a0g0 = cbind(a0g0, a0) }
    }
    bnamelast = bname 
  }
  #---------------------------------------------------------------------------
    
  if (time==1) { time_last = difftime2(time_start, time_last) }

  rownames(betag0) = xs_tmpnames 
  colnames(betag0) = c(1:lambda_n)
  dim( betag0 )
  betag0[is.na(betag0)] = 0 
  if (family != "cox") {
    a0g0 = as.vector(a0g0)
    names(a0g0) = names(a0g1) 
  }
#   colSums((betag0 != 0)*1)

  returnlist = list( lambda=lambda, gamma=gamma, betag0=betag0, betag1=betag1, a0g0=a0g0, a0g1=a0g1, cv.glmnet.fit=cv.glmnet.fit, family=family )
  class(returnlist) = "glmnetr"  
  return(returnlist)
}

###############################################################################################################
## get log likelihood for various model fits on the relaxed (1-gamma)*betag0 + gamma*betag ####################
## object=glmnetr.fit.train ; xs_new=testxs ; start_new=teststart ; y_new=testy_ ; event_new=testevent; lambda_n=NULL; lambda=lambda; gamma=gamma ; family="gaussian" ; 

#' Evaluation fit of leave out fold
#'
#' @description
#' This function derives the log likelihood for a leave out based upon the fit
#' of the input object.
#' 
#' @param object an output object from _cv.glmnetr_  
#' @param xs_new A new predictor matrix 
#' @param start_new A new vector of start times or the Cox model.  May be NULL.
#' @param y_new  a new outcome vector.
#' @param event_new event vector in case of the Cox model.  May be NULL for other models.  
#' @param family Model family, one of "cox", "gaussian" or "binomial".  
#' @param lambda_n length of the lambda vector.  
#' @param gamma  The gamma vector.  
#'
#' @return Returns the log likelihood of object fit using new data.  
#' 
#' @importFrom stats logLik glm residuals deviance  
#' @importFrom survival coxph coxph.control  
#'
glmnetrll_1fold = function(object, xs_new, start_new, y_new, event_new, family="cox", 
                           lambda_n=NULL, gamma=c(0,0.25,0.50,0.75,1)) {

  nobs = dim(xs_new)[1]
  
#  dim(object$betag0) ; dim(object$betag1) ; dim(xs_new)
  #  lambda_n = length(lambda) 
  lambda_n_local = dim(object$betag0)[2] 
  if (is.null(lambda_n)) { lambda_n = lambda_n_local }                          ## lambda_n gotten from main model 
  gamma_n = length(gamma)

  cvloglik_r    = matrix(rep(0,gamma_n*lambda_n), nrow=gamma_n, ncol=lambda_n)
  cvdevian_r    = matrix(rep(0,gamma_n*lambda_n), nrow=gamma_n, ncol=lambda_n)
  
  one = rep(1, length(y_new)) 
  if (family == "cox") {
    if (is.null(start_new)) {
      fitnull = coxph(Surv(y_new, event_new)~one, init = c(0), control=coxph.control(iter.max=0)) 
    } else {
      fitnull = coxph(Surv(start_new, y_new, event_new)~one, init = c(0), control=coxph.control(iter.max=0)) 
    }
    loglik_null = fitnull$loglik[1]
    loglik_sat = cox.sat.dev(y_new, event_new) 
  } else if (family=="binomial") { 
    p_ = sum(y_new)/length(y_new)
    loglik_null = sum( log(p_)*y_new + log(1-p_)*(1-y_new) )
    loglik_sat  = 0 
  } else if (family == "gaussian") { 
    fitnull = glm(y_new ~ one, family=family) 
    loglik_null = -fitnull$null.deviance/2
    loglik_sat  = 0 
  } else {
    fitnull = glm(y_new ~ one, family=family) 
    resid = residuals(fitnull, type="deviance")
    fit_sat = glm(y_new ~ resid, family=family) 
    loglik_null = -fitnull$null.deviance/2
    loglik_sat  = -fit_sat$deviance/2
  }
    
#  loglik_sat =  loglik_null + 0.5 * nulldev 
  
  if (family == "cox") { nevent = fitnull$nevent 
  } else if (family == "binomial") { nevent = sum(y_new)
  } else { nevent = NULL }
  
  if (lambda_n_local < lambda_n) { cat(paste0("lambda_n_local=", lambda_n_local, " < lambda_n=", lambda_n, "\n" )) }
#  print(" ")
#  print("Here in glmnetrll_1fold")
#  print(lambda_n)
#  print(dim(object$betag0))
#  print(dim(object$betag1))
  for (ja_ in 1:lambda_n_local) {
    for (ia_ in 1:gamma_n) {
      gam = gamma[ia_]
      beta_   = (1-gam) * object$betag0[,ja_] + gam * object$betag1[,ja_]
#      length(beta_)  ; class (beta_) ; 
#      beta_[beta_!=0]      
      score_new = xs_new %*% beta_ 
      if (family != "cox") { score_new = score_new + (1-gam)*object$a0g0[ja_] + gam*object$a0g1[ja_] }
      summary(score_new)
      if (family=="cox") { 
        if (is.null(start_new)) {
          fit1 = coxph(Surv(y_new, event_new)~score_new    , init = c(1), control=coxph.control(iter.max=20)) 
        } else {
          fit1 = coxph(Surv(start_new, y_new, event_new)~score_new    , init = c(1), control=coxph.control(iter.max=20)) 
        }
#        cvdevian_r[ia_,ja_] = coxnet.deviance(score_new,Surv(y_new,event_new))                
        cvdevian_r[ia_,ja_] = 2*(loglik_sat[1] - fit1$loglik[1]) 
      } else if (family=="gaussian") {
        cvdevian_r[ia_,ja_] = sum((y_new-score_new)^2) ; # var(y_new) ; var(score_new) ; 
      } else if (family=="binomial") {
        p_ = 1/(1+exp(-score_new)) 
        cvdevian_r[ia_,ja_] = -2*( t(log(p_))%*%y_new + t(log(1-p_))%*%(1-y_new) )
      } else {
#        fit1 = glm( y_new ~ score_new , family=family, start=c(0,1) ) 
        fit1 = glm( y_new ~ score_new , family=family ) 
        fit1$null.deviance
        fit1$deviance
        cvloglik_r[ia_,ja_] = logLik(fit1)[1] 
        cvdevian_r[ia_,ja_] = deviance(fit1)  
      }
    }
  }
  
  if (lambda_n_local < lambda_n) {
    for (ja_ in c((lambda_n_local+1):lambda_n)) {
      for (ia_ in 1:gamma_n) {
        cvloglik_r   [ia_,ja_] = max(cvloglik_r   [ia_,])
      }
    }
  }
  
  returnlist = list(cvloglik_r=cvloglik_r, cvdevian_r=cvdevian_r,  
                    cvloglik_sat=loglik_sat, cvlogliknull=loglik_null, nevent=nevent, nobs=nobs) 
  return( returnlist )
}

###################################################################################################
###################################################################################################
## object=cv.glmnet.fit.f ; object2=glmnetr.fit ; xs_new = xs ; start_new = start ; y_new = y_ ; event_new = event ; 

#' Get Deviance ratio.  
#' 
#' @description fit models to derive the devaince ratios.  
#'
#' @param object  - A _glmnet_ output object with relax=FALSE, i.e model fit for gamma=1. 
#' @param object2 - A _glmnetr_ output object with relaxed fits, i.e model fit for gamma=0. 
#' @param xs_new  - predictor matrix 
#' @param start_new - start times in case of usage in Cox model.  May be NULL.  
#' @param y_new   - outcome vector. 
#' @param event_new - event indicator in case of Cox model.  
#' @param family  - Model family, one of "cox", "gaussian" or "binomial".   
#'
#' @return - Deviance ratios.  
#' 
# @importFrom survival coxph coxph.control residuals.coxph 
# @importFrom stats glm logLik residuals.glm 
#' @importFrom stats glm logLik residuals 
#' @importFrom survival coxph coxph.control  
#' 
# object=cv.glmnet.fit.f ; object2=glmnetr.fit ; xs_new=xs ; start_new=start ; y_new=y_ ; event_new=event ; family=family ; 
glmnetr_devratio = function( object, object2, xs_new, start_new, y_new, event_new, family) {

  one = rep(1, length(y_new))
  
  if (family == "cox") {
    if (is.null(start_new)) {
      fitnull = coxph(Surv(y_new, event_new)~one, init = c(0), control=coxph.control(iter.max=0)) 
    } else {
      fitnull = coxph(Surv(start_new, y_new, event_new)~one, init = c(0), control=coxph.control(iter.max=0)) 
    } 
    loglik_null = fitnull$loglik[1]
    loglik_sat = cox.sat.dev(y_new, event_new) 
    nulldev     = 2*(loglik_sat - loglik_null)
  } else if (family == "binomial") {
    fitnull = glm(y_new ~ one ) 
    loglik_null = logLik(fitnull)[1]
    loglik_sat  = 0 
    nulldev = fitnull$null.deviance 
  }else if (family == "gaussian") {
    fitnull = glm(y_new ~ one ) 
    loglik_null = -sum((y_new-mean(y_new))^2)/2                                 ## not actually the log likelihood
    loglik_sat  = 0 
    nulldev = fitnull$null.deviance 
  } else {
    fitnull = glm(y_new ~ one ) 
    resid   = residuals(fitnull, type="deviance")
    fit_sat  = glm(y_new ~ resid) 
    loglik_null = logLik(fitnull)[1]
    loglik_sat  = logLik(fit_sat)[1]    
    nulldev = fit_sat$null.deviance 
  }

##  loglik_sat = 0                         ## for no ties data only Cox model  
##  nulldev = -2*loglik_null               ## for no ties data only Cox model
##  dev = 2*(loglik_sat - loglik)          ## general definition of model deviance

  lambda_n  = dim(object2$betag0)[2] 
  dev       = matrix(rep(0,2*lambda_n), nrow=lambda_n, ncol=2)
  devratio  = matrix(rep(0,2*lambda_n), nrow=lambda_n, ncol=2)
  
  for (gam in 0:1) {
    for (ja_ in 1:lambda_n) {
      beta_ = (1-gam) * object2$betag0[,ja_] + gam * object2$betag1[,ja_]
      score = xs_new %*% beta_ 
      if (family != "cox") { score = score + (1-gam)*object2$a0g0[ja_] + gam*object2$a0g1[ja_] }
      if        ( ( family=="cox") &   is.null(start_new)  ) { 
        fit1 = coxph(Surv(y_new, event_new)~score, init = c(1),control=coxph.control(iter.max=0)) 
      } else if ( ( family=="cox") & (!is.null(start_new)) ) { 
        fit1 = coxph(Surv(start_new, y_new, event_new)~score, init = c(1),control=coxph.control(iter.max=0)) 
      } else if (family=="binomial") {
        p_ = 1/(1+exp(-score)) 
        loglik = ( t(log(p_))%*%y_new + t(log(1-p_))%*%(1-y_new) )
      } else if (family=="gaussian") {
        loglik = -sum((y_new-score)^2)/2                                        ## not actually the log likelihood
      } else {
        fit1 = glm( y_new ~ score )    
      } 
      if (family=="cox") { 
        dev[ja_, gam+1] = 2*(loglik_sat[1] - fit1$loglik[1])
      } else if (family %in% c("binomial","gaussian")) {
        dev[ja_, gam+1] = - 2*loglik 
      } else {  
        dev[ja_, gam+1] = fit1$deviance 
      } 
    }
  }
  devratio = 100*(1 -dev/nulldev) 
  devratio = round(devratio, digits=2) 
  colnames(devratio) = c("DR g:0","DR g:1") 
  returnlist= list(devratio=devratio, loglik_sat=loglik_sat, loglik_null=loglik_null, nulldev=nulldev) 
  return(returnlist) 
}

# table = cbind(object$nzero[1:lambda_n],devratiop[,2],devratiop[,1], lambda)
# rownames(table) = c(1:lambda_n)
# colnames(table) = c("Df", "% Dev", "% Dev R", "Lambda")
# table

# cv.glmnet.fit$glmnet.fit$dev.ratio

# object$glmnet.fit
# object$relaxed$glmnet.fit
# class(object$relaxed$glmnet.fit)
# typeof(object$relaxed$glmnet.fit)

###############################################################################################################
###############################################################################################################

###############################################################################################################
###############################################################################################################
#' Generate example data 
#' 
#' @description
#' This function generates an example data set with specified number of observations, 
#' and predictors.  The first column in teh dsing matrix is identically equal to 1
#' for an intercept.  Columns 2 to 5 are for the 4 levels of a character variable, 
#' 6 to 11 for the 6 levels of a character variable.  Columns 12 to 17 are for 3
#' binomial predictors, again over paramtereized.  Such over paramtereization
#' can cause difficulties with the _glmnet()_ of the glmnet package.  
#'
#' @param nrows Sample size (>=100) for simulated data, default=1000.
#' @param ncols Number of columns (>=17) in design matrix, i.e. predictors, default=100. 
#' @param beta Vector of length <= ncols for "left most" coefficients.  If beta has 
#' length < ncols, then the values at length(beta)+1 to ncols are set to 0.
#' Default=NULL, where a beta of length 25 is assigned standard normal values.
#'
#' @return A list with elements xs for desing matrix, y_ for a quantitative outcome, 
#' yt for a survival time, event for an indicator of event (1) or censoring (0), 
#' in the Cox proportional hazards survival model setting, yb for 
#' yes/no (binomial) outcome data, and beta the beta used in random number generation. 
#' @export
#' 
#' @importFrom stats rnorm runif rbeta 
#' @importFrom Matrix rankMatrix
#' 
#' @seealso
#'   \code{\link{glmnetr}} , \code{\link{cv.glmnetr}} , \code{\link{nested.glmnetr}}  
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
glmnetr.simdata = function(nrows=1000, ncols=100, beta=NULL) {
  if (nrows<100) {nrows=100}
  if (ncols<17) {ncols=17}
  if (is.null(beta)) {
    beta = c(rep(0,ncols)) 
    beta[1:25]=rnorm(25)
  } else {
    beta_ = c(rep(0,ncols)) 
    beta_[1:length(beta)]=beta
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
  y_ = as.vector( xs %*% beta ) + rnorm(nrows) 
  ct = as.vector( -log(runif(nrows)) )
  yt = exp(y_) 
  event = as.vector( (yt < ct)*1 ) 
  yt    = as.vector( event*yt + (1-event)*ct ) 
  yb    = as.vector( (y_ > 1)*1 ) 
  runis = (runif(nrows) < 0.3)*1 
  start = as.vector( rbeta(nrows,2,4) * yt * runis )
  
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
  return(list(xs=xs,start=start,y_=y_, yt=yt, event=event,yb=yb, beta=beta))
}

###############################################################################################################
###############################################################################################################

#' Calculate the CoxPH saturated log-likelihood
#'
#' @description This function calculates the saturated log-likelihood for the Cox 
#' model using both the Efron and Breslow approximations for the case where all ties
#' at a common event time have the same weights (exp(X%*%B)).  For the simple case 
#' without ties the saturated log-likelihood is 0 as the contribution to the log-likelihood 
#' at that time point can be made arbitrarily close to 1 by assigning a large weight 
#' to the record corresponding to an event.  Similarly, in the 
#' case of ties one can assign a much larger weight to be associated with one of the
#' event times such that the associated record contributes a 1 to the likelihood.  Next one can 
#' assigns a very large weight to a second tie, but smaller than the first tie considered,
#' and this too will contribute a 1 to the likelihood.  Continuing in this way for this
#' and all time points with ties, 
#' the partial log-likelihood is 0, just like for the no-ties case.  Note, this is 
#' the same argument with which we derive the log-likelihood of 0 for the no ties
#' case. Still, to be consistent with others we derive the saturated log-likelihood
#' with ties under the constraint that all ties carry the same weights.
#'
#' @param y_ Time variable for a survival analysis, whether or not there is a start time
#' @param e_ Event indicator with 1 for event 0 otherwise. 
#'
#' @return Saturated log likelihood for the Efron and Breslow approximations.
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