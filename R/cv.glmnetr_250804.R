################################################################################
##### cv.glmnetr_yymmdd.R ######################################################
################################################################################
## do cross validation to choose tuning parameter, i.e. lambda and gamma #######
#' Get a cross validation informed relaxed lasso model fit. Available to the 
#' user but intended to be call from nested.glmnetr(). 
#'
#' @description
#' Derive an elastic net (including a relaxed lasso) model and identify 
#' hyperparameters, i.e. alpha, gamma and lambda, which give the best fit 
#' based upon cross validation.  It is analogous to (and uses) the cv.glmnet() 
#' function of the 'glmnet' package, but also tunes on alpha. 
#' 
#' @param trainxs predictor matrix 
#' @param trainy__ outcome vector
#' @param family  model family, "cox", "binomial" or "gaussian" (default) 
#' @param lambda the lambda vector.  May be NULL.  
#' @param gamma the gamma vector.  Default is c(0,0.25,0.50,0.75,1). 
#' @param alpha A vector for alpha values considetered when tuning, for example
#' c(0,0.2,0.4,0.6,0.8,1). Default is c(1) to fit the 
#' lasso model involving only the L1 penalty. The value for alpha of 0 corresponds to the 
#' ridge model involving only the L2 penalty.
#' @param folds_n number of folds for cross validation.  Default and generally recommended is 10.
#' @param fine use a finer step in determining lambda.  Of little value unless one 
#' repeats the cross validation many times to more finely tune the hyperparameters.  
#' See the 'glmnet' package documentation.  
#' @param path The path option from cv.glmnet(). 0 for FALSE reducing 
#' computation time when the numerics are stable, 1 to avoid cases where 
#' the path = 0 option might get very slow. 0 by default.
#' @param foldid a vector of integers to associate each record to a fold.  The integers should be between 1 and folds_n.
#' @param track indicate whether or not to update progress in the console.  Default of
#' 0 suppresses these updates.  The option of 1 provides these updates.  In fitting 
#' clinical data with non full rank design matrix we have found some R-packages to
#' take a vary long time or seemingly be caught in infinite loops.  Therefore we allow
#' the user to track the program progress and judge whether things are moving forward or 
#' if the process should be stopped. 
#' @param ... Additional arguments that can be passed to cv.glmnet() 
#' 
#' @details
#' This is the main program for model derivation.  As currently implemented the 
#' package requires the data to be input as vectors and matrices with no missing values 
#' (NA).  All data vectors and matrices must be numerical.  For factors (categorical variables) one
#' should first construct corresponding numerical variables to represent the factor 
#' levels.  To take advantage of the lasso model, one can use one hot coding
#' assigning an indicator for each level of each categorical variable, or creating 
#' as well other contrasts variables suggested by the subject matter.
#'  
#' @return A cross validation informed relaxed lasso model fit.
#' 
#' @seealso
#'   \code{\link{summary.cv.glmnetr}} , \code{\link{predict.cv.glmnetr}} , \code{\link{nested.glmnetr}} 
#'   
#' @author Walter Kremers (kremers.walter@mayo.edu)
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
#' sim.data=glmnetr.simdata(nrows=200, ncols=100, beta=NULL)
#' xs=sim.data$xs 
#' y_=sim.data$y_ 
#' event=sim.data$event
#' # for this example we use a small number for folds_n to shorten run time 
#' cv.glmnetr.fit = nested.glmnetr(xs, y_=y_, family="gaussian", folds_n=4, resample=0) 
#' plot(cv.glmnetr.fit)
#' plot(cv.glmnetr.fit, coefs=1)
#' summary(cv.glmnetr.fit)
#' 
cv.glmnetr = function(trainxs, trainy__, family, alpha=1, gamma=c(0,0.25,0.5,0.75,1), lambda=NULL, foldid=NULL, folds_n=NULL, fine=0, path=0, track=0, ... ) { 
  if (track >= 1) { cat( "  ... begin cv.glmnetr   ") }
  if (track >= 2) { time_start = diff_time(NULL, NULL) ; time_last=time_start ; 
  } else if (track == 1) { cat("\n") } 
  lambda_null = 0 
  if (is.null(lambda)) {
    if (track >= 1) { cat( "  ... search for lambda vector   ") }
    if (track == 1) { cat( "\n" ) }
    lambda_null = 1 
    
    temp_fit_1       = cv.glmnet( trainxs, trainy__, family=family, alpha=1, foldid=foldid, nfolds=folds_n , ... )   ## L1 penalty 
    cv_glmnet_fit_a0 = cv.glmnet( trainxs, trainy__, family=family, alpha=0, foldid=foldid, nfolds=folds_n , ... )   ## L2 penalty 
    if (fine == 1) { 
      lambda  = sort( unique( c(temp_fit_1$lambda, cv_glmnet_fit_a0$lambda) ) , decreasing = TRUE )
    } else if (fine == 2) { 
      lambda  = sort( temp_fit_1$lambda , decreasing = TRUE )
    } else {
      lambda2 = sort( cv_glmnet_fit_a0$lambda[ ( (cv_glmnet_fit_a0$lambda < min(temp_fit_1$lambda)) |
                                                   (cv_glmnet_fit_a0$lambda > max(temp_fit_1$lambda)) ) ] )
      lambda  = sort( unique( c(temp_fit_1$lambda, lambda2) ) , decreasing = TRUE) 
    } 
    
    if (track >= 2) { time_last = diff_time(time_start, time_last) }
  }
  
  cv_glmnet_fit_ = list()
  lalpha = length(alpha)
  for (j_ in c(1:lalpha)) {
    if (track >= 1) { cat( "  fitting model for alpha =", alpha[j_], "   ") }
    if (track == 1) { cat( "\n" ) }
    #      print(c(j_,alpha[j_]))
    alpha_ = alpha[j_]
    cv_glmnet_fit_[[j_]] = cv.glmnet( trainxs, trainy__, family=family, alpha=alpha_, lambda=lambda, gamma=gamma, relax=TRUE, foldid=foldid, nfolds=folds_n, path=path, ... ) 
    class(cv_glmnet_fit_[[j_]]) = "cv.glmnetr" 
    cv_glmnet_fit_[[j_]]$alpha  = alpha_
    cv_glmnet_fit_[[j_]]$family = family
    
    cv_glmnet_fit_[[j_]]$nzero.min   = cv_glmnet_fit_[[j_]]$nzero[cv_glmnet_fit_[[j_]]$index[1]]      ## CHECK ##
    cv_glmnet_fit_[[j_]]$measure.min = cv_glmnet_fit_[[j_]]$cvm[  cv_glmnet_fit_[[j_]]$index[1]]      ## CHECK ##
    
    index.min.g0 = which.min( cv_glmnet_fit_[[j_]]$relaxed$statlist[[ 1 ]]$cvm )
    cv_glmnet_fit_[[j_]]$relaxed$index.min.g0   = index.min.g0
    cv_glmnet_fit_[[j_]]$relaxed$nzero.min.g0   = cv_glmnet_fit_[[j_]]$nzero[ index.min.g0 ] 
    cv_glmnet_fit_[[j_]]$relaxed$lambda.min.g0  = cv_glmnet_fit_[[j_]]$lambda[ index.min.g0 ]
    cv_glmnet_fit_[[j_]]$relaxed$measure.min.g0 = cv_glmnet_fit_[[j_]]$relaxed$statlist[[1]]$cvm[ index.min.g0 ]
    
    indx = cv_glmnet_fit_[[j_]]$relaxed$index
    cvm.min_ =  cv_glmnet_fit_[[j_]]$relaxed$statlist[[ indx[1,2] ]]$cvm[ indx[1,1] ]
    if (j_ == 1) { 
      a.index.min = 1 
      cvm.min = cvm.min_ 
    } else if (cvm.min_ <= cvm.min ) {
      a.index.min = j_ 
      cvm.min = cvm.min_
    }
    if (track >= 2) { time_last = diff_time(time_start, time_last) } 
  }
  
  class(cv_glmnet_fit_) = "cv.glmnetr.list" 
  cv_glmnet_fit_$family = family 
  cv_glmnet_fit_$alpha  = alpha  
  
  cv_elastic_fit = cv_glmnet_fit_[[a.index.min]]
  class(cv_elastic_fit) = "cv.glmnetr.el"
  indx = cv_elastic_fit$relaxed$index
  index.elastic = c(      a.index.min ,       indx[1,2] ,                          indx[1,1]  ) 
  vals.elastic  = c(alpha[a.index.min], gamma[indx[1,2]], cv_elastic_fit$lambda[ indx[1,1] ], cvm.min , cv_elastic_fit$relaxed$nzero.min)
  names(index.elastic) = c("alpha", "gamma", "lambda")
  names(vals.elastic) = c("alpha", "gamma", "lambda", "measure", "NZero")
  #    print(index.elastic)
  cv_elastic_fit$index.elastic = index.elastic
  cv_elastic_fit$vals.elastic = vals.elastic
  
  cv_lasso_fit = cv_glmnet_fit_[[ lalpha ]] 
  
  cv_ridge_fit = cv_glmnet_fit_a0
  cv_ridge_fit$family = family 
  
  #    cv_elastic_fit
  #    cv_lasso_fit
  #    cv_ridge_fit
  
  if (track == 1) { print(vals.elastic) ; cat( "  ... finish cv.glmnetr \n") ; } 
  if (track >= 2) { print(vals.elastic) ; time_last = diff_time(time_start, time_last) ; cat( "  ... finish cv.glmnetr \n") ; } 
  
  return( list( cv_glmnet_fit_ = cv_glmnet_fit_ , cv_lasso_fit   = cv_lasso_fit   , 
                cv_ridge_fit   = cv_ridge_fit   , cv_elastic_fit = cv_elastic_fit , 
                lambda=lambda, lambda_null=lambda_null ) ) 
}
