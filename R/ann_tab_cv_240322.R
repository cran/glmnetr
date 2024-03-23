################################################################################
##### ann_tab_cv_yymmdd.R ######################################################
################################################################################
#' calculate cross-entry for multinomial outcomes
#'
#' @param xx the sigmoid of the link, i.e, the estimated probabilities, i.e. xx = 1/(1+exp(-xb)) 
#' @param yy the observed data as 0's and 1's 
#'
#' @return the cross-entropy on a per observation basis 
#' @export
#'
calceloss = function (xx,yy) {
  celoss = 1  
  if (is.vector(xx)) {
    if ((min(xx) < 0) | (max(xx) > 1)) {   ## this means xx is the xbeta and not the prob 
      xx = 1/(1+exp(-xx)) 
      xx[xx<1e-6] = 1e-6
      xx[xx > (1-1e-6)] = (1-1e-6) 
    }
#    celoss = - (t(log(xx)) %*% yy + t(log(1 - xx)) %*% (1 - yy) ) / length(xx)
    celoss = - ( sum(log(xx) * yy) + sum(log(1 - xx) * (1 - yy)) ) / length(xx)
  } else {
  yy = as.numeric(yy)
  ywide = matrix(rep(0,dim(xx)[1]*dim(xx)[2]),nrow=dim(xx)[1],ncol=dim(xx)[2])
  for (i_ in c(1:dim(xx)[1])) { ywide[i_,yy[i_]] = 1 }
  ywide
  exptemp = exp(xx)
  rowsum = exptemp %*% c(rep(1,dim(xx)[2]))
  ptemp = exptemp 
  for (i_ in c(1:dim(xx)[2])) {  ptemp[,i_] = exptemp[,i_]/rowsum }
  ywide*log(ptemp)
  celoss = - sum( ywide*log(ptemp) ) / dim(xx)[1] 
  }
  return( celoss ) 
}

################################################################################
################################################################################
#' 
#' Construct the weights for going from the observed data with an offset in column 1
#' to the first hidden layer 
#'
#' @param tnsr an input tensor which is to be modified to mimic the linear term of a 
#' generalized linear model, e.g a Cox or logistic regression model 
#' @param lasso 1 if the first column is the linear estimate from a linear model, 
#' often a lasso model 
#' @param lscale Scale used to allow ReLU to exend +/- lscale before capping the 
#' inputted linear estimated 
#' @param scale Scale used to transform the inital random paramter assingments by 
#' dividing by scale
#' @param trnspose 1 to transpose the matrix before returning, 0 to not. 
#' @param rreturn 1 (default) to return an R (numeric) vector, 0 to return a torch tensor 
#'
#' @return a weight matrix in tensor format
#'
wtzero = function(tnsr, lasso=0, lscale=5, scale=1, rreturn=1, trnspose=0 ) { 
  weight_r = as.matrix(tnsr)
  if (scale <= 0) { 
    weight_r = matrix( rep(0 ,dim(weight_r)[1] * dim(weight_r)[2]), nrow=dim(weight_r)[1], ncol=dim(weight_r)[2] ) 
  } else if (scale > 1) {
    weight_r = weight_r / scale 
  } 
  if (lasso != 0) {
    dmw = dim(weight_r) 
    row0 = c(rep(0,dmw[2]))  
    col0 = c(rep(0,dmw[1]))  
    weight_r[1,] = row0
    weight_r[2,] = row0
    weight_r[,1] = col0
    weight_r[1,1] =  1/lscale  
    weight_r[2,1] = -1/lscale  
    weight_r  
  }
  if (trnspose==1) { weight_r = t(weight_r) }
  if (rreturn ==1) { 
    weight = weight_r 
  } else {
    weight = torch_tensor(weight_r, dtype=torch_float(), requires_grad=TRUE)
  }
  return(weight)
}

################################################################################
#'
#' Construct the weights for going between two hidden layers, carrying forward an
#' offset term to mimic a linear model 
#' 
#' @param tnsr an input tensor wihc is to be modified to mimic the linear term of a 
#' generalized linear model, e.g a Cox or logistic regression model 
#' @param lasso 1 if the first column is the lienar estimate from a linear model, 
#' often a lasso model 
#' @param rreturn 1 (default) to return an R (numeric) vector, 0 to return a torch tensor 
#' @param trnspose 1 to transpose the matrix before returning, 0 to not. 
#'
#' @return a weight matrix in tensor format
#'
wtmiddle = function(tnsr, lasso=0, rreturn=1, trnspose=0 ) { 
  weight_r = as.matrix(tnsr)
  if (lasso < 0) { 
    weight_r = matrix( rep(0 ,dim(weight_r)[1] * dim(weight_r)[2]), nrow=dim(weight_r)[1], ncol=dim(weight_r)[2] ) 
  } else if (lasso > 1) {
    weight_r = weight_r / lasso 
  } 
  if (lasso != 0) {
    dmw = dim(weight_r) 
    row0 = c(rep(0,dmw[2]))  
    col0 = c(rep(0,dmw[1]))  
    weight_r[1,] = row0
    weight_r[2,] = row0 
    weight_r[,1] = col0
    weight_r[,2] = col0
    weight_r[1,1] = 1 
    weight_r[2,2] = 1 
    weight_r  
  }
  if (trnspose==1) { weight_r = t(weight_r) }
  if (rreturn ==1) { 
    weight = weight_r 
  } else {
    weight = torch_tensor(weight_r, dtype=torch_float(), requires_grad=TRUE)
  }
  return(weight)
}

################################################################################
#'
#' Construct the weights for going from the last hidden layer to the last layer of 
#' the model, not counting any activation, to carry forward an offset to mimic a linear model 
#'
#' @param tnsr an input tensor which is to be modified to mimic the linear term of a 
#' generalized linear model, e.g a Cox or logistic regression model 
#' @param lasso 1 if the first column is the linear estimate from a linear model, 
#' often a lasso model 
#' @param lscale Scale used to allow ReLU to exend +/- lscale before capping the 
#' inputted linear estimated 
#' @param scale Scale used to transform the inital random paramter assingments by 
#' dividing by scale
#' @param rreturn 1 (default) to return an R (numeric) vector, 0 to return a torch tensor 
#' @param trnspose 1 to transpose the matrix before returning, 0 to not. 
#'
#' @return a weight matrix in tensor format
#'
wtlast = function(tnsr, lasso=0, lscale=5, scale=1, rreturn=1, trnspose=0 ) { 
  weight_r = as.matrix(tnsr)
  if (scale <= 0) { 
    weight_r = matrix( rep(0 ,dim(weight_r)[1] * dim(weight_r)[2]), nrow=dim(weight_r)[1], ncol=dim(weight_r)[2] ) 
  } else if (scale > 1) {
    weight_r = weight_r / scale 
  } 
  if (lasso != 0) {
    dmw = dim(weight_r) 
    col1 = c(rep(1,dmw[1]))  
    weight_r[,1] =  col1 * lscale  
    weight_r[,2] = -col1 * lscale 
    weight_r  
  }
  if (trnspose==1) { weight_r = t(weight_r) }
  if (rreturn ==1) { 
    weight = weight_r 
  } else {
    weight = torch_tensor(weight_r, dtype=torch_float(), requires_grad=TRUE)
  }
  
  return(weight)
}

################################################################################
#'
#' Construct the bias terms for going from model layer to layer to carry forward 
#' an offset to mimic a linear model 
#'
#' @param tnsr an input tensor which is to be modified to mimic the linear term of a 
#' generalized linear model, e.g a Cox or logistic regression model 
#' @param lasso 1 if the first column is the linear estimate from a linear model, 
#' often a lasso model 
#' @param rreturn 1 (default) to return an R (numeric) vector, 0 to return a torch tensor 
#'
#' @return a weight matrix in tensor format

bsint = function(tnsr, lasso=0, rreturn=1) { 
  bias_r = as.numeric(tnsr)
  lng = length(bias_r)
  if ((lasso >= 1) | (lasso < 0)) {
    if (lng >= 2) {
      bias_r[(1:2)] = c(0,0) 
    }
  }
  if (lasso < 0) {                                                              ## for testing
    bias_r = matrix( rep(0 , lng), nrow=lng, ncol=1)                      
  } else if (lasso > 1) {
    if (lng > 2) {
      bias_r[(3:lng)] = bias_r[(3:lng)] / lasso 
    }
  } 
  if (rreturn == 1) { 
    bias = bias_r  
  } else {
    bias = torch_tensor(bias_r, dtype=torch_float(), requires_grad=TRUE)
  }
  return(bias)
}

################################################################################
################################################################################
#' Standardize a data set 
#'
#' @param datain The data matrix set to be standardized 
#' @param lasso 1 to not standardize the first column, 0 (default) to not 
#'
#' @return a standardized data matrix 
#'
dtstndrz = function(datain, lasso=0) {
  if (lasso >= 1) {
    lassocol =  datain[,1] 
  }
  myxs_1 = datain
  myxs_means = colMeans(myxs_1)
  myxs_2     = sweep(myxs_1, 2, myxs_means, FUN = "-") 
  myxs_sds  = sqrt( colSums(myxs_2^2) / (dim(myxs_1)[1]-1) ) 
  myxs_3 = myxs_2[,(myxs_sds > 1e-6)]
  myxs_sds = myxs_sds[(myxs_sds > 1e-6)]
  myxs_4    = sweep(myxs_3, 2, myxs_sds, FUN = "/") 
  if (lasso >= 1) {
    myxs_4[,1] = lassocol  
  }
  dataout = myxs_4 
  return(dataout)
}

################################################################################
################################################################################
#'
#' predicted values from an ann_tab_cv output object based upon the model and its
#' lasso model used for generating an offset
#'
#' @param lassomod a lasso model from a glmnetr() call used to genrte an offset
#' @param nnmodel a ann_tab_cv() output object for 
#' @param datain new data 
#' @param lasso 1 if an offset is to be added as column 1 for calculations (default), 
#' 0 to subset to the terms significant in a lasso model without adding the offset
#'
#' @return predictions from an nerual network model accounting from a lasso model 
#' 
prednn_tl = function (lassomod, nnmodel, datain, lasso=1) {
  lassopred = predict(lassomod, datain)
  lassobeta = predict(lassomod)
  family = lassomod$sample[1]
  datain_1 = datain[,(lassobeta[[1]] != 0)[2:length(lassobeta[[1]])]]           ## pick up non zero features, remove intercept from this list 
  if (family == "cox") {
    datain_1 = datain[,(lassobeta[[1]] != 0)]                                   ## pick up non zero feature betas
  } else {
    datain_1 = datain[,(lassobeta[[1]] != 0)[2:length(lassobeta[[1]])]]         ## pick up non zero features beta, remove intercept from this list
  }
  if (lasso == 1) {
    datain_2 = cbind(lasso=lassopred,datain_1)                                  ## add lasso prediction as first column 
  } else {
    datain_2 = datain_1                                                         ## add lasso prediction as first column 
  }
  preds = nnmodel$model(datain_2)   ## on probability scale  
  preds = as.numeric(preds)
  return(preds)
}

################################################################################
#### torch nn_sum nn_cumsum nn_mul nn_add 
################################################################################
# 
## myxs=simdata$xs ;  myy=simdata$yt ;  myevent=simdata$event ;  start=NULL ; family="cox" ;  fold_n=5 ;  epochs=200 ;  eppr=40 ;  wd=0 ;  lenz1=8 ;  lenz2=5 ;  mylr=0.01 ;  lasso = 0 ;
#' 
#' Fit an Artificial Neural Network model on "tabular" provided as a matrix, optionally 
#' allowing for an offset term 
#' 
#' @description Fit an Artificial Neural Network model for analysis of "tabular" data.  The 
#' model has two hidden layers where the number of terms in each layer is configurable 
#' by the user.  The activation function can also be switched between relu() (default) 
#' gelu() or sigmoid().  Optionally an offset term may be included. Model "family" may 
#' be "cox" to fit a generalization of the Cox proportional hazards model, "binomial" 
#' to fit a generalization of the logistic regression model and "gaussian" to fit a 
#' generalization of linear regression model for a quantitative response.  See the 
#' corresponding vignette for examples.     
#' 
#' @param myxs     predictor input - an n by p matrix, where n (rows) is sample size, and p (columns) 
#' the number of predictors.  Must be in matrix form for complete data, no NA's, no Inf's, etc.,
#' and not a data frame. 
#' @param mystart  an optional vector of start times in case of a Cox model. Class numeric of length same as number of patients (n)
#' @param myy      dependent variable as a vector: time, or stop time for Cox model, Y_ 0 or 1 for binomial (logistic), numeric for gaussian. 
#' Must be a vector of length same as number of sample size. 
#' @param myevent  event indicator, 1 for event, 0 for census, Cox model only.
#' Must be a numeric vector of length same as sample size.
#' @param myoffset   an offset term to be used when fitting the ANN.  Not yet implemented 
#' in its pure form.  Functionally an offset can be included in the first column
#' of the predictor or feature matrix myxs and indicated as such using the lasso option.   
#' @param family   model family, "cox", "binomial" or "gaussian" (default) 
#' @param fold_n   number of folds for each level of cross validation
#' @param epochs   number of epochs to run when tuning on number of epochs for fitting 
#' final model number of epochs informed by cross validation 
#' @param eppr     for EPoch PRint.  print summary info every eppr epochs. 0 will 
#' print first and last epochs, 0 for first and last epoch, -1 for minimal and -2 for none.  
#' @param lenz1    length of the first hidden layer in the neural network, default 16 
#' @param lenz2    length of the second hidden layer in the neural network, default 16 
#' @param actv     for ACTiVation function.  Activation function between layers, 
#' 1 for relu, 2 for gelu, 3 for sigmoid. 
#' @param drpot    fraction of weights to randomly zero out.  NOT YET implemented. 
#' @param mylr     learning rate for the optimization step in the neural network model fit
#' @param wd       a possible weight decay for the model fit, default 0 for not considered  
#' @param l1       a possible L1 penalty weight for the model fit, default 0 for not considered  
#' @param lasso    1 to indicate the first column of the input matrix is an offset 
#' term, often derived from a lasso model, else 0 (default)  
#' @param lscale Scale used to allow ReLU to exend +/- lscale before capping the 
#' inputted linear estimated 
#' @param scale Scale used to transform the inital random paramter assingments by 
#' dividing by scale
#' @param resetlw  1 as default to re-adjust weights to account for the offset every 
#' epoch.  This is only used in case lasso is set to 1.  
#' @param minloss default of 1 for minimizing loss, else maximizing agreement (concordance 
#' for Cox and Binomial, R-square for Gaussian), as function of epochs by cross validaition  
#' @param gotoend fit to the end of epochs.  Good for plotting and exploration 
#' @param seed an optional a numerical/integer vector of length 2, for R and torch 
#' random generators, default NULL to generate these.  Integers should be positive 
#' and not more than 2147483647.
#' @param foldid  a vector of integers to associate each record to a fold.  Should 
#' be integers from 1 and fold_n.
#' 
#' @return an artificial neural network model fit 
#' 
#### note, self is a "direct to" and not a function.  It can be found in ~tabnet/R/tab-network.R but not in ~tabnet/NAMESPACE 
#' @importFrom torch nn_relu nn_gelu nn_sigmoid nn_identity nn_softmax nn_module 
#' @importFrom torch nn_linear nnf_linear nn_parameter nn_sequential nn_dropout torch_tensor torch_float optim_adam 
#' @importFrom torch nnf_binary_cross_entropy nn_bce_loss nn_mse_loss nnf_mse_loss nn_cross_entropy_loss 
#' @importFrom survival Surv coxph coxph.control concordance 
#' 
#' @seealso
#'   \code{\link{ann_tab_cv_best}} , \code{\link{predict_ann_tab}}, \code{\link{nested.glmnetr}} 
#' 
#' @author Walter Kremers (kremers.walter@mayo.edu)
#' 
#' @export
#'
ann_tab_cv = function(myxs, mystart=NULL, myy, myevent=NULL, myoffset=NULL, family="binomial", fold_n=5, 
                      epochs=200, eppr=40, lenz1=16, lenz2=8, actv=1, drpot=0, mylr=0.005, wd=0, l1=0, 
                      lasso=0, lscale=5, scale=1, resetlw=1, minloss=1, gotoend=0, seed=NULL, foldid=NULL) { 
stratified = 1 

if (!(family %in% c("cox", "binomial", "gaussian"))) { 
  cat( "\n  ***** ann_cd_lin() is only set up for famlies cox, binomial or gaussian *****\n")
  cat( paste("  ***** not ", family, " *****\n")) 
  family = "gaussian" 
#  stop() 
}  
  
if (drpot > 0) {
#  drpot = 0 
  cat(paste("  drpot (for drop out) is not yet implemented and drpot is set to 0 \n"))
}  
  
if (!is.null(mystart)) { 
  start = NULL 
  cat(paste("  start time for Cox model not yet implemented\n   ann_tab_cv will stop "))
}

if (fold_n == 1) { 
  fold_n = 3
  cat(paste(" --- fold_n cannot be 1 and is so set to 3 ---"))
}
  
if (l1 > 0) {
  cat(paste(" l1 penalty not yet implemented... \n"))
}
  
#  myxs=xs_z0 ;  mystart=NULL ;  myy=y_ ;  myevent=NULL ;  myoffset=NULL ;  family="binomial" ;  fold_n=5 ;  
#  epochs=200 ;  eppr=40 ;  lenz1=32 ;  lenz2=8 ;  actv=1 ;  mylr=0.005 ;  drpot=0 ;  wd=0 ;  lasso=0 ;  resetlw=0 ;  
  
nobs  = dim(myxs)[1]
nfeat = dim(myxs)[2]
  
if (family == "cox") {
#  myy     = round(myy , digits = 8)
#  cat(paste("dim(myy)=",length(myy)),"\n") 
#  cat(paste("dim(myxs)= (",dim(myxs)[1], ",", dim(myxs)[2], ")\n")) 
  myorder = order(myy , decreasing = TRUE) 
#  cat(paste("dim(myorder)=",length(myorder)),"\n")
  myxs = as.matrix(myxs[myorder,])
  myy  = myy [myorder]
  myevent = myevent[myorder]
  if (!is.null(start)) {
    mystart = mystart[myorder]
  }
  if (!is.null(myoffset)) {
    myoffset = myoffset[myorder]
  }
} 
  
model = NULL 
mymodel = NULL
modelinit = NULL 

# n_obs = nrow(myxs)[1] 
if (family == "binomial") {
  nclass = length(levels(myy) )
  nclass = length(names(table(myy))) 
} else if (family == "cox") {
  nclass = length(table(myevent))
} else {
  nclass = 1 ; 
}
c(nobs, nfeat, nclass)

if (nclass == 2) { nclass = 1 }
#if (family == "gaussian") {nclass = 1 }
#if (family == "cox") {nclass = 1 }

cvloss     = matrix(rep(10,(fold_n*epochs)), nrow=fold_n, ncol=epochs) 
cvaccuracy = matrix(rep( 0,(fold_n*epochs)), nrow=fold_n, ncol=epochs) 
cvagree    = matrix(rep( 0,(fold_n*epochs)), nrow=fold_n, ncol=epochs) 
losstrain    = matrix(rep(10,(fold_n*epochs)), nrow=fold_n, ncol=epochs) 


if (is.null(foldid)) {
  if (is.null(seed)) {
    seedr = round(runif(1)*1e9) 
    seedt = round(runif(1)*1e9)
  } else {
    seedr = seed[1]
    seedt = seed[2]
  }
  ##----------------------------------------------------
  if (is.null(seedr)) { seedr = round(runif(1)*1e9) 
  } else if (is.na(seedr))  { seedr = round(runif(1)*1e9) } 
  ##----------------------------------------------------
  if (is.null(seedt)) { seedt = round(runif(1)*1e9) 
  } else if (is.na(seedt))  { seedt = round(runif(1)*1e9) } 
  ##----------------------------------------------------
  set.seed(seedr)                                                              ##<<=========
  foldid = get.foldid(myy, myevent, family, fold_n, stratified) 
#  print(table(foldid))
  ##----------------------------------------------------
} else {
  if (is.null(seed)) {
    seedr = NA 
    seedt = round(runif(1)*1e9)
  } else {
    seedr = seed[1]
    seedt = seed[2]
  }
  if (is.null(seedt)) { 
    seedr = 0 
    seedt = round(runif(1)*1e9) 
  } else if (is.na(seedt))  { 
    seedr = 0 
    seedt = round(runif(1)*1e9) 
  } 
}
seed = c(seedr=seedr, seedt=seedt) 
torch_manual_seed( seedt ) 

if        (actv == 1) { act_fn = nn_relu()      ; actvc = "relu" 
} else if (actv == 2) { act_fn = nn_gelu()      ; actvc = "gelu" 
} else if (actv == 3) { act_fn = nn_sigmoid()   ; actvc = "sigmoid" 
} else if (actv == 4) { act_fn = nn_sigmoid()$mul(2)$add(-1) ; actvc = "arctan" 
} else if (actv ==-1) { act_fn = nn_identity()  ; actvc = "identity" 
} else                { act_fn = nn_relu()      ; actvc = "relu" 
}

if        (family == "binomial"  ) { lastact_fn = nn_sigmoid() 
} else if (family == "multiclass") { lastact_fn = nn_softmax()
} else if (family == "gaussian"  ) { lastact_fn = nn_identity()
} else if (family == "cox"       ) { lastact_fn = nn_identity() 
} else { lastact_fn = nn_identity() } 

#-------------------------------------------------------------------------------
#--  this model piece (module) will allow input of adjusted weights  -----------

my_nn_linear0 <- nn_module(
#  classname = "my_nn_linear0",
  initialize = function(weight, bias) {
#    self$'weight' <- nn_parameter(weight)
#    self$'bias' <- nn_parameter(bias)
    self$weight <- nn_parameter(weight)
    self$bias <- nn_parameter(bias)
  },
  forward = function(x) {
    nnf_linear(x, self$weight, self$bias)
  }
)

#-------------------------------------------------------------------------------
#--  construct base model to get initial random weights and biases  ------------

  model0 = nn_sequential(
    nn_linear(nfeat, lenz1),
    act_fn , 
    nn_dropout(p = drpot),
    nn_linear(lenz1, lenz2),
    act_fn , 
    nn_dropout(p = drpot),
    nn_linear(lenz2,nclass),
    lastact_fn 
  )

##### extract weight and bias matrices and store in R data so are not updated automatically ####
# lasso = 0 
trnspose = 0 
weight0_r0 = wtzero(   model0$parameters$`0.weight` , lasso, lscale, scale, trnspose=trnspose)
weight3_r0 = wtmiddle( model0$parameters$`3.weight` , lasso, trnspose=trnspose)
weight6_r0 = wtlast(   model0$parameters$`6.weight` , lasso, lscale, scale, trnspose=trnspose)
bias0_r0 = bsint( model0$parameters$`0.bias` , lasso  )
bias3_r0 = bsint( model0$parameters$`3.bias` , lasso  )
bias6_r0 = bsint( model0$parameters$`6.bias` , lasso=0)

#-----------------------------------------------------

weight0 = torch_tensor(weight0_r0, dtype=torch_float(), requires_grad=TRUE)
weight3 = torch_tensor(weight3_r0, dtype=torch_float(), requires_grad=TRUE)
weight6 = torch_tensor(weight6_r0, dtype=torch_float(), requires_grad=TRUE)
bias0 = torch_tensor(bias0_r0, dtype=torch_float(), requires_grad=TRUE)
bias3 = torch_tensor(bias3_r0, dtype=torch_float(), requires_grad=TRUE)
bias6 = torch_tensor(bias6_r0, dtype=torch_float(), requires_grad=TRUE) 

#-------------------------------------------------------------------------------
#-----  adjust weights for lasso input and put back into model  ----------------

  model = nn_sequential(
    my_nn_linear0(weight0, bias0) ,
    act_fn , 
    nn_dropout(p = drpot) ,
    my_nn_linear0(weight3, bias3) ,
    act_fn , 
    nn_dropout(p = drpot),
    my_nn_linear0(weight6, bias6) ,
    lastact_fn 
  )

#-------------------------------------------------------------------------------
# model$parameters$'3.weight'
# model = resetlw_(lasso=0)
# model$parameters$'3.weight'
# model = resetlw_(lasso=1)
# model$parameters$'3.weight'
# weight3

# xs.t = torch_tensor(myxs, dtype = torch_float())
# y_.t = torch_tensor(myy , dtype = torch_float())
# y_pred = model(xs.t)

#   mymodel$parameters
#   myoptim$zero_grad()

#modelstore <<- mymodel 
#myclone = clone(mymodel)

#--  assign cost and optimizer  ------------------------------------------------
#--  no nn_ function for cox partial likelihood so this is written with in loops-
if (nclass > 2) { 
  loss_fn  = nn_cross_entropy_loss()  
  #    loss_fnf = nnf_cross_entropy()
} else if (family =="binomial") {
  loss_fn  = nn_bce_loss()  
} else if (family =="gaussian") {
  loss_fn  = nn_mse_loss()  
} else if (family =="cox") {
#  loss_fn  = nn_surv1_loss()
} 

## could use for binomial nn_binary_cross_entropy_with_logits loss function and 
## eliminate need for an activation function at last step 

fold = 1 
#####  do the folds  ###########################################################
for (fold in c(1:fold_n)) {
  x_train = as.matrix (myxs[foldid!=fold,])
  x_test  = as.matrix (myxs[foldid==fold,])
  y_train = as.numeric(myy [foldid!=fold ])
  y_test  = as.numeric(myy [foldid==fold ])
  
  if (family == "cox") {
    if (!is.null(mystart)) {
      s_train = mystart[foldid!=fold,]
      s_test  = mystart[foldid==fold,]
    }
    e_train = myevent[foldid!=fold]
    e_test  = myevent[foldid==fold]
  }
  
  if (!is.null(myoffset)) {
    o_train = myoffset[foldid!=fold ]
    o_test  = myoffset[foldid==fold ]
  }

  x_train_t     = torch_tensor(x_train, dtype = torch_float())
  x_test_t      = torch_tensor(x_test , dtype = torch_float())
  y_train_t     = torch_tensor(y_train, dtype = torch_float())
  y_test_t      = torch_tensor(y_test , dtype = torch_float())
  if (family == "cox") {
    e_train_t   = torch_tensor(e_train, dtype=torch_float())
    e_test_t    = torch_tensor(e_test , dtype=torch_float())
    if (!is.null(mystart)) {
      s_train_t = torch_tensor(s_train, dtype=torch_float())
      s_test_t  = torch_tensor(s_test , dtype=torch_float())      
    }
  }
  
  if (!is.null(myoffset)) {
    o_train_t = torch_tensor(o_train, dtype=torch_float())
    o_test_t  = torch_tensor(o_test , dtype=torch_float())
  }
  
  #####  reset model to same starting model for CV  ############################
  weight0 = torch_tensor(weight0_r0, dtype=torch_float(), requires_grad=TRUE)
  weight3 = torch_tensor(weight3_r0, dtype=torch_float(), requires_grad=TRUE)
  weight6 = torch_tensor(weight6_r0, dtype=torch_float(), requires_grad=TRUE)
  bias0 = torch_tensor(bias0_r0, dtype=torch_float(), requires_grad=TRUE)
  bias3 = torch_tensor(bias3_r0, dtype=torch_float(), requires_grad=TRUE)
  bias6 = torch_tensor(bias6_r0, dtype=torch_float(), requires_grad=TRUE) 
  
    model = nn_sequential(
      my_nn_linear0(weight0, bias0) ,
      act_fn ,  
      nn_dropout(p = drpot) ,
      my_nn_linear0(weight3, bias3) ,
      act_fn ,  
      nn_dropout(p = drpot) ,
      my_nn_linear0(weight6, bias6) ,
      lastact_fn 
    ) 
  
  ##---  reset optimizer with every fold  --------------------------------------
  myoptim = optim_adam(model$parameters, lr = mylr, weight_decay = wd)
  
  ##---  check model on TEST data before fitting (epochs) begins  --------------
#    if (!is.null(myoffset)) {  
#      y_pred_t = model(x_train_t)$squeeze(2)$add(o_train_t)
#      y_predtest_t  = model(x_test_t)$squeeze(2)$add(o_test_t)
#    }  

    if (family== "multiclass") {
      corrects = y_pred_t$argmax(dim=2) # +1
      accuracy0 = as.numeric( sum(corrects == y_test_t) ) / length(y_test) 
    } else if (family == "binomial") {
      y_pred_t = model(x_train_t)$squeeze(2)
      y_predtest_t  = model(x_test_t)$squeeze(2) 
      y_pred_r = as.numeric(y_pred_t)
      y_predtest_r  = as.numeric(y_predtest_t)
      losstrain0 = as.numeric(nnf_binary_cross_entropy(y_pred_t,y_train_t))
      corrects = as.numeric( (y_predtest_r >= 0.5)*(y_test > 0) + (y_predtest_r < 0.5)*(y_test < 1) ) 
      accuracy0 = sum(corrects) / length(corrects) 
      agree0    = concordance(y_test ~ y_predtest_r)[[1]]   
    } else if (family == "gaussian") { 
      y_pred_t = model(x_train_t)$squeeze(2)
      y_predtest_t  = model(x_test_t)$squeeze(2)
      y_pred_r = as.numeric(y_pred_t)
      y_predtest_r  = as.numeric(y_predtest_t)
      losstrain0 = as.numeric(nnf_mse_loss(y_predtest_t,y_test_t)) ; ## mean((y_predtest_r - y_test)^2) ## same as ... 
      agree0     = cor(y_test, y_predtest_r)  # ^2   
    } else if (family == "cox") { 
      if (!is.null(myoffset)) {
        y_pred_t = model(x_train_t)$squeeze(2)$add(o_test_t)      
        y_predtest_t = model(x_test_t)$squeeze(2)$add(o_test_t)      
      } else {
        y_pred_t = model(x_train_t)$squeeze(2)
        y_predtest_t = model(x_test_t)$squeeze(2)
      }
      y_pred_r = as.numeric(y_pred_t) 
      fit0 = coxph( Surv(y_train,e_train) ~ y_pred_r, init=c(1), control=coxph.control(iter.max=0)) 
      losstrain0 = - fit0$loglik[1] / sum(e_train)
      agreetrain0     =  fit0$concordance[6]
      
      y_predtest_r = as.numeric(y_predtest_t) 
      fit1 = coxph( Surv(y_test,e_test) ~ y_predtest_r, init=c(1), control=coxph.control(iter.max=0)) 
      losstest0 = - fit1$loglik[1] / sum(e_test)
      agreetest0     =  fit1$concordance[6]
    } 

  if (eppr >= -1) {    
    if (family == "binomial") {
      cat(" Fold: ", fold, " Epoch: 0 Train Loss: ", losstrain0, "Train Concordance:", agree0,"\n")
    } else if (family == "gaussian") {
      cat(" Fold: ", fold, " Epoch: 0 Train Loss: ", losstrain0, "Train R-square:", agree0,"\n")
    } else if (family == "multiclass") {
      cat(" Fold: ", fold, " Epoch: 0 Train Loss: ", losstrain0, "Train Accuracy:", accuracy0,"\n")   
    } else if (family == "cox") {
      cat(" Fold: ", fold, " Epoch: 0 Train Loss: ", losstrain0, "Train Concordance:", agreetrain0,"\n")   
    }
  }

  l1_t = torch_tensor(l1,dtype=torch_float())
  nobs_train_t = torch_tensor(length(y_train_t), dtype=torch_float())
      
  i_ = 1 
  ## train the net #############################################################
  for (i_ in 1:epochs) {
    if (family %in% c("binomial", "gaussian")) {
               
      if (family == "binomial") {
        y_pred <- model(x_train_t)  
        if (l1 > 0) {
          loss = loss_fn(y_pred, y_train_t)                                       
        } else {
          loss = loss_fn(y_pred, y_train_t) + weight0$abs()$sum()$mul(l1_t) + weight3$abs()$sum()$mul(l1_t) + weight6$abs()$sum()$mul(l1_t) 
        }
      } else if (family == "gaussian") {
        y_pred <- model(x_train_t)$squeeze()                     
        if (l1 > 0) {
#          if (i_ ==1) { print( "HERE 1") }
          loss = y_pred$sub(y_train_t)$pow(2)$sum()$div(nobs_train_t) + weight0$abs()$sum()$mul(l1_t) + weight3$abs()$sum()$mul(l1_t) + weight6$abs()$sum()$mul(l1_t) 
        } else {
          loss = loss_fn(y_pred, y_train_t)                                     ##<<---------------------------
        }
      } 
      myoptim$zero_grad()                                                       
      loss$backward() 
      myoptim$step()   
    }
    
    if (family == "cox") {
      if (!is.null(myoffset)) {
        y_pred <- model(x_train_t)$squeeze(2)$add(o_train_t) 
      } else {
        y_pred <- model(x_train_t)$squeeze(2)
      }
      myoptim$zero_grad()       ## myoptim_1 ?? 
#       loss = - e_train_t$mul(y_pred)$sum()$div(e_train_t$sum()) + y_pred$exp()$cumsum(1)$log()$mul(e_train_t)$sum()$div(e_train_t$sum())         
##     myoptim$zero_grad()     ## or myoptim_2 placement ??  These provide the same model fit. 
      if (l1 >  0) {
        loss = - e_train_t$mul(y_pred)$sum()$div(e_train_t$sum()) + y_pred$exp()$cumsum(1)$log()$mul(e_train_t)$sum()$div(e_train_t$sum())  
               + weight0$abs()$sum()$mul(l1_t) + weight3$abs()$sum()$mul(l1_t) + weight6$abs()$sum()$mul(l1_t) 
      } else {
        loss = - e_train_t$mul(y_pred)$sum()$div(e_train_t$sum()) + y_pred$exp()$cumsum(1)$log()$mul(e_train_t)$sum()$div(e_train_t$sum())  
      }
      loss$backward() 
      myoptim$step()                                                            ## why do I get a null here ?? 
    }
    
    model$parameters$`0.weight`[1:5,1:5]
    if ((resetlw ==1) & ((lasso < 0) | (lasso >=1))) { 
      weight0_r = wtzero(   model$parameters$`0.weight` , lasso, lscale, scale)
      weight3_r = wtmiddle( model$parameters$`3.weight` , lasso)
      weight6_r = wtlast(   model$parameters$`6.weight` , lasso, lscale, scale)
      bias0_r = bsint( model$parameters$`0.bias` , lasso )
      bias3_r = bsint( model$parameters$`3.bias` , lasso )
      bias6_r = bsint( model$parameters$`6.bias` , lasso ) 
      
      weight0 = torch_tensor(weight0_r, dtype=torch_float(), requires_grad=TRUE)
      weight3 = torch_tensor(weight3_r, dtype=torch_float(), requires_grad=TRUE)
      weight6 = torch_tensor(weight6_r, dtype=torch_float(), requires_grad=TRUE)
      bias0 = torch_tensor(bias0_r, dtype=torch_float(), requires_grad=TRUE)
      bias3 = torch_tensor(bias3_r, dtype=torch_float(), requires_grad=TRUE)
      bias6 = torch_tensor(bias6_r, dtype=torch_float(), requires_grad=TRUE) 
      
      model = nn_sequential(
        my_nn_linear0(weight0, bias0) ,
        act_fn ,  
        nn_dropout(p = drpot) ,
        my_nn_linear0(weight3, bias3) ,
        act_fn ,  
        nn_dropout(p = drpot) ,
        my_nn_linear0(weight6, bias6) ,
        lastact_fn 
      )    
    }
    model$parameters$`0.weight`[1:5,1:5]
    
#   check model on TRAIN data 
    ##-- figure out loss$item as this doesnt always agree with loss_fn ---------
    if (family == "cox") {
      losstrain[fold, i_] = loss$item() ## / sum(e_train)                          ## already divided by sum(e_train)
    } else {
      losstrain[fold, i_] = loss$item() 
    }
    
#   check model on TEST data 
    if (family == "multiclass") {
      y_predtest_t = model(x_test_t)$squeeze(2)
      corrects = y_predtest_t$argmax(dim=2) # +1
      cvaccuracy[fold,i_] = as.numeric( sum(corrects == y_test_t) ) / length(y_test) 
      
    } else if (family == "binomial") {
      if (!is.null(myoffset)) { 
        y_predtest_t = model(x_test_t)$squeeze(2)$add(o_test_t)
      } else { 
        y_predtest_t = model(x_test_t)$squeeze(2) 
      } 
      y_predtest_r = as.numeric( y_predtest_t )
      cvloss[fold, i_] = as.numeric(nnf_binary_cross_entropy(y_predtest_t,y_test))
      corrects = as.numeric( (y_predtest_r >= 0.5)*(y_test > 0) + (y_predtest_r < 0.5)*(y_test < 1) ) 
      # corrects = (y_predtest <= 0)*(y_test > 0) + (y_predtest > 0)*(y_test < 1)
      cvaccuracy[fold,i_] = sum(corrects) / length(corrects) 
      cvagree   [fold,i_] = concordance(y_test ~ y_predtest_r)[[1]]   
      
    } else if (family == "gaussian") { 
        y_predtest_t = model(x_test_t)$squeeze(2)
        y_predtest_r = as.numeric( y_predtest_t )
        cvloss [fold, i_] = as.numeric(nnf_mse_loss(y_predtest_t,y_test_t)) ; ## mean((y_predtest_r - y_test)^2) ## same as ... 
        cvagree[fold,i_] = cor(y_test, y_predtest_r)  # ^2   
        
    } else if (family == "cox") { 
      if (!is.null(myoffset)){
        y_predtest_t = model(x_test_t)$squeeze(2)$add(o_test_t)
      } else {
        y_predtest_t = model(x_test_t)$squeeze(2)
      }
      y_predtest_r = as.numeric(y_predtest_t)
      if (sum(is.na(y_predtest_r))>0) { y_predtest_r = rep(1,length(y_predtest_r)) }
      fit0 = coxph( Surv(y_test,e_test) ~ y_predtest_r, init=c(1), control=coxph.control(iter.max=0)) 
      cvloss [fold, i_] = - fit0$loglik[1] / sum(e_test) 
      cvagree[fold, i_] =  fit0$concordance[6]
    } 
    
    docat = 0 
    if ((eppr >= 0) & (i_ == 1)) { docat = 1 }
    if ((eppr >= 0) & (i_ == epochs)) { docat = 1 }
    if (eppr >  0) { if (i_ %% eppr == 0) { docat = 1 } }
    if (docat == 1) {
      if (family == "binomial") {
        cat(" Epoch:", i_,"Train Loss: ", loss$item(), "Test Loss: ", cvloss[fold, i_], "Test Concordance:", cvagree[fold,i_],"\n")
      } else if (family == "gaussian") {
        cat(" Epoch:", i_,"Train Loss: ", loss$item(), "Test Loss: ", cvloss[fold, i_], "Test R-square:", cvagree[fold,i_]^2,"\n")
      } else if (family == "multiclass") {
        cat(" Epoch:", i_,"Train Loss: ", loss$item(), "Test Loss: ", cvloss[fold, i_], "Test Accuracy:", cvaccuracy[fold,i_],"\n")        
      } else if (family == "cox") {
        cat(" Epoch:", i_,"Train Loss: ", loss$item(), "Test Loss: ", cvloss[fold, i_], "Test Concordance:", cvagree[fold,i_],"\n") 
      }
    } 
  } ##  end of epoch  ##########################################################

#  as.matrix(weight0_r[1:5,1:5])  
#  as.matrix(model$parameters$'0.weight'[1:5,1:5])  

  if (eppr >= 0) {
    cat(" fold ", fold, " (of " , fold_n,  ") minimizing validation loss at epoch ", which.min(cvloss[fold,]), "\n\n" )  
  }
}  ##  end of fold  ############################################################

cvlossmn     = colMeans(cvloss)
cvaccuracymn = colMeans(cvaccuracy)
cvagreemn    = colMeans(cvagree   )
which_loss   = which.min(cvlossmn)
which_agree  = which.max(cvagreemn^2)
if (length(which_agree)==0) { which_agree = which_loss }
which_accu   = which.max(cvaccuracymn)
cv_loss      = cvlossmn[which_loss] 
cv_agree     = cvagreemn[which_loss] 
cv_accuracy  = cvaccuracymn[which_loss] 

if (eppr >= 0) {
  if (family %in% c("cox","binomial")) {
    cat( " epoch minimizing CV loss = ", which_loss , " CV loss = ", cv_loss,  
         " CV concordance = ", cv_agree, "\n" )
  } else if (family %in% c("gaussian")) {
    cat( " epoch minimizing CV loss = ", which_loss , " CV loss = ", cv_loss,  
         " CV R-square = ", cv_agree^2, "\n" )
  } else if (nclass > 2) {
    cat( " epoch minimizing CV loss = ", which_loss , " CV loss = ", cv_loss,  
         " CV accuracy = ", cv_accuracy, "\n" )
  }
}

#-------------------------------------------------------------------------------
#--  set up data to train on the whole data set --------------------------------

if (minloss == 1) { whichone = which_loss 
} else { whichone = which_agree  }

# set up arrays to store model checks 
lossf    = c(rep(10,epochs)) 
lossf2   = c(rep(10,epochs)) 
agreef    = c(rep(0,epochs)) 
accuracyf = c(rep(0,epochs)) 

#-- set up data for training and evaluation ------------------------------------
x_train = as.matrix (myxs)
x_test  = as.matrix (myxs)
y_train = as.numeric(myy )
y_test  = as.numeric(myy )

if (family == "cox") {
  if (!is.null(mystart)) {
    s_train = mystart
    s_test  = mystart
  }
  e_train = myevent
  e_test  = myevent
}

if (!is.null(myoffset)) {
  o_train = myoffset 
  o_test  = myoffset 
}

x_train_t     = torch_tensor(x_train, dtype = torch_float())
x_test_t      = torch_tensor(x_test , dtype = torch_float())
y_train_t     = torch_tensor(y_train, dtype = torch_float())
y_test_t      = torch_tensor(y_test , dtype = torch_float())
if (family == "cox") {
  e_train_t   = torch_tensor(e_train, dtype=torch_float())
  e_test_t    = torch_tensor(e_test , dtype=torch_float())
  if (!is.null(mystart)) {
    s_train_t = torch_tensor(s_train, dtype=torch_float())
    s_test_t  = torch_tensor(s_test , dtype=torch_float())      
  }
}

if (!is.null(myoffset)) {
  o_train_t = torch_tensor(o_train, dtype=torch_float())
  o_test_t  = torch_tensor(o_test , dtype=torch_float())
}

#-- reset the starting NN model to the common starting model  ------------------

weight0 = torch_tensor(weight0_r0, dtype=torch_float(), requires_grad=TRUE)
weight3 = torch_tensor(weight3_r0, dtype=torch_float(), requires_grad=TRUE)
weight6 = torch_tensor(weight6_r0, dtype=torch_float(), requires_grad=TRUE)
bias0 = torch_tensor(bias0_r0, dtype=torch_float(), requires_grad=TRUE)
bias3 = torch_tensor(bias3_r0, dtype=torch_float(), requires_grad=TRUE)
bias6 = torch_tensor(bias6_r0, dtype=torch_float(), requires_grad=TRUE) 

  model = nn_sequential(
    my_nn_linear0(weight0, bias0) ,
    act_fn ,  
    nn_dropout(p = drpot) ,
    my_nn_linear0(weight3, bias3) ,
    act_fn ,  
    nn_dropout(p = drpot) ,
    my_nn_linear0(weight6, bias6) ,
    lastact_fn 
  ) 

## assign cost and optimizer ##### MUST follow the model definition ############

if (nclass>2) {
  loss_fn = nn_cross_entropy_loss()  
} else if (family == "binomial") {
  loss_fn = nn_bce_loss() 
} else if (family == "gaussian") {
  loss_fn = nn_mse_loss()  
}

myoptim = optim_adam(model$parameters, lr = mylr, weight_decay = wd)

## myoptim = optim_adam(model$parameters, lr = mylr)  

#--  check model on TRAIN/TEST (same) data  ------------------------------------  

  if (family== "multiclass") {
    y_pred_t = model(x_train_t)$squeeze(2)
    y_pred_r = as.numeric(y_pred_t)
    corrects = y_pred_t$argmax(dim=2) # +1
    accuracy0 = as.numeric( sum(corrects == y_train_t) ) / length(y_train)
    agree0 = accuracy0 
    
  } else if (family == "binomial") {
    if (!is.null(myoffset)) {  
      y_pred_t = model(x_train_t)$squeeze(2)$add(o_train_t)
    } else {
      y_pred_t = model(x_train_t)$squeeze(2)
    }
    y_pred_r = as.numeric(y_pred_t)
    losstrain0 = as.numeric(nnf_binary_cross_entropy(y_pred_t,y_train_t))
    corrects = as.numeric( (y_pred_r >= 0.5)*(y_train > 0) + (y_pred_r < 0.5)*(y_train < 1) ) 
    accuracy0 = sum(corrects) / length(corrects) 
    agree0    = concordance(y_train ~ y_pred_r)[[1]]   
    
  } else if (family == "gaussian") { 
    if (!is.null(myoffset)) {  
      y_pred_t = model(x_train_t)$squeeze(2)$add(o_train_t)
    } else {
      y_pred_t = model(x_train_t)$squeeze(2)
    }
    y_pred_r = as.numeric(y_pred_t)
    losstrain0 = as.numeric(nnf_mse_loss(y_pred_t,y_train_t)) ; ## mean((y_predtest_r - y_test)^2) ## same as ... 
    agree0     = cor(y_train, y_pred_r)  # ^2   
  } else if (family == "cox") { 
    if (!is.null(myoffset)) {
      y_pred_t = model(x_train_t)$squeeze(2)$add(o_train_t)      
    } else {
      y_pred_t = model(x_train_t)$squeeze(2)
    }
    y_pred_r = as.numeric(y_pred_t) 
    if (sum(is.na(y_pred_r))>0) { y_pred_r = rep(1,length(y_pred_r)) }
    fit0 = coxph( Surv(y_train,e_train) ~ y_pred_r, init=c(1), control=coxph.control(iter.max=0)) 
    losstrain0 = - fit0$loglik[1] / sum(e_train)
    agree0     =  fit0$concordance[6]
  } 

if (eppr >= -2) {  
  if (family == "binomial") {
    cat(" Epoch:   0 Full data Loss: ", losstrain0, "Train Concordance:", agree0,"\n")
  } else if (family == "gaussian") {
    cat(" Epoch:   0 Full data Loss: ", losstrain0, "Train R-square:", agree0 ^2,"\n")
  } else if (family == "multiclass") {
    cat(" Epoch:   0 Full data Loss: ", losstrain0, "Train Accuracy:", accuracy0,"\n")   
  } else if (family == "cox") {
    cat(" Epoch:   0 Full data Loss: ", losstrain0, "Train Concordance:", agree0,"\n")   
  }
}

## train the net on the whole data form the hyperparameter above ###############

nobs_train_t = torch_tensor(length(y_train_t), dtype=torch_float())

if (minloss == 1) { imax = which_loss
} else { imax = which_agree }
if (gotoend == 1) { imax = epochs  }

if (eppr >= -2) { cat(paste("  imax=", imax, "  minloss=", minloss, "  which_loss=", which_loss, "  which_agree=", which_agree, "  gotoend=", gotoend, "  l1=", l1, "\n")) }

i_ = 1 
for (i_ in 1:imax) {
  
  if (family %in% c("binomial", "gaussian")) {
    if (family == "binomial") {
      y_pred <- model(x_train_t)$squeeze()
      if (l1 > 0) {
#        if (i_ == 1) { print( "HERE 2") } 
        loss = loss_fn(y_pred, y_train_t) + weight0$abs()$sum()$mul(l1_t) + weight3$abs()$sum()$mul(l1_t) + weight6$abs()$sum()$mul(l1_t) 
      } else {
        loss = loss_fn(y_pred, y_train_t) 
      }
    } else if (family == "gaussian") {
      y_pred <- model(x_train_t)$squeeze()                     
      if (l1 > 0) {
#       loss = y_pred$sub(y_train_t)$pow(2)$sum()$div(nobs_train_t) #+weight0$abs()$sum()$mul(l1_t) + weight3$abs()$sum()$mul(l1_t) + weight6$abs()$sum()$mul(l1_t) 
        loss = y_pred$sub(y_train_t)$pow(2)$sum()$div(nobs_train_t) + weight0$abs()$sum()$mul(l1_t) + weight3$abs()$sum()$mul(l1_t) + weight6$abs()$sum()$mul(l1_t) 
      } else {
        loss = loss_fn(y_pred, y_train_t)                                       ##<<---------------------------
      }
      loss = loss_fn(y_pred, y_train_t)                                         ##<<---------------------------
    } 
    myoptim$zero_grad()                                                         ## MOD 230415 moved from before y_pred <- model() 
    loss$backward() 
    myoptim$step()   
  }
  
  if (family == "cox") {
    myoptim$zero_grad()  
    if (!is.null(myoffset)) {
      y_pred = model(x_train_t)$squeeze(2)$add(o_train_t)
#      loss = - e_train_t$mul(y_pred$add(o_train_t))$sum()$div(e_train_t$sum()) + y_pred$add(o_train_t)$exp()$cumsum(1)$log()$mul(e_train_t)$sum()$div(e_train_t$sum()) 
    } else {
      y_pred = model(x_train_t)$squeeze(2)
    }
    if (l1 >  0) {
      loss = - e_train_t$mul(y_pred)$sum()$div(e_train_t$sum()) + y_pred$exp()$cumsum(1)$log()$mul(e_train_t)$sum()$div(e_train_t$sum())         
             + weight0$abs()$sum()$mul(l1_t) + weight3$abs()$sum()$mul(l1_t) + weight6$abs()$sum()$mul(l1_t) 
    } else {
      loss = - e_train_t$mul(y_pred)$sum()$div(e_train_t$sum()) + y_pred$exp()$cumsum(1)$log()$mul(e_train_t)$sum()$div(e_train_t$sum())  
    }
    loss$backward() 
    myoptim$step()      
  }

  if ((resetlw ==1) & ((lasso < 0) | (lasso >=1))) { 
    weight0_r = wtzero(   model$parameters$`0.weight` , lasso, lscale, scale)
    weight3_r = wtmiddle( model$parameters$`3.weight` , lasso)
    weight6_r = wtlast(   model$parameters$`6.weight` , lasso, lscale, scale)
    bias0_r = bsint( model$parameters$`0.bias` , lasso )
    bias3_r = bsint( model$parameters$`3.bias` , lasso )
    bias6_r = bsint( model$parameters$`6.bias` , lasso ) 
    
    weight0 = torch_tensor(weight0_r, dtype=torch_float(), requires_grad=TRUE)
    weight3 = torch_tensor(weight3_r, dtype=torch_float(), requires_grad=TRUE)
    weight6 = torch_tensor(weight6_r, dtype=torch_float(), requires_grad=TRUE)
    bias0 = torch_tensor(bias0_r, dtype=torch_float(), requires_grad=TRUE)
    bias3 = torch_tensor(bias3_r, dtype=torch_float(), requires_grad=TRUE)
    bias6 = torch_tensor(bias6_r, dtype=torch_float(), requires_grad=TRUE) 
    
    model = nn_sequential(
      my_nn_linear0(weight0, bias0) ,
      act_fn ,  
      nn_dropout(p = drpot) ,
      my_nn_linear0(weight3, bias3) ,
      act_fn ,  
      nn_dropout(p = drpot) ,
      my_nn_linear0(weight6, bias6) ,
      lastact_fn 
    )    
  }
  
  y_pred_r = as.numeric(y_pred) 

#   check model on full data  
  lossf[i_] = loss$item() 
  if (family == "multiclass") {
    corrects = y_pred_t$argmax(dim=2) 
    accuracyf[i_] = as.numeric( sum(corrects == y_train_t) ) / length(y_train) 
    agreef   [i_] = accuracyf[i_] 
    lossf [i_] = loss$item() 
  } else if (family == "binomial") {
    corrects = (y_pred_r >=0.5)*(y_train > 0) + (y_pred_r <= 0.5)*(y_train < 1)
    accuracyf[i_] = sum(corrects) / length(corrects) 
    lossf [i_] = loss$item() 
    agreef   [i_] = concordance(y_train ~ y_pred_r)[[1]]   
  } else if (family == "gaussian") {
#    if (i_ == 1) { print("    HERE 3") }
    lossf [i_] = as.numeric(nnf_mse_loss(y_pred,y_train_t)) ;
    agreef[i_] = cor(y_train, y_pred_r) # ^2     
  } else if (family == "cox") { 
    if (sum(is.na(y_pred_r))>0) { y_pred_r = rep(1,length(y_pred_r)) }
    fit0 = coxph( Surv(y_train, e_train) ~ y_pred_r, init=c(1), control=coxph.control(iter.max=0)) 
    lossf [i_]  = - fit0$loglik[1] / sum(e_train) 
    agreef [i_]  =  fit0$concordance[6]
  }
  
  docat = 0 
  if ((eppr >=  0) & ((i_ <= 5) | (i_ == 10))) { docat = 1 }
  if ((eppr >= -2) & (i_ %in% c(which_loss, which_agree, epochs)))   { docat = 1 }
  if (eppr  >   0) { if (i_ %% eppr == 0) { docat = 1 } }
  if (docat == 1) {
    if (family == "cox") {
      cat(" Epoch:", i_,"Train Loss: ", lossf[i_], "Train Concordance:"  , agreef[i_], "\n")
    } else if (family == "binomial") {
      cat(" Epoch:", i_,"Train Loss: ", loss$item(), "Train Concordance:", agreef[i_], "\n")
    }else if (family == "gaussian") {
      cat(" Epoch:", i_,"Train Loss: ", loss$item(), "Train R-square:",  agreef[i_]^2, "\n")
    } else if (family == "multiclass") {
      cat(" Epoch:", i_,"Train Loss: ", loss$item(), "Train Accuracy:", accuracyf[i_], "\n")        
    }
  }
}

losslog = list(cvloss=cvloss, cvaccuracy= cvaccuracy, cvagree=cvagree, lossn=lossf, agreen=agreef) ;  

modelsumc = c(family, actvc) 
names(modelsumc) = c("family", "activation")  
modelsum = c(fold_n, epochs, lenz1, lenz2, actv, drpot, mylr, wd, l1, lasso, lscale, scale, 
             which_loss, which_agree, cv_loss, cv_agree, cv_accuracy, 
             loss$item(), agreef[i_], accuracyf[i_], agree0)
names(modelsum) = c("n folds", "epochs", "length Z1", "length Z2", "actv", "drpot", "mylr", "wd", "l1",
                    "lasso", "lscale", "scale", 
                    "which loss", "which agree", "CV loss", "CV Agree", "CV accuracy", 
                    "naive loss", "naive agree", "naive accuracy", "agree i_=0" )

parminit = list(weight0_r = weight0_r0, weight3_r = weight3_r0, weight6_r = weight6_r0, 
                bias0_r=bias0_r0, bias3_r=bias3_r0, bias6_r=bias6_r0)

weight0_r = as.matrix( model$parameters$`0.weight` )
weight3_r = as.matrix( model$parameters$`3.weight` )
weight6_r = as.matrix( model$parameters$`6.weight` )
bias0_r  = as.numeric( model$parameters$`0.bias` )
bias3_r  = as.numeric( model$parameters$`3.bias` )
bias6_r  = as.numeric( model$parameters$`6.bias` )

parmfinal= list(weight0_r = weight0_r, weight3_r = weight3_r, weight6_r = weight6_r, 
                bias0_r=bias0_r, bias3_r=bias3_r, bias6_r=bias6_r)

rlist = list(model=model, modelsum=modelsum, modelsumc=modelsumc, parminit=parminit, 
             parmfinal=parmfinal, losslog=losslog, seed=seed, foldid=foldid )

class(rlist) <- c("ann_tab_cv")

if (eppr > 0) { print(modelsum) } 

if ( docat==1 ) { cat("\n") } 

return( rlist ) 
}

################################################################################
################################################################################
#' Fit multiple Artificial Neural Network models on "tabular" provided as a matrix, and 
#' keep the best one.
#' 
#' @description Fit an multiple Artificial Neural Network models for analysis of "tabular" 
#' data using ann_tab_cv() and select the best fitting model according to cross
#' validaiton.   
#'
#' @param myxs     predictor input - an n by p matrix, where n (rows) is sample size, and p (columns) 
#' the number of predictors.  Must be in matrix form for complete data, no NA's, no Inf's, etc.,
#' and not a data frame. 
#' @param mystart  an optional vector of start times in case of a Cox model. Class numeric of length same as number of patients (n)
#' @param myy      dependent variable as a vector: time, or stop time for Cox model, Y_ 0 or 1 for binomial (logistic), numeric for gaussian. 
#' Must be a vector of length same as number of sample size. 
#' @param myevent  event indicator, 1 for event, 0 for census, Cox model only.
#' Must be a numeric vector of length same as sample size.
#' @param myoffset   an offset term to be ues when fitting the ANN.  Not yet implemented.   
#' @param family   model family, "cox", "binomial" or "gaussian" (default) 
#' @param fold_n   number of folds for each level of cross validation
#' @param epochs   number of epochs to run when tuning on number of epochs for fitting 
#' final model number of epochs informed by cross validation 
#' @param eppr     for EPoch PRint.  print summry info every eppr epochs. 0 will 
#' print first and last epochs, -1 nothing.  
#' @param lenz1    length of the first hidden layer in the neural network, default 16 
#' @param lenz2    length of the second hidden layer in the neural network, default 16 
#' @param actv     for ACTiVation function.  Activation function between layers, 
#' 1 for relu, 2 for gelu, 3 for sigmoid. 
#' @param drpot    fraction of weights to randomly zero out.  NOT YET implemented. 
#' @param mylr     learning rate for the optimization step in teh neural network model fit
#' @param wd       weight decay for the model fit.  
#' @param l1       a possible L1 penalty weight for the model fit, default 0 for not considered  
#' @param lasso    1 to indicate the first column of the input matrix is an offset 
#' term, often derived from a lasso model 
#' @param lscale Scale used to allow ReLU to extend +/- lscale before capping the 
#' inputted linear estimated 
#' @param scale Scale used to transform the initial random parameter assingments by 
#' dividing by scale
#' @param resetlw  1 as default to re-adjust weights to account for the offset every 
#' epoch.  This is only used in case lasso is set to 1 
#' @param minloss default of 1 for minimizing loss, else maximizing agreement (concordance 
#' for Cox and Binomial, R-square for Gaussian), as function of epochs by cross validation  
#' @param gotoend fit to the end of epochs.  Good for plotting and exploration 
#' @param bestof   how many models to run, from which the best fitting model will be selected.
#' @param seed an optional a numerical/integer vector of length 2, for R and torch 
#' random generators, default NULL to generate these.  Integers should be positive 
#' and not more than 2147483647.
#' @param foldid  a vector of integers to associate each record to a fold.  Should 
#' be integers from 1 and fold_n.  
#' 
#' @return an artificial neural network model fit 
#' 
#' @seealso
#'   \code{\link{ann_tab_cv}} , \code{\link{predict_ann_tab}}, \code{\link{nested.glmnetr}}
#' 
#' @author Walter Kremers (kremers.walter@mayo.edu)
#' 
#' @export
#'
ann_tab_cv_best = function(myxs, mystart=NULL, myy, myevent=NULL, myoffset=NULL, family="binomial", fold_n=5, 
                           epochs=200, eppr=40, lenz1=32, lenz2=8, actv=1, drpot=0, mylr=0.005, wd=0, l1=0, 
                           lasso=0, lscale=5, scale=1, resetlw=1, minloss=1, gotoend=0, bestof=10, seed=NULL, foldid=NULL) {
  stratified = 1 
  bestof = max(1, bestof)
  
#  seedr = NA 
  if (is.null(foldid)) {  ## foldid NOT SPECIFIED 
    if (is.null(seed)) { 
      seedr = round(runif(1)*1e9) 
    } else {
      seedr = seed[1]
      if (is.null(seedr)) { seedr = round(runif(1)*1e9)
      } else if (is.na(seedr)) { seedr = round(runif(1)*1e9) } 
    }
  } else {  ## foldid SPECIFIED  
    seedr = seed[1] 
    if (is.null(seedr)) { seedr = NA } 
  }
  
  if (is.null(seed)) { 
    seedt = round(runif(1)*1e9) 
  } else {
    seedt = seed[2]
    if (is.null(seedt)) { seedt = round(runif(1)*1e9) 
    } else if ( is.na(seedt) ) { seedt = round(runif(1)*1e9) } 
#    torch_manual_seed( seedt ) 
  }
  
  seed0 = c(seedr=seedr, seedt=seedt) 
  seed0
  
  if (is.null(foldid)) { 
    set.seed(seedr) 
    foldid = get.foldid(myy, myevent, family, fold_n, stratified) 
  } 
  
  if (eppr >= 0) { cat("\n Starting fit 1, for best of ", bestof ," model fits\n") }
  modeltemp = ann_tab_cv(myxs, mystart, myy, myevent, myoffset, family, fold_n, 
                        epochs, eppr, lenz1, lenz2, actv, drpot, mylr, wd, l1, 
                        lasso, lscale, scale, resetlw, minloss, gotoend, seed=c(NA,seed0[2]), foldid=foldid)
  cvloss  = modeltemp$modelsum[15]
  modelbesttemp = modeltemp 
  minloss = modeltemp$modelsum[15]
  concor0 = rep(0,bestof) ; 
  concor0[1] = modeltemp$modelsum[16] 
  seed1 = modeltemp$seed 
  seedts = rep(0,bestof)
  seedts[1] = modeltemp$seed[2] 
#  print(modeltemp$modelsum[c(15:20)])
  
  for (i_ in c(2:bestof)) {
    if (eppr >= 0) { cat(" Starting fit ", i_, ", for best of ", bestof ," model fits\n") }
    modeltemp = ann_tab_cv(myxs, mystart, myy, myevent, myoffset, family, fold_n, 
                           epochs, eppr, lenz1, lenz2, actv, drpot, mylr, wd, l1, 
                           lasso, lscale, scale, resetlw, minloss, gotoend, seed=c(NA,NA), foldid=foldid) 
    cvloss = modeltemp$modelsum[15]
    concor0[i_] = modeltemp$modelsum[16]
    seedts[i_] = modeltemp$seed[2] 
#    print(modeltemp$modelsum[c(15:20)])
    if (cvloss < minloss) {
      modelbesttemp = modeltemp 
      minloss = cvloss 
    }
  }
  
  seedbest = c(seed0[1], modelbesttemp$seed[2])
  names(seedbest) = c("seedr", "seedt")

  modelsum = c(modelbesttemp$modelsum,bestof=bestof)
  rlist = list( model=modelbesttemp$model, modelsum=modelsum, modelsumc=modelbesttemp$modelsumc, 
                parminit=modelbesttemp$parminit, parmfinal=modelbesttemp$parmfinal, 
                losslog = modelbesttemp$losslog, seed=seedbest, concor0=concor0, seed0=seed0, seedts=seedts ) 
  
#   list(model=model, modelsum=modelsum, modelsumc=modelsumc, parminit=parminit, 
#   parmfinal=parmfinal, losslog=losslog, seed=seed )  
  
  return(rlist)
}

################################################################################
################################################################################
