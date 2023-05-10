#' Get predicteds for an Artificial Neural Network model fit in nested.glmnetr()
#' 
#' All but one of the Artificial Neural Network (ANNs) fit by nested.glmnetr() are based 
#' upon a neural network model and input from a lasso model.  Thus a simple model(xs)
#' statement will not give the proper predicted values.  This function process information 
#' form the lasso and ANN model fits to give the correct predicteds.  Whereas the ann_tab_cv()
#' function ca be used to fit a model based upon an input data set it does not fit a 
#' lasso model to allow an informed starting point for the ANN fit.  The pieces fo this
#' are in nested.glmnetr().  To fit a cross validation (CV) informed ANN model fit
#' one can run nested.glmnetr() with folds_n = 0 to derive the full data models without 
#' doing a cross validation.   
#' 
#' @param object a output object from the nested.glmnetr() function  
#' @param xs new data of the same form used as input to nested.glmnetr() 
#' @param modl ANN model entry an integer from 1 to 5 indicating which "lasso informed" 
#' ANN is to be used for calculations.  The number corresponds to the position of the 
#' ensemble input from the nested.glmnetr() call.  The model must already be fit 
#' to calculate predicteds:   
#' 1 for ensemble[1] = 1, for model based upon raw data ; 
#' 2 for ensemble[2] = 1, terms indicated by relaxed lasso model plus lasso predictions as a predictor variable (features) ;
#' 3 for ensemble[3] = 1, terms indicated by relaxed lasso model with initial weights corresponding to offset but allowed to update. ; 
#' 4 for ensemble[4] = 1, terms indicated by relaxed lasso model ;   
#' 5 for ensemble[5] = 1, terms indicated by relaxed lasso model with initial weights corresponding to offset but not allowed to update. ; 
#' 
#' @return a vector of predicteds
#' 
#' @export
#'
predict_ann_tab = function(object, xs, modl=3) {
  
  modl_ = modl
  
  if (modl != 0) {  
    ensemble = object$ensemble
    if (ensemble[modl] == 0) {
      cat(paste(" model", modl, " not run in nested CV.  Check ensemble specification \n   "))
      cat(paste(ensemble)) 
    }
  }
  
  if (modl != 0) { modl_  = c(1,3,4,2,5)[modl] } 
  
  family = object$sample[1]
  calbeta = object$ann_zb$calbeta 
  lassopred = predict(object$cv_glmnet_fit, xs, lam="lambda.min" , gam="gamma.min" )
  lassobeta = predict(object$cv_glmnet_fit,     lam="lambda.min" , gam="gamma.min" )
  if (family == "cox") {
    lassopred = calbeta * lassopred 
  } else if (family %in% c("binomial", "gaussian")) {
    lassopred = calbeta[1] + calbeta[2] * lassopred 
  }
  
  xs_means = object$ann_zb$xs_means 
  xs_sds_ = object$ann_zb$xs_sds_ 
  xs_c     = sweep(xs, 2, xs_means, FUN = "-") 
  xs_z     = sweep(xs_c, 2, xs_sds_, FUN = "/") 
  xs_z0    = xs_z[,(xs_sds_ > 1e-8)]
  if (family %in% c("binomial", "gaussian")) { 
    xs_z1 = xs_z[,(lassobeta[[1]] != 0)[2:length(lassobeta[[1]])] ]   
  } else if (family == "cox") {
    xs_z1 = xs_z[,(lassobeta[[1]] != 0)] 
  }

  xs_z2 = cbind(lasso=lassopred,xs_z1)    
  
  if        (modl_ == 1) {  
    xs_t = torch_tensor(xs_z0, dtype=torch_float()) ; 
    pred = as.numeric( object$ann_fit_1$model(xs_t) ) 
  } else if (modl_ == 2) {  
    xs_t = torch_tensor(xs_z1, dtype=torch_float()) ; 
    pred = as.numeric( object$ann_fit_2$model(xs_t) ) 
  } else if (modl_ == 3) {  
    xs_t = torch_tensor(xs_z2, dtype=torch_float()) ; 
    pred = as.numeric( object$ann_fit_3$model(xs_t) ) 
  } else if (modl_ == 4) {  
    xs_t = torch_tensor(xs_z2, dtype=torch_float()) ; 
    pred = as.numeric( object$ann_fit_4$model(xs_t) ) 
  } else if (modl_ == 5) {  
    xs_t = torch_tensor(xs_z2, dtype=torch_float()) ; 
    pred = as.numeric( object$ann_fit_5$model(xs_t) ) 
  } else if (modl_ == 0) { pred = lassopred }  

#  str(object$ann_fit_5$model$parameters)
  
  return(pred)
  
}

################################################################################

#object = nested_bin_fit_ex_p2
#object = nested_cox_fit_exp0

#preds_ann3 = predict_ann_tab( object , xs) 
#summary(preds_ann3)
 
#preds_ann1 = predict_ann_tab( object , xs, modl=1) 
#summary(preds_ann1)

#lassopred = predict_ann_tab( object , xs, modl=0) 
#summary(lassopred)

# preds_ann1 = log(preds_ann1/(1-preds_ann1)) ;  preds_ann3 = log(preds_ann3/(1-preds_ann3)) ;

#plot(lassopred ,preds_ann3)
#plot(lassopred ,preds_ann1)
#plot(preds_ann3,preds_ann1) 



  
  