################################################################################
##### rederive_yymmdd.R ########################################################
################################################################################
#' Rederive XGB models not kept in nested.glmnetr() output
#' 
#' @description Because the XGBoost models sometimes take large amounts of 
#' storage one may decide to set keep=0 with in the doxgb list passed to 
#' nested.glmnetr().  This function allows the user to rederive the XGBoost 
#' models without doing the search.  Note, the random forest fitting routine 
#' does not allow for (start,stop) times. 
#'
#' @param object A nested.glmnetr() output object
#' @param xs Same xs used as input to ntested.glmnetr() for input object.
#' @param y_ Same y_ used as input to ntested.glmnetr() for input object.
#' @param event Same event used as input to ntested.glmnetr() for input object.
#' @param type Same type used as input to ntested.glmnetr() for input object.
#' @param tuned 1 (default) to derive the tuned model like with xgb.tuned(), 0
#' to derive the basic models like with xgb.simple().  
#'
#' @return an output like nested.glmnetr()$xgb.simple.fitX or 
#' nested.glmnetr()$xgb.tuned.fitX for X in c("", "F", "O")
#' 
#' @seealso
#'   \code{\link{xgb.tuned}} , \code{\link{xgb.simple}} , \code{\link{nested.glmnetr}} 
#' 
#' @importFrom xgboost xgb.train xgb.DMatrix 
#'
#' @export
#' 
rederive_xgb = function(object, xs, y_, event=NULL, type="base", tuned=1) {
  
#  print(dim(xs), length(y_) , length(event)) 
  
  if (is.matrix(y_)) { y_ = as.numeric(y_) }
  
  if (is.null(type)) { type = "base" } 

  stop_ = 0 
  
  if ( object$fits[2] == 0) { 
    stop_ = 1 
    warning("  Gradient Boosting Machine model was not fit when generating input object")
  } 
  
  if ( object$doxgb$keep == 1) { 
#    stop_ = 1 
    cat("  Note, Gradient Boosting Machine was kept in original run\n")
  } 
  
  if (stop_ == 0) {
    family = object$sample[1]
    
    if (type %in% c("feat", "offs")) {
      predminR = predict(object, xs, comment=0)
      ofst = object$lasso.intcal.naive[4] + object$lasso.lincal.naive[4] * predminR
      if (family=="cox") { ofst  = ofst - mean(ofst) }
    }
    
    if (family=="cox") { Surv.xgb = ifelse( event == 1, y_, -y_) }
    
    if (type == "base") {
      if (tuned == 0) {
        if (object$doxgb$keep == 1) {
          param.final   = object$xgb.simple.fit$param.final
          nrounds.final = object$xgb.simple.fit$doxgb$nrounds.final
          doxgb = object$xgb.simple.fit$doxgb 
        } else {
          param.final = object$xgb.simple.param.final
          nrounds.final = object$doxgb_simple$nrounds.final
          doxgb = object$doxgb_simple 
        }
      } else {
        if (object$doxgb$keep == 1) {
          param.final   = object$xgb.tuned.fit$param.final
          nrounds.final = object$xgb.tuned.fit$doxgb$nrounds.final
          doxgb = object$xgb.tuned.fit$doxgb 
        } else {
          param.final = object$xgb.tuned.param.final
          nrounds.final = object$doxgb_tuned$nrounds.final
          doxgb = object$doxgb_tuned 
        }
      } 
      if (family == "cox") { xgb.dat <- xgb.DMatrix(data = xs, label = Surv.xgb)
      } else {               xgb.dat <- xgb.DMatrix(data = xs, label = y_)     }
    } else if (type == "feat") {
      if (tuned == 0) {
        if (object$doxgb$keep == 1) {
          param.final   = object$xgb.simple.fitF$param.final
          nrounds.final = object$xgb.simple.fitF$doxgb$nrounds.final
          doxgb = object$xgb.simple.fitF$doxgb
        } else {
          param.final = object$xgb.simple.param.finalF
          nrounds.final = object$doxgb_simpleF$nrounds.final
          doxgb = object$doxgb_simpleF 
        }
      } else {
        if (object$doxgb$keep == 1) {
          param.final   = object$xgb.tuned.fitF$param.final
          nrounds.final = object$xgb.tuned.fitF$doxgb$nrounds.final
          doxgb = object$xgb.tuned.fitF$doxgb 
        } else {
          param.final = object$xgb.tuned.param.finalF
          nrounds.final = object$doxgb_tunedF$nrounds.final
          doxgb = object$doxgb_tunedF
        }
      }        
      if (family == "cox") { xgb.dat <- xgb.DMatrix(data = cbind(xs,ofst), label = Surv.xgb)
      } else {               xgb.dat <- xgb.DMatrix(data = cbind(xs,ofst), label = y_)     }
    } else if (type == "offs") {
      if (tuned == 0) {
        if (object$doxgb$keep == 1) {
          param.final   = object$xgb.simple.fitO$param.final
          nrounds.final = object$xgb.simple.fitO$doxgb$nrounds.final
          doxgb = object$xgb.simple.fitO$doxgb
        } else {
          param.final = object$xgb.simple.param.finalO
          nrounds.final = object$doxgb_simpleO$nrounds.final
          doxgb = object$doxgb_simpleO
        }
      } else {
        if (object$doxgb$keep == 1) {
          param.final   = object$xgb.tuned.fitO$param.final
          nrounds.final = object$xgb.tuned.fitO$doxgb$nrounds.final
          doxgb = object$xgb.tuned.fitO$doxgb
        } else {
          param.final = object$xgb.tuned.param.finalO
          nrounds.final = object$doxgb_tunedO$nrounds.final
          doxgb = object$doxgb_tunedO
        }
      }
      if (family == "cox") { xgb.dat <- xgb.DMatrix(data = xs, label = Surv.xgb, base_margin=ofst)
      } else {               xgb.dat <- xgb.DMatrix(data = xs, label = y_, base_margin=ofst)     }
    }
    
    seed = object$seed$seedr[1] 
    set.seed(seed)    
    xgb.trained = xgb.train( params = param.final, data = xgb.dat, nrounds = nrounds.final )
    doxgb$nrounds.final = nrounds.final
    xgb.trained$doxgb = doxgb
    xgb.trained$param.final = param.final
##    xgb.trained$data=xgb.dat                                           ##<<<<<<<<<<<<<--------------------
##    xgb.trained$seed = seed 
    
    return(xgb.trained)
  }
  
}

########### rederive xgb complete ##############################################

########### rederive rf model ##################################################
#' Rederive Random Forest models not kept in nested.glmnetr() output
#' 
#' @description Because the random forest models sometimes take large amounts of 
#' storage one may decide to set keep=0 within the dorf list passed to 
#' nested.glmnetr().  This function allows the user to rederive the random forest 
#' models without doing the search. Note, the random forest fitting routine does
#' not allow for (start,stop) times. 
#'
#' @param object A nested.glmnetr() output object
#' @param xs Same xs used as input to ntested.glmnetr() for input object.
#' @param y_ Same y_ used as input to ntested.glmnetr() for input object.
#' @param event Same event used as input to ntested.glmnetr() for input object.
#' @param type Same type used as input to ntested.glmnetr() for input object.
#'
#' @return an output like nested.glmnetr()$rf_tuned_fitX for X in c("", "F", "O")
#' 
#' @seealso
#'   \code{\link{rf_tune}} , \code{\link{nested.glmnetr}} 
#' 
#' @importFrom randomForestSRC rfsrc 
#' 
#' @export
#'
rederive_rf = function(object, xs, y_, event=NULL, type=NULL) {
  
  if (is.matrix(y_)) { y_ = as.numeric(y_) }
  
  if (is.null(type)) { type = "base" }
  
  family = object$sample[1]
  
  if        (family == "binomial") { splitrule="entropy"  
  } else if (family == "gaussian") { splitrule="mse"
  } else if (family == "cox"     ) { splitrule="logrank" } 
  
  stop_ = 0 
  if ( object$fits[3] == 0) { 
    stop_ = 1 
    warning("  Random Forest model was not fit when generating input object")
  }
  
  if ( type == "ofst" ) { 
    stop_ = 1 
    warning("  Random Forest model does not support use of offset")
  }
  
  if (stop_ == 0) {
    dorf = object$dorf 
    if (dorf$keep == 1) {
      ## not necessary but for example and testing 
      if        (type == "base") { object1 = object$rf_tuned_fit 
      } else if (type == "feat") { object1 = object$rf_tuned_fitF 
      } else if (type == "offs") { object1 = object$rf_tuned_fitO }
      mtry   = object1$rf_tuned$mtry 
      ntree  = object1$rf_tuned$ntree 
      nsplit = object1$rf_tuned$nsplit
    } else {
      ## the use case ###########################
      if        (type == "base") { object1 = object$dorf_base
      } else if (type == "feat") { object1 = object$dorf_feat
      } else if (type == "offs") { object1 = object$dorf_offs } 
      mtry  = object1$mtry
      ntree  = object1$ntreec[2]
      nsplit = object1$nsplit
    }
#    c(mtry, ntree, nsplit)
    
    if (type == "feat") {
      predminR = predict(object, xs, comment=0)
      ofst = object$lasso.intcal.naive[4] + object$lasso.lincal.naive[4] * predminR
      if (family == "cox") { ofst  = ofst - mean(ofst) }
    }
  
    seed = object$seed$seedr[1] 
    set.seed(seed) 
    if (family == "cox") {
      if (type == "feat") { df = data.frame(cbind(xs, ofst, y_, event)) 
      } else { df = data.frame(cbind(xs, y_, event)) }
      rf_tuned = rfsrc( Surv(y_, event) ~ . , data=df, mtry=mtry, nsplit=nsplit, ntree=ntree, 
                        splitrule=splitrule, membership = TRUE, importance=TRUE) 
    } else { 
      if (family == "binomial") { y_ = factor( y_ ) }
      if (type == "feat") { df = data.frame(cbind(xs, ofst, y_)) 
      } else { df = data.frame(cbind(xs, y_)) }
      rf_tuned = rfsrc( y_ ~ . , data=df, mtry=mtry, nsplit=nsplit, ntree=ntree, 
                        splitrule=splitrule, membership = TRUE, importance=TRUE) 
    }
    
    rffit = list(rf_tuned=rf_tuned, err.ratev=object1$err.ratev, mtryc=object1$mtryc, ntreec=object1$ntreec, seed=seed )
    class(rffit) <- c("rf_tune") 
    return(rffit)
  }
}
      
########### rederive rf model complete #########################################

########### rederive orf model #################################################
#' Rederive Oblique Random Forest models not kept in nested.glmnetr() output
#' 
#' @description Because the oblique random forest models sometimes take large 
#' amounts of storage one may decide to set keep=0 within the doorf list passed 
#' to nested.glmnetr().  This function allows the user to rederive the oblique 
#' random forest models without doing the search. Note, the oblique random 
#' forest fitting for survival data routine does not allow for (start,stop) 
#' times. 
#'
#' @param object A nested.glmnetr() output object
#' @param xs Same xs used as input to ntested.glmnetr() for input object.
#' @param y_ Same y_ used as input to ntested.glmnetr() for input object.
#' @param event Same event used as input to ntested.glmnetr() for input object.
#' @param type Same type used as input to ntested.glmnetr() for input object.
#'
#' @return an output like nested.glmnetr()$rf_tuned_fitX for X in c("", "F", "O")
#' 
#' @seealso
#'   \code{\link{orf_tune}} , \code{\link{nested.glmnetr}} 
#' 
#' @importFrom randomForestSRC rfsrc 
#' 
#' @export
#'
rederive_orf = function(object, xs, y_, event=NULL, type=NULL) {
  
  xs = xs[,(diag(cov(xs)) != 0)]
  
  if (is.matrix(y_)) { y_ = as.numeric(y_) }
  
  if (is.null(type)) { type = "base" }
  
  family = object$sample[1]
  
#  if        (family == "binomial") { splitrule="entropy" ; splitrule="mse" ; 
#  } else if (family == "gaussian") { splitrule="mse"
#  } else if (family == "cox"     ) { splitrule="logrank" } 
  
  stop_ = 0 
  if ( object$fits[8] == 0) { 
    stop_ = 1 
    warning("  Random Forest model was not fit when generating input object")
  }
  
  if ( type == "ofst" ) { 
    stop_ = 1 
    warning("  Random Forest model does not support use of offset")
  }
  
  if (stop_ == 0) {
    doorf = object$doorf 
    if (doorf$keep == 1) {
      ## not necessary but for example and testing 
      if        (type == "base") { object1 = object$orf_tuned_fit 
      } else if (type == "feat") { object1 = object$orf_tuned_fitF 
      } else if (type == "offs") { object1 = object$orf_tuned_fitO }
      mtry   = object1$orf_tuned$mtry 
      ntree  = object1$orf_tuned$ntree 
      nsplit = object1$orf_tuned$nsplit
    } else {
      ## the use case ###########################
      if        (type == "base") { object1 = object$doorf_base
      } else if (type == "feat") { object1 = object$doorf_feat
      } else if (type == "offs") { object1 = object$doorf_offs } 
      mtry  = object1$mtry
      ntree  = object1$ntreec[2]
      
      nsplit = object1$nsplit
    }
    #    c(mtry, ntree, nsplit)
    
    if (type == "feat") {
      predminR = predict(object, xs, comment=0)
      ofst = object$lasso.intcal.naive[4] + object$lasso.lincal.naive[4] * predminR
      if (family == "cox") { ofst  = ofst - mean(ofst) }
    }
    
    seed = object$seed$seedr[1] 
    set.seed(seed) 
    
    if (family == "cox") {
      if (type == "feat") { df = data.frame(cbind(xs, ofst, y_, event)) 
      } else { df = data.frame(cbind(xs, y_, event)) }
      orf_tuned = orsf( Surv(y_, event) ~ . , data=df, mtry=mtry, n_split=nsplit, n_tree=ntree) # , split_rule=split_rule) 
    } else { 
      if (family == "binomial") { y_ = factor( y_ ) }
      if (type == "feat") { df = data.frame(cbind(xs, ofst, y_)) 
      } else { df = data.frame(cbind(xs, y_)) }
      orf_tuned = orsf( y_ ~ . , data=df, mtry=mtry, n_split=nsplit, n_tree=ntree) # , split_rule=split_rule) 
    }
    
    orffit = list(orf_tuned=orf_tuned, oobag_funv=object1$oobag_funv, mtryc=object1$mtryc, ntreec=object1$ntreec, seed=seed )
    class(orffit) <- c("orf_tune") 
    return(orffit)
  }
}

########### rederive orf model complete ########################################
