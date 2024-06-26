% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rederive_240508.R
\name{rederive_xgb}
\alias{rederive_xgb}
\title{Rederive XGB models not kept in nested.glmnetr() output}
\usage{
rederive_xgb(object, xs, y_, event = NULL, type = "base", tuned = 1)
}
\arguments{
\item{object}{A nested.glmnetr() output object}

\item{xs}{Same xs used as input to ntested.glmnetr() for input object.}

\item{y_}{Same y_ used as input to ntested.glmnetr() for input object.}

\item{event}{Same event used as input to ntested.glmnetr() for input object.}

\item{type}{Same type used as input to ntested.glmnetr() for input object.}

\item{tuned}{1 (default) to derive the tuned model like with xgb.tuned(), 0
to derive the basic models like with xgb.simple().}
}
\value{
an output like nested.glmnetr()$xgb.simple.fitX or 
nested.glmnetr()$xgb.tuned.fitX for X in c("", "F", "O")
}
\description{
Because the XGBoost models sometimes take large amounts of 
storage one may decide to set keep=0 with in the doxgb list passed to 
nested.glmnetr().  This function allows the user to rederive the XGBoost 
models without doing the search.  Note, the random forest fitting routine 
does not allow for (start,stop) times.
}
\seealso{
\code{\link{xgb.tuned}} , \code{\link{xgb.simple}} , \code{\link{nested.glmnetr}}
}
