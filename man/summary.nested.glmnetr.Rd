% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.nested.glmnetr_221210.R
\name{summary.nested.glmnetr}
\alias{summary.nested.glmnetr}
\title{Summary of the fit of a cross validated tuned relaxed lasso model fit, inferred by nested cross validation.  .}
\usage{
\method{summary}{nested.glmnetr}(object, cvfit = FALSE, printg1 = FALSE, ...)
}
\arguments{
\item{object}{A nested.glmnetr output object.}

\item{cvfit}{Default is FALSE to describe the evaluation of the cross validated relaxed lasso model.
Option of TRUE will describe the cross validation tuned relaxed lasso model itself.}

\item{printg1}{TRUE to also print out the fully penalized lasso beta, else to suppress.
Only applies to cvfit=TRUE.}

\item{...}{Additional arguments passed to the summary function.}
}
\value{
- A fit summary, or a model summary.
}
\description{
Summarize the model fit from a nested.glmnetr output object, i.e. a nested, 
cross-validated relaxed lasso model.  Else summarize the cross-validated 
model fit.
}
\examples{
\donttest{
sim.data=glmnetr.simdata(nrows=1000, ncols=100, beta=NULL)
xs=sim.data$xs 
y_=sim.data$yt
event=sim.data$event
# for this example we use a small number for folds_n to shorten run time 
fit3 = nested.glmnetr(xs, NULL, y_, event, family="cox", folds_n=3)  
summary(fit3)
}

}
\seealso{
\code{\link{glmnetr.compcv}} , \code{\link{summary.cv.stepreg}} , \code{\link{nested.glmnetr}}
}
