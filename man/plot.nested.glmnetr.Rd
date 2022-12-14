% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.cv.glmnetr_221210.R
\name{plot.nested.glmnetr}
\alias{plot.nested.glmnetr}
\title{Plot the cross validated relaxed lasso deviances or coefficients from a nested.glmnetr call.}
\usage{
\method{plot}{nested.glmnetr}(
  x,
  gam = NULL,
  lambda.lo = NULL,
  title = NULL,
  plup = 0,
  coefs = FALSE,
  comment = TRUE,
  ...
)
}
\arguments{
\item{x}{A nested.glmnetr output object}

\item{gam}{A specific level of gamma for plotting.  By default gamma.min will be used.}

\item{lambda.lo}{A lower limit of lambda when plotting.}

\item{title}{A title}

\item{plup}{Plot upper 95 percent two-sided confidence intervals for the deviance plots.}

\item{coefs}{Default is FALSE to plot deviances.  Option of TRUE to plot coefficients.}

\item{comment}{Default of TRUE to write to console information on lam and gam selected for output.
FALSE will suppress this write to console.}

\item{...}{Additional arguments passed to the plot function.}
}
\value{
This program returns a plot to the graphics window, and may provide 
some numerical information to the R Console.
}
\description{
Plot the cross validated relaxed lasso deviances or coefficients from a nested.glmnetr call.
}
\examples{
\donttest{
sim.data=glmnetr.simdata(nrows=1000, ncols=100, beta=NULL)
xs=sim.data$xs 
y_=sim.data$yt
event=sim.data$event
# for this example we use a small number for folds_n to shorten run time 
fit3 = nested.glmnetr(xs, NULL, y_, event, family="cox", folds_n=3) 
plot(fit3)
plot(fit3, coefs=TRUE)
}

}
\seealso{
\code{\link{plot.glmnetr}} , \code{\link{plot.cv.glmnetr}} , \code{\link{nested.glmnetr}}
}
\author{
Walter Kremers (kremers.walter@mayo.edu)
}
