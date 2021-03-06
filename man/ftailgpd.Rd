% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gpd_mixture.R
\name{ftailgpd}
\alias{ftailgpd}
\title{Fit tail GPD to threshold exceedances where scale is smoothly varying as 
a function of periodic splines}
\usage{
ftailgpd(y, dates = NULL, df.lf = NULL, df.hf = NULL, Xb = NULL,
  b.formula = ~b.lf + b.hf + b.lf:b.hf, winter.ind = NULL,
  use.winter.shape = FALSE, return.basis = FALSE)
}
\arguments{
\item{dates}{vector of dates of each observation}

\item{df.lf}{natural spline degrees of freedom on date}

\item{df.hf}{periodic basis spline degrees of freedom on the 
day of the year}

\item{Xb}{basis design matrix for tails (log(scale) of GPD). Each basis 
function is in a  separate column.}

\item{b.formula}{formula for natural spline basis and periodic basis 
splines. Must be written in terms of b.lf and b.hf}

\item{winter.ind}{should a separate shape paramter be estimated for the 
winter months?}

\item{return.basis}{(logical) should the basis matrix Xb be returned?}
}
\value{
list 
           par: estimated parameters (last column corresponds to shape)
                first several correspond to log(scale) basis terms
           Xb: basis matrix
}
\description{
Fit tail GPD to exceedances where the log(shape) parameter 
varies as a linear combination of basis functions.
}
