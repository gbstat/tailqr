#' Make time basis matrix
#' @description makes a matrix of polynomial (in date) and natural spline
#' (in day of year) basis functions + interaction terms for set of dates.
#' Dates are scaled to be between (0,1) before defining splines.
#' @param dates vector of dates of each observation
#' @param df.lf natural spline degrees of freedom on date
#' @param df.hf periodic basis spline degrees of freedom on the 
#'              day of the year
#' @param b.formula formula for natural spline basis and periodic basis 
#' splines. Must be written in terms of b.lf and b.hf
#'
#' @return matrix of basis functions with all linear and interaction terms 
#'         evaluated at every date passed to function.
#' @export         
make_time_basis <- function(dates,
                            df.lf,
                            df.hf,
                            b.formula = ~ b.lf + b.hf + b.lf:b.hf) {
  my.env <- new.env()
  doy <- as.numeric(strftime(dates, format = "%j"))
  dt <- scale(dates, center = T)
  b.lf <- splines::ns(dt, df = df.lf)
  b.hf <- pbs::pbs(doy, df  = df.hf)
  # b.formula <- formula(b.formula)
  environment(b.formula) <- environment()
  # b.formula <- as.formula(as.character(b.formula), env = environment())
  Xb <- model.matrix(b.formula)
  return(Xb)
}