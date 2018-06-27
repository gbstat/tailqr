

#' Weight for single observation
#' @description Calculate weights for triangular weight functions, which 
#' increase from 0 to 1 between st and md, and decrease from 1 to 0 from md to
#' ed.
#' @param x input point to evaluate weights at
#' @param st starting point (see description)
#' @param md midpoint (see description)
#' @param ed end point (see description)
#'
#'@export
weight_one <- function(x, st, md, ed) {
  # Interior GPD type
  if ((st > -Inf) & (ed < Inf)) {
    if (x <= st) {
      return(0)
    }
    else if (x <= md) {
      return((x - st) / (md - st))
    }
    else if (x <= ed) {
      return(1 - (x - md) / (ed - md))
    }
    else{
      return(0)
    }
  }
  # QR type
  if (st == -Inf) {
    if (x <= md) {
      return(1)
    }
    else if (x <= ed) {
      return(1 - (x - md) / (ed - md))
    }
    else{
      return(0)
    }
  }
  # Last GPD type
  if (ed == Inf) {
    if (x <= st) {
      return(0)
    }
    else if (x <= md) {
      return((x - st) / (md - st))
    }
    else{
      return(1)
    }
  }
}

#' Calculate triangular weights over a grid
#' @param x grid of points on which to evaluate weights
#'
#' @inheritParams weight_one
#'
#' @return estimates of triangular weights at \code{x}
#' @export
calc_weight_over_x <- function(x, st, md, ed) {
  sapply(x, weight_one, st, md, ed)
}

#' Calculate tail weights
#' @description Calculate tail weights for weighted GPDs
#' @param x grid of points at which to evaluate weights
#' @param knots thresholds above which GPDs are defined
#' @return weight matrix for each GPD model along grid \code{x}
#' @export
#' @examples
#' xs <- seq(0, 10, l = 10000)
#' w = calc_tail_weights(xs, knots= 1:5)
#' matplot(xs,w, type = "l")
calc_tail_weights <- function(x, knots) {
  knots <- sort(knots)
  nk <- length(knots)
  sts <- c(-Inf, knots[-nk])
  mds <- knots
  eds <- c(knots[-1], Inf)
  mapply(
    calc_weight_over_x,
    st = sts,
    md = mds,
    ed = eds,
    MoreArgs = list(x = x)
  )
}

#' Calculate weighted quantile
#' @description Estimate a weighted quantile given several estimates of the
#' same quantile from different QR and GPD models
#' @param qs quantile estimates from different QR and GPD models
#' @param taus quantiles to estimate
#' @param thresh.taus vector of quantiles, above which the GPDs are defined
#' @param end.tau Last quantile, above which the outer most GPD gets all 
#' the weight and the remaining have weight zero
#' @param type One of "upper" or "lower" depending on if upper or lower 
#' quantiles are being estimated
#' @return weighted estimate of a quantile
#' @export
calc_weighted_quantile <- function(qs,
                                   taus,
                                   thresh.taus,
                                   end.tau,
                                   type = c("upper", "lower")) {
  type <- match.arg(type)
  if (type == "upper")
    thresh.taus = c(thresh.taus, end.tau)
  else
    thresh.taus = c(end.tau, thresh.taus)
  qs[is.na(qs)] <- 0
  W <- calc_tail_weights(taus, thresh.taus)
  qhat <- t(apply(qs, 1, function(x)
    rowSums(x * W)))
  return(qhat)
}

#' @export
calc_u_loss <- function(qhat, qemp, target.taus, return.all = FALSE){
  dpenalty <- dbeta(target.taus, shape1 = 1/2, shape2 = 1/2)
  ss <- (qemp - qhat)^2
  loss <- ss%*%dpenalty
  if(!return.all){
    return(sum(loss))
  }
  else{
    weighted_ss = t(t(ss)*dpenalty)
    return(list(loss = loss, ss  = ss, 
                weighted_ss = weighted_ss, loss_by_tau = colSums(weighted_ss)))
  }
}

#' Make winter indicator
#' @description Create an indicator for whether a given date is between 
#' December and February
#' @param dates vector of dates to create indicator for
#' @return indicator of which dates are in winter
#' @export
make_winter_ind <- function(dates) {
  (format(dates, "%m") %in% c('12', '01', '02')) * 1L
}
