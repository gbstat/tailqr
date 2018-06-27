#' Quantile regression for bulk of distribution
#' @description Fit quantile regression to bulk of data and extract exceedances
#' @param y response vector
#' @param taus quantiles to estimate
#' @param Xb (optional) matrix of basis functions. Each basis function is in a 
#'        separate column.
#' @param b.formula basis formula. Should be function of b.lf and b.hf        
#' @param return.basis (logical) should the basis matrix be returned?
#' @param return.yhat (logical) should matrix of predicted quantiles (see taus)
#'        be returned?
#' @param return.excd (logical) should vector of exceedances of highest quantile
#'        be returned?
#' @param type one of "upper" or "lower" if tail tau quantiles should be handled
#'        separately.
#' @param threshold.tau quantile above or below which taus are to be interpreted 
#' conditionally for distributions exceeding threshold.tau (either below or
#' above).
#' @inheritParams make_time_basis
#' @return A list containing the basis matrix (Xb), fitted quantiles (yhat), QR 
#' basis coefficients (betas), lower and upper exceedances, dates, indices, and 
#' thresholds (lexcd, lexcd.dates, lexcd.id, lthresh) and (uexcd, uexcd.dates, 
#' uexcd.id, uthresh)
#' @export
rqbulk <- function(y,
                   dates,
                   df.lf,
                   df.hf,
                   taus,
                   Xb = NULL,
                   b.formula = ~ b.lf + b.hf + b.lf:b.hf,
                   threshold.tau = NULL,
                   return.basis = TRUE,
                   return.yhat = TRUE,
                   return.excd = FALSE,
                   type = c("upper", "lower")) {
  type <- match.arg(type)
  if (!is.null(threshold.tau)) {
    if (type == "upper") {
      taus <- (taus - threshold.tau) / (1 - threshold.tau)
    }
    else{
      taus <- -(taus - threshold.tau) / threshold.tau
    }
  }
  if (is.null(Xb)) {
    Xb <- make_time_basis(dates, df.lf, df.hf, b.formula)
  }
  qrf <- quantreg::rq(y ~ -1 + Xb, tau = taus)
  betas <- qrf$coefficients
  if (return.yhat | return.excd) {
    yhat <- as.matrix(Xb %*% betas, drop = FALSE)
  }
  if (!return.basis) {
    Xb <- NULL
  }
  if (return.excd) {
    ntau <- length(taus)
    uexcd.id <- which(y > yhat[, ntau])
    uexcd.dates <- dates[uexcd.id]
    uexcd <- y[uexcd.id] - yhat[uexcd.id, ntau]
    uthresh <- yhat[, ntau]
    
    lexcd.id <- which(y < yhat[, 1])
    lexcd.dates <- dates[lexcd.id]
    lexcd <- y[lexcd.id] - yhat[lexcd.id, 1]
    lthresh <- yhat[, 1]
  }
  else{
    lthresh <- lexcd.dates <- lexcd <- lexcd.id <-
      uthresh <-  uexcd.dates <- uexcd <- uexcd.id <- NULL
  }
  if (!return.yhat) {
    yhat <- NULL
  }
  obj <- list(
    Xb = Xb,
    yhat = yhat,
    betas = betas,
    lexcd.dates = lexcd.dates,
    lexcd = lexcd,
    lexcd.id = lexcd.id,
    lthresh = lthresh,
    uexcd.dates = uexcd.dates,
    uexcd = uexcd,
    uexcd.id = uexcd.id,
    uthresh = uthresh
  )
  obj[sapply(obj, is.null)] <-
    NULL               # drop null elements
  return(obj)
}

#' Combined quantile regression
#' @description Fit different quantile regression models to the bulk quantiles
#' and the tail quantiles. This allows for a higher order basis design in
#' the bulk quantiles and a lower order basis design in the tails.
#' @inheritParams rqbulk
#' @param b.df.lf df for bulk low freq basis
#' @param b.df.hf df for bulk high freq basis
#' @param b.b.formula bulk formula involving b.lf and b.hf
#' @param t.df.lf df for tail low freq basis
#' @param t.df.hf df for tail high freq basis
#' @param t.b.formula tail formula (simpler) involving b.lf and b.hf
#' @param taus vector of intermediate (bulk) quantiles to estimate 
#' (between 0 and 1)
#' @param l.taus vector of lower tail quantiles to estimate (between 0 and 1)
#' @param u.taus vector of higher tail quantiles to estimate (between 0 and 1)
#' @param pred.dates (optional) vector of dates to make predictions on 
#' (if different from dates)
#' @return List of lower, mid, and upper quantile estimates
#' @export
rqcomb <-
  function(y,
           dates,
           b.df.lf,
           b.df.hf,
           b.b.formula,
           t.df.lf,
           t.df.hf,
           t.b.formula,
           taus,
           l.taus,
           u.taus,
           pred.dates = NULL) {
    if (is.null(pred.dates)) {
      pred.same <- TRUE
      pred.dates <- dates
    }
    else{
      pred.same <- FALSE
    }
    ntau <- length(taus)
    # Bulk
    fit <- rqbulk(
      y,
      dates,
      df.lf = b.df.lf,
      df.hf = b.df.hf,
      b.formula = b.b.formula,
      taus = taus,
      return.yhat = ifelse(pred.same, TRUE, FALSE),
      return.basis = FALSE,
      return.excd = T
    )
    if (!pred.same) {
      Xb.pred.dates <- make_time_basis(
        pred.dates,
        df.lf = b.df.lf,
        df.hf = b.df.hf,
        b.formula = b.b.formula
      )
      fit$yhat <-
        qr_predict(pred.dates, fit$betas, Xb = Xb.pred.dates)
      Xb.lowdim.pred.dates <- make_time_basis(
        pred.dates,
        df.lf = t.df.lf,
        df.hf = t.df.hf,
        b.formula =  t.b.formula
      )
      Xb.lowdim <-
        make_time_basis(dates,
                        df.lf = t.df.lf,
                        df.hf = t.df.hf,
                        b.formula =  t.b.formula)
    }
    else{
      Xb.lowdim.pred.dates <- Xb.lowdim <-
        make_time_basis(dates,
                        df.lf = t.df.lf,
                        df.hf = t.df.hf,
                        b.formula =  t.b.formula)
    }
    # Tails
    # Upper
    ufit <-
      rqbulk(
        fit$uexcd,
        fit$uexcd.dates,
        Xb = Xb.lowdim[fit$uexcd.id, ],
        return.basis = FALSE,
        taus = u.taus,
        threshold.tau = taus[ntau],
        type = "upper"
      )
    uexcd.hat <- Xb.lowdim.pred.dates %*% ufit$betas
    uyhat <- uexcd.hat + fit$yhat[, ntau]
    # Lower
    lfit <-
      rqbulk(
        -fit$lexcd,
        fit$lexcd.dates,
        Xb = Xb.lowdim[fit$lexcd.id, ],
        return.basis = FALSE,
        taus = l.taus,
        threshold.tau = taus[1],
        type = "lower"
      )
    lexcd.hat <- Xb.lowdim.pred.dates %*% lfit$betas
    lyhat <- -lexcd.hat + fit$yhat[, 1]
    
    return(list(yhat = cbind(lyhat, fit$yhat, uyhat)))
  }