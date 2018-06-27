#' Extract compounding threshold exceedances above outer QR quantiles and 
#' GPD quantiles
#' @param y response vector
#' @param dates vector of dates of each observation
#' @param df.lf natural spline degrees of freedom on date
#' @param df.hf periodic basis spline degrees of freedom on the 
#'              day of the year
#' @param high.taus Upper quantile thresholds for fitting GPDs and extracting 
#' exceedances
#' @param low.taus Lower quantile thresholds for fitting GPDs and extracting 
#' exceedances
#' @param bulk.Xb Basis function matrix for bulk of distribution
#' @param tail.Xb Basis function matrix for tails of distribution
#' @param tail.df.lf natural spline degrees of freedom on date for tail
#' @param tail.df.hf periodic basis spline degrees of freedom on the 
#'              day of the year for tail
#' @param bulk.formula R formula for quantile regression on bulk
#' @param tail.formula R formula for log(scale) GPD parameters in tails
#' @param use.winter.shape (logical) Should a separate shape parameter be fit
#' for the winter and rest?
#'
#' @return List containing exceedances above and below the high.taus and 
#' low.taus quantiles of the temperature distributions as estimated by QR and 
#' GPD models
#' @export
get_exceedances <- function(y, 
                            dates, 
                            df.lf, 
                            df.hf, 
                            high.taus = NULL, 
                            low.taus = NULL,
                            bulk.Xb = NULL, 
                            tail.Xb = NULL,
                            tail.df.lf = df.lf,
                            tail.df.hf = df.hf,
                            bulk.formula = ~ b.lf + b.hf + b.lf:b.hf,
                            tail.formula = ~ b.lf + b.hf + b.lf:b.hf,
                            use.winter.shape = FALSE
){
  if(any(low.taus !=sort(low.taus))|any((high.taus != sort(high.taus)))){
    stop("taus should be in order")
  }
  if(is.null(bulk.Xb)){
    bulk.Xb <- make_time_basis(dates, df.lf, df.hf, bulk.formula)
  }
  if(is.null(tail.Xb)){
    tail.Xb <- make_time_basis(dates, tail.df.lf, tail.df.hf, tail.formula)
  }
  n.high <- length(high.taus)
  n.low <- length(low.taus)
  qrf <- quantreg::rq(y ~ -1 + bulk.Xb, tau = c(max(low.taus), min(high.taus)))
  betas <- qrf$coefficients
  yhat <- as.matrix(bulk.Xb%*%betas, drop = FALSE)
  bulk.Xb <- NULL
  if(use.winter.shape){
    winter.ind <- make_winter_ind(dates)
  }
  else{
    winter.ind = NULL
  }
  # Upper
  if(!is.null(high.taus)){
    upper.list = list()
    upmat <- matrix(NA, nrow = ncol(tail.Xb) + ifelse(use.winter.shape, 2,1),
                    ncol = n.high)
    thresh <- c(yhat[,2])
    for(i in 1:n.high){
      if(i != 1)  
        thresh <- c(nxt_thresh)
      uexcd.id <- which(y > thresh)
      uexcd.dates <- dates[uexcd.id]
      uexcd <- y[uexcd.id] - thresh[uexcd.id]
      upmat[,i] <- ftailgpd(uexcd, 
                            uexcd.dates, 
                            Xb = tail.Xb[uexcd.id,], 
                            use.winter.shape = use.winter.shape)$par
      if(i != n.high){
        gp_pars <- get_gp_pars(upmat[,i], tail.Xb, winter.ind = winter.ind)
        nxt_thresh <- qtailgpd(high.taus[i + 1], gp_pars$scale,gp_pars$shape, 
                               threshold = thresh, high.taus[i],type = "upper")
      }
      upper.list[[i]] <- list(excd.id = uexcd.id, 
                              excd.dates = uexcd.dates,
                              excd = uexcd, 
                              thresh = thresh,
                              tau = high.taus[i])
    }
  }
  # Lower
  if(!is.null(low.taus)){
    lower.list = list()
    lpmat <- matrix(NA, nrow = ncol(tail.Xb) + + ifelse(use.winter.shape, 2,1),
                    ncol = n.low)
    thresh <- c(yhat[,1])
    for(i in n.low:1){
      if(i != n.low)  
        thresh <- c(nxt_thresh)
      lexcd.id <- which(y < thresh)
      lexcd.dates <- dates[lexcd.id]
      lexcd <- y[lexcd.id] - thresh[lexcd.id]
      lpmat[,i] <- ftailgpd(-lexcd, 
                            lexcd.dates, 
                            Xb = tail.Xb[lexcd.id,],
                            use.winter.shape = use.winter.shape)$par
      if(i != 1){
        gp_pars <- get_gp_pars(lpmat[,i], tail.Xb, winter.ind = winter.ind)
        nxt_thresh <- qtailgpd(low.taus[i - 1], gp_pars$scale, gp_pars$shape, 
                               threshold = thresh, low.taus[i], type = "lower")
      }
      lower.list[[i]] <- list(excd.id = lexcd.id,
                              excd.dates = lexcd.dates,
                              excd = lexcd,
                              thresh = thresh, 
                              tau = low.taus[i])
    }
  }
  return(list(upper = upper.list, lower = lower.list, 
              upmat = upmat, lpmat = lpmat, tail.Xb = tail.Xb))
}

#' Fit tail GPD to threshold exceedances where scale is smoothly varying as 
#' a function of periodic splines
#' @description Fit tail GPD to exceedances where the log(shape) parameter 
#' varies as a linear combination of basis functions.
#' @inheritParams make_time_basis
#' @param Xb basis design matrix for tails (log(scale) of GPD). Each basis 
#' function is in a  separate column.
#' @param winter.ind should a separate shape paramter be estimated for the 
#' winter months?
#' @param return.basis (logical) should the basis matrix Xb be returned?
#' @return list 
#'            par: estimated parameters (last column corresponds to shape)
#'                 first several correspond to log(scale) basis terms
#'            Xb: basis matrix
#' @export            
ftailgpd <-
  function(y,
           dates = NULL,
           df.lf = NULL,
           df.hf = NULL,
           Xb = NULL,
           b.formula = ~ b.lf + b.hf + b.lf:b.hf,
           winter.ind = NULL,
           use.winter.shape = FALSE,
           return.basis = FALSE) {
    if (is.null(Xb)) {
      Xb <- make_time_basis(dates, df.lf, df.hf, b.formula)
    }
    if (use.winter.shape) {
      if (is.null(winter.ind)) {
        winter.ind <- make_winter_ind(dates)
      }
      out <- extRemes::fevd(
        y,
        data = data.frame(Xb),
        threshold = 0,
        scale.fun = ~ Xb - 1,
        shape.fun = ~ winter.ind,
        use.phi = TRUE,
        type = "GP"
      )
    }
    else{
      out <- extRemes::fevd(
        y,
        data = data.frame(Xb),
        threshold = 0,
        scale.fun = ~ Xb - 1,
        use.phi = TRUE,
        type = "GP"
      )
    }
    if (return.basis) {
      return(list(par = out$results$par, Xb = Xb))
    }
    else{
      return(list(par = out$results$par))
    }
  }


#' Fit collection of GPD models to compounding exceedances as defined in 
#' \code{get_exceedances}
#'
#' @param excd.list list of exceedances created by \code{get_exceedances}
#' @param Xb basis design matrix for tails (log(scale) of GPD). Each basis 
#' function is in a  separate column.
#' @param type one of "upper" or "lower" if fitting upper or lower tails
#' @param return.qs (logical) should quantile estimates be returned
#' @param pred.taus vector of quantiles to predict
#' @param use.winter.shape (logical) should a separate shape parameter be 
#' estimated for the winter months?
#' @param all.winter.ind indicator for which dates correspond to winter months
#' @param Xb.pred.time Basis matrix for prediction dates.
#'
#' @return
#' List of fitted parameter matrix (one column for each GPD order: log(scale), 
#' shape) and estimated quantiles (qs, of dimension time, taus, model).
#' @export
fgpdset <- function(excd.list,
                    Xb,
                    type = c("upper", "lower"),
                    return.qs = FALSE,
                    pred.taus = NULL,
                    use.winter.shape = FALSE,
                    all.winter.ind = NULL,
                    Xb.pred.time = NULL) {
  if (is.null(Xb.pred.time)) {
    Xb.pred.time <- Xb
  }
  type <- match.arg(type)
  n.gpd <- length(excd.list)
  pmat <-
    matrix(NA,
           nrow = ncol(Xb) + ifelse(use.winter.shape, 2, 1),
           ncol = n.gpd)
  for (i in 1:n.gpd) {
    if (type == "upper") {
      pmat[, i] <- ftailgpd(
        excd.list[[i]]$excd,
        excd.list[[i]]$excd.dates,
        Xb = Xb[excd.list[[i]]$excd.id,],
        use.winter.shape = use.winter.shape
      )$par
    }
    else{
      pmat[, i] <- ftailgpd(
        -excd.list[[i]]$excd,
        excd.list[[i]]$excd.dates,
        Xb = Xb[excd.list[[i]]$excd.id,],
        use.winter.shape = use.winter.shape
      )$par
    }
  }
  if (!return.qs) {
    return(list(pmat = pmat))
  }
  else{
    np.taus <- length(pred.taus)
    ndates <- nrow(Xb.pred.time)
    qs <- array(NA, dim = c(ndates, np.taus, n.gpd))
    for (i in 1:n.gpd) {
      pall <-
        get_gp_pars(pmat[, i], Xb.pred.time, winter.ind = all.winter.ind)
      qs[, , i] <- qtailgpd(
        pred.taus,
        pall$scale,
        pall$shape,
        threshold = excd.list[[i]]$thresh,
        thresh_tau = excd.list[[i]]$tau,
        type = type
      )
    }
    return(list(pmat = pmat, qs = qs))
  }
}

#' Get GPD parameters
#' @description Get GPD parameters from matrix of basis coefficients and basis function
#' design matrix.
#' @param bpars basis coefficients in same order as columns of Xb (basis matrix)
#' @param Xb basis function design matrix for prediction time points
#' @param winter.ind if a separater shape should be used for the winter months,
#' this is an indicator for which time points correspond to winter months.
#' @return list of fitted GPD parameters at each time point (scale and shape)
#' @export
get_gp_pars <- function(bpars, Xb, winter.ind = NULL){
  nscl.par <- ncol(Xb)
  scl <- exp(Xb%*%bpars[1:nscl.par]) # assumes terms columns are log-scale basis
  if(is.null(winter.ind)){
    shp <- rep(bpars[nscl.par + 1], length(scl))  # last term is shape
  }
  else{
    Xb.shp <- model.matrix(~winter.ind)
    shp <- Xb.shp %*% bpars[(nscl.par + 1):(nscl.par + 2)]
  }
  return(list(scale = scl, shape = shp))
}


#' Calculate tail quantiles from GPD
#' @param taus vector of quantiles to estimate (refers to cumulative probs of
#'             entire distribution, not just distribution of exceedances)
#' @param scale Vector of fitted scale parameters, one for each date.
#' @param shape Vector of fitted shape parameters, one for each date.
#' @param threshold Vector of thresholds above which the GPD was estimated
#' @param thresh_tau Corresponding estimate of CDF at threshold parameters
#' @param type (string) either "upper" for upper tail or "lower" for lower tail.
#'
#' @return matrix of quantile estimates (rows: dates, columns: taus)
#' @export
qtailgpd <- function(taus,
                     scale,
                     shape,
                     threshold = NULL,
                     thresh_tau = NULL,
                     type = c("upper", "lower")) {
  type <- match.arg(type)
  if (!is.null(thresh_tau)) {
    if (is.null(threshold)) {
      stop("threshold must be supplied")
    }
    if (type == "lower") {
      thresh_tau <- 1 - thresh_tau
      taus <- 1 - taus
      threshold <- -threshold
    }
    gpd_taus = (taus - thresh_tau) / (1 - thresh_tau)
  }
  else{
    gpd_taus <- taus
  }
  n <- length(scale)
  ntau <- length(taus)
  qs <- matrix(NA_real_, nrow = n, ncol = ntau)
  taul0 <- (gpd_taus <= 0)
  if (any(taul0)) {
    qs[, taul0] <- NA
  }
  for (i in 1:n) {
    qs[i, !taul0] <-  extRemes::qevd(
      gpd_taus[!taul0],
      threshold = 0,
      scale = scale[i],
      shape = shape[i],
      type = "GP"
    )
  }
  if (!is.null(threshold)) {
    qs <- qs + threshold
  }
  if (type == "lower") {
    qs <- -qs
  }
  return(qs)
}


#' Fit a weighted quantile regression and gpd model
#' @description Estimate quantiles using mixture of quantile regression and 
#' GPDs. Use triangular weighting functions to tradeoff between QR fit and GPD
#' fits.
#' @inheritParams rqcomb
#' @param y response vector samples
#' @param dates dates corresponding to response observations
#' @param thresh.l.taus lower threshold quantiles for defining collection of 
#' GPDs
#' @param thresh.u.taus upper threshold quantiles for defining collection of 
#' GPDs
#' @param use.winter.shape (logical) should a separate shape parameter be 
#' estimated for winter months?
#'
#' @export
fit_weighted_rqgpd <- function(y,
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
                               thresh.l.taus,
                               thresh.u.taus,
                               use.winter.shape = FALSE,
                               pred.dates = NULL) {
  if (is.null(pred.dates)) {
    pred.dates <- dates
  }
  
  qrf <- rqbulk(
    y,
    dates,
    df.lf = b.df.lf,
    df.hf = b.df.hf,
    b.formula = b.b.formula,
    taus = c(l.taus, taus, u.taus),
    return.yhat = F,
    return.excd = F,
    return.basis = F
  )
  
  Xb.pred <-
    make_time_basis(pred.dates,
                    df.lf = b.df.lf,
                    df.hf = b.df.hf,
                    b.formula = b.b.formula)
  qrf$yhat <- qr_predict(pred.dates, qrf$betas, Xb = Xb.pred)
  
  e <- get_exceedances(
    y,
    dates,
    df.lf = b.df.lf,
    df.hf = b.df.hf,
    high.taus = thresh.u.taus,
    low.taus = thresh.l.taus,
    tail.df.lf = t.df.lf,
    tail.df.hf = t.df.hf,
    bulk.formula = b.b.formula,
    tail.formula = t.b.formula,
    use.winter.shape = use.winter.shape
  )
  
  Xb.lowdim <-
    make_time_basis(dates,
                    df.lf = t.df.lf,
                    df.hf = t.df.hf,
                    b.formula = t.b.formula)
  Xb.pred.lowdim <-
    make_time_basis(pred.dates,
                    df.lf = t.df.lf,
                    df.hf = t.df.hf,
                    b.formula = t.b.formula)
  if (use.winter.shape) {
    all.winter.ind <- make_winter_ind(pred.dates)
  }
  else{
    all.winter.ind <- NULL
  }
  
  low.seq <- 1:length(l.taus)
  mid.seq <- (max(low.seq) + 1):(max(low.seq) + length(taus))
  up.seq <- (max(mid.seq) + 1):(max(mid.seq) + length(u.taus))
  # Lower
  lqs <- fgpdset(
    e$lower,
    Xb.lowdim,
    type = "lower",
    return.qs = T,
    pred.taus = l.taus,
    use.winter.shape = use.winter.shape,
    all.winter.ind = all.winter.ind,
    Xb.pred.time = Xb.pred.lowdim
  )$qs
  lqs = abind::abind(lqs, qrf$yhat[, low.seq], along = 3)
  # Upper
  uqs <- fgpdset(
    e$upper,
    Xb.lowdim,
    type = "upper",
    return.qs = T,
    pred.taus = u.taus,
    use.winter.shape = use.winter.shape,
    all.winter.ind = all.winter.ind,
    Xb.pred.time = Xb.pred.lowdim
  )$qs
  uqs = abind::abind(qrf$yhat[, up.seq],
              uqs, along = 3)
  # Calculate weighted quantile estimates
  lqs <-
    calc_weighted_quantile(lqs,
                           l.taus,
                           thresh.l.taus,
                           end.tau = 0,
                           type = "lower")
  uqs <-
    calc_weighted_quantile(uqs,
                           u.taus,
                           thresh.u.taus,
                           end.tau = 1,
                           type = "upper")
  return(cbind(lqs, qrf$yhat[, mid.seq], uqs))
}