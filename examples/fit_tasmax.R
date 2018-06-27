library(quantreg)
library(splines)
library(pbs) 
library(lubridate)
library(PCICt)
library(extRemes)
library(fields)
library(classInt)
library(abind)
library(dplyr)
library(tidyr)
library(abind)
library(tailqr)

# Load tasmax data from 40 ensembles --------------------------------------
data.dir <- "./data"
load(file.path(data.dir, "seattle.Rdata"))
  
# Parameters --------------------------------------------------------------
blk_size <- 10                        # Number of years per block
nblk <- 2                             # Number of blocks
# Thresholds  
l.thresh.taus <- c(0.05, 0.075, 0.1)  # Lower threshold quantiles to define GPDs
u.thresh.taus <- c(0.9, 0.925, 0.95)  # Upper threshold quantiles to define GPDs
# Quantiles to predict
l.taus <- seq(0.001, 0.1, l = 5)      # Lower quantiles to predict
u.taus <- seq(0.9, 0.999, l = 5)      # Upper quantiles to predict
n.l.tau <- length(l.taus)             # Number of lower quantiles to predict
n.u.tau <- length(u.taus)             # Number of upper quantiles to predict
all.taus <- c(l.taus, u.taus)         # Combine the quantiles to predict
sub <- 1:(365*blk_size*nblk)          # Define subset of data to predict
sdates <- date[sub]                   # Subset the dates of observations
sy <- ymat[sub, 1]                    # First ensemble member, subset of obs

# Fit QR-GPD --------------------------------------------------------------
# Bulk 
qrf <- rqbulk(
  sy,
  sdates,
  df.lf = 3,
  df.hf = 3,
  b.formula = ~ b.lf + b.hf + b.lf:b.hf ,
  taus = c(l.taus, u.taus),
  return.yhat = T,
  return.excd = T,
  return.basis = T
)
# Extract compounded exceedances
e <- get_exceedances(
  sy,
  sdates,
  df.lf = 3,
  df.hf = 3,
  high.taus = u.thresh.taus,
  low.taus = l.thresh.taus,
  tail.df.lf = 3,
  tail.df.hf = 3,
  bulk.formula = ~ b.lf + b.hf + b.lf:b.hf,
  tail.formula = ~ b.lf + b.hf,
  use.winter.shape = F
)
# Use low-dimensional spline basis for tail
Xb.lowdim <- make_time_basis(sdates, df.lf = 3, df.hf = 3, 
                             b.formula = ~ b.lf + b.hf)
# Lower tail
lqs <- fgpdset(e$lower,
               Xb.lowdim, type = "lower", return.qs = T, 
               pred.taus = l.taus, use.winter.shape = F)$qs
lqs = abind(lqs, qrf$yhat[,1:n.l.tau], along = 3)
# Upper tail
uqs <- fgpdset(e$upper,
               Xb.lowdim, type = "upper", return.qs = T, 
               pred.taus = u.taus, use.winter.shape = F)$qs
uqs = abind(qrf$yhat[,(n.l.tau + 1):(n.l.tau + n.u.tau)], uqs, along = 3)
lqs <- calc_weighted_quantile(lqs, l.taus, l.thresh.taus, end.tau = 0.04, 
                              type = "lower")
uqs <- calc_weighted_quantile(uqs, u.taus, u.thresh.taus, end.tau = 0.96, 
                              type = "upper")
gpd.qs <- cbind(lqs, uqs)

# Plot fit ----------------------------------------------------------------
# Define color pallet
extra <- 3
pal <- tim.colors(length(all.taus) + extra)
pal <- pal[-((length(l.taus)+ 1):(length(l.taus) + extra))]
pdf("./fig/seattle_tail_quantiles.pdf")
par(mar=c(5, 4, 4, 2) + 0.1)
par(xpd = T, mar = par()$mar + c(0,0,2,0))
plot(sdates,sy, cex = 0.1, pch = 21, ylab = "Max Daily Temperature (K)", xlab = "Date", 
     main = "Upper and Lower Quantile Estimates", ylim = range(ymat))
matplot(sdates, gpd.qs[sub, ], type = "l",add = T, col = pal, lty = 1)
legend(min(sdates),max(sy)+ 15,legend = round(all.taus,3) , col = pal,
       lty = 1, lwd = 2, ncol = length(pal), cex = 0.5, bty = "n")
dev.off()
