# Quantile Regression - GPD Mixture Estimation

This R package can be used for estimating periodic time series data using quantile regression for the bulk of a distribution, and a weighted mixture of quantile regression estimates and parametric Generalized Pareto Distribution (GPD) estimates for the tails of a distribution. 

Installation can be done using the `devtools` package as follows:

1. First, install the `devtools` package in R
```R
install.packages("devtools")
```

2. Then, load the devtools package.
```R
library(devtools)
```

3. Install the pacakge from GitHub with
```R
install_github("gbstat/tailqr")
```

For an example of how the `tailqr` package can be used to estimate the tail quantiles of daily maximum temperature (tasmax) data, see the `fit_model.R` file in the `/examples` directory.
