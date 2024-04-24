
<!-- README.md is generated from README.Rmd. Please edit that file -->

# SKMSD

<!-- badges: start -->

[![R-CMD-check](https://github.com/tsnm1/SKMSD/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/tsnm1/SKMSD/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/tsnm1/SKMSD/branch/main/graph/badge.svg)](https://app.codecov.io/gh/tsnm1/SKMSD?branch=main)
<!-- badges: end -->

Simultaneous Knockoff for Multi-Source Count Data (SKMSD) is a novel
method proposed for enhancing multi-source count data. This method is
primarily based on the Gaussian copula framework of the zero-inflated
Poisson–Gamma (ZIPG) distribution for knockoff generation. It is
combined with the Simultaneous Knockoff method to address variable
selection in multi-source data.

## Installation

You can install the development version of SKMSD from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("tsnm1/SKMSD")
```

## Example

Here is a basic example of how to solve a common problem in variable
selection during the implementation of SKMSD:

``` r
# library(SKMSD)
# library(pscl)

## basic example code

data(data_pg_copula)
i = 1 # or 2
data_K_j <- data_pg_copula[[1]][[1]][[i]][[1]][[1]] 
count_K_j = data_K_j    #[data_K_j[,1]==c,]
n_x1 <- 3
data_x <- as.data.frame(count_K_j[,c(1:n_x1+1)])
W <- as.data.frame(count_K_j[,-c(1:(n_x1+2))])
M <- count_K_j[,n_x1+2]
n_w <- dim(W)[1]
n_data = 2
y <- rep(rep(c(1,2),n_data),rep(n_w/2/n_data,n_data*2))
class_K <- rep(c(1:n_data),rep(n_w/n_data,n_data))
T_var = 1:(dim(W)[2]/10)
name_data <- names(table(class_K))
fdr = 0.2

# ZIPG_DE_S3 <- SKMSD(W = W, class_K = class_K, data_x = data_x, M = M, y = y, T_var = T_var,
#                    fdr = fdr, test_statistic = 'DE', filter_statistics = 3)
# ZIPG_GLM_S3 <- SKMSD(W = W, class_K = class_K, data_x = data_x, M = M, y = y, T_var = T_var,
#                      fdr = fdr, test_statistic = 'GLM', filter_statistics = 3)
# ZIPG_RF_S3 <- SKMSD(W = W, class_K = class_K, data_x = data_x, M = M, y = y, T_var = T_var,
#                      fdr = fdr, test_statistic = 'RF', filter_statistics = 3)
```

## Simulated data generation

In general, our package also includes a function (simulate_datasets_3)
that can simulate the generation of multi-source count data. An example
is provided below:

``` r
# data_pg_copula <- simulate_datasets_3(
#   times = 1, n_data = 2, n_1_all = c(400, 400), n_p_all = c(400, 800), n_p = c(40, 80),
#   diff = 0, prob_max = 0.5, marginal1 = "pg", copula1 = TRUE
# )
# i = 1 # i = 1,2;
# W = data_pg_copula[[1]][[1]][[i]][[1]][[1]] # i = 1,2;
```

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this. -->
<!-- You can also embed plots, for example: -->
<!-- ```{r pressure, echo = FALSE} -->
<!-- # plot(pressure) -->
<!-- ``` -->
<!-- In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN. -->
