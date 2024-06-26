% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/others.R
\name{simulate_datasets_3}
\alias{simulate_datasets_3}
\title{Generation of simulated count data.}
\usage{
simulate_datasets_3(
  times = 1,
  n_data = 2,
  n_1_all = c(400, 400),
  n_p_all = c(200, 400, 600, 800),
  n_p = c(20, 40, 60, 80),
  diff = seq(10, 60, 10),
  prob_max = seq(0.1, 0.9, 0.2),
  rho_p = 0.5,
  n_x1 = 3,
  rho_x = 0.1,
  M_max = 100,
  fc1 = 1,
  marginal1 = "pg",
  copula1 = TRUE
)
}
\arguments{
\item{times}{The number of simulated data generated, default is 1;}

\item{n_data}{The number of multi-source data, default is 2;}

\item{n_1_all}{The sample size of data from different sources, default is c(400,400);}

\item{n_p_all}{The number of variables, value or vector of count data, default is c(200, 400, 600, 800);}

\item{n_p}{The default number of significant variables was n_p = c(20, 40, 60, 80).}

\item{diff}{The number of different significant variables in multi-source data sets was different;}

\item{prob_max}{The maximum zero proportion of different count data;}

\item{rho_p}{Calculation parameters of covariance matrix based on copula mode;}

\item{n_x1}{Number of covariates;}

\item{rho_x}{Parameters of the way in which the covariates were generated.}

\item{M_max}{The mean value of count data;}

\item{fc1}{Mean differences of significant variables in different distributions;}

\item{marginal1}{The marginal distribution of the count data variable, c('nb','pg'), is actually a ZINB and ZIPG distribution with zero inflation.}

\item{copula1}{Boolean value, whether to generate simulated data based on copula, defaults to TRUE.}
}
\value{
data_list: Generated simulation data, different levels represent different numbers, sample sizes, significant variable differences and different zero proportions.
res_prob_0: Zero proportional actual value.
}
\description{
Generation of simulated count data.
}
\examples{
data_pg_copula <- simulate_datasets_3(
  times = 1, n_data = 2, n_1_all = c(400, 400), n_p_all = c(50, 100, 400), n_p = c(10, 20, 40),
  diff = 0, prob_max = 0.5, marginal1 = "pg", copula1 = TRUE, fc1 = 1,
)
data_K_j <- data_pg_copula[[1]][[1]][[i]][[1]][[1]] # i = 1,2;

}
