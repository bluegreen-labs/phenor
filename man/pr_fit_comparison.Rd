% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pr_fit_comparison.R
\name{pr_fit_comparison}
\alias{pr_fit_comparison}
\title{Fit multiple models for comparison}
\usage{
pr_fit_comparison(
  random_seeds = c(1, 12, 40),
  models = c("LIN", "TT", "TTs", "PTT", "PTTs", "M1", "M1s", "AT", "SQ", "SQb", "SM1",
    "SM1b", "PA", "PAb", "PM1", "PM1b", "UN", "UM1", "SGSI", "AGSI"),
  data = phenor::phenocam_DB,
  method = "GenSA",
  control = list(max.call = 10, temperature = 10000),
  par_ranges = system.file("extdata", "parameter_ranges.csv", package = "phenor",
    mustWork = TRUE),
  ncores = 1
)
}
\arguments{
\item{random_seeds}{a vector with random seeds for cross validation}

\item{models}{list of models to compare}

\item{data}{which standard or custom dataset to use}

\item{method}{optimization method to use (default = GenSA)
- GenSA :  Generalized Simulated Annealing algorithm
- genoud : GENetic Optimization Using Derivatives
- BayesianTools: various bayesian based optimization tools}

\item{control}{additional optimization control parameters
(default = list(max.call = 5000, temperature = 10000))}

\item{par_ranges}{location of the parameter ranges of the models}

\item{ncores}{number of cores to use to calculate model comparisons,
system specific and defaults to 1 threat (default = 1)}
}
\description{
Model comparison routine to facilitate model development
and quick comparisons of the skill of various models.
}
\examples{

# estimate will return the best estimated parameter set given the
# validation data
\dontrun{
my_comparison <- pr_fit_comparison(random_seeds = c(38,1),
 models = c("TT","PTT"),
 dataset = "phenocam_DB",
 par_ranges = "parameter_ranges.csv")
}
}
\keyword{comparison}
\keyword{model,}
\keyword{phenology,}
\keyword{validation,}
