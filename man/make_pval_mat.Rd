% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/make_pval_mat.R
\name{make_pval_mat}
\alias{make_pval_mat}
\title{make_pval_mat}
\usage{
make_pval_mat(solution)
}
\arguments{
\item{solution}{Gene blockclustering solution from inside run_ICLite}
}
\value{
Matrix of p-values from Spearman's rho calculation between gene module scoring and immune cell log ratios
}
\description{
Calculated matrix of p-values for correlation between gene modules and cells
}
