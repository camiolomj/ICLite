% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mod_score.R
\name{mod_score}
\alias{mod_score}
\title{mod_score}
\usage{
mod_score(module)
}
\arguments{
\item{module}{A character string of genes matched to global object "scaled_expression_data" which is generated within the run_ICLite wrapper}
}
\value{
Numeric value for module gene expression
}
\description{
Returns the mean scaled expression value for an input set of genes
}
