% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_clustering.R
\name{test_clustering}
\alias{test_clustering}
\title{test_clustering}
\usage{
test_clustering(test_mat, num_clust)
}
\arguments{
\item{test_mat}{A binary matrix generated as part of "run_ICLite"}

\item{num_clust}{The assumed number of clusters to be used for blockclustering}
}
\value{
A blockcluster solution from the input binary matrix using the assumed number of clusters
}
\description{
Creates a blockclustering solution from a binary matrix
}
