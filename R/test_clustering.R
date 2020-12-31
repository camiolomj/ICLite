#' @title test_clustering
#' @description Creates a blockclustering solution from a binary matrix
#' @param test_mat A binary matrix generated as part of "run_ICLite"
#' @param num_clust The assumed number of clusters to be used for blockclustering
#' @return A blockcluster solution from the input binary matrix using the assumed number of clusters
#' @export

test_clustering<-function(test_mat, num_clust){

  cluster_solution<-blockcluster::coclusterBinary (test_mat, semisupervised = FALSE,
                                     nbcocluster = c(num_clust, num_clust),
                                     strategy = blockcluster::coclusterStrategy (),
                                     nbCore = 1, model = "pi_rho_epsilonkl")
  if(cluster_solution@successful==T){
    return(cluster_solution)
  }else(return(NA))


}
