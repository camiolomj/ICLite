#' @title solution_by_num_clust
#' @description Runs blockclustering for binary matrices created inside run_ICLite using varied numbers of clusters
#' @param number_of_clusters Input vector for number of clusters
#' @return A list of clustering solutions
#' @export

solution_by_num_clust<-function(number_of_clusters){

  return(mapply(test_clustering, test_mat_list, number_of_clusters))

}
