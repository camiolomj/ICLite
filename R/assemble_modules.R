#' @title assemble_modules
#' @description Creates modules from gene clustering results
#' @param clustered_genes Gene clustering results from blockclustering solution inside run_ICLite
#' @return Gene modules

assemble_modules<-function(clustered_genes){

  return(split(x = names(clustered_genes), f = clustered_genes))

}
