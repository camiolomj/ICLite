#' @title quantify_successful_connections
#' @description Tallies significant positive correlations
#' @param pval_mat p-values for gene to cell interations
#' @param cor_mat rho values for gene to cell interations
#' @return Integer value of positive correlations
#' @export

quantify_successful_connections<-function(pval_mat, cor_mat){

  return(length(which(pval_mat<.05&cor_mat>0)))

}
