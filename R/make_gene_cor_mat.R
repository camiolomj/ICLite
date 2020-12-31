#' @title make_gene_cor_mat
#' @description Contructs a symmetric correlation matrix from the input gene expression data
#' @param expr_data A matrix of normalized gene expression data where columns represent individuals and rows represent features (e.g. genes)
#' @return A matrix of Spearman rho values
#' @export

make_gene_cor_mat<-function(expr_data){

  return(cor(t(expr_data), method = "spearman"))

}
