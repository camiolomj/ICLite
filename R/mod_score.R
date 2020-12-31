#' @title mod_score
#' @description Returns the mean scaled expression value for an input set of genes
#' @param module A character string of genes matched to global object "scaled_expression_data" which is generated within the run_ICLite wrapper
#' @return Numeric value for module gene expression
#' @export

mod_score<-function(module){

  return(rowMeans(scaled_expression_data[,match(module, colnames(scaled_expression_data))]))

}
