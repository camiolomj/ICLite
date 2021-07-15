#' @title mod_score
#' @description Returns the mean scaled expression value for an input set of genes
#' @param module A character string of genes matched to global object "scaled_expression_data" which is generated within the run_ICLite wrapper
#' @return Numeric value for module gene expression
#' @export

mod_score<-function(module){

  if(length(module)>1){
    return(rowMeans(scaled_expression_data[,match(module, colnames(scaled_expression_data))]))
  }else{
    return(scaled_expression_data[,match(module, colnames(scaled_expression_data))])
  }
}
