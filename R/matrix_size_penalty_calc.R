#' @title matrix_size_penalty_calc
#' @description Calculated matrix size penalty for optimal solution identification
#' @param solution Gene blockclustering solution from inside run_ICLite
#' @return Unscaled penalty factor
#' @export

matrix_size_penalty_calc<-function(solution){
  if(is.na(solution)==F){
    return(nrow(solution@data))
  }else{return(NA)}
}
