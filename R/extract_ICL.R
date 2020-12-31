#' @title extract_ICL
#' @description Extracts ICL value from blockclustering solution inside run_ICLite
#' @param solution Gene blockclustering solution from inside run_ICLite
#' @return ICL value
#' @export

extract_ICL<-function(solution){

  if(is.na(solution)==F){
    return(solution@ICLvalue)
  }else{return(NA)}
}
