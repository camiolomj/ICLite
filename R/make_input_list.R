#' @title make_input_list
#' @description Contructs list of input parameters
#' @param input_connectivities
#' @param input_rho
#' @return List of input combinations

make_input_list<-function(input_connectivities, input_rho){

  return(as.list(crossing(input_connectivities, input_rho)))

}
