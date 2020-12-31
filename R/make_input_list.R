#' @title make_input_list
#' @description Contructs list of input parameters
#' @param input_connectivities vector of input connectivity requirements.  Larger numbers create smaller modules
#' @param input_rho Vector of input rho exclusion values ranging from 0.1 to 0.9.  Anything below this cutoff is converted to a 0 on binary scale.  All interactions above are entered as 1 in the binary matrix.  Larger numbers create smaller modules.
#' @return List of input combinations
#' @export
#'
make_input_list<-function(input_connectivities, input_rho){

  return(as.list(crossing(input_connectivities, input_rho)))

}
