#' @title create_binary_mat
#' @description Creates a binary matrix based off rho cutoff and connectivity values input into run_ICLite
#' @param rho_cutoff A vector of rho exclusion values ranging from 0.1 to 0.9.
#' @param min_connectivity A vector of minimum connectivity values.  Higher cutoffs will result in smaller gene modules
#' @return A binary matrix based off input cutoff parameters

create_binary_mat<-function(rho_cutoff, min_connectivity){
  filtered_cor_data<-ifelse(bal.cor.matrix>rho_cutoff, 1,0)
  col.sums.bal<-colSums(filtered_cor_data)
  genes_that_pass<-colnames(filtered_cor_data[,which(col.sums.bal>=min_connectivity)])
  return(filtered_cor_data[match(genes_that_pass, rownames(filtered_cor_data)),match(genes_that_pass, colnames(filtered_cor_data))])
}
