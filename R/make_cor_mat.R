#' @title make_cor_mat
#' @description Calculated matrix of rho values for correlation between gene modules and cells
#' @param solution Gene blockclustering solution from inside run_ICLite
#' @return Matrix of rho estimates from Spearman's rho calculation between gene module scoring and immune cell log ratios

make_cor_mat<-function(solution){

  gene_clusters<-solution@rowclass
  names(gene_clusters)<-rownames(solution@data)
  gene_modules<-assemble_modules(gene_clusters)
  mod_score_mats<-mod_score_matrix_assembly(gene_modules)
  pt_score_mat<-mod_score_matrix_assembly(gene_modules)
  return(cell_to_mod_cor(pt_score_mat))

}
