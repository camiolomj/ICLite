#' @title test_for_cell_connections
#' @description Quantifies number of positive correlations between gene modules and cells for an input gene clustering solution
#' @param solution A blockclustering solution
#' @return An integer value of total positive correlations between cells and gene modules

test_for_cell_connections<-function(solution){

  if(is.na(solution)==F){
    gene_clusters<-solution@rowclass
    names(gene_clusters)<-rownames(solution@data)
    gene_modules<-assemble_modules(gene_clusters)
    pt_score_mat<-mod_score_matrix_assembly(gene_modules)
    pmat<-cell_to_mod_pval(pt_score_mat)
    cmat<-cell_to_mod_cor(pt_score_mat)
    return(quantify_successful_connections(pmat, cmat))
  }else{return(0)}
}
