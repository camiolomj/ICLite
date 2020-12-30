#' @title mod_score_matrix_assembly
#' @description Wrapper for scoring module expression across individuals
#' @param mod_list A list of character strings indicating gene module membership
#' @return A matrix of module gene expression across individuals

mod_score_matrix_assembly<-function(mod_list){

  total_modules<-length(mod_list)
  score_mat<-matrix(nrow = nrow(scaled_expression_data), ncol = total_modules)

  for(r in 1:total_modules){
    score_mat[,r]<-mod_score(mod_list[[r]])
  }

  colnames(score_mat)<-paste0("Module ", 1:total_modules)
  rownames(score_mat)<-rownames(scaled_expression_data)
  return(score_mat)

}
