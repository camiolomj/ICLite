#' @title run_ICLite: Wrapper function for ICLite package
#'
#' @description Runs the ICLite algorithm based on input parameters provided
#' by user.  Creates output graphical files based on gene module and cell interactions,
#' CSV files of gene modules and global environment object "gene_module_lists" that
#' may be used for downstream Gene Ontology analysis.  Example data can be loaded using
#' load_IMSA_data()
#'
#' Finding the optimal solution for a data set may require running multiple iterations
#' of ICLite with varying input parameters.  Running separate tests for correlation
#' between genes in a data set is recommended to determine minimum values for input
#' rho exclusion.  For general practice, we recommend using no lower than 0.4.  The
#' number of assumed gene clusters should be considered in relation to the total size
#' of the transcriptional data set.  Though ICLite does penalize for over-clustering,
#' it will only consider solutions from the input vector.  Therefore, initial runs
#' may benefit from a broad array of values that may be narrowed on successive iterations.
#'
#' Users should include a vector of at least 2 input parameters for minimum connectivity,
#' rho cutoff and number of clusters.
#'
#' @param gene_expression_data A matrix of normalized gene expression data where columns represent individuals and rows represent features (e.g. genes)
#' @param immune_cell_logratios A matrix of corresponding cell log ratios where rows represent individuals and columns represent features (e.g. cells).  Transformation of percentage values using the "compositions" package is recommended prior to use
#' @param input_connectivities A vector of minimum connectivity values.  Higher cutoffs will result in smaller gene modules
#' @param input_rho A vector of rho exclusion values ranging from 0.3 to 0.9.  Gene correlations below this value are converted to 0 in binary space while those above are converted to 1.  Higher cutoffs will result in smaller gene modules
#' @param number_of_clusters A vector of assumed number of clusters to be used for blockclustering
#' @return Module to gene connectivity solution and gene module memberships to working directory and global environment
#' @export

run_ICLite<-function(gene_expression_data, immune_cell_logratios, input_connectivities,
                     input_rho, number_of_clusters){

    number_of_clusters<<-as.list(number_of_clusters)

    if(length(input_rho)<2){
      print("Please include at least 2 possible rho cutoffs")
    }else if(length(input_connectivities)<2){
      print("Please include at least 2 possible connectivity cutoffs")
    }else if(length(number_of_clusters)<2){
      print("Please include at least 2 possible gene cluster number inputs")
    }else if(max(input_rho)>0.9){
      print("Please use input rho values between 0.3 and 0.9")
    }else if(min(input_rho)<0.3){
      print("Please use input rho values between 0.3 and 0.9")
    }else{

      ##Make correlation matrix from transcriptional data
      bal.cor.matrix<<-make_gene_cor_mat(gene_expression_data)

      ##Make scaled expression data to be used in module scoring
      scaled_expression_data<<-scale(t(gene_expression_data))

      ##Make list of input parameter combinations from ranges of connectivity cutoffs and rho exclusions
      mat_input_list<<-make_input_list(input_connectivities, input_rho)

      ##Create binary matrices from gene expression correlation values based on input parameters
      test_mat_list<<-mapply(create_binary_mat, mat_input_list$input_rho, mat_input_list$input_connectivities)

      ##Generate gene clustering solutions using blockCluster
      test_solutions<<-lapply(number_of_clusters, solution_by_num_clust)
      unpacked_solutions<<-unlist(test_solutions)

      if(length(which(is.na(unpacked_solutions)==T))==length(unpacked_solutions)){
        print("Gene clustering failed.  Please vary your input parameters")
      }else{
        ##Extract ICL values for evaluation of clustering solutions
        ICL_values<<-lapply(unpacked_solutions, extract_ICL)

        ##Create matrix size score to penalize small solutions
        mat_sizes<<-lapply(unpacked_solutions, matrix_size_penalty_calc)

        ##Calculate number of positive correlations between gene modules and cell log ratios
        connection_values<<-lapply(unpacked_solutions, test_for_cell_connections)

        ##Index of solution input parameters
        choice_mat<<-as.data.frame(tidyr::crossing(number_of_clusters, tidyr::crossing(input_connectivities, input_rho)))

        ##Weighted scoring of solutions
        if(length(ICL_values)-length(which(is.na(ICL_values)==T))>1){
                tradeoff_score<<-scale(unlist(ICL_values))+ 1.6*scale(unlist(mat_sizes))+
          1.2*(scale(unlist(connection_values)/(as.numeric(choice_mat$number_of_clusters))))
        }else{tradeoff_score<-ifelse(!is.na(ICL_values), 1,0)}

        ##Identify optimal conditions
        chosen_num_clust<<-as.numeric(choice_mat[which(tradeoff_score==max(tradeoff_score, na.rm = T)),1])
        chosen_connectivity<<-as.numeric(choice_mat[which(tradeoff_score==max(tradeoff_score, na.rm = T)),2])
        chosen_rho<<-as.numeric(choice_mat[which(tradeoff_score==max(tradeoff_score, na.rm = T)),3])
        accepted_solution<<-unlist(unpacked_solutions[[which(tradeoff_score==max(tradeoff_score, na.rm = T))]])

        ##Construct final modules and export them as csv files
        all_mods<<-unique(accepted_solution@rowclass)
        gene_mod_membership<<-factor(accepted_solution@rowclass)
        gene_module_lists<<-list()

        for(g in 1:length(all_mods)){

          write.csv(rownames(accepted_solution@data[gene_mod_membership==all_mods[g],]),
                  paste0("Gene Module ", all_mods[g]+1, " genes.csv"))
        }

        ##Create cell vs module correlation plots for significant interactions
        cell_type_short<<-colnames(immune_cell_logratios)
        cell_type_short<<-gsub("[.]", "_", cell_type_short)

        for(g in 1:length(all_mods)){

          mod<-rownames(accepted_solution@data[gene_mod_membership==all_mods[g],])
          gene_module_lists[[g]]<<-mod
          score_for_mod<-mod_score(mod)

          for(c in 1:ncol(immune_cell_logratios)){

            test<-cor.test(score_for_mod, immune_cell_logratios[,c], method = "spearman")
            if(test$p.value<0.1&test$estimate>0){

              png(paste0(cell_type_short[c], " vs Gene Mod ", (all_mods[g]+1), ".png"), height = 800,
                  pointsize = 12, width = 800, res = 300)
              par(mar =c(3.9, 3.9, 2.5, 1.6))
              plot(score_for_mod~immune_cell_logratios[,c], pch = 16, cex.main = 1.2, cex.axis = 1.2,
                cex.lab = ifelse(nchar(cell_type_short[c])>10, .95, 1.1),
                ylab = paste0("Module ", (all_mods[g]+1), " Score"),
                xlab = paste0(cell_type_short[c]),
                main = paste(paste0("rho = ", round(test$estimate, 3)),
                           paste0("p-value = ", round(test$p.value, 4)), sep="\n"))
              abline(lm(score_for_mod~immune_cell_logratios[,c]), lty = 2, col = "red")
              dev.off()
            }
          }
        }

        names(gene_module_lists)<-paste0("Module ", (all_mods[g]+1))

        ##Create a corrplot for modules versus cells
        my_palette <<- colorRampPalette(c("purple", "black", "yellow"))(n = 799)
        cell.v.gene.cor<<-make_cor_mat(accepted_solution)
        cell.v.gene.pval<<-make_pval_mat(accepted_solution)
        colnames(cell.v.gene.cor)<-paste0("Module ", 1:chosen_num_clust)
        rownames(cell.v.gene.cor)<-cell_type_short

        png("mod solution corrplot.png", height = 5*500, width = 5*500, res = 300)
        par(mar = c( 5, 3, 3,8))
        corrplot::corrplot(cell.v.gene.cor, method="circle", is.corr=FALSE,
                 col=my_palette, tl.col="black",
                 p.mat = cell.v.gene.pval, sig.level = 0.05, insig = "blank")
        dev.off()
      }
    }
  }
