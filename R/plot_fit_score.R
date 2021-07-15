#' @title plot_fit_score
#'
#' @description Can be used after completed ICLite run to plot output fit scoring vs input
#' rho, connectivity and number of assumed clusters.  Size of dots are related to scaled fit
#' scoring.  Coloration indicates number of clusters in gene clustering solution.
#'
#' @return Plot of ICLite fit scoring by input parameters
#' @export

plot_fit_score<-function(){
  converted_tradeoff<-tradeoff_score+min(tradeoff_score, na.rm = T)
  converted_tradeoff[is.na(converted_tradeoff)]<-0
  Fit_Score<-20^converted_tradeoff
  combo_mat<-cbind(choice_mat, Fit_Score)
  combo_mat$number_of_clusters<-factor(as.numeric(choice_mat$number_of_clusters))

  b<-ggplot2::ggplot(combo_mat, aes(x=input_connectivities, y=input_rho)) +
    geom_point(aes(colour = number_of_clusters, size = Fit_Score)) +
    ylim(0.3, (max(combo_mat$input_rho)+0.025))+
    xlim((min(combo_mat$input_connectivities)-50), (max(combo_mat$input_connectivities)+50))

  png("ICLite fit scoring.png", height = 1200, width = 1200, res = 300)

  plot(b)

  dev.off()

  print(paste0("The optimal conditions are: rho = ", chosen_rho))
  print(paste0("Connectivity = ", chosen_connectivity))
  print(paste0("Number of clusters = ", chosen_num_clust))
}
