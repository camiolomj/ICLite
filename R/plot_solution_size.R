#' @title plot_solution_size
#'
#' @description Can be used after completed ICLite run to plot output correlation matrix
#' size according to input parameters.  Uses objects from the global environment generated
#' by the run_ICLite wrapper to populate graph data.
#'
#' @return Plot of ICLite gene matrix according to input rho and connectivity values
#' @export

plot_solution_size<-function(){

  png("ICLite rho vs genes included.png", height = 1500, width = 1500, res = 300)

    plot(c(na.omit(unlist(mat_sizes)))~
         choice_mat$input_connectivities[which(is.na(unlist(mat_sizes))==F)]+
        choice_mat$input_rho[which(is.na(unlist(mat_sizes))==F)], lty = 1,
       pch = 19, cex = 2,
       ylab = "Genes Included in Module Assembly",
       xlab = "Rho Cutoff Value",
       col = RColorBrewer::brewer.pal(n = 9, name = "Set1")[1:length(unique(choice_mat$input_connectivities))]
            [match(choice_mat$input_connectivities[which(is.na(unlist(mat_sizes))==F)],
            c(unique(choice_mat$input_connectivities)))]
    )

    legend("topright", legend = paste0("Connectivity = ", unique(choice_mat$input_connectivities)),
           fill = RColorBrewer::brewer.pal(n = 9, name = "Set1")[1:length(unique(choice_mat$input_connectivities))],
           cex = 0.8)

  dev.off()

}
