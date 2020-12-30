#' @title load_IMSA_data
#' @description Loads gene expression and immune log ratios from the IMSA data set
#' @return loads immune_cell_logratios and gene_expression_data objects into global environment


load_IMSA_data<-function(){

  immune_cell_logratios<<-read.csv(system.file("extdata", "IMSA_cell_logratios.csv", package = "ICLite"),
                                header = T, row.names = 1)

  gene_expression_data<<-as.matrix(read.csv(system.file("extdata", "IMSA_BAL_data.csv", package = "ICLite"),
                                         header = T, row.names = 1))

}
