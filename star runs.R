
###########################################Install Packages###################################################
library(devtools)
install_github("camiolomj/ICLite")
install_github("cran/blockcluster")
install_github("cran/ape")
install.packages("compositions")##Optional; for transformation of cell percentages if necessary
BiocManager::install("org.Hs.eg.db")
BiocManager::install("topGO")
BiocManager::install("GOSemSim")

#######################################Unique Functions#######################################################

map2color<-function(x,pal,limits=NULL){
  if(is.null(limits)) limits=range(x)
  pal[findInterval(x,seq(limits[1],limits[2],length.out=length(pal)+1), all.inside=TRUE)]
}


###########################################Package Loading####################################################
library(ICLite)
library(ggplot2)
library('org.Hs.eg.db')
library(topGO)
library(GOSemSim)
library(ape)
library(RColorBrewer)

##########################################Generate Initial ICLite Results####################################

##Create input parameters for the initial run
input_connectivities<-c(150,300)
input_rho<-c(0.4, 0.5, 0.6)
number_of_clusters<-as.list(c(12,24,36))

##Load the IMSA dataset
load_IMSA_data()

##Rename the immune cells
cell_type_short<-c("Other Monocytic",
                   "Alveolar Macrophage",
                   "CD206+ CD11b+ CCR4+ Mac",
                   "CD206+ FceRI+ CCR4+ Mac",
                   "CD206+ FceRI+ CCR4- Mac",
                   "CD14+ Monocytes",
                   "FceRI+ CD127+ CCR4+ Innate",
                   "FceRI+ CD127+ CCR4- Innate",
                   "FceRI+ CD127- Innate",
                   "Naive CD4 T cell",
                   "CD4 TRM",
                   "CD4 CCR5+ CD161+",
                   "CD4 EM",
                   "CD4 CM",
                   "CD4 Treg",
                   "Naive CD8 T cells",
                   "CD8 TRM",
                   "CD8 EM",
                   "CD8 CM",
                   "CD8 Treg",
                   "gdT cell",
                   "Other CD3+",
                   "B cells",
                   "NK or NK T cells")

colnames(immune_cell_logratios)<-cell_type_short

##Initialize the first run of ICLite
set.seed(95)
run_ICLite(
  gene_expression_data = gene_expression_data,
  immune_cell_logratios = immune_cell_logratios,
  input_connectivities = input_connectivities,
  input_rho = input_rho,
  number_of_clusters = number_of_clusters
)

##save your workspace
save.image(file='full toy run ICLite.RData')

##plot gene matrix size vs input parameters
plot_solution_size()

##plot fit scoring
plot_fit_score()

##########################################Generate Focused ICLite Results####################################

##Create more focused input parameters based off previous results
input_connectivities<-c(75,150,225)
input_rho<-c(0.4, 0.425, 0.45)
number_of_clusters<-as.list(c(22,24,26))

##Initialize the second run of ICLite
set.seed(95)
run_ICLite(
  gene_expression_data = gene_expression_data,
  immune_cell_logratios = immune_cell_logratios,
  input_connectivities = input_connectivities,
  input_rho = input_rho,
  number_of_clusters = number_of_clusters
)

save.image(file='focused toy run ICLite.RData')

##plot gene matrix size vs input parameters for new run
plot_solution_size()

##plot fit scoring for new run
plot_fit_score()

###############################Perform Gene ontology enrichment analysis###########################
names(gene_module_lists)<-paste0("Module ",unique(accepted_solution@rowclass)+1)

##Create a color palette for bar plots
mypal <- colorRampPalette( c( "white", "yellow", "red" ) )( 200 )
mod_GOresults<-list()

##Set GO terms to display
nodes_to_display<-4

##Set node number for topGO
node_size<-15

##Create background gene list by converting all gene names from transcriptional data set
background_names<-unique(mapIds(org.Hs.eg.db, rownames(gene_expression_data), 'ENTREZID', 'SYMBOL'))

##Cycle GO process for all modules from ICLite

for(n in 1:length(gene_module_lists)){
  
  if(length(gene_module_lists[[n]]>1)){
    
    ##Generate a GOI list
    GOI<-unique(mapIds(org.Hs.eg.db, gene_module_lists[[n]], 'ENTREZID', 'SYMBOL'))
    geneList <- factor(as.integer(background_names %in% GOI))
    names(geneList) <- background_names
    
    ##Run GO analysis with topGO
    GOdata <- new("topGOdata", allGenes = geneList,
                  nodeSize = node_size, ontology = "BP", annot = annFUN.org, mapping = "org.Hs.eg.db")
    resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
    
    ##Generate table
    upRes <- GenTable(GOdata, classicFisher = resultFisher, ranksOf = "classicFisher", topNodes = nodes_to_display, numChar = 40)
    upRes$classicFisher[which(upRes$classicFisher=="<1e-30"|upRes$classicFisher=="< 1e-30")]<-1*10^-30
    upRes_order<-upRes[order(-log10(as.numeric(upRes$classicFisher)), decreasing = F),]
    allRes<-rbind(upRes_order)
    
    ##Make an output barplot
    plot_data<- -log10(as.numeric(allRes$classicFisher))
    png(paste0("GO_barplot_", names(gene_module_lists)[n], "_module.png"), height = 200*nodes_to_display, width = 800*5, res = 300, pointsize = 15)
    par(mar=c(5.1,40.1,1.6,2.1))
    barplot(plot_data, horiz = T, axes = T, col = map2color(plot_data, mypal),
           names.arg = allRes$Term, axisnames = T, las = 1, cex.axis = 1.25, cex.main = 1.15,
           xlab = "-log10(p-value)",
           main = paste0(names(gene_module_lists[n]), " GO enrichment"))
    abline(v = 0, lty = 1, col = "black")
    abline(v = -log10(.05), lty = 2, col = "grey50")
    dev.off()
    
    ##Store GO results for semantic clustering
    mod_GOresults[[n]]<-resultFisher
    
  }
}

###########################Perform Semantic Similarity Clustering of Modules##########################

##Filter for significant GO enrichment
go_terms_list<-list()

for(g in 1:length(mod_GOresults)){
  
  go_terms_list[[g]]<-mod_GOresults[[g]]@score[which(mod_GOresults[[g]]@score<0.05)]
  
}

##Create GO library
hsGO <- godata('org.Hs.eg.db', ont="BP")

##Create a semantic similarity score matrix
GO_semantic_mat<-matrix(0, ncol = length(gene_module_lists), nrow = length(gene_module_lists))


for(c in 1:ncol(GO_semantic_mat)){
  
  gs1<-rownames(as.matrix(go_terms_list[[c]]))
  
  for(r in 1:nrow(GO_semantic_mat)){
    
    gs2<-rownames(as.matrix(go_terms_list[[r]]))
    GO_semantic_mat[r,c]<-mgoSim(gs1, gs2, semData=hsGO, measure="Wang", combine="BMA")
    
  }
}

rownames(GO_semantic_mat)<-unique(accepted_solution@rowclass)+1
colnames(GO_semantic_mat)<-unique(accepted_solution@rowclass)+1

##Create a phylogram for visualization

GO_hclust = hclust(dist(GO_semantic_mat), method = "ward.D2")
GO_dendro = as.dendrogram(GO_hclust)

tree_cut<-6
GO_semantic_clusters<- cutree(GO_hclust, k = tree_cut)
GO_phylo<-as.phylo(GO_hclust)

cols <- sample(brewer.pal(8, "Dark2"), 7)

png("ICLite module phylogram.png", height = 2400, width = 2400, res = 300)
plot(GO_phylo, type = "unrooted", cex = 1.2, 
     no.margin = F, label.offset = 0.03, tip.color=cols[GO_semantic_clusters])
dev.off()

############################Using Module Scoring for Model Building##################

##Create mod score values for cohort
individual_mod_scores<-do.call(cbind, lapply(gene_module_lists, mod_score))
colnames(individual_mod_scores)<-paste0("Mod_", unique(accepted_solution@rowclass)+1)

#########################Examples from Troubleshooting###############################


##Create input parameters that are too stringent 

input_connectivities<-c(1000, 2000)
input_rho<-c(0.85, 0.95)
number_of_clusters<-as.list(c(24,36))

##Initialize run of ICLite
set.seed(95)
run_ICLite(
  gene_expression_data = gene_expression_data,
  immune_cell_logratios = immune_cell_logratios,
  input_connectivities = input_connectivities,
  input_rho = input_rho,
  number_of_clusters = number_of_clusters
)

##Slightly better parameters
input_connectivities<-c(1000, 2000)
input_rho<-c(0.45, 0.55)
number_of_clusters<-as.list(c(24,36))

set.seed(95)
run_ICLite(
  gene_expression_data = gene_expression_data,
  immune_cell_logratios = immune_cell_logratios,
  input_connectivities = input_connectivities,
  input_rho = input_rho,
  number_of_clusters = number_of_clusters
)

##Output test matrix size:
for(m in 1:length(test_mat_list)){
  
  print(attributes(test_mat_list[[m]])$dim)
  
}

##Working parameters
input_connectivities<-c(750, 1000)
input_rho<-c(0.45, 0.55)
number_of_clusters<-as.list(c(24,36))

set.seed(95)
run_ICLite(
  gene_expression_data = gene_expression_data,
  immune_cell_logratios = immune_cell_logratios,
  input_connectivities = input_connectivities,
  input_rho = input_rho,
  number_of_clusters = number_of_clusters
)