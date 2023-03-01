
####################
#### Constantes ####
####################

# Inspired of OkabeIto color-blind friendly palette
color_code = c(
  "CTRL" = "#009E73",
  "CTRL2" = "#009E73",
  "CTRL_CTRL2" = "#009E73",
  "DON" = "#0072B2",
  "2DG" = "#E69F00",
  "AOA" = "#D30707",
  "CTRLaK" = "#66FF99",
  "DONaK" = "#56B4E9",
  "2DGaK" = "#F0E442",
  "AOAaK" = "#FF8989",
  "VPA" = "#FF59B5"
)

# Lists of interesting genes to study
list_genes_TH_hemato = c('SMAD6', 'GATA1', 'GATA2', 'RUNX1', 'TAL1', 'HHEX', 'ZFPM1', 'FLI1', 'CBFA2T3', 'SPI1', 'ERG')
list_genes_metabo_tca = c('ACLY', 'ACSS1', 'IDH1', ' MDH1', 'MDH2', 'PDHA1', 'PDHA2', 'PDK1')
list_genes_metabo_glycolysis = c('G6PD', ' HK1', 'HK2', 'HK3', 'LDHA', 'LDHB', 'PKM', 'GLUT1', 'SLC1A5')
list_genes_metabo_glut = c('GLS', 'SLC1A5', 'ASCT2')
list_genes_metabo = c('ACACA', 'BCAT1', 'FASN',  'GOT1', 'GPAT1', 'MTOR')



###################
#### Functions ####
###################


# Finds last subfolder generated in the parent folder
pic_last_dir = function(parent_folder){
  dir = list.dirs(path = parent_folder, recursive = F, full.names = F)
  dir = dir[length(dir)]
  return(paste0(parent_folder,dir))
}



# Loads an RData file and returns it
loadRData <- function(fileName){
  load(fileName)
  get(ls()[ls() != "fileName"])
}

##########################
#### Function for GO ####
##########################

# https://ycl6.github.io/GO-Enrichment-Analysis-Demo/
# Ontology terms : BP -> Biological Process, MF -> Molecular Function, CC -> Cellular Component
# https://www.researchgate.net/figure/The-GO-terms-of-the-BP-CC-and-MF-categories-enrichment-of-the-118-differentially_fig1_336970414

### Perform gene ontology with clusterProfiler package (https://doi.org/10.1016/j.xinn.2021.100141)
# https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html
# https://yulab-smu.top/biomedical-knowledge-mining-book/clusterprofiler-go.html
GO_clusterProfiler_fun = function(tab_corr, list_gene, gene_universe, title_plot) {
  # gene_universe : list of background genes 
  # tab_corr : tab with the correspondence between gene symbol and Ensembl_ID
  
   # Convert Gene Name (gene) into Ensembl ID (ensembl_gene_id)
   ############################################################
    if(colnames(list_gene) != "gene") { colnames(list_gene) = "gene"}
    list_gene = left_join(x = list_gene,
                          y = tab_corr,
                          by = "gene")

    # Extract gene classification
    #############################
    ggo_BP = groupGO(gene = list_gene$ensembl_gene_id,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  level = 3,
                  keyType = "ENSEMBL",
                  readable = TRUE)
    ggo_CC = groupGO(gene = list_gene$ensembl_gene_id,
                     OrgDb = org.Hs.eg.db,
                     ont = "CC",
                     level = 3,
                     keyType = "ENSEMBL",
                     readable = TRUE)
    ggo_MF = groupGO(gene = list_gene$ensembl_gene_id,
                     OrgDb = org.Hs.eg.db,
                     ont = "MF",
                     level = 3,
                     keyType = "ENSEMBL",
                     readable = TRUE)
    ggo_BP = as.data.frame(ggo_BP@result) %>% dplyr::mutate(ontology = "BP")
    ggo_CC = as.data.frame(ggo_CC@result) %>% dplyr::mutate(ontology = "CC")
    ggo_MF = as.data.frame(ggo_MF@result) %>% dplyr::mutate(ontology = "MF")
    df_ggo = rbind(ggo_BP, ggo_CC, ggo_MF)
    print("ggo calculated")
    
    # GO over-representation analysis
    #################################
    ego = enrichGO(gene = list_gene$ensembl_gene_id,
                   OrgDb = org.Hs.eg.db,
                   universe = gene_universe,
                   ont = "ALL", # toutes les catégories
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   readable = TRUE,
                   keyType = "ENSEMBL")
   
    # Reduce term redundancy
    if (nrow(ego@result) != 0) { ego = simplify(ego,
                                                cutoff = 0.7,
                                                by = "p.adjust",
                                                select_fun = min,
                                                measure = "Wang", 
                                                semData = NULL) }
    
    df_ego = as.data.frame(ego@result)
    print("ego calculated")

    # Plot enrichment
    #################
    if (nrow(df_ego) != 0) {
      
      plot_dot = dotplot(ego, showCategory = 25) +
        ggtitle(title_plot) +
        facet_grid(ONTOLOGY ~ ., scales="free")
      
      plot_bar = barplot(ego, showCategory = 25) +
        ggtitle(title_plot) +
        facet_grid(ONTOLOGY ~ ., scales="free")
     
      output = list(ggo = df_ggo,
                    ego_result = df_ego,
                    dotplot_25 = plot_dot,
                    barplot_25 =  plot_bar)
      
    } else {
      
      print("No enriched category found.")
      
      output = list(ggo = df_ggo,
                    ego_result = df_ego)
      
    }
  
  return(output)

}
