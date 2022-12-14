
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

# Performs gene ontology
GO_fun = function(list_gene, output_suffix, directory_output) {
  
  # Chercher la correspondance entre les gènes de notre liste et la BDD
  name_correspondance = bitr(list_gene$gene, 
                             fromType = "SYMBOL", 
                             toType = "ENTREZID", 
                             OrgDb = "org.Hs.eg.db")
  colnames(name_correspondance) = c("gene","entrezgene_id")
  
  # Garder seulement les gènes dont on a trouvé la correspondance du nom dans la BDD
  list_gene = inner_join(list_gene, 
                         name_correspondance, 
                         by = "gene")
  
  # Rechercher les catégories de gene ontology
  ego = enrichGO(gene = list_gene$entrezgene_id,
                 OrgDb = org.Hs.eg.db,
                 ont = "ALL", # toutes les catégories
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2,
                 readable = TRUE)
  
  # Possibilité de réduire la redondance des termes GO en les groupant par un score de similarité
  ego_simplified = simplify(ego,
                            cutoff = 0.7,
                            by = "p.adjust",
                            select_fun = min)
  
  # Enregistrement des objets R
  save(ego, file = paste0(directory_output,"ego_", output_suffix, ".rda")) 
  save(ego_simplified, file = paste0(directory_output,"ego_", output_suffix,".rda"))  
  
  # Enregistrement des listes de termes GO
  write.csv2(ego@result,
             file = paste0(directory_output, 
                           "Summary_ego_", 
                           output_suffix, 
                           ".csv")) 
  
  write.csv2(ego_simplified@result, 
             file = paste0(directory_output, 
                           "Summary_ego_", 
                           output_suffix, 
                           ".csv"))  
  
  # Enregistrement des plots
  pdf(file = paste0(directory_output,"plot_ego_50_", output_suffix,".pdf"), 
      width = 12, height = 18)
  print(dotplot(ego,showCategory = 50)) 
  dev.off()
  
  pdf(file = paste0(directory_output,"plot_ego_10_", output_suffix,".pdf"), 
      width = 12, height = 18)
  print(dotplot(ego,showCategory = 10)) 
  dev.off()
  
  pdf(file = paste0(directory_output,"plot_ego_simplified_50_", output_suffix,".pdf"), 
      width = 12, height = 18)
  print(dotplot(ego_simplified,showCategory = 50)) 
  dev.off() 
  
  pdf(file = paste0(directory_output,"plot_ego_simplified_10_", output_suffix,".pdf"), 
      width = 12, height = 18)
  print(dotplot(ego_simplified,showCategory = 10)) 
  dev.off()
  
  # Sortie de fonction
  output = list(ego = ego, 
                ego_simplified = ego_simplified)
  return(output)
  
}
