
start.time <- Sys.time()

#*********************
# Link with Snakefile
#*********************

"Apply quality control filter on the dataset

Usage:
  scATACseq_SKMK_QualityControl.R [options] <cond> <seurat_obj> <output_seurat> <output_seurat_filtered> <output_plots_panel> <output_plots_list> <output_tab>
  scATACseq_SKMK_QualityControl.R -h | --help

Options:
  -h, --help               Show this screen
" -> doc

library(docopt)
arguments <- docopt(doc)

#*************
# Dependencies
#*************

library(Seurat)
library(Signac)
library(ggplot2)
library(dplyr)
library(gridExtra)
library(grid)
library(BRGenomics) # function tidyChromosomes()


#**************
# Input loading
#**************

seurat_obj = readRDS(file = arguments$seurat_obj)
name_cond = arguments$cond

  
#************************
# Initial number of cells
#************************

# Extract number of cells and peaks 
nbcells = ncol(seurat_obj)
nbpeaks = nrow(seurat_obj)

# Store data in a tab
df_summary = data.frame("Conditions" = name_cond,
                        "Initial_Cell_Number" = nbcells,
                        "Initial_Peak_Number" = nbpeaks)


#***************************
# Nucleosome banding pattern
#***************************

seurat_obj = NucleosomeSignal(object = seurat_obj)
seurat_obj$nucleosome_group = ifelse(seurat_obj$nucleosome_signal > "2", 'NS > 2', 'NS < 2')
  
qc_plot_nucl_fragments = FragmentHistogram(object = seurat_obj) +
  geom_vline(xintercept = 147, color = "red") +  # no nucleosome limit
  geom_vline(xintercept = 294, color = "red") +  # one nucleosome limit
  geom_vline(xintercept = 441, color = "red") +  # two nucleosomes limit
  geom_vline(xintercept = 588, color = "red") +  # three nucleosomes limit
  labs(title = "Nucleosome banding pattern") +
  theme(title = element_text(face = "bold"))
  
qc_plot_nucl_signal = VlnPlot(
  object = seurat_obj,
  features = 'nucleosome_signal',
  pt.size = 0.1) + 
  theme(legend.position = 'none',
        axis.title.x = element_blank(),
        axis.text.x = element_blank()) +  
  geom_hline(yintercept = 2, color = "red") 
  
qc_plot_nucl_quality = FragmentHistogram(
    object = seurat_obj, 
    group.by = 'nucleosome_group') +
    labs(title = "Nucleosome signal") +
    theme(title = element_text(face = "bold"))


#**********************************
# Number of fragments and blacklist
#**********************************

seurat_obj$pct_reads_in_peaks = seurat_obj$peak_region_fragments / 
    seurat_obj$passed_filters * 100
seurat_obj$blacklist_ratio = seurat_obj$blacklist_region_fragments /
    seurat_obj$peak_region_fragments
  
qc_plot_pct_reads_vln = VlnPlot(object = seurat_obj,
                                features = 'pct_reads_in_peaks',
                                pt.size = 0.1) + 
    theme(legend.position = 'none',
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
qc_plot_pct_reads_point = qplot(log10(seurat_obj$passed_filters),
                                  seurat_obj$pct_reads_in_peaks) +
    labs(x = "log10 # of fragments", y = "% reads in peaks") +
    geom_hline(yintercept = 60, color = "firebrick") +
    geom_vline(xintercept = log10(5000), color = "firebrick")
  
qc_plot_blacklist = VlnPlot(object = seurat_obj,
                            features = 'blacklist_ratio',
                            pt.size = 0.1) + 
    theme(legend.position = 'none',
          axis.title.x = element_blank(),
          axis.text.x = element_blank())
  
qc_plot_peak_region_fragments = VlnPlot(object = seurat_obj,
                                        features = 'peak_region_fragments',
                                        pt.size = 0.1) + 
    theme(legend.position = 'none',
          axis.title.x = element_blank(),
          axis.text.x = element_blank())


#*********************
# TSS enrichment score 
#*********************

seurat_obj = TSSEnrichment(object = seurat_obj, fast = FALSE)
seurat_obj$high.tss = ifelse(seurat_obj$TSS.enrichment > 2, 'High', 'Low')
seurat_obj$pct_reads_in_TSS = seurat_obj$TSS_fragments / seurat_obj$passed_filters * 100
  
qc_plot_tss = TSSPlot(seurat_obj, group.by = 'high.tss') + 
    NoLegend() +
    labs(title = "TSS enrichment score") +
    theme(title = element_text(face = "bold"))
  
qc_plot_tss_reads =  qplot(seurat_obj$pct_reads_in_TSS, seurat_obj$pct_reads_in_peaks) +
    labs(x = "% reads in TSS", y = "% reads in peaks") +
    geom_abline(intercept = 0, slope = 1, linetype = 2, color = "firebrick")
  
 
#*************************************
# Visualize QC metrics and save plots
#*************************************

list_graphs = list(qc_plot_nucl_fragments, qc_plot_nucl_signal, qc_plot_nucl_quality,
              qc_plot_pct_reads_vln, qc_plot_pct_reads_point, qc_plot_blacklist,
              qc_plot_peak_region_fragments, qc_plot_tss, qc_plot_tss_reads)
ggsave(plot = arrangeGrob(grobs = list_graphs, 
                          ncol = 3, 
                          nrow = 3,
                          top = textGrob(paste("QC metrics for", name_cond))),
         filename = arguments$output_plots_panel,
         width = 40, height = 40)

saveRDS(object = list_graphs, file = arguments$output_plots_list)
saveRDS(object = seurat_obj, file = arguments$output_seurat)


#***********************************
# Number of cells after QC filtering  
#***********************************

# Initialize summary df with initial numbers of cells and peaks
vec_cell = c(paste0(name_cond, "_NbCells"))
vec_peak = c(paste0(name_cond, "_NbPeaks"))
df_filter = data.frame(
  "Condition" = c(vec_cell, vec_peak),
  "CreateChromatinAssay_filter" = c(df_summary$Initial_Cell_Number, 
                                    df_summary$Initial_Peak_Number))

# Apply first cell filter on nucleosome signal 
seurat_obj = subset(x = seurat_obj, subset = nucleosome_signal < 2)
df_filter$Nucleosome_filter = c(ncol(seurat_obj), nrow(seurat_obj))

# Apply cell filter on number of fragments in cell
seurat_obj =  subset(x = seurat_obj, 
                     subset = peak_region_fragments > 3000 & 
                       peak_region_fragments < 50000)
df_filter$Fragments_filter =  c(ncol(seurat_obj), nrow(seurat_obj))

# Apply cell filter on percentage of fragments in peaks
seurat_obj =  subset(x = seurat_obj, subset = pct_reads_in_peaks > 50)
df_filter$PctReads_filter =  c(ncol(seurat_obj), nrow(seurat_obj))

# Apply cell filter on blacklist fragments
seurat_obj =  subset(x = seurat_obj, subset = blacklist_ratio < 0.05)
df_filter$Blacklist_filter =  c(ncol(seurat_obj), nrow(seurat_obj))

# Apply cell filter on TSS enrichment score
seurat_obj =  subset(x = seurat_obj, subset = TSS.enrichment > 2)
df_filter$TSS_filter =  c(ncol(seurat_obj), nrow(seurat_obj))


#***********************************
# Number of peaks after QC filtering  
#***********************************

## Clean chromosome list - Keep only standard chromosomes
# Filter in the Grange object
seurat_obj@assays$peaks@ranges = tidyChromosomes(gr = seurat_obj@assays$peaks@ranges,
                                                 keep.X = TRUE,
                                                 keep.Y = TRUE,
                                                 keep.M = FALSE,
                                                 keep.nonstandard = FALSE,
                                                 genome = "hg38")
# Subset in the Count matrix as well
list_peak_clean = unlist(as.vector(seurat_obj@assays$peaks@ranges$peak_name))
table(rownames(seurat_obj) %in% list_peak_clean) # we need to remove the FALSE peaks
seurat_obj = subset(seurat_obj, features = list_peak_clean) # works only because we have only one assay
table(rownames(seurat_obj) %in% list_peak_clean)
df_filter$Chromosome_filter = c(ncol(seurat_obj), nrow(seurat_obj))


#**********************
# Save filtered output
#**********************

saveRDS(object = seurat_obj, file = arguments$output_seurat_filtered)
write.table(x = df_filter, file = arguments$output_tab, row.names = FALSE)


#**********
# Rsession
#**********

# Show running time 
end.time <- Sys.time()
time.taken <- round(end.time - start.time,2)
cat("Total execution time :", time.taken)

# Clean working space 
rm(list = ls())
gc()

# Show package version
sessionInfo()
