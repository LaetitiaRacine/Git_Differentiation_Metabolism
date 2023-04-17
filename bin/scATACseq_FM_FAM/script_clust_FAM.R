library(Seurat)
library(Signac)
library(knitr)
library(stringr)
library(ggplot2)
library(dplyr)
library(kableExtra)
# library(future)
# options(future.globals.maxSize=2*2^30)
# plan(strategy = "multicore", workers = 6)

print("script_clust_FAM")
directory_output = "/shared/projects/humancd34_diff_rna_atacseq/scATACseq/exp/"
seurat_obj = readRDS("/shared/projects/humancd34_diff_rna_atacseq/scATACseq/data/merged_seurat_qc_filter_annot_normreduc.rds")
current_date = format(Sys.time(), "%Y%m%d")
DefaultAssay(seurat_obj) <- 'peaks'
Idents(seurat_obj) = "peaks_snn_res.0.25"

clust_markers = FindAllMarkers(object = seurat_obj,
                      assay = "peaks",
                      features = NULL, #use all genes
                      logfcthreshold = 0.25,
                      test.use = "LR", # not default
                      slot = "data",
                      min.pct = 0.05, # not default
                      min.diff.pct = -Inf,
                      node = NULL,
                      verbose = TRUE,
                      only.pos = TRUE,
                      max.cells.per.ident = Inf,
                      random.seed = 42, # not default
                      latent.vars = NULL,
                      min.cells.feature = 3,
                      min.cells.group = 3,
                      mean.fxn = NULL,
                      fc.name = NULL,
                      base = 2,
                      return.thresh = 2,
                      densify = FALSE)

saveRDS(object = clust_markers, paste0(directory_output,
                              "clust_markers_FAM_",
                              current_date,
                              ".rds"))

