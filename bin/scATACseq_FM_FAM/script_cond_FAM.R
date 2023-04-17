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

print("script_cond_FAM")
directory_output = "/shared/projects/humancd34_diff_rna_atacseq/scATACseq/exp/"
seurat_merged = readRDS("/shared/projects/humancd34_diff_rna_atacseq/scATACseq/data/merged_seurat_qc_filter_annot_normreduc.rds")
current_date = format(Sys.time(), "%Y%m%d")
DefaultAssay(seurat_merged) <- 'peaks'
Idents(seurat_merged) = "orig.ident"

# https://satijalab.org/seurat/reference/findallmarkers

print("Start of the function")
test = FindAllMarkers(object = seurat_merged,
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
                      densify = FALSE
)

saveRDS(object = test, paste0(directory_output, 
                              "cond_markers_FAM_", 
                              current_date, 
                              ".rds"))

