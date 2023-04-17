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

print("script_cond_FM")
directory_output = "/shared/projects/humancd34_diff_rna_atacseq/scATACseq/exp/"
seurat_merged = readRDS("/shared/projects/humancd34_diff_rna_atacseq/scATACseq/data/merged_seurat_qc_filter_annot_normreduc.rds")
current_date = format(Sys.time(), "%Y%m%d")
DefaultAssay(seurat_merged) <- 'peaks'
Idents(seurat_merged) = "orig.ident"

# https://satijalab.org/seurat/reference/findmarkers

cond_markers_CTRL = FindMarkers(object = seurat_merged,
                                ident.1 = "CTRL",
                                ident.2 = c("DON", "2DG", "AOA"),
                                group.by = NULL,
                                subset.ident = NULL,
                                assay = "peaks",
                                slot = "data",
                                reduction = NULL,
                                features = NULL, #use all genes
                                logfc.threshold = 0.25,
                                test.use = "LR",
                                min.pct = 0.05,
                                min.diff.pct = -Inf,
                                verbose = TRUE,
                                only.pos = TRUE,
                                max.cells.per.ident = Inf,
                                random.seed = 42,
                                latent.vars = NULL,
                                min.cells.feature = 3,
                                min.cells.group = 3,
                                mean.fxn = NULL,
                                fc.name = NULL,
                                base = 2,
                                densify = FALSE)

cond_markers_CTRL$condition = "CTRL"
saveRDS(cond_markers_CTRL, paste0(directory_output, 
                                  "cond_markers_CTRL_",
                                  current_date, 
                                  ".rds"))

cond_markers_DON = FindMarkers(object = seurat_merged,
                               ident.1 = "DON",
                               ident.2 = c("CTRL", "2DG", "AOA"),
                               group.by = NULL,
                               subset.ident = NULL,
                               assay = "peaks",
                               slot = "data",
                               reduction = NULL,
                               features = NULL, #use all genes
                               logfc.threshold = 0.25,
                               test.use = "LR",
                               min.pct = 0.05,
                               min.diff.pct = -Inf,
                               verbose = TRUE,
                               only.pos = TRUE,
                               max.cells.per.ident = Inf,
                               random.seed = 42,
                               latent.vars = NULL,
                               min.cells.feature = 3,
                               min.cells.group = 3,
                               mean.fxn = NULL,
                               fc.name = NULL,
                               base = 2,
                               densify = FALSE)
cond_markers_DON$condition = "DON"
saveRDS(cond_markers_DON, paste0(directory_output, 
                                 "cond_markers_DON_",
                                 current_date, ".rds"))

cond_markers_2DG = FindMarkers(object = seurat_merged,
                               ident.1 = "2DG",
                               ident.2 = c("CTRL", "DON", "AOA"),
                               group.by = NULL,
                               subset.ident = NULL,
                               assay = "peaks",
                               slot = "data",
                               reduction = NULL,
                               features = NULL, #use all genes
                               logfc.threshold = 0.25,
                               test.use = "LR",
                               min.pct = 0.05,
                               min.diff.pct = -Inf,
                               verbose = TRUE,
                               only.pos = TRUE,
                               max.cells.per.ident = Inf,
                               random.seed = 42,
                               latent.vars = NULL,
                               min.cells.feature = 3,
                               min.cells.group = 3,
                               mean.fxn = NULL,
                               fc.name = NULL,
                               base = 2,
                               densify = FALSE)

cond_markers_2DG$condition = "2DG"
saveRDS(cond_markers_2DG, paste0(directory_output, 
                                 "cond_markers_2DG_",
                                 current_date, ".rds"))

cond_markers_AOA = FindMarkers(object = seurat_merged,
                               ident.1 = "AOA",
                               ident.2 = c("CTRL", "2DG", "DON"),
                               group.by = NULL,
                               subset.ident = NULL,
                               assay = "peaks",
                               slot = "data",
                               reduction = NULL,
                               features = NULL, #use all genes
                               logfc.threshold = 0.25,
                               test.use = "LR",
                               min.pct = 0.05,
                               min.diff.pct = -Inf,
                               verbose = TRUE,
                               only.pos = TRUE,
                               max.cells.per.ident = Inf,
                               random.seed = 42,
                               latent.vars = NULL,
                               min.cells.feature = 3,
                               min.cells.group = 3,
                               mean.fxn = NULL,
                               fc.name = NULL,
                               base = 2,
                               densify = FALSE)

cond_markers_AOA$condition = "AOA"
saveRDS(cond_markers_AOA, paste0(directory_output, 
                                 "cond_markers_AOA_",
                                 current_date, ".rds"))
