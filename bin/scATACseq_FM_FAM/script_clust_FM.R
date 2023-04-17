library(Seurat)
library(Signac)
library(knitr)
library(dplyr)
# library(future)
# options(future.globals.maxSize=2*2^30)
# plan(strategy = "multicore", workers = 6)


"
Usage:
script_clust_FM.R <var_ident1>
" -> doc

library(docopt)
arguments <- docopt(doc)
var_ident1 = arguments$var_ident1
print(var_ident1)
print(setdiff(seq(0,15), var_ident1))
print("script_clust_FM")

directory_output = "/shared/projects/humancd34_diff_rna_atacseq/scATACseq/exp/"
seurat_obj = readRDS("/shared/projects/humancd34_diff_rna_atacseq/scATACseq/data/merged_seurat_qc_filter_annot_normreduc.rds")
current_date = format(Sys.time(), "%Y%m%d")
DefaultAssay(seurat_obj) <- 'peaks'
Idents(seurat_obj) = "peaks_snn_res.0.25"


clust_markers = FindMarkers(object = seurat_obj,
                            ident.1 = var_ident1,
                            ident.2 = setdiff(seq(0,15), var_ident1),
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

saveRDS(clust_markers, paste0(directory_output, "clust_markers_", 
                              var_ident1, "_", current_date, ".rds"))