
start.time <- Sys.time()

#*********************
# Link with Snakefile
#*********************

"Create Seurat Object from Cell Ranger outputs

Usage:
  scATACseq_SKMK_CreateSeuratObject.R [options] <matrixh5_file> <fragment_file> <metadata_file> <output>
  scATACseq_SKMK_CreateSeuratObject.R -h | --help

Options:
  -h, --help               Show this screen
" -> doc

library(docopt)
arguments <- docopt(doc)


#*************
# Dependencies
#*************

library(Signac)
library(Seurat)
library(stringr)
library(EnsDb.Hsapiens.v86)   # v.86 = hg38 | v.75 = hg19
library(GenomeInfoDb)


#**************
# Input loading
#**************

matrixh5 = arguments$matrixh5_file
fragment = arguments$fragment_file
singlecell = arguments$metadata_file

name_cond = str_extract(string = singlecell, pattern = "[:alnum:]+(?=_singlecell.csv)")


#*********************
# Create Seurat Object
#*********************

counts <- Seurat::Read10X_h5(filename = matrixh5)

metadata <- read.csv(
  file = singlecell,
  header = TRUE,
  row.names = 1
)

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = fragment,
  min.cells = 100,
  min.features = 200
)

seurat_obj <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata)


#*****************
# Add annotations
#*****************

# Extract annotations from EnsDb database
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86, standard.chromosomes = TRUE)
seqlevelsStyle(annotations) = "UCSC"
print(table(annotations$type))

# Add annotation database into seurat object (needed for downstream analysis (QC functions))
Annotation(seurat_obj) = annotations

# Visualize first object of the list to see the annotations
seurat_obj@assays$peaks@annotation


#*****************
# Update orig.ident
#*****************

Idents(object = seurat_obj) = name_cond
seurat_obj$orig.ident = name_cond

#*****************
# Add peak_name as grange metadata
#*****************

peak_name = data.frame(paste0(seurat_obj@assays$peaks@ranges@seqnames, "-", seurat_obj@assays$peaks@ranges@ranges))
colnames(peak_name) = "peak_name"
mcols(seurat_obj@assays$peaks@ranges)$peak_name = peak_name

#************
# Save output
#************

saveRDS(object = seurat_obj, file = arguments$output)


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

