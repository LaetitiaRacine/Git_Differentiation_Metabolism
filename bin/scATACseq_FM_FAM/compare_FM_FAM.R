
library(dplyr)
library(generics)

directory = "/shared/projects/humancd34_diff_rna_atacseq/scATACseq/"

# Load FindMarkers list
cond_FM_2DG = readRDS(paste0(directory, "cond_markers_2DG_20230416.rds"))
cond_FM_AOA = readRDS(paste0(directory, "cond_markers_AOA_20230416.rds"))
cond_FM_DON = readRDS(paste0(directory, "cond_markers_DON_20230416.rds"))
cond_FM_CTRL = readRDS(paste0(directory, "cond_markers_CTRL_20230416.rds"))

# Load FindAllMarkers list
cond_FAM = readRDS(paste0(directory, "cond_markers_FAM_20230416.rds"))
cond_FAM_2DG = cond_FAM %>% dplyr::filter(cluster == "2DG") %>% dplyr::rename("condition" = "cluster") %>% dplyr::select(-gene)
cond_FAM_DON = cond_FAM %>% dplyr::filter(cluster == "DON") %>% dplyr::rename("condition" = "cluster") %>% dplyr::select(-gene)
cond_FAM_AOA = cond_FAM %>% dplyr::filter(cluster == "AOA") %>% dplyr::rename("condition" = "cluster") %>% dplyr::select(-gene)
cond_FAM_CTRL = cond_FAM %>% dplyr::filter(cluster == "CTRL") %>% dplyr::rename("condition" = "cluster") %>% dplyr::select(-gene)


# Compare
common_2DG <- generics::intersect(cond_FM_2DG, cond_FAM_2DG)
common_AOA <- generics::intersect(cond_FM_AOA, cond_FAM_AOA)
common_DON <- generics::intersect(cond_FM_DON, cond_FAM_DON)
common_CTRL <- generics::intersect(cond_FM_CTRL, cond_FAM_CTRL)

tab = data.frame(condition = c("2DG", "AOA", "DON", "CTRL"),
                 nbpeaks_FAM = c(nrow(cond_FAM_2DG), nrow(cond_FAM_AOA), nrow(cond_FAM_DON), nrow(cond_FAM_CTRL)),
                 nbpeaks_FM = c(nrow(cond_FM_2DG), nrow(cond_FM_AOA), nrow(cond_FM_DON), nrow(cond_FM_CTRL)),
                 nb_common_row = c(nrow(common_2DG), nrow(common_AOA), nrow(common_DON), nrow(common_CTRL)))
tab
