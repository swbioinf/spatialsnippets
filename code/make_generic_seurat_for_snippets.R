library(Seurat)
library(tidyverse)


# Data paths
data_dir               <- file.path("~/projects/spatialsnippets/datasets/GSE234713_IBDcosmx_GarridoTrigo2023/processed_data")
seurat_file_01_loaded  <- file.path(data_dir, "GSE234713_CosMx_IBD_seurat_01_loaded.RDS")
seurat_file_01_generic <- file.path(data_dir, "GSE234713_CosMx_IBD_seurat_01_generic4snippettests.RDS")


# Load
so <- readRDS(seurat_file_01_loaded)



# Make names 'generic'
so$sample     <- as.factor(so$tissue_sample)
so$individual <- as.factor(so$sample) # this is not a paried design, might never be useful.
so$slide      <- as.factor(so$orig.ident)
so$cluster    <- as.factor(so$celltype_subset)
so$celltype   <- as.factor(so$celltype_subset)
# so fov, exists.
so$orig.group <- so$group
generic_group_names <- c(CD='GroupA', UC='GroupB', HC='Control')
so$group     <- factor(as.character(generic_group_names[as.character(so$orig.group)]), levels=generic_group_names)

saveRDS(so, seurat_file_01_generic)
