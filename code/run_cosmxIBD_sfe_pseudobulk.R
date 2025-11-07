
library(alabaster.sfe)
library(SpatialExperiment)
library(scuttle)
library(scater)
library(tidyverse)
library(limma)
library(DT)
library(edgeR)
library(BiocParallel)


data_dir              <- file.path("~/projects/spatialsnippets/datasets/GSE234713_IBDcosmx_GarridoTrigo2023/processed_data")
sfe_01_loaded <- file.path(data_dir, "GSE234713_CosMx_IBD_sfe_01_loaded")
se_pseudobulk_file.fov    <- file.path(data_dir, "GSE234713_CosMx_IBD_pseudobulk_fov_se.RDS")
se_pseudobulk_file.sample <- file.path(data_dir, "GSE234713_CosMx_IBD_pseudobulk_sample_se2.RDS") #### TEST _ Updated


sfe <- readObject(sfe_01_loaded)



min_reads_per_cell <- 200


sfe <- sfe[,sfe$total_count>= min_reads_per_cell]

# # at the FOV level
# sfe$pdb_sample <- paste0(sfe$tissue_sample, '_', sfe$fov, '_', sfe$celltype_subset)
# se.pdb <- aggregateAcrossCells(sfe, ids=sfe$pdb_sample,
#                                BPPARAM = MulticoreParam(workers=8) # Use 8 cores, requires BiocParallel
# )
# saveRDS(se.pdb, se_pseudobulk_file.fov)

# At the sample elvel
sfe$pdb_sample <- paste0(sfe$tissue_sample, '_', sfe$celltype_subset)
se.pdb <- aggregateAcrossCells(sfe, ids=sfe$pdb_sample,
                               BPPARAM = MulticoreParam(workers=8) # Use 8 cores, requires BiocParallel
)
saveRDS(se.pdb, se_pseudobulk_file.sample)

print("Done.")
