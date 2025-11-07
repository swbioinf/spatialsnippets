library(alabaster.sfe)
library(SpatialExperiment)
library(SpatialFeatureExperiment)
library(tidyverse)
library(Voyager)


data_dir               <- file.path("~/projects/spatialsnippets/datasets/GSE234713_IBDcosmx_GarridoTrigo2023/processed_data")
sfe_01_loaded          <- file.path(data_dir, "GSE234713_CosMx_IBD_sfe_01_loaded")
sfe_file_01_generic    <- file.path(data_dir, "GSE234713_CosMx_IBD_sfe_01_generic4snippettests")


sfe <- readObject(sfe_01_loaded)


# Make names 'generic'
sfe$sample     <- as.factor(sfe$tissue_sample)
sfe$individual <- as.factor(sfe$sample) # this is not a paried design, might never be useful.
sfe$slide      <- as.factor(sfe$slide_name)
sfe$cluster    <- as.factor(sfe$celltype_subset)
sfe$celltype   <- as.factor(sfe$celltype_subset)

sfe$orig.group <- sfe$group
generic_group_names <- c(CD='GroupA', UC='GroupB', HC='Control')
sfe$group     <- factor(as.character(generic_group_names[as.character(sfe$orig.group)]), levels=generic_group_names)

saveObject(sfe, sfe_file_01_generic)


################################################################################

# Also, make a generic version of the pdb, since that takes so long
se_pseudobulk_file.sample <- file.path(data_dir, "GSE234713_CosMx_IBD_pseudobulk_sample_se2.RDS") #### TEST _ Updated
sfe_file_01_generic_min200_pdb <- file.path(data_dir, "GSE234713_CosMx_IBD_sfe_01_generic4snippettests_min200_pdb.RDS")


se.pdb <- readRDS(se_pseudobulk_file.sample)


# Make names 'generic'
se.pdb$sample     <- as.factor(se.pdb$tissue_sample)
se.pdb$individual <- as.factor(se.pdb$sample) # this is not a paried design, might never be useful.
se.pdb$slide      <- as.factor(se.pdb$slide_name)
se.pdb$cluster    <- as.factor(se.pdb$celltype_subset)
se.pdb$celltype   <- as.factor(se.pdb$celltype_subset)

se.pdb$orig.group <- se.pdb$group
generic_group_names <- c(CD='GroupA', UC='GroupB', HC='Control')
se.pdb$group     <- factor(as.character(generic_group_names[as.character(se.pdb$orig.group)]), levels=generic_group_names)

saveRDS(se.pdb, file = sfe_file_01_generic_min200_pdb)



