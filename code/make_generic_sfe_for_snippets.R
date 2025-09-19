library(alabaster.sfe)
library(SpatialExperiment)
library(SpatialFeatureExperiment)
library(tidyverse)
library(Voyager)


data_dir               <- file.path("~/projects/spatialsnippets/datasets/GSE234713_IBDcosmx_GarridoTrigo2023/processed_data")
sfe_01_loaded          <- file.path(data_dir, "GSE234713_CosMx_IBD_sfe_01_loaded")
sfe_file_01_generic <- file.path(data_dir, "GSE234713_CosMx_IBD_sfe_01_generic4snippettests.RDS")

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



