
# Build another subset object- this one to have ALL of the samples (so example will work)
# But subsampled to only 5 fovs of each.
# THis means we can still run the samples, but the RAM requirements are much better.
# ( Full dataset will not process with 16G ram. )
# This shoudl be a good example for self directed learning.

# There is already a 2-condition only subset for use in training - that is 6 samples with 4 fovs each.
# See: https://github.com/swbioinf/spatial-taster-workshop/blob/main/code/subsample_sfe_for_workshop.R


################################################################################
# Libraries

library(SpatialFeatureExperiment)
library(alabaster.sfe)
library(tidyverse)
library(Voyager)

################################################################################
# Config
data_dir              <- file.path("~/projects/spatialsnippets/datasets/GSE234713_IBDcosmx_GarridoTrigo2023/processed_data")
sfe_03_banksyniche   <- file.path(data_dir, "GSE234713_CosMx_IBD_sfe_03_banksyniche")


sfe_midi   <- file.path(data_dir, "GSE234713_CosMx_IBD_midi3group03_sfe")


sfe <- readObject(sfe_03_banksyniche)


# Replace cell_ID-like names with cell
# Cell ID is not unique across samples
colnames(sfe) <- sfe$cell

# first 5 fovs of all samples.
sfe <- sfe[,sfe$fov   <= 5]
sfe$tissue_sample <- droplevels(sfe$tissue_sample)

saveObject(sfe, sfe_midi )


