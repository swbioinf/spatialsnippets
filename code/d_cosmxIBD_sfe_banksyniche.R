################################################################################
# Libraries

library(SpatialFeatureExperiment)
library(alabaster.sfe)  # Bioconductor 3.21 and above.
library(tidyverse)
library(Banksy)
library(Voyager) #plotting


################################################################################
# Config
data_dir              <- file.path("~/projects/spatialsnippets/datasets/GSE234713_IBDcosmx_GarridoTrigo2023/processed_data")
sfe_02_celltypeanno   <- file.path(data_dir, "GSE234713_CosMx_IBD_sfe_01_celltypeanno")
sfe_03_banksyniche   <- file.path(data_dir, "GSE234713_CosMx_IBD_sfe_03_banksyniche")

spe_banksy_res         <- file.path(data_dir, "GSE234713_CosMx_IBD_banksy_res.RDS")


# lambda 	A numeric vector in ∈[0,1]∈[0,1] specifying a spatial weighting parameter. Larger values (e.g. 0.8) incorporate more spatial neighborhood and find spatial domains, while smaller values (e.g. 0.2) perform spatial cell-typing.
# 0.2  suggested for 'celltypeing' (not using)
# 0.8for domain finding.
lamdas <- c( 0.6)

# PCs to use (50 seems like alot?)
# ElbowPlot(se) => 20 is plenty enough
npcs <- 20

# number of nearest spatial neighbours:
k_geom <- 20 # 30?

# clustering resolutions
resolutions <- c(0.5, 0.3)

################################################################################
# 1) Process

sfe <- readObject(sfe_02_celltypeanno)

# TEMP subset for testing
#sfe <- sfe[,sfe$group %in% c('HC','CD')]

# try to limit to just hvgs
hvgs <- rownames(rowData(sfe))[rowData(sfe)$hvg]
sfe.hvg <- sfe[hvgs,]
sfe.hvg$tissue_sample<- droplevels(sfe.hvg$tissue_sample)

#split by sample
get_one_sample <- function(the_sample){sfe.hvg[,sfe.hvg$tissue_sample == the_sample]}
sfe.list <- lapply(FUN=get_one_sample, X=levels(sfe.hvg$tissue_sample))

# Get the underlaying banky computation
sfe.list <- lapply(sfe.list, computeBanksy, assay_name = 'logcounts', k_geom = k_geom, parallel=TRUE, num_cores=12)
sfe.merged <- do.call(cbind, sfe.list)

rm(sfe.hvg, sfe.list)
sfe.merged <- Banksy::runBanksyPCA(sfe.merged, lambda = lamdas, npcs = npcs, group='tissue_sample', seed=12)
sfe.merged <- Banksy::clusterBanksy(sfe.merged, lambda = lamdas, npcs = npcs, resolution = resolutions, seed=12)


# Save interium object, in case the order gets messed up again.
print("save sfemerged")
saveRDS(colData(sfe.merged),  spe_banksy_res )
print("merged saved")


################################################################################
# 2) Store it.
# Cell_ID is not unique
# cell is unique

niches_table <- readRDS(spe_banksy_res     )

niche_lookup <- setNames(niches_table$clust_M0_lam0.6_k50_res0.5, nm=as.character( niches_table$cell ))
niche_lookup2 <- setNames(niches_table$clust_M0_lam0.6_k50_res0.3, nm=as.character( niches_table$cell ))


sfe$clust_M0_lam0.6_k50_res0.5 <- niche_lookup[as.character(sfe$cell)]
sfe$clust_M0_lam0.6_k50_res0.3 <- niche_lookup2[as.character(sfe$cell)]

sfe$niche <- factor(paste0('n', as.character(sfe$clust_M0_lam0.6_k50_res0.3)),
                    levels=paste0('n',levels(sfe$clust_M0_lam0.6_k50_res0.3)))

saveObject(sfe, path = sfe_03_banksyniche)
