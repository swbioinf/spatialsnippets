
library(OSTA.data)
library(DropletUtils)
library(SpatialExperimentIO)
library(Voyager)
library(BiocParallel)
library(HDF5Array)
library(alabaster.sfe)
library(scran)
library(scater)
library(igraph)
library(scrapper)

out_dir  <- file.path("~/projects/spatialsnippets/datasets/OSTA.data/")
saved.sfe_path <- file.path(out_dir, "CosMx1k_MouseBrain1_loadedwithclusters")
num_pcs=15
BPPARAM  = MulticoreParam(12)

#https://new.bioconductor.org/packages/release/bioc/vignettes/OSTA.data/inst/doc/OSTA.data.html


# retrieval
id <- "CosMx1k_MouseBrain1"
pa <- OSTA.data_load(id, mol=FALSE)
dir.create(td <- tempfile())
unzip(pa, exdir=td)
# importing
(spe <- readCosmxSXE(td, addTx=FALSE))
#
#
# library(SpaceTrooper)
# # prepare data for 'SpaceTrooper'
# cos <- updateCosmxSPE(spe, td) #, sampleName="sample_01")
# cos <- readAndAddPolygonsToSPE(cos)
# cos$in_tissue <- TRUE
# cos
#
#

sfe <- toSpatialFeatureExperiment(spe)
# Centroids only, haven't figured out how to convert the polygons over.
#plotSpatialFeature(sfe, "sample_id",colGeometryName = "centroids") +
#  theme(legend.position = "none")

print('Normalise...')
#spe <- normalizeRnaCounts (sfe)

# rmeove zero/tiny counts
min_count <- 20
sfe <- sfe[,colSums(counts(sfe))>=min_count]

# Tnen norm.
sfe <- logNormCounts(sfe, BPPARAM=BPPARAM)



gene_variances <- modelGeneVar(sfe, BPPARAM=BPPARAM)
hvg <- getTopHVGs(gene_variances , prop=0.3)
rowData(sfe)$hvg <- rownames(sfe) %in% hvg

print("run pca...")
sfe <- fixedPCA(sfe, subset.row = hvg)
print("run UMAP...")
sfe <- runUMAP(sfe, pca=num_pcs, BPPARAM=BPPARAM )

print("cluster...")
set.seed(12) # set seed for consistant clustering.

g <- buildSNNGraph(sfe, k=20,  use.dimred = 'PCA', BPPARAM=BPPARAM)
print('built SNN graph')
saveRDS(g, file.path("~/snn_graph_object.rds"))
sfe$snn.cluster  <- igraph::cluster_louvain(g)$membership
print('clustered SNN graph')

sfe$snn.cluster <- as.factor(sfe$snn.cluster)
sfe$cluster_code <- factor(paste0("c",sfe$snn.cluster), levels=paste0("c",levels(sfe$snn.cluster)))

saveObject(sfe, saved.sfe_path )

