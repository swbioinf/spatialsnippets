


## Load libraries and data
library(alabaster.sfe)
library(SpatialExperiment)
library(tidyverse)
library(Voyager)
library(scater)
library(scran) # scoremarkers
library(SingleR)
library(schard) #https://github.com/cellgeni/schard

data_dir              <- file.path("~/projects/spatialsnippets/datasets/GSE234713_IBDcosmx_GarridoTrigo2023/processed_data")

# Data object.
sfe_01_loaded         <- file.path(data_dir, "GSE234713_CosMx_IBD_sfe_01_loaded")
sfe_02_celltypeanno   <- file.path(data_dir, "GSE234713_CosMx_IBD_sfe_01_celltypeanno")
#se_pseudobulk_file    <- file.path(data_dir, "GSE234713_CosMx_IBD_pseudobulk_sample_se2.RDS")


# Preprocessed singleR results
predictions_broad_file    <- file.path(data_dir,"singleR_tabula_sapiens_predictions", "predictions_with_ts_intestine_broad_cell_class.RDS")
predictions_detailed_file <-  file.path(data_dir,"singleR_tabula_sapiens_predictions","predictions_with_ts_intestine_detailed_celltype.RDS")

sfe <- readObject(sfe_01_loaded)



ts_largeintestine_h5ad <- file.path("~/projects/spatialsnippets/reference/tabula_sapiens/tabula_sapiens_large_intestine_82e3b450-6704-43de-8036-af1838daa7df.h5ad")
sce.ts_intestine = schard::h5ad2sce(ts_largeintestine_h5ad)




# Not needed, but first filter down to matched genes in our panel
sce.ts_intestine.genename <- sce.ts_intestine[rowData(sce.ts_intestine)$feature_name %in% rownames(sfe),]

# Are there any duplicates (we'd need to handle them, but there aren't)
# takes the count of each feature, then checks that there aren't any >1
stopifnot(sum(table(rowData(sce.ts_intestine.genename)$feature_name) != 1 ) == 0)

# just rename the genes to the gene names
rownames(sce.ts_intestine.genename) <-  rowData(sce.ts_intestine.genename)$feature_name

# Pull out the normalised matrix.
# Quirk of this coming from the python world, the normalised assay is called 'X'
ref_norm_matrix <- assay(sce.ts_intestine.genename, 'X')



#-------------------------------------------------------------------------------
library(BiocParallel) # allow parallelisation with MulticoreParam().

#norm_matrix <- GetAssayData(so, assay = 'RNA', layer = 'data')

predictions_broad <- SingleR::SingleR(test = sfe,
                                      ref   = ref_norm_matrix,
                                      labels = sce.ts_intestine.genename$broad_cell_class,
                                      aggr.ref = TRUE, # builds a pseudobulk reference , speedier processing
                                      BPPARAM = MulticoreParam(workers=12)
)
# Save results to disk
saveRDS(predictions_broad, predictions_broad_file)


predictions_detailed <- SingleR::SingleR(test = sfe,
                                         ref   = ref_norm_matrix,
                                         labels = sce.ts_intestine.genename$cell_type,
                                         aggr.ref = TRUE, # builds a pseudobulk reference , speedier processing
                                         BPPARAM = MulticoreParam(workers=12)
)

# Save results to disk
saveRDS(predictions_detailed, predictions_detailed_file)
