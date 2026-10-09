library(Seurat)
#library(SeuratObject)
library(tidyverse)
library(progressr) # used within LoadNanostring functions
library(data.table)
# revision of loading the data to use LoadNanostring function from Seurat now
# have to load metadata manually after, but that's easier than a whole other function!
#
# This should make data more consistant with typical case (ie. don't put negatives in their own assay anymore.)
# Will require some changes to downstream analysis, but worth it for consistancy.
#
# Neeed to have older sp version atm: https://github.com/satijalab/seurat/issues/10527



dataset_dir      <- '~/projects/spatialsnippets/datasets/'
project_data_dir <- file.path(dataset_dir,'GSE234713_IBDcosmx_GarridoTrigo2023')
sample_dir            <- file.path(project_data_dir, "raw_data_for_sfe/")
annotation_file       <- file.path(project_data_dir,"GSE234713_CosMx_annotation.csv.gz")
data_dir              <- file.path(project_data_dir, "processed_data/")
tempdir               <- '~/temp/' # working on a server with a very small default tmp, which causes silent incomplete fread data.table reads on very large data (e.g. tx)


seurat_file_00_raw    <- file.path(data_dir, "GSE234713_CosMx_IBD_seurat_00_raw_v2.RDS")
seurat_file_01_loaded <- file.path(data_dir, "GSE234713_CosMx_IBD_seurat_01_loaded_v2.RDS")

# config
min_count_per_cell <- 100
max_pc_negs        <- 1.5
max_avg_neg        <- 0.5

sample_codes <- c(HC="Healthy controls",UC="Ulcerative colitis",CD="Crohn's disease")



#the_sample <- 'GSM7473682_HC_a'

load_sample_into_seurat <- function(the_sample){

  so <- LoadNanostring(file.path(sample_dir, the_sample),
                       fov = the_sample,
                       assay='RNA')

  # sample info
  # This is easy beecause there is only one sample per slide.
  so$individual_code <- factor(substr(the_sample,12,16))
  so$tissue_sample   <- factor(substr(the_sample,12,16))
  so$group     <- factor(substr(the_sample, 12, 13), levels=names(sample_codes))
  so$condition <- factor(as.character(sample_codes[so$group]), levels=sample_codes)


  #Seurat ignores the metadata (or can only load certain (not all) hardcoded columns if specified)
  # So manually add all of the metadata!
  metadata.file = file.path(sample_dir,the_sample,  paste0(the_sample,'_metadata_file.csv'))
  so <- add_metadata_file_content_to_seurat_metadata(so, metadata.file)


  # Add fov name tool (fov came with metadata)
  so$fov_name        <- paste0(so$individual_code,"_", str_pad(so$fov, 3, 'left',pad='0'))

  # Previously put neg probes into their own assay, but going to stop doing this, as it doesn't seem standard.
  # Will,however, populate the probe annotation
  # only have negative probes in ths one,not system controls,
  neg_probes <- rownames(so)[grepl(x=rownames(so), pattern="NegPrb")]
  so@assays[['RNA']]@meta.data$target <- rownames(so@assays[['RNA']])
  so@assays[['RNA']]@meta.data$CodeClass <- ifelse(so@assays[['RNA']]@meta.data$target %in% neg_probes , "Negative","RNA")


  neg_matrix          <- GetAssayData(so,assay = 'RNA', layer = 'counts')[neg_probes,]
  so$nCount_negprobes <-   colSums(neg_matrix)
  so$avg_neg          <-  colMeans(neg_matrix)

  return(so)

}






add_metadata_file_content_to_seurat_metadata <- function(so, metadata.file) {
  # Pull in full metadata

  md <- data.table::fread(metadata.file,  sep = ',',  data.table = FALSE,  verbose = FALSE, tmpdir  = tempdir)

  # Seurat cell identifiers constructed within ReadNanostring (matrix rownames) as:
  # https://github.com/satijalab/seurat/blob/5f67c701ecc6738e9e47dc21dbd202d335df5ec3/R/preprocessing.R#L2251
  # Combination of Cell ID (for non-zero cell_IDs) and FOV are assumed to be unique. Used to create barcodes / rownames.
  #bcs <- paste0(as.character(tx$cell_ID), "_", tx$fov)
  # Checking
  ## expect a small number of 'cells' to have been filtered out.
  ##table(md$so_cell_id %in% colnames(so))
  # ANd everything from seurat should be in metadata
  #table(colnames(so) %in% md$so_cell_id)
  # MD is not a data.table here (and its only ever going to be one slide so that's fine)!
  md$so_cell_id <- paste0(as.character(md$cell_ID), "_", md$fov)
  coords <- match(rownames(so@meta.data), md$so_cell_id)
  md <- md[coords,]
  stopifnot(all(md$so_cell_id == colnames(so)))
  md$so_cell_id <- NULL # remove col, don't need it.

  so <- AddMetaData(so, md)

}





samples <- c('GSM7473682_HC_a','GSM7473683_HC_b','GSM7473684_HC_c',
             'GSM7473685_UC_a','GSM7473686_UC_b','GSM7473687_UC_c',
             'GSM7473688_CD_a','GSM7473689_CD_b','GSM7473690_CD_c')
sample_prefix <- paste0(substr(samples, 12,15))

# Allow skipping of this step
LOAD_RAW=TRUE
if (LOAD_RAW) {

so.list <- lapply(FUN=load_sample_into_seurat, X=samples)

#NB: merge is in SeuratObject packages, but must be called without ::
options(future.globals.maxSize= 10*1024^3) # 10G.
so.raw <- merge(so.list[[1]], y=so.list[2:length(so.list)], add.cell.ids=sample_prefix)
#Error in getGlobalsAndPackages(expr, envir = envir, globals = globals) :
#  The total size of the 8 globals exported for future expression (‘FUN()’) is 5.90 GiB.. This exceeds the maximum #allowed size of 500.00 MiB (option 'future.globals.maxSize'). The three largest globals are ‘FUN’ (5.90 GiB of #class ‘function’), ‘p’ (25.22 KiB of class ‘function’) and ‘slot<-’ (1.36 KiB of class ‘function’)
rm(so.list)

# save
saveRDS(so.raw, seurat_file_00_raw)

}





# load previously saved.
so.raw <- readRDS(seurat_file_00_raw)

# Negative probe handling
so.raw$pc_neg <-  ( so.raw$nCount_negprobes / so.raw$nCount_RNA ) * 100

#so.raw[["negprobes"]] <- JoinLayers(so.raw[["negprobes"]]) # For caluclating these, need to have the negprobes merged
#so.raw$avg_neg <-  colMeans(so.raw[["negprobes"]])   # only defined firsts sample.



# Pull in annotation

anno_table <- read_csv(annotation_file)

anno_table <- as.data.frame(anno_table)
rownames(anno_table) <- anno_table$id

head(so.raw@meta.data)
head(anno_table)

so.raw$full_cell_id      <- as.character(rownames(so.raw@meta.data))
so.raw$celltype_subset   <- factor(anno_table[so.raw$full_cell_id,]$subset)
so.raw$celltype_SingleR2 <- factor(anno_table[so.raw$full_cell_id,]$SingleR2)
so.raw$fov_name          <- factor(so.raw$fov_name)
so.raw$group          <- factor(so.raw$group, levels=c("CD","UC","HC"))
so.raw$condition      <- factor(so.raw$condition, levels=c("Crohn's disease",  'Ulcerative colitis', 'Healthy controls'))
table(is.na(so.raw$celltype_subset))




so <- so.raw[ ,so.raw$nCount_RNA >= min_count_per_cell &
                so.raw$avg_neg <= max_avg_neg &
                !(is.na(so.raw$celltype_subset) )]






