library(Seurat)
#library(SeuratObject)
library(tidyverse)
library(progressr) # used within LoadNanostring functions


dataset_dir      <- '~/projects/spatialsnippets/datasets/'
project_data_dir <- file.path(dataset_dir,'GSE234713_IBDcosmx_GarridoTrigo2023')
sample_dir            <- file.path(project_data_dir, "raw_data_for_sfe/")
annotation_file       <- file.path(project_data_dir,"GSE234713_CosMx_annotation.csv.gz")
data_dir              <- file.path(project_data_dir, "processed_data/")


seurat_file_00_raw    <- file.path(data_dir, "GSE234713_CosMx_IBD_seurat_00_raw_v2.RDS")
seurat_file_01_loaded <- file.path(data_dir, "GSE234713_CosMx_IBD_seurat_01_loaded_v2.RDS")

# config
min_count_per_cell <- 100
max_pc_negs        <- 1.5
max_avg_neg        <- 0.5

sample_codes <- c(HC="Healthy controls",UC="Ulcerative colitis",CD="Crohn's disease")


# Load LoadNanostring.X for keeping metadata e.t.c
# As of September 2026, readnanostring explicitly names feilds that can be loaded from metadata (and there are more than named there)
# other aspesct of this function may be outdated
source('code/LoadNanostring.R')



load_sample_into_seurat <- function(the_sample){

  # Using a modified version of the LoadNanostring function, which keeps metadata, and loads a litte more effiently for large datasets (which is not this one.)
  # ..../datasets//GSE234713_IBDcosmx_GarridoTrigo2023/raw_data_for_sfe//GSM7473682_HC_a
  so <- LoadNanostring.X(file.path(sample_dir, the_sample),
                         project = the_sample,
                         assay='RNA',
                         fov=the_sample, tempdir= '~/tmp')

  # sample info
  so$individual_code <- factor(substr(so$orig.ident,12,16))
  so$tissue_sample   <- factor(substr(so$orig.ident,12,16))
  so$group     <- factor(substr(the_sample, 12, 13), levels=names(sample_codes))
  so$condition <- factor(as.character(sample_codes[so$group]), levels=sample_codes)

  so$fov_name        <- paste0(so$individual_code,"_", str_pad(so$fov, 3, 'left',pad='0'))


  # Put neg probes into their own assay.
  neg_probes <- rownames(so)[grepl(x=rownames(so), pattern="NegPrb")]
  neg_matrix         <- GetAssayData(so,assay = 'RNA', layer = 'counts')[neg_probes,]
  #so[["negprobes"]] <- CreateAssayObject(counts = neg_matrix)
  so[["negprobes"]] <- CreateAssay5Object(counts = neg_matrix)

  ## and remove from the main one
  rna_probes  <- rownames(so)[(! rownames(so) %in% neg_probes)]
  so[['RNA']] <- subset( so[['RNA']],features = rna_probes)


  return(so)

}

samples <- c('GSM7473682_HC_a','GSM7473683_HC_b','GSM7473684_HC_c',
             'GSM7473685_UC_a','GSM7473686_UC_b','GSM7473687_UC_c',
             'GSM7473688_CD_a','GSM7473689_CD_b','GSM7473690_CD_c')[1:2]
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

so.raw[["negprobes"]] <- JoinLayers(so.raw[["negprobes"]]) # For caluclating these, need to have the negprobes merged
so.raw$avg_neg <-  colMeans(so.raw[["negprobes"]])   # only defined firsts sample.



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






