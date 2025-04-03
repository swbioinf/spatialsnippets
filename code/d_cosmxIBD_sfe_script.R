################################################################################
## LIBRARIES

library(tidyverse)
library(patchwork)
library(SpatialFeatureExperiment)
library(alabaster.base)
library(alabaster.sfe) # SaveObject
library(alabaster.bumpy)
library(scran)
library(scater)
library(scuttle)
library(bluster) # clustering parameters, needed for passing htreads to NNGraphParam
library(BiocParallel)
library(data.table) # fread fast file reading
library(BumpyMatrix) # internal represeation of molecules

################################################################################
## CONFIG

dataset_dir      <- '~/projects/spatialsnippets/datasets/'
project_data_dir <- file.path(dataset_dir,'GSE234713_IBDcosmx_GarridoTrigo2023')
sample_dir            <- file.path(project_data_dir, "raw_data_for_sfe/")
annotation_file       <- file.path(project_data_dir,"GSE234713_CosMx_annotation.csv.gz")
data_dir              <- file.path(project_data_dir, "processed_data/")


sfe_01_loaded <- file.path(data_dir, "GSE234713_CosMx_IBD_sfe_01_loaded")

# config
min_count_per_cell <- 50 # previusly 100
min_detected_genes_per_cell <- 20 # from paper.
max_pc_negs        <- 1.5
max_avg_neg        <- 0.5

sample_codes <- c(HC="Healthy controls",UC="Ulcerative colitis",CD="Crohn's disease")

################################################################################
## FUNCTIONS

##5 has annot geom issue
#sfe_data_dir = sample_dirs[5]
#the_sample = sample_names[5]
load_sfe_with_molecules <- function (sfe_data_dir, the_sample) {
  # This function loads SFE without molecules
  # then loads the molecules and adds them in (which avoids an error, maybe ram related)
  # It also checks for weird missing polygons - miss one cell and your cell outlines are just feature outlines.


  # Where add_molecules=TRUE.
  #Error: Capacity error: array cannot contain more than 2147483646 bytes, have 2148177769
  # Pulls cell ids from cell_id
  #"1_1" "2_1" "3_1" "4_1" "5_1" "6_1"
  # "165_349" "166_349" "167_349" "168_349" "169_349" "170_349"
  # <cellnum>_<fov>
  # Updated to sfe 1.9.3 (for alabaster.sfe support) - sample_id is now recorded, no need to edit.
  print("Read SFE without molecules")
  sfe  <- readCosMX(data_dir=sfe_data_dir,
                    sample_id = the_sample,
                    add_molecules = FALSE # True yeilds errors above.
  )

  # the tx file within the data dir
  tx_file <- file.path(sfe_data_dir, paste0(the_sample,"_tx_file.csv"))
  # Read in transcript coordinates, but only keep cellular tx (strips out alot)
  #fov cell_ID      cell x_local_px y_local_px x_global_px y_global_px     z target  CellComp
  #1       1   c_1_1_1       4255        140   10366.964    126795.1     2   Cd74 Cytoplasm
  #1       1   c_1_1_1       4255        140   10366.944    126795.5     7 Ifitm3 Cytoplasm
  #1       1   c_1_1_1       4252        139   10363.784    126796.5     6  Acta2 Cytoplasm
  #1       1   c_1_1_1       4250        150   10362.394    126785.8     3  Ifit1 Cytoplasm
  mol_table <- fread(tx_file,
                     #nrow=100000,
                     select=c('cell_ID', 'fov', 'target', 'x_global_px', 'y_global_px', 'CellComp'))
  #dim(mol_table)
  #It is possible (at least with cyto2 segmenation) to get a 'None' CellComp, and a actuall Cell_ID
  # (those cell calls are dodgy looking though - probably filter them out Later )
  #mol_table <- mol_table[CellComp != 'None', ] # Assigned to cells only
  # So insead filter the no cellID ones!
  print("Reading molecules")
  mol_table <- mol_table[cell_ID != 0, ]
  mol_table <- mol_table[,cell_id:=paste0(cell_ID,"_",fov)]
  mol_table <- mol_table[,!c('CellComp', "fov","cell_ID")]
  #dim(mol_table)


  # construct BumpyMatrix
  print("formatting molecules")
  mol_bm <- splitAsBumpyMatrix(
    mol_table[, c("x_global_px", "y_global_px")],
    row = mol_table$target, col = mol_table$cell_id,
    sparse =TRUE) # sparse=TRUE is important, and not default !! 10fold less ram.
  rm(mol_table)




  # there may be more 'cells' from the tx file than in the object, but these
  # are quietly ignored. But everything in sft should be in mol_bm!
  stopifnot(all(rownames(sfe) %in% rownames(mol_bm) ))

  #stopifnot(all( colnames(sfe) %in% colnames(mol_bm) ))
  # Cells with no tx transcripts? - Possible!
  # Silently drop them.
  #missing_cells <- colnames(sfe)[!colnames(sfe) %in% colnames(mol_bm) ]
  #colData(sfe)[missing_cells,] %>% as.data.frame() %>% View()
  #table(colSums(counts(sfe)[, missing_cells]))
  #genuinely, nope, they have no transcirpts.
  common_cells <- intersect(colnames(mol_bm),colnames(sfe))
  sfe <- sfe[,common_cells]
  mol_bm <- mol_bm[rownames(sfe),common_cells] # match order and drop cells not in sfe, ditto genes

  # now save the molecules in the object
  print("Storing molecules")
  assay(sfe, 'molecules') <- mol_bm


  # Where there is missing polygon data for sfe,
  # remove the corresponding cell
  # Otherwise you end up with spatial annotations instead of segmentations,
  # which means files cant be joined (and probably other issues!)
  # Affects one (1!) cell in a dataset - I do not care to see why.
  # Should do nothing otherwise.
  print("Checking for polygonless cells")
  sfe<-check_and_rm_cells_without_polygons(sfe)


  # Move negatives into altExpr
  probes <- rownames(rowData(sfe))
  rowData(sfe)$target    <- probes
  rowData(sfe)$CodeClass <- factor(ifelse(grepl("^SystemControl",probes), "FalseCode",
                                          ifelse(grepl("^NegPrb",probes), "NegProbe","RNA" )),
                                   levels=c("RNA","NegProbe","FalseCode"))
  rowData(sfe)$CodeClass <- droplevels(rowData(sfe)$CodeClass )
  # An attempt at storing negative probes in AltExp, but is causing issues with saveing.
  # unsure if sfe altexp is supported?
  # leaving in main assay for now.
  # No Falsecodes, only negatives
  #sfe.neg <- sfe[rowData(sfe)$CodeClass != "RNA",] # you put your neg probes in
  #sfe     <- sfe[rowData(sfe)$CodeClass == "RNA",] # you take your neg probes out
  #altExp(sfe,"NegPrb") <- sfe.neg                  # you put your neg probes in
  ## and you shake the data out

  return(sfe)
}






load_one_sample_as_sfe <- function(sfe_data_dir, the_sample) {

  print(the_sample)
  sfe <- load_sfe_with_molecules(sfe_data_dir, the_sample)

  # sample code nice and short, and may be inferred
  sample_code <- substr(the_sample,12,16)
  sfe$sample_id  <- the_sample
  sfe$slide_name <- the_sample

  # and some experimental info
  sfe$individual_code <- substr(the_sample,12,16)
  sfe$tissue_sample   <- substr(the_sample,12,16)
  sfe$group           <- factor(substr(the_sample, 12,13), levels=names(sample_codes))
  sfe$condition       <- factor(as.character(sample_codes[sfe$group]), levels=sample_codes)
  sfe$fov_name        <- paste0(sfe$individual_code,"_", str_pad(sfe$fov, 3, 'left',pad='0'))

  # cell labels: need to match that of the annotation file:
  #                id   subset             SingleR2
  # HC_c_2_1 HC_c_2_1   stroma            Pericytes
  # HC_c_3_1 HC_c_3_1   stroma          Endothelium
  # HC_c_4_1 HC_c_4_1   stroma            Pericytes


  #Have this:

  #      fov   cell_ID      Area AspectRatio CenterX_local_px
  #1         1         1      6153        0.67             2119
  #2         1        10      5971        0.63             1265
  #3         1       100      1785        1.11             2954
  #4         1      1000      3946        0.59             2959

  #     individual_code tissue_sample    group        condition    fov_name
  #                HC_a          HC_a       HC Healthy controls    HC_a_001
  #                HC_a          HC_a       HC Healthy controls    HC_a_001
  #                HC_a          HC_a       HC Healthy controls    HC_a_001

  #> table(table(sfe$cell_ID))
  #   1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19
  # 699  110   56   78   79  110  100  165   21    4   24   38   24    5   66  265   36  190 1451


  # HC_c_4_1
  # <>>sample code>_<cell_ID>
  # Where cell_ID includes <number>_<fovnum>: "3518_20"
  sfe$cell <- paste0(sfe$tissue_sample, "_",sfe$cell_ID)

  sfe$total_count <- colSums(counts(sfe))
  sfe$distinct_genes <- colSums(counts(sfe)!=0)
  sfe$total_count_log10 <- log10(sfe$total_count)


  #Negative probes

  sfe$neg_count   <- colSums(counts(sfe[rowData(sfe)$CodeClass == "NegProbe", ]))
  sfe$avg_neg     <- colMeans(counts(sfe[rowData(sfe)$CodeClass == "NegProbe", ]))
  sfe$pc_neg<- sfe$neg_count / (sfe$neg_count + sfe$total_count) * 100


  return(sfe)
}



check_and_rm_cells_without_polygons <-function(sfe) {
  # Theory:
  # A cell polygon is empty ()
  # Yeilding folling error during parquet construction:
  #">>> Constructing cell polygons
  # >>> ..removing 1 empty polygons"
  #
  # This means that theres one more 'cell' than 'cell segmentation genometriy'
  #
  # cellSeg() function will look for a match in those counts and store
  # geometries in
  # * colGeometries if cell counts match
  # * annotGeometires if they dont.
  #
  # Consequently, you can't join / cbind sfe object unless their colgeometries match.
  #
  #    Error in value[[3L]](cond) :
  #    failed to combine 'int_colData' in 'cbind(<SpatialFeatureExperiment>)':
  #    failed to rbind column 'colGeometries' across DataFrame objects: the DFrame objects to combine must have
  #    the same column names
  #
  #
  # Seems to be rare, not ever datasset, and can be just one cell..
  #
  # This function will turn annotGeom back into colGeom, buy remooving the missing cells
  #
  if ('cellSeg' %in% names(annotGeometries(sfe))) {

    print('Found cellSeg in annotGeometries instead of colGeometries')

    geomet_cells <- rownames(annotGeometries(sfe)$cellSeg)
    sfe_cells    <- colnames(sfe)
    common_cells <- intersect(sfe_cells, geomet_cells)

    # Can find cells with geometrey, but not in object.
    # possibly the ones with no counts?
    #stopifnot(all(geomet_cells %in% sfe_cells))
    geomet_cells_without_sfe <- setdiff( geomet_cells, sfe_cells)
    print(paste("dropping",length(geomet_cells_without_sfe), "cells from cell seg geometry (absent in count matrix)"))


    # Warn on number of cells dropped
    cells_without_geomet <- setdiff(sfe_cells, geomet_cells)
    print(paste("dropping",length(cells_without_geomet), "cells from counts matrix (no segmentation outline)"))

    # Drop cells
    sfe <- sfe[,common_cells]
    #annotGeometries(sfe)$cellSeg <- annotGeometries(sfe)$cellSeg[common_cells,]
    #annogeo <- annotGeometries(sfe)$cellSeg
    #annogeo <- annotGeometries(sfe)$cellSeg[common_cells,]

    # Put it back in using the accessor, all proper like
    cellSeg(sfe) <- annotGeometries(sfe)$cellSeg[common_cells,]

    # Remove and
    # for consistantcy, unname the empty list
    annotGeometries(sfe)$cellSeg <- NULL
    annotGeometries(sfe) <- unname(annotGeometries(sfe))

    # Dropping the sampleID, which got added for some reason.
    keep_cols <- colnames(cellSeg(sfe))
    keep_cols <- keep_cols[keep_cols != 'sample_id']
    cellSeg(sfe) <- cellSeg(sfe)[,keep_cols]

  }
  return(sfe)
}




do_basic_preprocessing <- function(sfe,
                                   hvg_prop = 0.3,
                                   num_pcs  = 20,
                                   sampleblock = NULL,
                                   BPPARAM  = MulticoreParam(16) # and everything else others
) {



  # Normalization.
  print('Normalise...')
  sfe <- logNormCounts(sfe, BPPARAM=BPPARAM)

  # Highly variable genes, but ignoring individual level data.
  # Attempts to model technival vs bio variation
  # Choosing 30% of genes, since we've only got 1000, and they are selected for cell type discernment.
  # 30% is arbitrary!
  print("Model gene variance...")

  if (is.null(sampleblock)) {
    gene_variances <- modelGeneVar(sfe, BPPARAM=BPPARAM)
  } else {
    to_block <- as.factor(colData(sfe)[,sampleblock])
    gene_variances <- modelGeneVar(sfe, BPPARAM=BPPARAM, block=to_block)
  }

  print("get top hvg...")
  # only consider RNA probes
  rna_probes             <- rowData(sfe)$target[rowData(sfe)$CodeClass == "RNA"]
  gene_vairances.RNAonly <- gene_variances[rna_probes, ]
  hvg <- getTopHVGs(gene_vairances.RNAonly, prop=hvg_prop)
  rowData(sfe)$hvg <- rowData(sfe)$target %in% hvg

  print("run pca...")
  sfe <- fixedPCA(sfe, subset.row = hvg)
  print("run UMAP...")
  sfe <- runUMAP(sfe, pca=num_pcs, BPPARAM=BPPARAM )

  print("cluster...")
  set.seed(12) # set seed for consistant clustering.
  #sfe$nn.cluster <- clusterCells(sfe, use.dimred="PCA", BLUSPARAM = bluster::NNGraphParam(num.threads=16))
  # try snn
  ## Haven't seen evidence of this using multiple thresads.
  # Split the graph building and clustering steps, as gettting some untraceable issues.
  g <- buildSNNGraph(sfe, k=20,  use.dimred = 'PCA', BPPARAM=BPPARAM)
  print('built SNN graph')
  saveRDS(g, file.path("~/snn_graph_object.rds"))

  sfe$snn.cluster  <- igraph::cluster_louvain(g)$membership

  print('clustered SNN graph')


  sfe$snn.cluster <- as.factor(sfe$snn.cluster)
  sfe$cluster_code <- factor(paste0("c",sfe$snn.cluster), levels=paste0("c",levels(sfe$snn.cluster)))

  print("Done")
  return(sfe)

}

################################################################################
## Go

#Load an annotate all samples in the directory full of samples (each as flatfiles)

if ( FALSE ) { # ALREADY RUN

  # Get list of samples
  sample_names <- list.files(sample_dir)
  sample_dirs  <- file.path(sample_dir, sample_names)

  # Use multiple apply to run each with corresponding path.
  sfe_list <- mapply(FUN=load_one_sample_as_sfe ,
                     sfe_data_dir=sample_dirs,
                     the_sample= sample_names)

  # This is ineffient, but it works - cbind can join pairs of sfe object.
  sfe <- do.call(cbind, sfe_list)

  print("Merged. Now apply annotation.")


  # Add the cell annotation from the paper.
  anno_table <- read_csv(annotation_file)
  anno_table <- as.data.frame(anno_table)
  rownames(anno_table) <- anno_table$id

  head(colData(sfe))
  head(anno_table)


  sfe$celltype_subset   <- factor(anno_table[sfe$cell,]$subset)
  sfe$celltype_SingleR2 <- factor(anno_table[sfe$cell,]$SingleR2)

  # and foactorise a few things now we've got the full table
  sfe$fov_name          <- factor(sfe$fov_name)
  sfe$individual_code   <- factor(sfe$individual_code)
  sfe$tissue_sample     <- factor(sfe$tissue_sample)

  # There are a few unannotted, will remove
  table(is.na(sfe$celltype_subset))

  ncol(sfe)
  sfe <- sfe[ ,sfe$distinct_genes >= min_detected_genes_per_cell &
                sfe$avg_neg <= max_avg_neg &
                !(is.na(sfe$celltype_subset) )]
  ncol(sfe)


}
## Preprocessing
#hvg_prop = 0.3
#num_pcs  = 20
#sampleblock = NULL
#num_threads = 16
#BPPARAM  = MulticoreParam(num_threads)
#saveObject(sfe, "~/testX2")
sfe <- readObject("~/testX2")


#sfe <- sfe[,1:400]
sfe <- do_basic_preprocessing(sfe, num_pcs = 15)

## Save
saveObject(sfe, sfe_01_loaded)





