#library(progressr)



#' @name nanostring-helpers
#' @rdname nanostring-helpers
#'
#' @return data frame containing counts for cells based on a single class of segmentation (eg Nuclear)
#'
#' @keywords internal
#'
#' @noRd
#'
#'

#library(magrittr)
#library(dplyr)
#library(Matrix)
build.cellcomp.matrix.X <- function(mols.df, class=NULL) {
  if (!is.null(class)) {
    if (!(class %in% c("Nuclear", "Membrane", "Cytoplasm"))) {
      stop(paste("Cannot subset matrix based on segmentation:", class))
    }
    mols.df <- mols.df[mols.df$CellComp == class,]  # subset based on cell class
  }
  
  # Remove counts outside of cells
  mols.df <- mols.df[CellComp != "None", ]
  
  # This is a big table
  # Fread uses datatable, keep using the more effienct datatable functions
  mols.df[, `:=`(bc= paste0(cell_ID, "_", mols.df$fov) )  ] 
  
  print("make counts:")
  counts <- mols.df[, .(n = .N), by = .(target, bc)]
  targets   <- unique(counts$target)
  barcodes  <- unique(counts$bc)
  
  print("make matrix:")
  # Convert cell, target and count into a sparse count matrix, using indicies
  mtx <- Matrix::sparseMatrix(
    i = match(counts$bc, barcodes),
    j = match(counts$target, targets),
    x = counts$n,
    dimnames = list(barcodes, targets)
  )
  

  print ("made count matrix")
  return(mtx)
}









#' #### TODO: email nanostring and ask for .mtx!!!!
#' 
#' #' Load a nanostring flatfile .mtx formatted file
#' #'
#' #' The SLIDE_exprMat_file.csv.gz files exported from aotmx 
#' #' have format <fov><cell_id><gene ...> , Stored in a dense format.
#' #' This funtion will read that in chunks and build a sparse matrix.
#' #'
#' #'
#' #' @param file exprsMat file from nanostring flat file format
#' #' @param chunk_size How many lines to read at once? Larger values need more RAM,
#' #' but may run quicker.
#' #' 
#' #' @importFrom readr read_csv_chunked
#' #'
#' ReadExprMatFile <- function(file, chunk_size=20000) {
#'   
#' 
#'   ## Check headers.
#'   
#'   ## If these aren't the first two columns, we're working with a different format.
#'   first_cols     <- c('fov','cell_ID')
#'   old_first_cols <- c('fov','cell_ID','cell')
#'   start_col <- 3 
#'   headers <- colnames(read.delim(file, sep =",", nrows=1))
#'   if (! all(headers[1:2] == first_cols)){
#'     stop(paste0("Expected first two headers in mtx file to be",
#'          paste(first_cols, collapse=','),
#'          ", instead saw: ", paste(headers[1:10], collapse=,',' )))
#'   }
#'   if (all(headers[1:3] == old_first_cols)){
#'     warning("Detected old version of mtx file (flat file export < v1.1.1). Ignroing cell column")
#'     start_col <- 4
#'   }
#'   
#'   
#'   # Local funtion to reformat a small dense table
#'   # into a sparse matrix with rownames.
#'   make_chunk_sparse <- function(chunk, pos){
#' 
#'     # Remove any counts with a cell ID of 0
#'     chunk <- chunk[!chunk$cell_ID == 0,]
#'     
#'     # Expect format to begin fov, cell_ID
#'     # fov cell_ID  A1BG   A2M  AAAS  AAK1 
#'     # Or more rarely, old format (ignore cell)
#'     # fov cell_ID cell A1BG   A2M  AAAS  AAK1 
#'     counts  <- as.sparse(chunk[start_col:ncol(chunk)])
#'     
#'     # Combination of Cell ID (for non-zero cell_IDs) 
#'     # and FOV are assumed to be unique. Used to create barcodes / rownames.
#'     bcs <- paste(chunk$cell_ID, chunk$fov, sep = "_")
#'     rownames(counts) <- bcs
#'     
#'     return(counts)
#'     
#'   }
#'   
#'   print(paste0(c("XXXX ", chunk_size)))
#'   
#'   # Read the file in chunks
#'   tx <- read_csv_chunked(
#'     file=file, 
#'     col_names = TRUE, 
#'     chunk_size = chunk_size, 
#'     show_col_types = FALSE,
#'     callback = ListCallback$new(make_chunk_sparse)
#'   )
#' 
#'   # do.call is not particularly effiecient, but it works.
#'   # Would be better to use a bind_rows-like function that works on 
#'   # a Sparse Matrix, if one is available.
#'   tx <- do.call(what = rbind, tx) 
#'   return(tx)
#' }
#' 



#' Read and Load Nanostring SMI data
#'
#' @param data.dir Directory containing all Nanostring SMI files with
#' default filenames
#' @param mtx.file Path to Nanostring cell x gene matrix CSV
#' @param metadata.file Contains metadata including cell center, area,
#' and stain intensities
#' @param molecules.file Path to molecules file
#' @param segmentations.file Path to segmentations CSV
#' @param type Type of cell spatial coordinate matrices to read; choose one
#' or more of:
#' \itemize{
#'  \item \dQuote{centroids}: cell centroids in pixel coordinate space
#'  \item \dQuote{segmentations}: cell segmentations in pixel coordinate space
#' }
#' @param mol.type Type of molecule spatial coordinate matrices to read;
#' choose one or more of:
#' \itemize{
#'  \item \dQuote{pixels}: molecule coordinates in pixel space
#' }
#' @param metadata List of metadata columns to include. Omit to read all 
#' availble metadata.
#' @param mols.filter Filter molecules that match provided string
#' @param genes.filter Filter genes from cell x gene matrix that match
#' provided string
#' @param fov.filter Only load in select FOVs. Nanostring SMI data contains
#' 30 total FOVs.
#' @param subset.counts.matrix If the counts matrix should be built from
#' molecule coordinates for a specific segmentation; One of:
#' \itemize{
#'  \item \dQuote{Nuclear}: nuclear segmentations
#'  \item \dQuote{Cytoplasm}: cell cytoplasm segmentations
#'  \item \dQuote{Membrane}: cell membrane segmentations
#' }
#' @param cell.mols.only If TRUE, only load molecules within a cell
#'
#' @return \code{ReadNanostring}: A list with some combination of the
#' following values:
#' \itemize{
#'  \item \dQuote{\code{matrix}}: a
#'  \link[Matrix:dgCMatrix-class]{sparse matrix} with expression data; cells
#'   are columns and features are rows
#'  \item \dQuote{\code{centroids}}: a data frame with cell centroid
#'   coordinates in three columns: \dQuote{x}, \dQuote{y}, and \dQuote{cell}
#'  \item \dQuote{\code{pixels}}: a data frame with molecule pixel coordinates
#'   in three columns: \dQuote{x}, \dQuote{y}, and \dQuote{gene}
#' }
#'
#' @importFrom future.apply future_lapply
#'
#' @export
#'
#' @order 1
#'
#' @concept preprocessing
#'
#' @template section-progressr
#' @template section-future
#'
#' @templateVar pkg data.table
#' @template note-reqdpkg
#'
#'



# mtx.file = NULL
# metadata.file = NULL
# molecules.file = NULL
# segmentations.file = NULL
# type = 'centroids'
# mol.type = 'pixels'
# metadata = NULL
# mols.filter = NA_character_
# genes.filter = NA_character_
# fov.filter = NULL
# subset.counts.matrix = NULL
# cell.mols.only = TRUE


ReadNanostring.X <- function(
    data.dir,
    mtx.file = NULL,
    metadata.file = NULL,
    molecules.file = NULL,
    segmentations.file = NULL,
    type = 'centroids',
    mol.type = 'pixels',
    metadata = NULL,
    mols.filter = NA_character_,
    genes.filter = NA_character_,
    fov.filter = NULL,
    subset.counts.matrix = NULL,
    cell.mols.only = TRUE,
    tempdir = NULL ####NEW
) {
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Please install 'data.table' for this function")
  }
  
  # Argument checking
  type <- match.arg(
    arg = type,
    choices = c('centroids', 'segmentations'),
    several.ok = TRUE
  )
  mol.type <- match.arg(
    arg = mol.type,
    choices = c('pixels'),
    several.ok = TRUE
  )
  
  use.dir <- all(vapply(
    X = c(mtx.file, metadata.file, molecules.file),
    FUN = function(x) {
      return(is.null(x = x) || is.na(x = x))
    },
    FUN.VALUE = logical(length = 1L)
  ))
  
  if (use.dir && !dir.exists(paths = data.dir)) {
    stop("Cannot find Nanostring directory ", data.dir)
  }
  
  ####--------------------------------------------------------------------------
  # Identify input files
  files <- c(
    #### matrix = mtx.file %||% '[_a-zA-Z0-9]*_exprMat_file.csv',
    metadata.file = metadata.file %||% '[_a-zA-Z0-9]*_metadata_file.csv',
    molecules.file = molecules.file %||% '[_a-zA-Z0-9]*_tx_file.csv',
    segmentations.file = segmentations.file %||% '[_a-zA-Z0-9]*-polygons.csv'
  )
  
  # use custom tempdir if defined.
  tempdir = tempdir %||% tempdir()
  
  
  files <- vapply(
    X = files,
    FUN = function(x) {
      x <- as.character(x = x)
      if (isTRUE(x = dirname(path = x) == '.')) {
        fnames <- list.files(
          path = data.dir,
          pattern = x,
          recursive = FALSE,
          full.names = TRUE
        )
        return(sort(x = fnames, decreasing = TRUE)[1L])
      } else {
        return(x)
      }
    },
    FUN.VALUE = character(length = 1L),
    USE.NAMES = TRUE
  )
  files[!file.exists(files)] <- NA_character_
  
  if (all(is.na(x = files))) {
    stop("Cannot find Nanostring input files in ", data.dir)
  }
  
  ####--------------------------------------------------------------------------
  ## Read Data from files
  
  # Metadata
  #### Set metadata feilds to all metadata, or allow a subset
  #### If metadata is null, read everything.
  #### Change in behaviour, this would previously read nothing.
  names <-colnames(read.delim(files[['metadata.file']], sep =",", nrows=1))

  if (is.null(metadata)) {
    metadata <- names
  } else {
    metadata <- match.arg(
      arg = metadata,
      choices = names,
      several.ok = TRUE  
    )
  }
    
  # Checking for loading spatial coordinates
  if (!is.na(x = files[['metadata.file']])) {
    pprecoord <- progressor()
    pprecoord(
      message = "Preloading cell spatial coordinates",
      class = 'sticky',
      amount = 0
    )
    md <- data.table::fread(
      file = files[['metadata.file']],
      sep = ',',
      data.table = FALSE,
      verbose = FALSE,
      tmpdir  = tempdir
    )
    #md <- data.table(read.csv(
    #  file = files[['metadata.file']],
    #  stringsAsFactors = FALSE
    #))
    
    print(paste0("MD read: ", nrow(md)))
    
    # filter metadata file by FOVs
    if (!is.null(x = fov.filter)) {
      md <- md[md$fov %in% fov.filter,]
    }
    pprecoord(type = 'finish')
  }
  
  
  ## Segmentations
  if (!is.na(x = files[['segmentations.file']])) {
    ppresegs <- progressor()
    ppresegs(
      message = "Preloading cell segmentation vertices",
      class = 'sticky',
      amount = 0
    )
    segs <- data.table::fread(
      file = files[['segmentations.file']],
      sep = ',',
      data.table = FALSE,
      verbose = FALSE,
      tmpdir = tempdir
    )
    #segs <- data.table(read.csv(files[['segmentations.file']], stringsAsFactors = FALSE))
    
    
    print(paste0("segs read: ", nrow(segs)))
    
    
    # filter metadata file by FOVs
    if (!is.null(x = fov.filter)) {
      segs <- segs[segs$fov %in% fov.filter,]
    }
    ppresegs(type = 'finish')
  }
  
  ## Molecule coordinates - AND matrix?
  
  # Check for loading of molecule coordinates
  if (!is.na(x = files[['molecules.file']])) {
    ppremol <- progressor()
    ppremol(
      message = "Preloading molecule coordinates",
      class = 'sticky',
      amount = 0
    )
    mx <- data.table::fread(
      file = files[['molecules.file']],
      sep = ',',
      verbose = FALSE,
      tmpdir = tempdir
    )
    #mx <- data.table( read.csv(file = files[['molecules.file']], stringsAsFactors = FALSE))
    
    
    
    
    print(paste0("mx read: ", nrow(mx)))
    print(paste0("mx unique cells: ", length(unique(mx$cell))  ))
    
    
    # filter molecules file by FOVs
    if (!is.null(x = fov.filter)) {
      mx <- mx[mx$fov %in% fov.filter,]
    }
    
    
    
    # Molecules outside of a cell have a cell_ID of 0
    if (cell.mols.only) {
      mx <- mx[mx$cell_ID != 0,]
    }
    
    print(paste0("mx filtere cellular: ", nrow(mx)))

    
    if (!is.na(x = mols.filter)) {
      ppremol(
        message = paste("Filtering molecules with pattern", mols.filter),
        class = 'sticky',
        amount = 0
      )
      mx <- mx[!grepl(pattern = mols.filter, x = mx$target), , drop = FALSE]
    }
    ppremol(type = 'finish')
    mols <- rep_len(x = files[['molecules.file']], length.out = length(x = mol.type))
    names(x = mols) <- mol.type
    files <- c(files, mols)
    files <- files[setdiff(x = names(x = files), y = 'molecules.file')]
  }
  files <- files[!is.na(x = files)]
  
  outs <- list("matrix"=NULL, "pixels"=NULL, "centroids"=NULL)
  if (!is.null(metadata)) {
    outs <- append(outs, list("metadata" = NULL))
  }
  if ("segmentations" %in% type) {
    outs <- append(outs, list("segmentations" = NULL))
  }
  
  
  ####--------------------------------------------------------------------------
  # Process to relevant outputs
  # Matrix, centroids and segmentations
  
  
  for (otype in names(x = outs)) {
    outs[[otype]] <- switch(
      EXPR = otype,
      'matrix' = {
        
        ptx <- progressor()
        ptx(message = 'Building counts matrix', class = 'sticky', amount = 0)
        
        # build a cell matrix from just the molecules file, 
        # There is no need to read from the 'mtx' file at all 
        # (the 'mtx' file is not MM format, it is a dense tsv with extra cols, 
        # And is difficult to read efficently, and has no extra information)
        tx <- t(build.cellcomp.matrix.X(mols.df=mx, class=subset.counts.matrix))
        
        ptx(type = 'finish')
        
        if (!is.na(x = genes.filter)) {
          ptx(
            message = paste("Filtering genes with pattern", genes.filter),
            class = 'sticky',
            amount = 0
          )
          tx <- tx[!grepl(pattern = genes.filter, x = rownames(x = tx)), , drop = FALSE]
          ptx(type = 'finish')
        }
        # only keep cells with counts greater than 0
        tx <- tx[, which(colSums(tx) != 0)]

        tx        
      },
        
        # ptx <- progressor()
        # ptx(message = 'Reading counts matrix', class = 'sticky', amount = 0)
        # if (!is.null(subset.counts.matrix)) {
        #   
        #   #TODO : Implement subset counts matrix
        #   
        #   tx <- ReadExprMatFile(files[[otype]], chunk_size = chunk_size)
        #   ###tx <- build.cellcomp.matrix(mols.df=mx, class=subset.counts.matrix)
        # } else {
        #   
        #   tx <- ReadExprMatFile(files[[otype]], chunk_size = chunk_size)
        #   
        #   tx <- data.table::fread(
        #     file = files[[otype]],
        #     sep = ',',
        #     data.table = FALSE,
        #     verbose = FALSE
        #   )
        #   
        #   #TODO: Subset fov
        #   # TODO : cellIDs of zero?
        #   
        #   # # Combination of Cell ID (for non-zero cell_IDs) and FOV are assumed to be unique. Used to create barcodes / rownames.
        #   # bcs <- paste0(as.character(tx$cell_ID), "_", tx$fov)
        #   # rownames(x = tx) <- bcs
        #   # # remove all rows which represent counts of mols not assigned to a cell for each FOV
        #   # tx <- tx[!tx$cell_ID == 0,]
        #   # # filter fovs from counts matrix
        #   # if (!is.null(x = fov.filter)) {
        #   #   tx <- tx[tx$fov %in% fov.filter,]
        #   # }
        #   # tx <- subset(tx, select = -c(fov, cell_ID))
        # }
        
        # Now flip cells to cols
        #tx <- t(tx)

        #tx
        # },
      'centroids' = {
        pcents <- progressor()
        pcents(
          message = 'Creating centroid coordinates',
          class = 'sticky',
          amount = 0
        )
        pcents(type = 'finish')
        data.frame(
          x = md$CenterX_global_px,
          y = md$CenterY_global_px,
          cell = paste0(as.character(md$cell_ID), "_", md$fov),
          stringsAsFactors = FALSE
        )
      },
      'segmentations' = {
        pcents <- progressor()
        pcents(
          message = 'Creating segmentation coordinates',
          class = 'sticky',
          amount = 0
        )
        pcents(type = 'finish')
        data.frame(
          x = segs$x_global_px,
          y = segs$y_global_px,
          cell = paste0(as.character(segs$cellID), "_", segs$fov),  # cell_ID column in this file doesn't have an underscore
          stringsAsFactors = FALSE
        )
      },
      'metadata' = {
        pmeta <- progressor()
        pmeta(
          message = 'Loading metadata',
          class = 'sticky',
          amount = 0
        )
        pmeta(type = 'finish')
        df <- md[,metadata] 
        df$cell <- paste0(as.character(md$cell_ID), "_", md$fov)
        rownames(df) <- df$cell
        
        print(paste0("unique md df:", nrow(df)))
        
        df
      },
      'pixels' = {
        ppixels <- progressor()
        ppixels(
          message = 'Creating pixel-level molecule coordinates',
          class = 'sticky',
          amount = 0
        )
        df <- data.frame(
          x = mx$x_global_px,
          y = mx$y_global_px,
          gene = mx$target,
          stringsAsFactors = FALSE
        )
        ppixels(type = 'finish')
        df
      },
      # 'microns' = {
      #   pmicrons <- progressor()
      #   pmicrons(
      #     message = "Creating micron-level molecule coordinates",
      #     class = 'sticky',
      #     amount = 0
      #   )
      #   df <- data.frame(
      #     x = mx$global_x,
      #     y = mx$global_y,
      #     gene = mx$gene,
      #     stringsAsFactors = FALSE
      #   )
      #   pmicrons(type = 'finish')
      #   df
      # },
      stop("Unknown Nanostring input type: ", outs[[otype]])
    )
  }
  return(outs)
}












########################################################################
#' @inheritParams ReadAkoya
#' @param data.dir Path to folder containing Nanostring SMI outputs
#'
#' @return \code{LoadNanostring}: A \code{\link[SeuratObject]{Seurat}} object
#'
#' @importFrom SeuratObject Cells CreateCentroids CreateFOV
#' CreateSegmentation CreateSeuratObject
#'
#' @export
#'
#' @rdname ReadNanostring
#'
#'
#' #### TODO: Should fov.filter be included?
#' 
#' 
#' 
#' 
 
#data.dir=slide_1_file
#fov="test"
#assay = 'Nanostring'
#fov.filter = NULL
LoadNanostring.X <- function(data.dir, fov, assay = 'Nanostring',
                             fov.filter = NULL,
                             tempdir    = NULL,
                             ...) {
  data <- ReadNanostring.X(
    data.dir   = data.dir,
    type       = c("centroids", "segmentations"),
    fov.filter = fov.filter,
    tempdir    = tempdir
  )
  
  segs  <- CreateSegmentation(data$segmentations)
  cents <- CreateCentroids(data$centroids)
  segmentations.data <- list(
    "centroids"    = cents,
    "segmentation" = segs
  )
  coords <- CreateFOV(
    coords    = segmentations.data,
    #type      = c("centroids", "segmentations"),
    molecules = data$pixels,
    assay     = assay
  )
  
  # subset annotation and matrix to cells shared by both
  # TODO: check why - seems like a tiny handful of very low count cells, 
  # without metadata.
  # maybe just controls probes?
  cells <- intersect(colnames(data$matrix), data$metadata$cell)
  

  obj <- CreateSeuratObject(counts    = data$matrix[,cells], 
                            assay     = assay,
                            meta.data = data$metadata,
                            ...)

  # subset both object and coords based on the cells shared by both
  cells <- intersect(
    Cells(x = coords, boundary = "segmentation"),
    Cells(x = coords, boundary = "centroids")
  )
  cells <- intersect(Cells(obj), cells)
  coords <- subset(x = coords, cells = cells)
  obj[[fov]] <- coords
  return(obj)
}





