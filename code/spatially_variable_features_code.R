



FindSpatiallyVariableFeatures.Seurat_EDITED <- function (object, assay = NULL, slot = "scale.data", features = NULL,
          image = NULL, selection.method = c("markvariogram", "moransi"),
          r.metric = 5, x.cuts = NULL, y.cuts = NULL, nfeatures = 2000,
          verbose = TRUE, ...)
{
  assay <- assay %||% DefaultAssay(object = object)
  features <- features %||% rownames(x = object[[assay]])
  image <- image %||% Seurat:::DefaultImage(object = object) # ADDING ::: for workaround
  tc <- GetTissueCoordinates(object = object[[image]])


  # MODIFICED FOR LOCAL COPY
  object[[assay]] <- FindSpatiallyVariableFeatures.StdAssay_EDITED(object = object[[assay]],
                                                   slot = slot, features = features, spatial.location = tc,
                                                   selection.method = selection.method, r.metric = r.metric,
                                                   x.cuts = x.cuts, y.cuts = y.cuts, nfeatures = nfeatures,
                                                   verbose = verbose, ...)
  object <- SeuratObject:::LogSeuratCommand(object = object)
}






FindSpatiallyVariableFeatures.StdAssay_EDITED <- function (object, layer = "scale.data", spatial.location, selection.method = c("markvariogram",
                                                                                                                                  "moransi"), features = NULL, r.metric = 5, x.cuts = NULL,
                                                             y.cuts = NULL, nfeatures = nfeatures, verbose = TRUE, ...)
{
  features <- features %||% rownames(x = object)
  if (selection.method == "markvariogram" && "markvariogram" %in%
      names(x = Misc(object = object))) {
    features.computed <- names(x = Misc(object = object,
                                        slot = "markvariogram"))
    features <- features[!features %in% features.computed]
  }
  data <- GetAssayData(object = object, layer = layer)
  data <- as.matrix(x = data[features, ])
  data <- data[Seurat:::RowVar(x = data) > 0, ]
  if (nrow(x = data) != 0) {
    svf.info <- FindSpatiallyVariableFeatures(object = data,
                                              spatial.location = spatial.location, selection.method = selection.method,
                                              r.metric = r.metric, x.cuts = x.cuts, y.cuts = y.cuts,
                                              verbose = verbose, ...)
  }
  else {
    svf.info <- c()
  }
  if (selection.method == "markvariogram") {
    if ("markvariogram" %in% names(x = Misc(object = object))) {
      svf.info <- c(svf.info, Misc(object = object, slot = "markvariogram"))
    }
    suppressWarnings(expr = Misc(object = object, slot = "markvariogram") <- svf.info)
    svf.info <- ComputeRMetric(mv = svf.info, r.metric)
    svf.info <- svf.info[order(svf.info[, 1]), , drop = FALSE]
  }
  if (selection.method == "moransi") {
    colnames(x = svf.info) <- paste0("MoransI_", colnames(x = svf.info))
    svf.info <- svf.info[order(svf.info[, 2], -abs(svf.info[,
                                                            1])), , drop = FALSE]
  }
  var.name <- paste0(selection.method, ".spatially.variable")
  var.name.rank <- paste0(var.name, ".rank")
  svf.info[[var.name]] <- FALSE
  svf.info[[var.name]][1:(min(nrow(x = svf.info), nfeatures))] <- TRUE
  svf.info[[var.name.rank]] <- 1:nrow(x = svf.info)

  ############################
  ## This line, attempting to store metadata,  produces error:

  print(">>>> USING EDITED FUNCTION!!!! <<<")

  # Error in `LayerData<-`:
  # ! 'layer' must be a single non-empty string
  #https://github.com/satijalab/seurat/issues/8226
  #object[names(x = svf.info)] <- svf.info

  ## Substituting:
  object <- AddMetaData(object, metadata = svf.info) # Works

  ############################


  return(object)
}
