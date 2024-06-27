# Modified version of Seurats BuildNicheAssay that runs acros all fovs at once.
BuildNicheAssay.using_all_fovs <- function(
    object,
    group.by,
    assay = "niche",
    neighbors.k = 20,
    niches.k = 4
) {

  sum.mtx.list <- list()

  # matched order
  group_by_cols <-unique(unlist(object[[group.by]]))


  for (fov in names(object@images)) {
    coords <- GetTissueCoordinates(object[[fov]], which = "centroids")
    cells <- coords$cell

    rownames(coords) <- cells
    coords <- as.matrix(coords[ , c("x", "y")])
    neighbors    <- FindNeighbors(coords, k.param = neighbors.k)
    neighbors$nn <- neighbors$nn[Cells(object[[fov]]), Cells(object[[fov]])]

    # Merge with previous
    #neighbors.all <- merge_neighbours(neighbors.all, neighbors)


    # Need a list of all those cells, used later.
    # put all cells into one list (being sure order matches network, might not match object order.)
    cells <- colnames(neighbors$nn)


    # Continuouing on the BuildNicheAssay function
    # build cell x cell type matrix
    ct.mtx <- matrix(
      data = 0,
      nrow = length(cells),
      ncol = length(unlist(unique(object[[group.by]])))
    )

    rownames(ct.mtx) <- cells
    colnames(ct.mtx) <- group_by_cols #unique(unlist(object[[group.by]]))

    cts <- object[[group.by]]
    for (i in 1:length(cells)) {
      ct <- as.character(cts[cells[[i]], ])
      ct.mtx[cells[[i]], ct] <- 1
    }

    # recorde neibhour cell type assay for this fov/image
    sum.mtx.list[[fov]] <- as.matrix(neighbors$nn %*% ct.mtx)

  }

  # Combine the niche assay matricies into one.
  # Cell type should be common across slides.
  sum.mtx <- bind_rows(lapply(sum.mtx.list, FUN=as.data.frame)) #avoid a direct tibble conversion

  # create niche assay
  niche.assay <- CreateAssayObject(counts = t(sum.mtx))
  object[[assay]] <- niche.assay
  DefaultAssay(object) <- assay

  # Not all cells in niche assay? Perhaps isolated?? Resulting in NA?

  # cluster niches assay
  object <- ScaleData(object)
  results <- kmeans(
    x = t(object[[assay]]@scale.data),
    centers = niches.k,
    nstart = 30
  )
  object$niches <- results[["cluster"]]


  return(object)


}
