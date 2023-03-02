#' Convert a distance matrix to a data frame
#' @export
DistToDF <- function(inDist) {
  if (class(inDist) != "dist") stop("wrong input type")
  A <- attr(inDist, "Size")
  B <- if (is.null(attr(inDist, "Labels"))) sequence(A) else attr(inDist, "Labels")
  if (isTRUE(attr(inDist, "Diag"))) attr(inDist, "Diag") <- FALSE
  if (isTRUE(attr(inDist, "Upper"))) attr(inDist, "Upper") <- FALSE
  data.frame(
    row = B[unlist(lapply(sequence(A)[-1], function(x) x:A))],
    col = rep(B[-length(B)], (length(B) - 1):1),
    value = as.vector(inDist)
  )
}


#' Create a pseudobulk of input matrix given cell.group
#' @importFrom stringr str_split_fixed
#' @importFrom SeuratObject CreateSeuratObject
#' @importFrom magrittr  %>%
#' @importFrom dplyr group_by pull bind_rows row_number
#' @importFrom tidyr drop_na nest
#' @importFrom tibble as_tibble
#' @importFrom Matrix sparse.model.matrix
#' @param counts A sparse matrix
#' @param cell.group A vector of assigned group/group (one per cell)
#' @param k.replicates Number of replicates
#' @param seed Integer seed
#' @param return.seurt bool Whether seurat object should be returned
#' @examples
#' data("pbmc_small")
#' counts <- GetAssayData(object = pbmc_small, assay = "RNA", slot = "counts")
#' cell.group <- pbmc_small$groups
#' k.replicates <- 3
#' counts.agg <- CreateKPseudoBulk(counts = counts, cell.group = cell.group, k.replicates = k.replicates)
#' @export
CreateKPseudoBulk <- function(counts,
                              cell.group,
                              k.replicates = 6,
                              seed = 42,
                              return.seurat = FALSE) {
  set.seed(seed = seed)
  metadata <- data.frame(
    barcode = colnames(counts),
    group = as.character(cell.group)
  )
  metadata_group <- list()

  for (group in setdiff(x = unique(x = metadata$group), y = NA)) {
    metadata.subset <- metadata[metadata$group == group, ] %>%
      drop_na() %>%
      as_tibble()
    # shuffle the rows just in case they were ordered somehow earlier
    metadata.subset <- metadata.subset[sample(nrow(metadata.subset)), ]
    partitioned_data <- metadata.subset %>%
      group_by((row_number() - 1) %/% (n() / k.replicates)) %>%
      nest() %>%
      pull(data)
    names(partitioned_data) <- paste0("replicate", names(partitioned_data))
    metadata_joined <- bind_rows(partitioned_data, .id = "replicate") %>% drop_na()
    metadata_group[[group]] <- metadata_joined
  }
  metadata_final <- bind_rows(metadata_group) %>% as.data.frame()
  rownames(metadata_final) <- metadata_final$barcode

  if (length(x = unique(x = metadata_final$group))==1) {
    mymodel.matrix <- sparse.model.matrix(
      object = ~ 0 + replicate,
      data = metadata_final
    )
  } else {
    mymodel.matrix <- sparse.model.matrix(
      object = ~ 0 + replicate:group,
      data = metadata_final
    )
  }
  colnames(x = mymodel.matrix) <- sapply(
    X = colnames(x = mymodel.matrix),
    FUN = function(name) {
      name <- gsub(pattern = "group", replacement = "", x = name)
      return(paste0(rev(x = unlist(x = strsplit(x = name, split = ":"))),
        collapse = "__"
      ))
    }
  )
  counts.use <- counts[, metadata_final$barcode]
  counts.agg <- as.matrix(x = (counts.use %*% mymodel.matrix))
  if (return.seurat) {
    object <- CreateSeuratObject(counts = counts.agg)
    metadata_split <- str_split_fixed(string = colnames(x = counts.agg),
                                      pattern = "__", n = 2)
    object$group <- metadata_split[, 1]
    object$replicate <- metadata_split[, 2]
    return(object)
  }
  return(counts.agg)
}

#' Calculate average distance given a set of replicate distances
#' @export
AverageReplicates <- function(dist.matrix) {
  if (class(dist.matrix)[1] == "dist") {
    dists.df <- DistToDF(inDist = dist.matrix)
  } else {
    dists.df <- dist.matrix %>%
      as.data.frame() %>%
      tibble::rownames_to_column() %>%
      tidyr::pivot_longer(-rowname)
    colnames(dists.df) <- c("row", "col", "value")
  }


  dists.df$row_replicate <- stringr::str_split_fixed(dists.df$row, pattern = "_", n = 3)[, 3]
  dists.df$col_replicate <- stringr::str_split_fixed(dists.df$col, pattern = "_", n = 3)[, 3]

  dists.df$row2 <- stringr::str_split_fixed(dists.df$row, pattern = "_replicate", n = 2)[, 1]
  dists.df$col2 <- stringr::str_split_fixed(dists.df$col, pattern = "_replicate", n = 2)[, 1]

  dists.df$row2 <- gsub("_", "-", dists.df$row2)
  dists.df$col2 <- gsub("_", "-", dists.df$col2)


  dists.df.mean <- dists.df %>%
    group_by(row2, col2) %>%
    summarise(value = mean(value))
  # return(dists.df.mean)

  dists.mean <- usedist::pivot_to_numeric_matrix(data = dists.df.mean, row2, col2, value)
  # dists.mean <- as.dist(dists.mean)
  return(dists.mean)
}
