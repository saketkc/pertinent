#' @export
getSeuObj <- function(counts) {
  total.counts <- colSums(counts)
  bc.rank <- DropletUtils::barcodeRanks(counts)

  bc.plot <- qplot(bc.rank$total, bc.rank$rank, geom = "line") + geom_vline(
    xintercept = metadata(bc.rank)$knee,
    color = "blue", linetype = 2
  ) + geom_vline(
    xintercept = metadata(bc.rank)$inflection,
    color = "green", linetype = 2
  ) + annotate("text", y = 1000, x = 1.5 * c(
    metadata(bc.rank)$knee,
    metadata(bc.rank)$inflection
  ), label = c("knee", "inflection"), color = c(
    "red",
    "green"
  )) + scale_x_log10() + scale_y_log10() + labs(
    y = "Barcode rank",
    x = "Total UMI count"
  )


  sampled.counts <- counts[, total.counts > metadata(bc.rank)$inflection]
  all.genes <- rownames(sampled.counts)

  seu.scaled <- CreateSeuratObject(counts = sampled.counts) %>%
    NormalizeData(verbose = FALSE) %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(verbose = FALSE, features = all.genes) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE) %>%
    FindNeighbors(
      dims = 1:30,
      verbose = FALSE
    ) %>%
    FindClusters(resolution = 0.5, verbose = FALSE)
  return(list(seu = seu.scaled, plot = bc.plot))
}

#' SCT vs STD workflow comparison
#' @param seu Seurat object
#' @param dims Number of dimensions for PCA/UMAP
#' @param nfeatures Number of variable features for STD workflow
#' @export
sct_vs_std <- function(seu, dims = 1:30, nfeatures = 3000, idents = c("seurat_clusters"),
                       outdir = NULL, overwrite = TRUE) {
  std <- NormalizeData(seu) %>%
    FindVariableFeatures(
      selection.method = "vst",
      nfeatures = nfeatures
    ) %>%
    ScaleData(verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = dims, verbose = FALSE) %>%
    FindNeighbors(dims = dims, verbose = FALSE) %>%
    FindClusters(verbose = FALSE)

  if (!is.null(outdir)) {
    dir.create(outdir, showWarnings = FALSE)
    seurat_clusters <- std$seurat_clusters
    seurat_clusters <- as.data.frame(seurat_clusters)
    embeddings <- as.data.frame(Embeddings(std[["umap"]]))
    write.table(embeddings, file.path(outdir, "STD_umap.tsv"), sep = "\t")
    write.table(seurat_clusters, file.path(outdir, "STD_seuratclusters.tsv"),
      sep = "\t"
    )
    saveRDS(std, file.path(outdir, paste0("STDObject", "_dims", paste(dims[1],
      dims[length(dims)],
      sep = ":"
    ), "_varfeatures", nfeatures, ".RDS")))
    DropletUtils::write10xCounts(x = std[["RNA"]]@counts, path = file.path(
      outdir,
      "STD_counts_10X"
    ), overwrite = overwrite)
    write.table(as.data.frame(seu@meta.data), file.path(
      outdir, "STD_counts_10X",
      "metadata.tsv"
    ), sep = "\t")
  }

  sct <- SCTransform(seu, ncells = 5000, variable.features.n = nfeatures, verbose = FALSE) %>%
    RunPCA(verbose = FALSE) %>%
    RunUMAP(dims = dims, verbose = FALSE) %>%
    FindNeighbors(
      dims = dims,
      verbose = FALSE
    ) %>%
    FindClusters(verbose = FALSE)

  if (!is.null(outdir)) {
    sct.obj <- sct[["SCT"]]
    vst.out <- sct.obj@misc$vst.out
    embeddings <- as.data.frame(Embeddings(sct[["umap"]]))
    write.table(embeddings, file.path(outdir, "SCT_umap.tsv"), sep = "\t")

    seurat_clusters <- sct$seurat_clusters
    seurat_clusters <- as.data.frame(seurat_clusters)
    write.table(seurat_clusters, file.path(outdir, "SCT_seuratclusters.tsv"),
      sep = "\t"
    )
    write.table(sct.obj@meta.features, file.path(outdir, "SCT_metafeatures.tsv"),
      sep = "\t"
    )
    write.table(vst.out$model_pars_fit, file.path(outdir, "SCT_model_pars_fit.tsv"),
      sep = "\t"
    )
    write.table(vst.out$model_pars, file.path(outdir, "SCT_model_pars.tsv"),
      sep = "\t"
    )
    write.table(vst.out$cell_attr, file.path(outdir, "SCT_cell_pars.tsv"), sep = "\t")
    write.table(vst.out$gene_attr, file.path(outdir, "SCT_gene_attrs.tsv"), sep = "\t")

    DropletUtils::write10xCounts(x = sct.obj@counts, path = file.path(
      outdir,
      "SCT_corrected_counts_10X"
    ), overwrite = overwrite)

    write.table(as.data.frame(seu@meta.data), file.path(
      outdir, "SCT_corrected_counts_10X",
      "metadata.tsv"
    ), sep = "\t")
    saveRDS(sct, file.path(outdir, paste0("SCTObject", "_dims", paste(dims[1],
      dims[length(dims)],
      sep = ":"
    ), "_varfeatures", nfeatures, ".RDS")))
  }
  to_return <- list(sct = sct, std = std)
  for (ident in idents) {
    if (ident %in% names(std@meta.data)) {
      Idents(std) <- std[[ident]]
      Idents(sct) <- sct[[ident]]
      p1 <- DimPlot(std, label = TRUE) + NoLegend() + ggtitle("Log normalized")
      p2 <- DimPlot(sct, label = TRUE) + NoLegend() + ggtitle("SCTransform")
      p3 <- p1 | p2
      p3
      ggsave(file.path(outdir, paste0("STD_vs_SCT_dimplot_ident_", ident, ".pdf")),
        width = 8, height = 4
      )
      ggsave(file.path(outdir, paste0("STD_vs_SCT_dimplot_ident_", ident, ".png")),
        width = 8, height = 4
      )
      to_return[[ident]] <- p3
    }
  }
  return(to_return)
}

#' Modified function to support differente filenames
#' @importFrom Matrix readMM
#' @export
Read10XCustom <- function(data.dir = NULL, gene.column = 2, unique.features = TRUE,
                          strip.suffix = FALSE, barcode_filename = "barcodes.tsv", gene_filename = "genes.tsv",
                          matrix_filename = "matrix.mtx") {
  full.data <- list()
  for (i in seq_along(along.with = data.dir)) {
    run <- data.dir[i]
    if (!dir.exists(paths = run)) {
      stop("Directory provided does not exist")
    }
    barcode.loc <- file.path(run, barcode_filename)
    gene.loc <- file.path(run, gene_filename)
    features.loc <- file.path(run, "features.tsv.gz")
    matrix.loc <- file.path(run, matrix_filename)
    # Flag to indicate if this data is from CellRanger >= 3.0
    pre_ver_3 <- file.exists(gene.loc)
    if (!pre_ver_3) {
      addgz <- function(s) {
        return(paste0(s, ".gz"))
      }
      barcode.loc <- addgz(s = barcode.loc)
      matrix.loc <- addgz(s = matrix.loc)
    }
    if (!file.exists(barcode.loc)) {
      stop("Barcode file missing. Expecting ", basename(path = barcode.loc))
    }
    if (!pre_ver_3 && !file.exists(features.loc)) {
      stop("Gene name or features file missing. Expecting ", basename(path = features.loc))
    }
    if (!file.exists(matrix.loc)) {
      stop("Expression matrix file missing. Expecting ", basename(path = matrix.loc))
    }
    data <- readMM(file = matrix.loc)
    cell.names <- readLines(barcode.loc)
    if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
      cell.names <- as.vector(x = as.character(x = sapply(
        X = cell.names, FUN = ExtractField,
        field = 1, delim = "-"
      )))
    }
    if (is.null(x = names(x = data.dir))) {
      if (i < 2) {
        colnames(x = data) <- cell.names
      } else {
        colnames(x = data) <- paste0(i, "_", cell.names)
      }
    } else {
      colnames(x = data) <- paste0(names(x = data.dir)[i], "_", cell.names)
    }
    feature.names <- read.delim(file = ifelse(test = pre_ver_3, yes = gene.loc,
      no = features.loc
    ), header = FALSE, stringsAsFactors = FALSE)
    if (any(is.na(x = feature.names[, gene.column]))) {
      warning("Some features names are NA. Replacing NA names with ID from the opposite column requested",
        call. = FALSE, immediate. = TRUE
      )
      na.features <- which(x = is.na(x = feature.names[, gene.column]))
      replacement.column <- ifelse(test = gene.column == 2, yes = 1, no = 2)
      feature.names[na.features, gene.column] <- feature.names[
        na.features,
        replacement.column
      ]
    }
    if (unique.features) {
      fcols <- ncol(x = feature.names)
      if (fcols < gene.column) {
        stop(paste0(
          "gene.column was set to ", gene.column, " but feature.tsv.gz (or genes.tsv) only has ",
          fcols, " columns.", " Try setting the gene.column argument to a value <= to ",
          fcols, "."
        ))
      }
      rownames(x = data) <- make.unique(names = feature.names[, gene.column])
    }
    # In cell ranger 3.0, a third column specifying the type of data was added and we
    # will return each type of data as a separate matrix
    if (ncol(x = feature.names) > 2) {
      data_types <- factor(x = feature.names$V3)
      lvls <- levels(x = data_types)
      if (length(x = lvls) > 1 && length(x = full.data) == 0) {
        message("10X data contains more than one type and is being returned as a list containing matrices of each type.")
      }
      expr_name <- "Gene Expression"
      if (expr_name %in% lvls) {
        # Return Gene Expression first
        lvls <- c(expr_name, lvls[-which(x = lvls == expr_name)])
      }
      data <- lapply(X = lvls, FUN = function(l) {
        return(data[data_types == l, , drop = FALSE])
      })
      names(x = data) <- lvls
    } else {
      data <- list(data)
    }
    full.data[[length(x = full.data) + 1]] <- data
  }
  # Combine all the data from different directories into one big matrix, note this
  # assumes that all data directories essentially have the same features files
  list_of_data <- list()
  for (j in 1:length(x = full.data[[1]])) {
    list_of_data[[j]] <- do.call(cbind, lapply(X = full.data, FUN = `[[`, j))
    # Fix for Issue #913
    list_of_data[[j]] <- as(object = list_of_data[[j]], Class = "dgCMatrix")
  }
  names(x = list_of_data) <- names(x = full.data[[1]])
  # If multiple features, will return a list, otherwise a matrix.
  if (length(x = list_of_data) == 1) {
    return(list_of_data[[1]])
  } else {
    return(list_of_data)
  }
}

#' convert column to factor
#' @export
Factorify <- function(values) {
  values.f <- factor(as.character(values), levels = sort(unique(as.character(values))))
  return(values.f)
}

#' Run SCT normalization
#' @export
DoSCT <- function(object) {
  object <- SCTransform(object, vst.flavor = "v2")
  object <- RunPCA(object, verbose = FALSE)
  object <- RunUMAP(object, dims = 1:30, verbose = FALSE)

  object <- FindNeighbors(object, dims = 1:30, verbose = FALSE)
  object <- FindClusters(object, verbose = FALSE)
}


#' Run Log normalization
#' @export
DoLogNorm <- function(object) {
  DefaultAssay(object) <- "RNA"
  object <- NormalizeData(object, scale.factor = median(object$nCount_RNA))
  object <- FindVariableFeatures(object)
  object <- ScaleData(object)

  object <- RunPCA(object, verbose = FALSE)
  object <- RunUMAP(object, dims = 1:30, verbose = FALSE)

  object <- FindNeighbors(object, dims = 1:30, verbose = FALSE)
  object <- FindClusters(object, verbose = FALSE)
  return(object)
}

#' Variable Feattures for SCT
#' @export
#'
VariableFeaturesSCTModel <- function(object, nfeatures = 3000, ...) {
  feature.attr <- SCTResults(object = object, slot = "feature.attributes")
  feature.variance <- feature.attr[, "residual_variance"]
  names(x = feature.variance) <- row.names(x = feature.attr)
  feature.variance <- sort(x = feature.variance, decreasing = TRUE)
  return(head(x = names(x = feature.variance), n = nfeatures))
}
