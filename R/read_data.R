#' Read raw output from kallisto/bustools
#'
#' @param barcodes_file Location of bustools output.
#' @return list containing spliced, unspliced and merged counts
#' @export
ReadKallistoOutput <- function(barcodes_file, genes_file, mtx_file) {
  mtx <- readMM(mtx_file)
  mtx <- t(mtx)
  mtx <- as(mtx, Class = "dgCMatrix")
  gene_ids <- read.csv(barcodes_file, header = F)$V1
  barcodes <- read.csv(genes_file, header = F)$V1
  rownames(mtx) <- make.unique(toupper(tx2g[gene_ids, "gene_name"]))
  colnames(mtx) <- barcodes
  return(mtx)
}

#' Read output from spliced mode of kallisto/bustools
#'
#' @param base.dir Location of bustools output.
#' @return list containing spliced, unspliced and merged counts
#' @export
ReadKBSpliced <- function(base.dir) {
  spliced.counts <- ReadKallistoOutput(
    barcodes_file = file.path(base.dir, "spliced.barcodes.txt"),
    genes_file = file.path(base.dir, "spliced.genes.txt"), mtx_file = file.path(
      base.dir,
      "spliced.mtx"
    )
  )
  unspliced.counts <- ReadKallistoOutput(
    barcodes_file = file.path(base.dir, "unspliced.barcodes.txt"),
    genes_file = file.path(base.dir, "unspliced.genes.txt"), mtx_file = file.path(
      base.dir,
      "unspliced.mtx"
    )
  )
  cells.common <- intersect(colnames(spliced.counts), colnames(unspliced.counts))

  spliced.counts <- spliced.counts[, colnames(spliced.counts) %in% cells.common]
  unspliced.counts <- unspliced.counts[, colnames(unspliced.counts) %in% cells.common]


  merged.counts <- spliced.counts + unspliced.counts
  return(list(merged.counts = merged.counts, spliced.counts = spliced.counts, unspliced.counts = unspliced.counts))
}

#' Fetch supplemetary files from GEO
#' @importFrom XML htmlParse xpathSApply
#' @importFrom httr GET parse_url
#' @export
FetchGEOFiles <- function(gse, download.dir = getwd(), download.files = FALSE, ...) {
  url.prefix <- "https://ftp.ncbi.nlm.nih.gov/geo/series/"
  gse_prefix <- paste0(substr(x = gse, start = 1, stop = nchar(gse) - 3), "nnn")
  url <- paste0(url.prefix, gse_prefix, "/", gse, "/", "suppl", "/")

  response <- GET(url = url)
  html_parsed <- htmlParse(file = response)
  links <- xpathSApply(doc = html_parsed, path = "//a/@href")
  suppl_files <- as.character(grep(pattern = "^G", x = links, value = TRUE))

  file.url <- paste0(url, suppl_files)
  file_list <- data.frame(filename = suppl_files, url = file.url)

  if (download.files) {
    names(file.url) <- suppl_files
    download_file <- function(url, filename, ...) {
      message(paste0("Downloading ", filename, " to ", download.dir))
      download.file(url = url, destfile = file.path(download.dir, filename), mode = "wb", ...)
      message("Done!")
    }
    lapply(seq_along(file.url), function(y, n, i) {
      download_file(y[[i]], n[[i]], ...)
    },
    y = file.url, n = names(file.url)
    )
  }

  return(file_list)
}


#' Load in data from remote or local mtx files
#' Adapted and inspired from Seurat's Read10X
#'
#' Enables easy loading of sparse data matrices
#'
#' @param mtx Name or remote URL of the mtx file
#' @param cells Name or remote URL of the cells/barcodes file
#' @param features Name or remote URL of the features/genes file
#' @param feature.column Specify which column of features files to use for feature/gene names; default is 2
#' @param cell.column Specify which column of cells file to use for cell names; default is 1
#' @param unique.features Make feature names unique (default TRUE)
#' @param strip.suffix Remove trailing "-1" if present in all cell barcodes.
#'
#' @return A sparse matrix containing the expression data.
#'
#' @importFrom Matrix readMM
#' @importFrom utils read.delim
#' @importFrom httr build_url parse_url
#' @importFrom tools file_ext
#'
#'
#' @export
#' @concept preprocessing
#'
#' @examples
#' \dontrun{
#' # For local files:
#'
#' expression_matrix <- ReadMtx(genes = "count_matrix.mtx.gz", features = "features.tsv.gz", cells = "barcodes.tsv.gz")
#' seurat_object <- CreateSeuratObject(counts = expression_matrix)
#'
#' # For remote files:
#'
#' expression_matrix <- ReadMtx(
#'   mtx = "http://localhost/matrix.mtx",
#'   cells = "http://localhost/barcodes.tsv",
#'   features = "http://localhost/genes.tsv"
#' )
#' seurat_object <- CreateSeuratObject(counts = data)
#' }
#'
ReadMtx <- function(mtx,
                    cells,
                    features,
                    cell.column = 1,
                    feature.column = 2,
                    unique.features = TRUE,
                    strip.suffix = FALSE) {
  mtx <- build_url(url = parse_url(url = mtx))
  cells <- build_url(url = parse_url(url = cells))
  features <- build_url(url = parse_url(url = features))
  all_files <- list("Expression matrix" = mtx, "Barcode" = cells, "Gene name" = features)

  check_file_exists <- function(filetype, filepath) {
    if (grepl(pattern = "^:///", x = filepath)) {
      filepath <- gsub(pattern = ":///", replacement = "", x = filepath)
      if (!file.exists(paths = filepath)) {
        stop(paste(filetype, "file missing. Expecting", filepath), call. = FALSE)
      }
    }
  }

  # check if all files exist
  lapply(seq_along(all_files), function(y, n, i) {
    check_file_exists(n[[i]], y[[i]])
  }, y = all_files, n = names(all_files))

  # convenience fucntion to read local or remote tab delimited files
  readTableUri <- function(uri) {
    if (grepl(pattern = "^:///", x = uri)) {
      uri <- gsub(pattern = ":///", replacement = "", x = uri)
      textcontent <- read.table(file = uri, header = FALSE, sep = "\t", row.names = NULL)
    } else {
      if (file_ext(uri) == "gz") {
        textcontent <- read.table(
          file = gzcon(url(uri), text = TRUE),
          header = FALSE, sep = "\t", row.names = NULL
        )
      } else {
        textcontent <- read.table(
          file = uri, header = FALSE,
          sep = "\t", row.names = NULL
        )
      }
    }
    return(textcontent)
  }

  # read barcodes
  cell.barcodes <- readTableUri(uri = cells)
  bcols <- ncol(x = cell.barcodes)
  if (bcols < cell.column) {
    stop(paste0(
      "cell.column was set to ", cell.column,
      " but ", cells, " only has ", bcols, " columns.",
      " Try setting the cell.column argument to a value <= to ", bcols, "."
    ))
  }
  cell.names <- cell.barcodes[, cell.column]

  if (all(grepl(pattern = "\\-1$", x = cell.names)) & strip.suffix) {
    cell.names <- as.vector(x = as.character(x = sapply(
      X = cell.names,
      FUN = ExtractField,
      field = 1,
      delim = "-"
    )))
  }

  # read features
  feature.names <- readTableUri(uri = features)
  fcols <- ncol(x = feature.names)
  if (fcols < feature.column) {
    stop(paste0(
      "feature.column was set to ", feature.column,
      " but ", features, " only has ", fcols, " column(s).",
      " Try setting the feature.column argument to a value <= to ", fcols, "."
    ))
  }
  if (any(is.na(x = feature.names[, feature.column]))) {
    na.features <- which(x = is.na(x = feature.names[, feature.column]))
    replacement.column <- ifelse(test = feature.column == 2, yes = 1, no = 2)
    if (replacement.column > fcols) {
      stop(
        paste0(
          "Some features names are NA in column ", feature.column,
          ". Try specifiying a different column."
        ),
        call. = FALSE,
        immediate. = TRUE
      )
    } else {
      warning(
        paste0(
          "Some features names are NA in column ", feature.column,
          ". Replacing NA names with ID from column ", replacement.column, "."
        ),
        call. = FALSE,
        immediate. = TRUE
      )
    }
    feature.names[na.features, feature.column] <- feature.names[na.features, replacement.column]
  }

  feature.names <- feature.names[, feature.column]
  if (unique.features) {
    feature.names <- make.unique(names = feature.names)
  }

  # read mtx
  if (grepl(pattern = "^:///", x = mtx)) {
    mtx <- gsub(pattern = ":///", replacement = "", x = mtx)
    data <- readMM(mtx)
  } else {
    if (file_ext(mtx) == "gz") {
      data <- readMM(gzcon(url(mtx)))
    } else {
      data <- readMM(mtx)
    }
  }

  colnames(x = data) <- cell.names
  rownames(x = data) <- feature.names

  return(data)
}

#' Convert parsebio's long dataframe to a sparse matrix
#' @importFrom tidyr pivot_wider
#' @importFrom magrittr %>%
#' @importFrom Matrix sparseMatrix
SparsifyParseDataframe <- function(df, all_genes, all_cells) {
  tmp <- pivot_wider(data = df, names_from = "bc_wells", values_from = "n", values_fill = 0) %>% as.data.frame()
  present_genes <- tmp$gene_name_id
  rownames(x = tmp) <- tmp$gene_name_id
  tmp$gene_name_id <- NULL
  tmp <- make.sparse(mat = tmp)
  present_cells <- colnames(x = tmp)

  missing_cells <- setdiff(x = all_cells, y = present_cells)

  # these genes are not in the matrix, but we want to include thse explicitly
  missing_genes <- setdiff(x = all_genes, y = present_genes)
  if (length(x = missing_genes) > 0) {
    missing_genes.cm <- sparseMatrix(i = {}, j = {}, dims = c(length(x = missing_genes), ncol(x = tmp)))
    rownames(x = missing_genes.cm) <- missing_genes
    tmp <- rbind(tmp, missing_genes.cm)
  }

  if (length(x = missing_cells) > 0) {
    missing_cells.cm <- sparseMatrix(i = {}, j = {}, dims = c(nrow(x = tmp), length(x = missing_cells)))
    colnames(x = missing_cells.cm) <- missing_cells
    tmp <- cbind(tmp, missing_cells.cm)
  }
  return(tmp)
}

#' Sort rownames and colnames of a matrix
#' @export
SortMatrixByName <- function(mat) {
  mat <- mat[sort(x = rownames(x = mat)), sort(x = colnames(x = mat))]
  return(mat)
}

#' Rename
#' @importFrom stringr str_split_fixed
RenameRows <- function(mat, pattern = "__", gene_field = 2) {
  rownames(x = mat) <- str_split_fixed(string = rownames(x = mat), pattern = pattern, n = Inf)[, gene_field]
  return(mat)
}

#' Read output from parsebio pipeline
#' @importFrom data.table fread
#' @importFrom Matrix readMM t
#' @importFrom stringr str_split_fixed
#' @importFrom SeuratObject CreateSeuratObject CreateAssayObject
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange filter group_by select tally
#' @importFrom sparseMatrixStats rowSums2
#' @export
ReadParsebioOutput <- function(path, add.hexR.assay = FALSE, add.polyT.assay = FALSE) {
  process.dir <- file.path(dirname(path = dirname(path = path)), "process")
  tscp_df <- fread(file = file.path(process.dir, "tscp_assignment.csv.gz"))

  tscp_df$gene_name[is.na(tscp_df$gene_name)] <- tscp_df$gene[is.na(tscp_df$gene_name)]
  tscp_df$gene_name[tscp_df$gene_name == ""] <- tscp_df$gene[tscp_df$gene_name == ""]
  tscp_df$gene_name_id <- paste0(tscp_df$gene, "__", tscp_df$gene_name)

  counts <- make.sparse(mat = t(x = readMM(file.path(path, "DGE.mtx"))))
  metadata <- read.csv(file = file.path(path, "cell_metadata.csv"), row.names = 1)

  genes <- read.csv(file = file.path(path, "all_genes.csv"))
  genes$gene_name[is.na(genes$gene_name)] <- genes$gene_id[is.na(genes$gene_name)]
  genes$gene_name[genes$gene_name == ""] <- genes$gene_id[genes$gene_name == ""]
  genes$gene_name_id <- paste0(genes$gene_id, "__", genes$gene_name)
  all_genes <- genes$gene_name_id

  colnames(x = counts) <- rownames(x = metadata)
  rownames(x = counts) <- all_genes

  counts <- RenameRows(mat = counts) %>% SortMatrixByName()

  seu <- CreateSeuratObject(counts = counts, meta.data = metadata, min.cells = 1, min.features = 1)
  cells <- Cells(seu)

  tscp_df_filtered <- tscp_df %>%
    filter(bc_wells %in% cells) %>%
    arrange(bc_wells, cell_barcode, polyN)
  tscp_df_filtered_summary_bygene <- tscp_df_filtered %>%
    group_by(bc_wells, gene_name_id) %>%
    tally()
  tscp_df_filtered_summary_bygene_byrt <- tscp_df_filtered %>%
    group_by(bc_wells, gene_name_id, rt_type) %>%
    tally()
  tscp_df_filtered_summary_bygene_byrt_R <- tscp_df_filtered_summary_bygene_byrt %>%
    filter(rt_type == "R") %>%
    select(-rt_type)
  tscp_df_filtered_summary_bygene_byrt_T <- tscp_df_filtered_summary_bygene_byrt %>%
    filter(rt_type == "T") %>%
    select(-rt_type)

  tmp_R <- SparsifyParseDataframe(tscp_df_filtered_summary_bygene_byrt_R, all_genes = all_genes, all_cells = cells) %>%
    RenameRows() %>%
    SortMatrixByName()
  tmp_T <- SparsifyParseDataframe(tscp_df_filtered_summary_bygene_byrt_T, all_genes = all_genes, all_cells = cells) %>%
    RenameRows() %>%
    SortMatrixByName()

  if (add.hexR.assay) {
    seu[["HexR"]] <- CreateAssayObject(counts = tmp_R, min.cells = 0, min.features = 0)
  }
  if (add.polyT.assay) {
    seu[["PolyT"]] <- CreateAssayObject(counts = tmp_T, min.cells = 0, min.features = 0)
  }

  R_sum <- rowSums2(x = tmp_R)
  T_sum <- rowSums2(x = tmp_T)
  all_sum <- rowSums2(x = counts)

  genewise_summary <- data.frame(gene_name_id = all_genes, R_reads = R_sum, T_reads = T_sum)
  genewise_summary$gene_id <- str_split_fixed(string = genewise_summary$gene_name_id, pattern = "__", n = 2)[, 1]
  genewise_summary$gene_name <- str_split_fixed(string = genewise_summary$gene_name_id, pattern = "__", n = 2)[, 2]
  genes_to_keep <- genewise_summary %>%
    filter(all_sum > 0) %>%
    pull(gene_name) %>%
    make.unique()
  seu <- subset(seu, features = intersect(rownames(seu), genes_to_keep))
  return(list(object = seu, summary = genewise_summary))
}
