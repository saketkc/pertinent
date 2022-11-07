#' Read raw output from kallisto/bustools
#'
#' @importFrom Matrix readMM t
#' @param barcodes_file Location of bustools output.
#' @return list containing spliced, unspliced and merged counts
#' @export
ReadKallistoOutput <- function(barcodes_file, genes_file, mtx_file) {
  mtx <- readMM(file = mtx_file)
  mtx <- t(x = mtx)
  mtx <- make.sparse(mat = mtx)
  gene_ids <- read.csv(file = barcodes_file, header = F)$V1
  barcodes <- read.csv(file = genes_file, header = F)$V1
  rownames(mtx) <- make.unique(names = toupper(x = tx2g[gene_ids, "gene_name"]))
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
  cells.common <- intersect(x = colnames(x = spliced.counts), y = colnames(x = unspliced.counts))

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
#' @importFrom Matrix readMM rowSums t
#' @importFrom stringr str_split_fixed
#' @importFrom SeuratObject CreateSeuratObject CreateAssayObject
#' @importFrom magrittr %>%
#' @importFrom dplyr arrange filter full_join group_by rename select summarise tally
#' @importFrom sparseMatrixStats rowSums2
#' @export
ReadParsebioOutput <- function(path, add.hexR.assay = TRUE, add.polyT.assay = TRUE, verbose = TRUE) {
  process.dir <- file.path(dirname(path = dirname(path = path)), "process")
  if (verbose) {
    message("Reading tscp file...")
  }
  tscp_df <- fread(file = file.path(process.dir, "tscp_assignment.csv.gz"), showProgress = FALSE)

  tscp_df$gene_name[is.na(x = tscp_df$gene_name)] <- tscp_df$gene[is.na(x = tscp_df$gene_name)]
  tscp_df$gene_name[tscp_df$gene_name == ""] <- tscp_df$gene[tscp_df$gene_name == ""]
  tscp_df$gene_name_id <- paste0(tscp_df$gene, "__", tscp_df$gene_name)

  counts <- make.sparse(mat = t(x = readMM(file = file.path(path, "DGE.mtx"))))
  metadata <- read.csv(file = file.path(path, "cell_metadata.csv"), row.names = 1)

  genes <- read.csv(file = file.path(path, "all_genes.csv"))
  genes$gene_name[is.na(x = genes$gene_name)] <- genes$gene_id[is.na(x = genes$gene_name)]
  genes$gene_name[genes$gene_name == ""] <- genes$gene_id[genes$gene_name == ""]
  genes$gene_name_id <- paste0(genes$gene_id, "__", genes$gene_name)
  all_genes <- genes$gene_name_id

  colnames(x = counts) <- rownames(x = metadata)
  rownames(x = counts) <- all_genes
  counts <- SortMatrixByName(mat = counts)

  if (verbose) {
    message("Creating object...")
  }
  seu <- CreateSeuratObject(counts = counts %>% RenameRows(), meta.data = metadata, min.cells = 1, min.features = 1)
  cells <- Cells(seu)

  if (verbose) {
    message("Filtering tscp file...")
  }
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

  if (verbose) {
    message("Creating hexR and polyA matrices ...")
  }
  tmp_R <- SparsifyParseDataframe(tscp_df_filtered_summary_bygene_byrt_R, all_genes = all_genes, all_cells = cells) %>%
    SortMatrixByName()

  tmp_T <- SparsifyParseDataframe(tscp_df_filtered_summary_bygene_byrt_T, all_genes = all_genes, all_cells = cells) %>%
    SortMatrixByName()

  gene_summary <- data.frame(gene_name_id = all_genes, T_reads = rowSums(tmp_R)[all_genes], R_reads = rowSums(tmp_T), total_reads = rowSums(x = counts)[all_genes])
  gene_id <- str_split_fixed(string = gene_summary$gene_name_id, pattern = "__", n = 2)[, 1]
  gene_name <- str_split_fixed(string = gene_summary$gene_name_id, pattern = "__", n = 2)[, 2]
  gene_summary$gene_id <- gene_id
  gene_summary$gene_name <- gene_name
  gene_summary <- gene_summary %>%
    select(gene_name_id, gene_id, gene_name, R_reads, T_reads, total_reads)

  if (add.hexR.assay) {
    seu[["HexR"]] <- CreateAssayObject(counts = tmp_R %>% RenameRows(), min.cells = -1, min.features = -1)
  }
  if (add.polyT.assay) {
    seu[["PolyT"]] <- CreateAssayObject(counts = tmp_T %>% RenameRows(), min.cells = -1, min.features = -1)
  }
  return(list(object = seu, summary = gene_summary))
}
