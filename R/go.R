#' Convert text string to ratio
#' @param x String
#' @return ratio Float
#'
#' @export
convert_to_ratio <- function(x) {
  return(eval(parse(text = x)))
}


#' Format GO output of EnrichR
#' @param go Dataframe as returned by running enrichR
#' @return df Formatted dataframe
#'
#' @export
FormatEnrichrOutput <- function(go) {
  go$GeneRatio <- sapply(go$Overlap, convert_to_ratio)
  go$GeneRatio <- as.numeric(go$GeneRatio)
  go$Counts <- sapply(strsplit(go$Overlap, "/"), `[`, 1)
  go$Counts <- as.numeric(go$Counts)
  go$p.adj <- go$Adjusted.P.value
  return(go)
}

#' Plot Gene ontology terms on a dotplot
#' @param go Formatted data frame
#' @param title Title of the plot
#' @export
plot_GO <- function(go, title) {
  p <- ggplot(go, aes(x = GeneRatio, y = Term)) +
    geom_point(aes(
      size = Counts,
      color = p.adj
    )) +
    theme_bw(base_size = 20) +
    scale_colour_gradient(limits = c(
      0,
      0.1
    ), low = "red") +
    ylab(NULL) +
    ggtitle(title)
  return(p)
}

#' Run clusterProfiler GLM on
#' @param df Dataframe with columns named gene and celltype
#' @param ont Ontology type BP/MF/CC
#' @importFrom dplyr left_join
#' @importFrom magrittr %>%
#' @importFrom stringr str_wrap
#' @importFrom ggplot2 theme_bw scale_x_discrete scale_y_discrete  scale_color_gradientn
#' @importFrom RColorBrewer brewer.pal
#' @export
GOClusterGenes <- function(df, ont = "BP") {
  library(org.Hs.eg.db)
  gene.df <- clusterProfiler::bitr(unique(df$gene),
    fromType = "SYMBOL",
    toType = c("ENTREZID"),
    OrgDb = org.Hs.eg.db
  )
  colnames(gene.df)[1] <- c("gene")
  df <- left_join(df, gene.df)
  df <- df[, c("ENTREZID", "cluster")] %>% unique()

  formula_res <- clusterProfiler::compareCluster(ENTREZID ~ cluster,
    data = df,
    fun = "enrichGO",
    OrgDb = org.Hs.eg.db,
    ont = ont
  )
  p <- dotplot(formula_res) + theme_bw(base_size = 11) + scale_x_discrete(guide = guide_axis(angle = 30)) + scale_y_discrete(labels = function(egoBP) str_wrap(egoBP, width = 80)) +
    scale_color_gradientn(colors = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)))

  return(list(formula = formula_res, plot = p))
}


#' Do enrichment analysis
#'
#' @export
DoEnrichment <- function(genes) {
  library(enrichR)
  dbs <- c("GO_Molecular_Function_2021", "GO_Cellular_Component_2021", "GO_Biological_Process_2021")
  enriched <- enrichR::enrichr(genes, dbs)
  return(enriched)
}


#' @export
DoEnrichmentCelltypeWise <- function(df, ont = "BP", universe = NULL) {
  library(clusterProfiler)
  library(org.Hs.eg.db)
  gene.df <- bitr(unique(df$gene),
    fromType = "SYMBOL",
    toType = c("ENTREZID"),
    OrgDb = org.Hs.eg.db
  )
  colnames(gene.df)[1] <- c("gene")
  df <- left_join(df, gene.df)
  df <- df[, c("ENTREZID", "celltype")] %>% unique()

  if (!is.null(universe)) {
    universe.df <- bitr(unique(universe),
      fromType = "SYMBOL",
      toType = c("ENTREZID"),
      OrgDb = org.Hs.eg.db
    )
    universe.entrez <- universe.df$ENTREZID
    formula_res <- compareCluster(ENTREZID ~ celltype,
      data = df,
      fun = "enrichGO",
      OrgDb = org.Hs.eg.db,
      ont = ont,
      universe = universe.entrez
    ) %>% filter(qvalue < 0.05)
  } else {
    formula_res <- compareCluster(ENTREZID ~ celltype,
      data = df,
      fun = "enrichGO",
      OrgDb = org.Hs.eg.db,
      ont = ont
    ) %>% filter(qvalue < 0.05)
  }


  rownames(gene.df) <- gene.df$ENTREZID
  enrichment_df <- formula_res@compareClusterResult
  enrichment_df$genes <- unlist(lapply(enrichment_df$geneID, FUN = function(x) {
    paste(unlist(gene.df[unlist(stringr::str_split(x, "/")), "gene"]), collapse = ",")
  }))
  return(list(formula_res = formula_res, enricment_df = enrichment_df))
}


#' @export
DoEnrichmentClusterprofiler <- function(genes, ont = "BP", universe = NULL) {
  library(clusterProfiler)
  library(org.Hs.eg.db)
  df <- data.frame(gene = genes)
  gene.df <- bitr(df$gene,
    fromType = "SYMBOL",
    toType = c("ENTREZID"),
    OrgDb = org.Hs.eg.db
  )
  colnames(gene.df)[1] <- c("gene")

  df <- left_join(df, gene.df) %>%
    unique() %>%
    drop_na() %>%
    as.data.frame()
  rownames(df) <- df$ENTREZID

  if (!is.null(universe)) {
    universe.df <- bitr(unique(universe),
      fromType = "SYMBOL",
      toType = c("ENTREZID"),
      OrgDb = org.Hs.eg.db
    )
    universe.entrez <- universe.df$ENTREZID

    enrichment <- enrichGO(
      gene = df$ENTREZID,
      universe = universe.entrez,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.1
    )
  } else {
    enrichment <- enrichGO(
      gene = df$ENTREZID,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = ont,
      pAdjustMethod = "BH",
      pvalueCutoff = 0.05,
      qvalueCutoff = 0.1
    )
  }
  if (dim(enrichment@result)[1] >= 1) {
    enrichment@result$gene_name <- unlist(lapply(enrichment@result$geneID, FUN = function(x) {
      paste(df[stringr::str_split_fixed(x, pattern = "\\/", n = Inf)[1, ], "gene"], collapse = ",")
    }))
  }

  to_return <- list(enrichment = enrichment, plot = barplot(enrichment))
  return(to_return)
}
