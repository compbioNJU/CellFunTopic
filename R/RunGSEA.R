
#' Perform Gene Set Enrichment Analysis (GSEA) on Seurat object
#'
#' Perform Gene Set Enrichment Analysis (GSEA) on Seurat object. GSEA is implemented using clusterProfiler package.
#'
#' @param SeuratObj Seurat object
#' @param by GO KEGG Reactome MSigDb WikiPathways DO NCG DGN. May be a character vector. Will be ignored if parameter "TERM2GENE" is not NULL.
#' @param TERM2GENE Customized terms. Annotation of TERM TO GENE mapping, a data.frame of 2 column with term and gene. NULL by default. If not NULL, parameter "by" will be ignored.
#' @param pvalueCutoff threshold of p-value in GSEA result
#' @param minpct minimum expression percent in each cluster, parameter used to filter Differential expression genes to run GSEA. 0 as default.
#' @param category MSigDB collection abbreviation, such as H, C1, C2, C3, C4, C5, C6, C7. If NULL, all gene sets.
#' @param subcategory MSigDB sub-collection abbreviation, such as CGP or BP.
#'
#' @importFrom clusterProfiler GSEA gseGO gseKEGG bitr read.gmt
#' @importFrom DOSE gseDO gseNCG gseDGN
#' @importFrom magrittr `%>%`
#' @importFrom foreach foreach `%do%`
#' @return Seurat object
#' @export
#'
#' @examples
#' \dontrun{
#' SeuratObj <- RunGSEA(SeuratObj, by = 'GO')
#' }
#'
RunGSEA <- function(SeuratObj,
                    by = 'GO',
                    TERM2GENE = NULL,
                    minpct = 0,
                    pvalueCutoff = 1,
                    category = NULL,
                    subcategory = NULL) {

  if (is.null(TERM2GENE) & (length(by) > 1)) {
    for (element in by) {
      SeuratObj <- RunGSEA_sub(SeuratObj,
                               by = element,
                               TERM2GENE = TERM2GENE,
                               minpct = minpct,
                               pvalueCutoff = pvalueCutoff,
                               category = category,
                               subcategory = subcategory)
    }
  } else {
    SeuratObj <- RunGSEA_sub(SeuratObj,
                             by = by,
                             TERM2GENE = TERM2GENE,
                             minpct = minpct,
                             pvalueCutoff = pvalueCutoff,
                             category = category,
                             subcategory = subcategory)
  }
  return(SeuratObj)
}


RunGSEA_sub <- function(SeuratObj,
                        by = 'GO',
                        TERM2GENE = NULL,
                        minpct = 0.25,
                        pvalueCutoff = 1,
                        category = NULL,
                        subcategory = NULL) {

  by <- match.arg(by, choices = c("GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"))
  if (!"species" %in% names(slot(object = SeuratObj, name = 'misc'))) {
    stop("Please run 'readData()' and specify the 'species' parameter.")
  }
  if (!"Allmarkers" %in% names(slot(object = SeuratObj, name = 'misc'))) {
    stop("Please run 'RunSeurat()' or 'Seurat::FindAllMarkers()', differentially expressed genes at
         SeuratObj@misc$Allmarkers are required to perform GSEA.")
  }
  if (!"featureData" %in% names(slot(object = SeuratObj, name = 'misc'))) {
    stop("Please run 'DetectGeneIDtype()' to determine the gene ID type of your data.")
  }
  species <- slot(object = SeuratObj, name = 'misc')[["species"]]
  Allmarkers <- slot(object = SeuratObj, name = 'misc')[["Allmarkers"]]
  Allmarkers$cluster <- as.character(Allmarkers$cluster)
  GeneIDtype <- slot(object = SeuratObj, name = 'misc')[["featureData"]][['GeneIDtype']]
  # gene ID transition
  if (GeneIDtype != "ENTREZID") {
    OrgDb <- getOrgDb(species)
    bitrDF <- bitr(Allmarkers$gene, fromType = GeneIDtype, toType = "ENTREZID", OrgDb = OrgDb, drop = TRUE)
    Allmarkers <- merge(Allmarkers, bitrDF, by.x = "gene", by.y = GeneIDtype)
  }


  if (!is.null(TERM2GENE)) {
    out <- foreach(cls=unique(Allmarkers$cluster)) %do% {
      submarkers <- Allmarkers %>% dplyr::filter(cluster==cls & pct.1 >= minpct)
      geneList <- sort(setNames(submarkers$avg_log2FC, submarkers$ENTREZID), decreasing = T)
      res <- GSEA(geneList     = geneList,
                  TERM2GENE    = TERM2GENE,
                  pvalueCutoff = pvalueCutoff,
                  verbose      = FALSE)
      df <- cbind(res@result, cluster=cls)
      df
    }
  } else {
    if (by == 'GO') {
      OrgDb <- getOrgDb(species)
      out <- foreach(cls=unique(Allmarkers$cluster)) %do% {
        submarkers <- Allmarkers %>% dplyr::filter(cluster==cls & pct.1 >= minpct)
        geneList <- sort(setNames(submarkers$avg_log2FC, submarkers$ENTREZID), decreasing = T)
        res <- gseGO(geneList     = geneList,
                     OrgDb        = OrgDb,
                     ont          = "ALL",
                     keyType      = "ENTREZID",
                     nPerm        = 1000,
                     minGSSize    = 10,
                     maxGSSize    = 500,
                     pvalueCutoff = pvalueCutoff,
                     verbose      = FALSE)
        df <- cbind(res@result, cluster=cls)
        df
      }
    }

    if (by == 'KEGG') {
      kegg_code <- get_kegg_code(species)
      out <- foreach(cls=unique(Allmarkers$cluster)) %do% {
        submarkers <- Allmarkers %>% dplyr::filter(cluster==cls & pct.1 >= minpct)
        geneList <- sort(setNames(submarkers$avg_log2FC, submarkers$ENTREZID), decreasing = T)
        res <- gseKEGG(geneList          = geneList,
                       organism          = kegg_code,
                       keyType           = 'ncbi-geneid',
                       nPerm             = 1000,
                       minGSSize         = 10,
                       maxGSSize         = 500,
                       pvalueCutoff      = pvalueCutoff,
                       pAdjustMethod     = "BH",
                       verbose           = FALSE,
                       use_internal_data = FALSE)
        df <- cbind(res@result, cluster=cls)
        df
      }
    }

    if (by == 'Reactome') {
      organism <- speToOrg(species)
      out <- foreach(cls=unique(Allmarkers$cluster)) %do% {
        submarkers <- Allmarkers %>% dplyr::filter(cluster==cls & pct.1 >= minpct)
        geneList <- sort(setNames(submarkers$avg_log2FC, submarkers$ENTREZID), decreasing = T)
        res <- ReactomePA::gsePathway(geneList     = geneList,
                          organism     = organism,
                          nPerm        = 1000,
                          minGSSize    = 10,
                          maxGSSize    = 500,
                          pvalueCutoff = pvalueCutoff,
                          verbose      = FALSE)
        df <- cbind(res@result, cluster=cls)
        df
      }
    }

    if (by == 'MSigDb') {
      if (!species %in% msigdbr::msigdbr_show_species()) {
        stop("This species is not supported by MSigDb.")
      }
      MSigDb_df <- msigdbr::msigdbr(species = species, category = category, subcategory = subcategory)
      TERM2GENE <- MSigDb_df %>% dplyr::select(gs_name, entrez_gene)
      if (!is.null(TERM2GENE)) {
        out <- foreach(cls=unique(Allmarkers$cluster)) %do% {
          submarkers <- Allmarkers %>% dplyr::filter(cluster==cls & pct.1 >= minpct)
          geneList <- sort(setNames(submarkers$avg_log2FC, submarkers$ENTREZID), decreasing = T)
          res <- GSEA(geneList     = geneList,
                      TERM2GENE    = TERM2GENE,
                      pvalueCutoff = pvalueCutoff,
                      verbose      = FALSE)
          df <- cbind(res@result, cluster=cls)
          df
        }
      }
    }

    if (by == 'WikiPathways') {
      Organisms <- rWikiPathways::listOrganisms()
      if (!species %in% Organisms) {
        stop("This species is not supported by WikiPathways.")
      }
      wpgmtfile <- rWikiPathways::downloadPathwayArchive(organism = species, format = "gmt")
      wp2gene <- read.gmt(wpgmtfile)
      wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
      wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
      wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
      out <- foreach(cls=unique(Allmarkers$cluster)) %do% {
        submarkers <- Allmarkers %>% dplyr::filter(cluster==cls & pct.1 >= minpct)
        geneList <- sort(setNames(submarkers$avg_log2FC, submarkers$ENTREZID), decreasing = T)
        res <- GSEA(geneList     = geneList,
                    TERM2GENE    = wpid2gene,
                    TERM2NAME    = wpid2name,
                    pvalueCutoff = pvalueCutoff,
                    verbose      = FALSE)
        df <- cbind(res@result, cluster=cls)
        df
      }
    }

    if (by == 'DO') {
      out <- foreach(cls=unique(Allmarkers$cluster)) %do% {
        submarkers <- Allmarkers %>% dplyr::filter(cluster==cls & pct.1 >= minpct)
        geneList <- sort(setNames(submarkers$avg_log2FC, submarkers$ENTREZID), decreasing = T)
        res <- gseDO(geneList,
                     nPerm         = 100,
                     minGSSize     = 120,
                     pvalueCutoff  = pvalueCutoff,
                     pAdjustMethod = "BH",
                     verbose       = FALSE)
        df <- cbind(res@result, cluster=cls)
        df
      }
    }

    if (by == 'NCG') {
      out <- foreach(cls=unique(Allmarkers$cluster)) %do% {
        submarkers <- Allmarkers %>% dplyr::filter(cluster==cls & pct.1 >= minpct)
        geneList <- sort(setNames(submarkers$avg_log2FC, submarkers$ENTREZID), decreasing = T)
        res <- gseNCG(geneList,
                      nPerm         = 100,
                      minGSSize     = 120,
                      pvalueCutoff  = pvalueCutoff,
                      pAdjustMethod = "BH",
                      verbose       = FALSE)
        df <- cbind(res@result, cluster=cls)
        df
      }
    }

    if (by == 'DGN') {
      out <- foreach(cls=unique(Allmarkers$cluster)) %do% {
        submarkers <- Allmarkers %>% dplyr::filter(cluster==cls & pct.1 >= minpct)
        geneList <- sort(setNames(submarkers$avg_log2FC, submarkers$ENTREZID), decreasing = T)
        res <- gseDGN(geneList,
                      nPerm         = 100,
                      minGSSize     = 120,
                      pvalueCutoff  = pvalueCutoff,
                      pAdjustMethod = "BH",
                      verbose       = FALSE)
        df <- cbind(res@result, cluster=cls)
        df
      }
    }
  }

  result <- do.call("rbind", out)
  if (!is.null(TERM2GENE)) {
    slot(object = SeuratObj, name = 'misc')[["GSEAresult_customized"]] <- result
  } else {
    slot(object = SeuratObj, name = 'misc')[[paste0("GSEAresult_", by)]] <- result
  }

  return(SeuratObj)
}
