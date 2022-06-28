#' Title
#' detect what gene ID type that rownames(counts) belongs to
#' @param SeuratObj Seurat object
#' @importFrom AnnotationDbi keys
#' @importFrom dplyr `%>%`
#' @return
#' @export
#'
#' @examples
#' #' \dontrun{
#' SeuratObj <- DetectGeneIDtype(SeuratObj)
#' GeneIDtype <- slot(object = SeuratObj, name = 'misc')[["featureData"]][['GeneIDtype']]
#' }
#'
DetectGeneIDtype <- function(SeuratObj) {
  if (!inherits(x = SeuratObj, what = 'Seurat')) {
    stop("Seurat object should be provided.")
  }
  species <- slot(object = SeuratObj, name = 'misc')[["species"]]
  OrgDb <- getOrgDb(species)
  featureData <- list()
  # featureData[['OrgDb']] <- OrgDb
  allGenelist <- list()
  for (ID in keytypes(OrgDb)) {
    allGenelist[[ID]] <- keys(OrgDb, keytype = ID)
  }
  # featureData[['allGene']] <- allGenelist
  # genes <- slot(object = SeuratObj, name = 'assays')[['RNA']] %>% slot(name = 'counts') %>% rownames()
  genes <- rownames(x = SeuratObj)
  GeneIDtype <- sapply(allGenelist, function(x){
    length(intersect(genes, x))
  }) %>% which.max() %>% names()
  featureData[['GeneIDtype']] <- GeneIDtype
  slot(object = SeuratObj, name = 'misc')[["featureData"]] <- featureData
  return(SeuratObj)
}



#' Calculate the proportion of transcripts mapping to mitochondrial genes
#'
#' @param SeuratObj Seurat object
#' @param by use_internal_data AnnotationHub biomaRt
#'
#' @importFrom AnnotationHub AnnotationHub query
#' @importFrom AnnotationDbi keys
#' @importFrom biomaRt useMart searchDatasets useDataset getBM
#' @importFrom dplyr `%>%`
#' @importFrom Seurat PercentageFeatureSet
#' @importFrom clusterProfiler bitr
#' @return
#' @export
#'
#' @examples
#' SeuratObj <- CalMTpercent(SeuratObj, by = 'use_internal_data')
#' SeuratObj <- CalMTpercent(SeuratObj, by = 'AnnotationHub')
#' SeuratObj <- CalMTpercent(SeuratObj, by = 'biomaRt')
#'
#'
CalMTpercent <- function(SeuratObj, by = 'use_internal_data') {
  by <- match.arg(by, choices = c('use_internal_data', 'AnnotationHub', 'biomaRt'))
  if (by == 'use_internal_data') {
    # Calculate nCount and nFeature
    if (!"nCount_RNA" %in% names(slot(object = SeuratObj, name = 'meta.data'))) {
      n.calc <- CalcN(object = slot(object = SeuratObj, name = 'assays')[['RNA']])
      if (!is.null(x = n.calc)) {
        names(x = n.calc) <- paste(names(x = n.calc), 'RNA', sep = '_')
        slot(object = SeuratObj, name = 'meta.data')[,names(x = n.calc)] <- n.calc
      }
    }
    SeuratObj <- DetectGeneIDtype(SeuratObj)
    genes <- rownames(x = SeuratObj)
    # find mitochondrial genes
    utils::data(list="txid_gene_MT", package="CellFunMap")
    GeneIDtype <- slot(object = SeuratObj, name = 'misc')[["featureData"]][['GeneIDtype']]
    species <- slot(object = SeuratObj, name = 'misc')[["species"]]
    if (GeneIDtype == "ENTREZID") {
      MTgenes <- genes[genes %in% txid_gene_MT$ENTREZID]
    } else {
      # gene ID transition
      OrgDb <- getOrgDb(species)
      bitrDF <- bitr(genes, fromType = GeneIDtype, toType = "ENTREZID", OrgDb = OrgDb, drop = TRUE)
      MTgenes <- bitrDF[bitrDF$ENTREZID %in% txid_gene_MT$ENTREZID, GeneIDtype]
    }
    rm(txid_gene_MT, envir = .GlobalEnv)
    spes <- c("Homo sapiens", "Mus musculus", "Arabidopsis thaliana")
    if (species %in% spes) {
      utils::data(list="MTGENE", package="CellFunMap")
      MTgenes <- c(MTgenes, unlist(sapply(MTGENE[[species]], function(x){genes[genes %in% x]})))
      rm(MTGENE, envir = .GlobalEnv)
    } else {
      utils::data(list="txid_species_gene_MT", package="CellFunMap")
      sss <- unique(txid_species_gene_MT[grep(species, txid_species_gene_MT$species, ignore.case = TRUE), c("ENTREZID", "SYMBOL")])
      MTgenes <- c(MTgenes, unlist(sapply(as.list(sss), function(x){genes[genes %in% x]})))
      rm(txid_species_gene_MT, envir = .GlobalEnv)
    }
    # Calculate the proportion of transcripts mapping to mitochondrial genes
    SeuratObj <- PercentageFeatureSet(SeuratObj, features = unique(MTgenes), col.name = "percent.mt")
    slot(object = SeuratObj, name = 'misc')[["featureData"]][["MTgenes"]] <- unique(MTgenes)
    return(SeuratObj)
  }

  if (by == 'AnnotationHub') {
    if (!inherits(x = SeuratObj, what = 'Seurat')) {
      stop("Seurat object should be provided.")
    }
    species <- slot(object = SeuratObj, name = 'misc')[['species']]
    EnsDb <- getEnsDb_AH(species)
    featureData <- list()
    # featureData[['EnsDb']] <- EnsDb

    # detect what gene ID type that rownames(counts) belongs to
    allGenelist <- list()
    for (ID in keytypes(EnsDb)) {
      allGenelist[[ID]] <- keys(EnsDb, keytype = ID)
    }
    # featureData[['allGene']] <- allGenelist
    # genes <- slot(object = SeuratObj, name = 'assays')[['RNA']] %>% slot(name = 'counts') %>% rownames()
    genes <- rownames(x = SeuratObj)
    GeneIDtype <- sapply(allGenelist, function(x){
      length(intersect(genes, x))
    }) %>% which.max() %>% names()
    featureData[['GeneIDtype']] <- ifelse(GeneIDtype == "GENEID", "ENSEMBL", GeneIDtype)
    # featureData[['GeneIDtype']] <- GeneIDtype

    geneTable <- select(EnsDb, keys = genes, columns = c("GENEID", "ENTREZID", "SYMBOL", "SEQNAME", "GENEBIOTYPE", "GENENAME"), keytype = GeneIDtype)
    featureData[['geneTable']] <- geneTable
    # Calculate nCount and nFeature
    if (!"nCount_RNA" %in% names(slot(object = SeuratObj, name = 'meta.data'))) {
      n.calc <- CalcN(object = slot(object = SeuratObj, name = 'assays')[['RNA']])
      if (!is.null(x = n.calc)) {
        names(x = n.calc) <- paste(names(x = n.calc), 'RNA', sep = '_')
        slot(object = SeuratObj, name = 'meta.data')[,names(x = n.calc)] <- n.calc
      }
    }
    # Calculate the proportion of transcripts mapping to mitochondrial genes
    MTgenes <- geneTable[geneTable$SEQNAME == "MT", GeneIDtype]
    SeuratObj <- PercentageFeatureSet(SeuratObj, features = unique(MTgenes), col.name = "percent.mt")
    slot(object = SeuratObj, name = 'misc')[["featureData"]] <- featureData
    return(SeuratObj)
  }

  if (by == 'biomaRt') {
    SeuratObj <- DetectGeneIDtype(SeuratObj)
    GeneIDtype <- slot(object = SeuratObj, name = 'misc')[["featureData"]][['GeneIDtype']]
    if (!requireNamespace('biomaRt', quietly = TRUE)) {
      stop("Please install biomaRt from Bioconductor at first")
    }
    species <- slot(object = SeuratObj, name = 'misc')[['species']]
    ensembl=useMart("ensembl")
    string <- strsplit(tolower(species), split = ' ')[[1]]
    pattern <- paste0(substr(string[1], 1, 1), string[2])
    searchResult <- searchDatasets(mart = ensembl, pattern = pattern)
    mart = useDataset(searchResult$dataset[1], mart=ensembl)
    # slot(object = SeuratObj, name = 'misc')[["featureData"]][['biomaRt']] <- mart
    filters <- switch(GeneIDtype,
                      ENSEMBL = "ensembl_gene_id",
                      ENTREZID = "entrezgene_id",
                      SYMBOL = "external_gene_name",
                      REFSEQ = "refseq_mrna")
    # genes <- slot(object = SeuratObj, name = 'assays')[['RNA']] %>% slot(name = 'counts') %>% rownames()
    genes <- rownames(x = SeuratObj)
    geneTable <- getBM(attributes=c('ensembl_gene_id', 'external_gene_name', 'entrezgene_id', 'refseq_mrna', "chromosome_name", "start_position","end_position", "band"),
                       filters= filters, values = genes, mart = mart)
    slot(object = SeuratObj, name = 'misc')[["featureData"]][['geneTable']] <- geneTable
    # Calculate nCount and nFeature
    if (!"nCount_RNA" %in% names(slot(object = SeuratObj, name = 'meta.data'))) {
      n.calc <- CalcN(object = slot(object = SeuratObj, name = 'assays')[['RNA']])
      if (!is.null(x = n.calc)) {
        names(x = n.calc) <- paste(names(x = n.calc), 'RNA', sep = '_')
        slot(object = SeuratObj, name = 'meta.data')[,names(x = n.calc)] <- n.calc
      }
    }
    # Calculate the proportion of transcripts mapping to mitochondrial genes
    MTgenes <- geneTable[geneTable$chromosome_name == "MT", filters]
    SeuratObj <- PercentageFeatureSet(SeuratObj, features = unique(MTgenes), col.name = "percent.mt")
    return(SeuratObj)
  }
}





#' Title
#' quality control according to 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'CellsPerGene'
#'
#' @param SeuratObj Seurat object
#' @param minCountsPerCell minCountsPerCell
#' @param maxCountsPerCell maxCountsPerCell
#' @param minFeaturesPerCell minFeaturesPerCell
#' @param maxFeaturesPerCell maxFeaturesPerCell
#' @param maxPercent.mt maxPercent.mt
#' @importFrom dplyr `%>%` filter
#' @importFrom Seurat VlnPlot FeatureScatter GetAssayData
#'
#' @return
#' @export
#'
#' @examples
#'
ScFunQC <- function(SeuratObj,
                    minCountsPerCell = NULL,
                    maxCountsPerCell = NULL,
                    minFeaturesPerCell = NULL,
                    maxFeaturesPerCell = NULL,
                    maxPercent.mt = NULL) {

  # Visualization before QC
  if (!dir.exists(paths = "./CellFunMap_output/plots")) {
    dir.create("./CellFunMap_output/plots", recursive = TRUE)
  }
  pdf(file = "./CellFunMap_output/plots/beforeQC.pdf")
  vp <- VlnPlot(SeuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  print(vp)
  plot1 <- FeatureScatter(SeuratObj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(SeuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  print(plot1 + plot2)
  dev.off()

  # Filter cells
  meta.data <- slot(object = SeuratObj, name = 'meta.data')
  if (is.null(minCountsPerCell)) {
    minCountsPerCell <- quantile(meta.data$nCount_RNA,probs=0.02)
  }
  if (is.null(maxCountsPerCell)) {
    maxCountsPerCell <- quantile(meta.data$nCount_RNA,probs=0.98)
  }
  if (is.null(minFeaturesPerCell)) {
    minFeaturesPerCell <- quantile(meta.data$nFeature_RNA,probs=0.02)
  }
  if (is.null(maxFeaturesPerCell)) {
    maxFeaturesPerCell <- quantile(meta.data$nFeature_RNA,probs=0.98)
  }
  if (is.null(maxPercent.mt)) {
    maxPercent.mt <- 20
  }
  cells <- meta.data %>% filter(nCount_RNA > minCountsPerCell, nCount_RNA < maxCountsPerCell,
                                nFeature_RNA > minFeaturesPerCell, nFeature_RNA < maxFeaturesPerCell,
                                percent.mt < maxPercent.mt) %>% rownames()
  # Filter features
  counts <- GetAssayData(object = SeuratObj)
  CellsPerGene <- Matrix::rowSums(counts > 0)
  features <- names(CellsPerGene)[CellsPerGene > ncol(counts)*0.001]

  if (length(x = cells) == 0) {
    stop("No cells found", call. = FALSE)
  }
  if (length(x = features) == 0) {
    stop("No features found", call. = FALSE)
  }
  # Subset SeuratObj according to cells and features
  SeuratObj <- subset(SeuratObj, cells = cells, features = features)

  # Visualization after QC
  pdf(file = "./CellFunMap_output/plots/afterQC.pdf")
  vp <- VlnPlot(SeuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  print(vp)
  dev.off()

  return(SeuratObj)
}














