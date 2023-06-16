#' Detect the gene ID type of Seurat object
#'
#' @param SeuratObj Seurat object
#' @importFrom magrittr `%>%`
#' @export
#'
#' @examples
#' \dontrun{
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
  for (ID in AnnotationDbi::keytypes(OrgDb)) {
    allGenelist[[ID]] <- AnnotationDbi::keys(OrgDb, keytype = ID)
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
#' @importFrom biomaRt useMart searchDatasets useDataset getBM
#' @importFrom magrittr `%>%`
#' @importFrom Seurat PercentageFeatureSet
#' @return Seurat object
#' @export
#'
#' @examples
#' \dontrun{
#' SeuratObj <- CalMTpercent(SeuratObj, by = 'use_internal_data')
#' SeuratObj <- CalMTpercent(SeuratObj, by = 'AnnotationHub')
#' SeuratObj <- CalMTpercent(SeuratObj, by = 'biomaRt')
#' }
#'
CalMTpercent <- function(SeuratObj, by = 'use_internal_data') {
  by <- match.arg(by, choices = c('use_internal_data', 'AnnotationHub', 'biomaRt'))
  if (by == 'use_internal_data') {
    # Calculate nCount and nFeature
    if (!"nCount_RNA" %in% names(slot(object = SeuratObj, name = 'meta.data'))) {
      n.calc <- caln(object = slot(object = SeuratObj, name = 'assays')[['RNA']])
      if (!is.null(x = n.calc)) {
        names(x = n.calc) <- paste(names(x = n.calc), 'RNA', sep = '_')
        slot(object = SeuratObj, name = 'meta.data')[,names(x = n.calc)] <- n.calc
      }
    }
    SeuratObj <- DetectGeneIDtype(SeuratObj)
    genes <- rownames(x = SeuratObj)
    # find mitochondrial genes
    GeneIDtype <- slot(object = SeuratObj, name = 'misc')[["featureData"]][['GeneIDtype']]
    species <- slot(object = SeuratObj, name = 'misc')[["species"]]
    if (GeneIDtype == "ENTREZID") {
      MTgenes <- genes[genes %in% txid_gene_MT$ENTREZID]
    } else {
      # gene ID transition
      OrgDb <- getOrgDb(species)
      bitrDF <- clusterProfiler::bitr(genes, fromType = GeneIDtype, toType = "ENTREZID", OrgDb = OrgDb, drop = TRUE)
      MTgenes <- bitrDF[bitrDF$ENTREZID %in% txid_gene_MT$ENTREZID, GeneIDtype]
    }
    spes <- c("Homo sapiens", "Mus musculus", "Arabidopsis thaliana")
    if (species %in% spes) {
      MTgenes <- c(MTgenes, unlist(sapply(MTGENE[[species]], function(x){genes[genes %in% x]})))
    } else {
      sss <- unique(txid_species_gene_MT[grep(species, txid_species_gene_MT$species, ignore.case = TRUE), c("ENTREZID", "SYMBOL")])
      MTgenes <- c(MTgenes, unlist(sapply(as.list(sss), function(x){genes[genes %in% x]})))
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
    for (ID in AnnotationDbi::keytypes(EnsDb)) {
      allGenelist[[ID]] <- AnnotationDbi::keys(EnsDb, keytype = ID)
    }
    # featureData[['allGene']] <- allGenelist
    # genes <- slot(object = SeuratObj, name = 'assays')[['RNA']] %>% slot(name = 'counts') %>% rownames()
    genes <- rownames(x = SeuratObj)
    GeneIDtype <- sapply(allGenelist, function(x){
      length(intersect(genes, x))
    }) %>% which.max() %>% names()
    featureData[['GeneIDtype']] <- ifelse(GeneIDtype == "GENEID", "ENSEMBL", GeneIDtype)
    # featureData[['GeneIDtype']] <- GeneIDtype

    geneTable <- AnnotationDbi::select(EnsDb, keys = genes, columns = c("GENEID", "ENTREZID", "SYMBOL", "SEQNAME", "GENEBIOTYPE", "GENENAME"), keytype = GeneIDtype)
    featureData[['geneTable']] <- geneTable
    # Calculate nCount and nFeature
    if (!"nCount_RNA" %in% names(slot(object = SeuratObj, name = 'meta.data'))) {
      n.calc <- caln(object = slot(object = SeuratObj, name = 'assays')[['RNA']])
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
    ensembl <- useMart("ensembl")
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
      n.calc <- caln(object = slot(object = SeuratObj, name = 'assays')[['RNA']])
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





#' quality control according to 'nCount_RNA', 'nFeature_RNA', 'percent.mt', 'CellsPerGene'
#'
#' @param SeuratObj Seurat object
#' @param minCountsPerCell minCountsPerCell
#' @param maxCountsPerCell maxCountsPerCell
#' @param minFeaturesPerCell minFeaturesPerCell
#' @param maxFeaturesPerCell maxFeaturesPerCell
#' @param maxPercent.mt maxPercent.mt
#' @param plot whether to plot quality control results
#'
#' @importFrom magrittr `%>%`
#' @importFrom Seurat VlnPlot FeatureScatter GetAssayData
#' @importFrom stats quantile
#'
#' @return Seurat object
#' @export
#'
#' @examples
#' \dontrun{
#' SeuratObj <- QCfun(SeuratObj)
#' }
#'
QCfun <- function(SeuratObj,
                  minCountsPerCell = NULL,
                  maxCountsPerCell = NULL,
                  minFeaturesPerCell = NULL,
                  maxFeaturesPerCell = NULL,
                  maxPercent.mt = NULL,
                  plot = FALSE) {

  # Visualization before QC
  if (plot & !dir.exists(paths = "./CellFunTopic_output/plots")) {
    dir.create("./CellFunTopic_output/plots", recursive = TRUE)
  }
  if (plot) {
    pdf(file = "./CellFunTopic_output/plots/beforeQC.pdf")
    vp <- VlnPlot(SeuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    print(vp)
    plot1 <- FeatureScatter(SeuratObj, feature1 = "nCount_RNA", feature2 = "percent.mt")
    plot2 <- FeatureScatter(SeuratObj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
    print(plot1 + plot2)
    dev.off()
  }

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
  if ('percent.mt' %in% colnames(meta.data)) {
    cells <- meta.data %>% dplyr::filter(nCount_RNA > minCountsPerCell, nCount_RNA < maxCountsPerCell,
                                         nFeature_RNA > minFeaturesPerCell, nFeature_RNA < maxFeaturesPerCell,
                                         percent.mt < maxPercent.mt) %>% rownames()
  } else {
    cells <- meta.data %>% dplyr::filter(nCount_RNA > minCountsPerCell, nCount_RNA < maxCountsPerCell,
                                         nFeature_RNA > minFeaturesPerCell, nFeature_RNA < maxFeaturesPerCell) %>% rownames()
  }

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
  if (plot) {
    pdf(file = "./CellFunTopic_output/plots/afterQC.pdf")
    vp <- VlnPlot(SeuratObj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
    print(vp)
    dev.off()
  }

  return(SeuratObj)
}














