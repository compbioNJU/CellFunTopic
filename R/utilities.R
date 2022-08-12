
#' Get OrgDb of a species
#'
#' @param species such as "Homo sapiens"
#'
#'
getOrgDb <- function(species) {
  spes <- c("Anopheles gambiae",
            "Arabidopsis thaliana",
            "Bos taurus",
            "Canis lupus familiaris",
            "Caenorhabditis elegans",
            "Gallus gallus",
            "Pan troglodytes",
            "Streptomyces coelicolor",
            "Escherichia coli strain K12",
            "Escherichia coli strain Sakai",
            "Drosophila melanogaster",
            "Toxoplasma gondii",
            "Homo sapiens",
            "Plasmodium falciparum",
            "Mus musculus",
            "Sus scrofa",
            "Rattus norvegicus",
            "Macaca mulatta",
            "Xenopus laevis",
            "Saccharomyces cerevisiae",
            "Danio rerio",
            "Myxococcus xanthus")

  if (species %in% spes) {
    OrgDb <- requireOrgDb(species)
  } else {
    OrgDb <- getOrgDb_AH(species)
  }
  return(OrgDb)
}

#' Get OrgDb of a species from package AnnotationHub
#'
#' @param species species name
#'
#' @importFrom AnnotationHub AnnotationHub query
#'
#'
getOrgDb_AH <- function(species) {

  if (!requireNamespace('AnnotationHub', quietly = TRUE)) {
    stop("Please install AnnotationHub from Bioconductor at first")
  }
  ah <- AnnotationHub()
  if (!(species %in% unique(ah$species))) {
    stop(sprintf("Sorry, the species %s is not supported by AnnotationHub package", species))
  }
  queryResult <- query(ah, pattern = c(species, "OrgDb"), ignore.case = TRUE)
  OrgDb <- ah[[names(queryResult@.db_uid)]]
  return(OrgDb)
}


#' Get EnsDb of a species from package AnnotationHub
#'
#' @param species species name
#'
#' @importFrom AnnotationHub AnnotationHub query
#'
#'
getEnsDb_AH <- function(species) {
  if (!requireNamespace('AnnotationHub', quietly = TRUE)) {
    stop("Please install AnnotationHub from Bioconductor at first")
  }
  ah <- AnnotationHub()
  if (!(species %in% unique(ah$species))) {
    stop(sprintf("Sorry, the species %s is not supported by AnnotationHub package", species))
  }
  queryResult <- query(ah, pattern = c(species, "EnsDb"), ignore.case = TRUE)
  accs <- names(queryResult@.db_uid)
  EnsDb <- ah[[accs[length(accs)]]]
  return(EnsDb)
}


#' require OrgDb Package of a species
#'
#' @param species species name
#'
#'
requireOrgDb <- function(species) {
  OrgDb <- org2db(species)
  suppressPackageStartupMessages(require(OrgDb, character.only = TRUE))
  OrgDb <- eval(parse(text=OrgDb))
  return(OrgDb)
}

#' get OrgDb name according to species
#'
#' @param species species name
#'
#'
org2db <- function(species) {
  OrgDb <- switch(species,
                  `Anopheles gambiae`             = "org.Ag.eg.db",
                  `Arabidopsis thaliana`          = "org.At.tair.db",
                  `Bos taurus`                    = "org.Bt.eg.db",
                  `Canis lupus familiaris`        = "org.Cf.eg.db",
                  `Caenorhabditis elegans`        = "org.Ce.eg.db",
                  `Gallus gallus`                 = "org.Gg.eg.db",
                  `Pan troglodytes`               = "org.Pt.eg.db",
                  `Streptomyces coelicolor`       = "org.Sco.eg.db",
                  `Escherichia coli strain K12`   = "org.EcK12.eg.db",
                  `Escherichia coli strain Sakai` = "org.EcSakai.eg.db",
                  `Drosophila melanogaster`       = "org.Dm.eg.db",
                  `Toxoplasma gondii`             = "org.Tgondii.eg.db",
                  `Homo sapiens`                  = "org.Hs.eg.db",
                  `Plasmodium falciparum`         = "org.Pf.plasmo.db",
                  `Mus musculus`                  = "org.Mm.eg.db",
                  `Sus scrofa`                    = "org.Ss.eg.db",
                  `Rattus norvegicus`             = "org.Rn.eg.db",
                  `Macaca mulatta`                = "org.Mmu.eg.db",
                  `Xenopus laevis`                = "org.Xl.eg.db",
                  `Saccharomyces cerevisiae`      = "org.Sc.sgd.db",
                  `Danio rerio`                   = "org.Dr.eg.db",
                  `Myxococcus xanthus`            = "org.Mxanthus.db"
  )
  return(OrgDb)
}

#' get 'organism' that ReactomePA::gsePathway needs.
#'
#' change latin name of species to trivial name
#' @param species species name
#'
#'
speToOrg <- function(species) {
  species <- gsub(" ", "_", species)
  organism <- switch(species,
                  Anopheles_gambiae             = "anopheles",
                  Arabidopsis_thaliana          = "arabidopsis",
                  Bos_taurus                    = "bovine",
                  Canis_lupus_familiaris        = "canine",
                  Caenorhabditis_elegans        = "celegans",
                  Gallus_gallus                 = "chicken",
                  Pan_troglodytes               = "chimp",
                  Streptomyces_coelicolor       = "coelicolor",
                  Escherichia_coli_strain_K12   = "ecolik12",
                  Escherichia_coli_strain_Sakai = "ecsakai",
                  Drosophila_melanogaster       = "fly",
                  Toxoplasma_gondii             = "gondii",
                  Homo_sapiens                  = "human",
                  Plasmodium_falciparum         = "malaria",
                  Mus_musculus                  = "mouse",
                  Sus_scrofa                    = "pig",
                  Rattus_norvegicus             = "rat",
                  Macaca_mulatta                = "rhesus",
                  Xenopus_laevis                = "xenopus",
                  Saccharomyces_cerevisiae      = "yeast",
                  Danio_rerio                   = "zebrafish"
  )
  return(organism)
}


#' Get scientific name abbreviate(kegg_code) used in the 'organism' parameter of gseKEGG function
#'
#' @param species species name
#'
#' @importFrom clusterProfiler search_kegg_organism
#'
#'
get_kegg_code <- function(species) {
  result <- search_kegg_organism(species, by='scientific_name')
  if (nrow(result) == 0) {
    stop("Can't find ", species, "'s kegg_code, please search on http://www.genome.jp/kegg/catalog/org_list.html")
  }
  kegg_code <- result[1,1]
  return(kegg_code)
}


#' Get pathway IDs and their corresponding description
#'
#' @param SeuratObj Object of class "Seurat"
#' @param by "GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ID2Description(SeuratObj, by = "GO")
#' ID2Description(SeuratObj, by = "KEGG")
#' }
#'
ID2Description <- function(SeuratObj, by = "GO") {
  by <- match.arg(by, choices = c("GO", "KEGG", "Reactome", "MSigDb", "WikiPathways", "DO", "NCG", "DGN"))
  GSEAresult <- slot(object = SeuratObj, name = 'misc')[[paste0("GSEAresult_", by)]]
  pws <- GSEAresult %>% dplyr::select(ID, Description) %>% unique %>% tibble::deframe()
  return(pws)
}


# Calculate nCount and nFeature
#
# @param object An Assay object
caln <- function(object) {
  return(list(
    nCount = Matrix::colSums(x = GetAssayData(object = object, slot = 'counts')),
    nFeature = Matrix::colSums(x = GetAssayData(object = object, slot = 'counts') > 0)
  ))
}











