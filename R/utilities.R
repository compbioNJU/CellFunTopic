
#' Get OrgDb of a species
#'
#' @param species such as "Homo sapiens"
#'
#' @return
#' @export
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
#' @param species
#'
#' @importFrom AnnotationHub AnnotationHub query
#'
#' @return
#' @export
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
#' @param species
#'
#' @importFrom AnnotationHub AnnotationHub query
#'
#' @return
#' @export
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
#' @param species
#'
#' @return
#' @export
#'
requireOrgDb <- function(species) {
  OrgDb <- org2db(species)
  require(OrgDb, character.only = TRUE)
  OrgDb <- eval(parse(text=OrgDb))
  return(OrgDb)
}

#' get OrgDb name according to species
#'
#' @param species
#'
#' @return
#' @export
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
#' change latin name of species to trivial name
#' @param species
#'
#' @return
#' @export
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
#' @param species
#'
#' @importFrom clusterProfiler search_kegg_organism
#'
#' @return
#' @export
#'
get_kegg_code <- function(species) {
  result <- search_kegg_organism(species, by='scientific_name')
  if (nrow(result) == 0) {
    stop("Can't find ", species, "'s kegg_code, please search on http://www.genome.jp/kegg/catalog/org_list.html")
  }
  kegg_code <- result[1,1]
  return(kegg_code)
}


