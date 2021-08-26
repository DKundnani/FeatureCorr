#' Gene-Tissue Expression (GTEx) Dataset v7
#'
#' A dataset containing the expression of >42000 genes in 54 human tissues
#'
#' @format A data frame with 42136 rows and 54 columns:
#' \describe{
#'   \item{Description}{Gene names(feature identifiers)}
#' }
#' @source \url{https://www.genome.gov/Funded-Programs-Projects/Genotype-Tissue-Expression-Project}
"GTEXv7"

#' The Cancer Genome Atlas, Pan Cancer Study Dataset for selected genes
#'
#' A dataset containing the expression of
#' 40 genes in 10071 human cancer samples(tissues or cells)
#'
#' @format A data frame with 40 rows and 10071 columns:
#' \describe{
#'   \item{rownames}{Gene names(feature identifiers)}
#' }
#' @source \url{https://www.cbioportal.org/}
"TCGA40"

#' Gene-Tissue Expression (GTEx) Dataset v7
#'
#' Log 2 transformed dataset containing the expression of
#' 40 genes in 10071 human cancer samples(tissues or cells)
#'
#' @format A data frame with 40 rows and 10071 columns:
#' \describe{
#'   \item{rownames}{Gene names(feature identifiers)}
#' }
#'
"logTCGA40"
