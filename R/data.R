#' Grantham distance matrix
#'
#' A distance matrix containing the Grantham distance between amino acids.
#' Grantham, R. (1974-09-06). "Amino acid difference formula to help explain protein evolution". Science.
#'
#' @format A data frame with 19 rows and 20 columns containg amino acid abreviations and distance between them.
#' @source \url{https://en.wikipedia.org/wiki/Amino_acid_replacement#Grantham's_distance}
#' @return NA
"grantham"

#' Sandberg distance matrix
#'
#' A distance matrix containing the Sandberg distance between amino acids.
#' Sandberg, Maria, et al. "New chemical descriptors relevant for the design of biologically active peptides. A multivariate characterization of 87 amino acids." Journal of medicinal chemistry 41.14 (1998): 2481-2491.
#'
#' @format A data frame with 19 rows and 20 columns containg amino acid abreviations and distance between them.
#' @return NA
"sandberg_dist"


#' Aligned HLA sequences
#'
#' A list containing all the aligned HLA sequences downloaded from the IMGT database downloaded 3rd October 2019.
#' Robinson, James, et al. "IMGT/HLA and IMGT/MHC: sequence databases for the study of the major histocompatibility complex." Nucleic acids research 31.1 (2003): 311-314.
#'
#' @format A list containing aligned sequences for all alleles for HLA class I and II genes
#' @source \url{ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/alignments/}
#' @return NA
"aligned_HLA_seq"

#' Chowell 2019 Clinical data
#'
#' Clinical data from the Chowell 2019 paper including HLA diversity as measured using the Grantham method.
#' Chowell, Diego, et al. "Evolutionary divergence of HLA class I genotype impacts efficacy of cancer immunotherapy." Nature medicine 25.11 (2019): 1715-1720.
#'
#' @format A data frame with 314 rows and 9 columns containg HLA genotypes for patiens and calculated mean Grantham scores
#' @source \url{https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-019-0639-4/MediaObjects/41591_2019_639_MOESM3_ESM.xlsx}
#' @return NA
"Chowell_2019_clinical"

