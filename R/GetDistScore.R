#' Method for calculating the distance between HLA alleles from a pre-defined HLA distance matrix
#'
#' Given two HLA alleles and a distance matrix output the distance between the two alleles
#' for the given HLA allele
#'
#' @param hla1 HLA allele for genes HLA-A, HLA-B or HLA-C, in digit format A:01:01
#' @param hla2 HLA allele for genes HLA-A, HLA-B or HLA-C, in digit format A:01:01
#' @param dist.mat distance matrix giving distance between HLA alleles
#' @return Distance between two HLA alleles
#' @name GetDistScore
#' @export

GetDistScore <- function(hla1,hla2,dist.mat) {
    if(length(hla1) != length(hla2)) stop('Different lenghts of HLA alleles')
    hla1 <- gsub('HLA-','',hla1)
    hla2 <- gsub('HLA-','',hla2)
    a1 <- sapply(hla1, FUN = function(x) which(row.names(dist.mat) == x))
    a2 <- sapply(hla2, FUN = function(x) which(row.names(dist.mat) == x))
    d1 <- sapply(seq_len(length(a1)), FUN = function(x) {
        d1.1 <- dist.mat[a1[[x]], a2[[x]]]
        ifelse(length(d1.1) > 0, d1.1, NA)
    })
    return(d1)
}
