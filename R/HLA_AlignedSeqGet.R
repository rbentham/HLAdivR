#' Method for extracting aligned protein sequences for HLA genes
#'
#' Upon the aligned protein sequences from IMGT produce the protein sequence
#' for the given HLA allele
#'
#' @param hla.allele HLA allele for genes HLA-A, HLA-B or HLA-C, in digit format A:01:01
#' @param full logical of whether to output full sequence or sequence in relation to default HLA allele
#' @return Aligned protein sequence of HLA allele
#' @name HLA_AlignedSeqGet
#' @export

HLA_AlignedSeqGet <- function(hla.allele, full=FALSE){

    hla.allele <- FixAllele(hla.allele)
    hla.class <- strsplit(hla.allele,'\\*')[[1]][1]

    #data("aligned_HLA_seq")
    allign.data <- aligned_HLA_seq[[hla.class]]

    loc1 <- which(allign.data$X1 == hla.allele)
    hla.seq <- as.character(t(as.matrix(allign.data[loc1,-1])))

    # remove nas
    hla.seq <- hla.seq[!is.na(hla.seq)]
    # join
    hla.seq <- paste(hla.seq, collapse = '')
    if(full){
        ref.allele <- as.character(aligned_HLA_seq[[hla.class]][1,1])
        ref.seq <- HLA_AlignedSeqGet(ref.allele)

        # replace all '-' with ref.seq
        hla.split <- strsplit(hla.seq,'')[[1]]
        ref.split <- strsplit(ref.seq,'')[[1]]
        hla.split[which(hla.split == '-')] <- ref.split[which(hla.split == '-')]
        hla.seq <- paste(hla.split,collapse = '')
    }

    return(hla.seq)
}
