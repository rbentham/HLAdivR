HLA_AllignedSeqGet <- function(hla.allele, full=FALSE){

    hla.allele <- FixAllele(hla.allele)
    hla.class <- strsplit(hla.allele,'\\*')[[1]][1]

    allign.data <- prot.list[[hla.class]]

    loc1 <- which(allign.data$X1 == hla.allele)
    hla.seq <- as.character(t(as.matrix(allign.data[loc1,-1])))

    # remove nas
    hla.seq <- hla.seq[!is.na(hla.seq)]
    # join
    hla.seq <- paste(hla.seq, collapse = '')
    if(full){
        ref.allele <- as.character(prot.list[[hla.class]][1,1])
        ref.seq <- HLA_AllignedSeqGet(ref.allele)

        # replace all '-' with ref.seq
        hla.split <- strsplit(hla.seq,'')[[1]]
        ref.split <- strsplit(ref.seq,'')[[1]]
        hla.split[which(hla.split == '-')] <- ref.split[which(hla.split == '-')]
        hla.seq <- paste(hla.split,collapse = '')
    }

    return(hla.seq)
}
