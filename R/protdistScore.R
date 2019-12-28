# path to protdist must be supplied
protdistScore <- function(hla.1, hla.2, hla.class1, exons23, method='PAM'){
    protdist.path <- '~/Documents/Programs/phylip-3.695/exe/protdist.app/Contents/MacOS'
    hla1.seq <- HLA_AllignedSeqGet(hla.1, full = TRUE)
    hla2.seq <- HLA_AllignedSeqGet(hla.2, full = TRUE)

    # If only looking at exons 2/3 trim the sequences
    if(exons23){
        hla1.seq <- ifExon23(hla1.seq, hla.class1)
        hla2.seq <- ifExon23(hla2.seq, hla.class1)
    }

    # Check if these are the same length - if not pad with *
    pad.list <- HLASeqPad(hla1.seq, hla2.seq)
    hla1.seq <- pad.list[[1]]
    hla2.seq <- pad.list[[2]]

    hla.mat <- rbind(strsplit(hla1.seq,'')[[1]],strsplit(hla2.seq,'')[[1]])
    hla.phy <- phyDat(hla.mat, type = 'AA')
    hla.ps <- as.proseq(hla.phy)
    hla.score <-as.numeric(Rprotdist(path = protdist.path,
                                     hla.ps,
                                     model = method,quiet = TRUE))
    return(hla.score)
}
