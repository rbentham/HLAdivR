HLADiversityScore <- function(hla.1, hla.2,
                                diversity.measure = 'grantham',
                                exons23 = TRUE){

    if(!diversity.measure %in% c('grantham','sandberg','pdist','PAM','JTT')){
        stop('diversity measure must be either "grantham", "sandberg" or "pdist"')
    }

    hla.1 <- FixAllele(hla.1)
    hla.2 <- FixAllele(hla.2)

    # Set up score functions
    if(diversity.measure == 'pdist'){
        score.fun <- pdistanceGetScore
        dist.mat <- NULL
    }else{
        score.fun <- getScore
        dist.mat <- dist.list[[diversity.measure]]
    }

    # Get class of allele
    hla.class1 <- strsplit(hla.1,'\\*')[[1]][1]
    hla.class2 <- strsplit(hla.2,'\\*')[[1]][1]
    if(hla.class1 != hla.class2){
        warning('No alignment currently for HLA alleles of different classes')
        return(NA)
    }

    if(diversity.measure %in% c('PAM','JTT')){
        score <- protdistScore(hla.1,hla.2,hla.class1,exons23,method = diversity.measure)
        return(score)
    }

    # Get sequences for ref and two alleles
    ref.allele <- as.character(prot.list[[hla.class1]][1,1])
    ref.seq <- HLA_AllignedSeqGet(ref.allele)

    hla1.seq <- HLA_AllignedSeqGet(hla.1)
    hla2.seq <- HLA_AllignedSeqGet(hla.2)

    # If only looking at exons 2/3 trim the sequences
    if(exons23){
        hla1.seq <- ifExon23(hla1.seq, hla.class1)
        hla2.seq <- ifExon23(hla2.seq, hla.class1)
        ref.seq <- ifExon23(ref.seq, hla.class1)
    }

    if(nchar(hla1.seq) == 0 & nchar(hla2.seq) == 0) return(0)

    # Check if these are the same length - if not pad with *
    pad.list <- HLASeqPad(hla1.seq, hla2.seq)
    hla1.seq <- pad.list[[1]]
    hla2.seq <- pad.list[[2]]

    # identify locations where hla.1/2 differ from the reference
    a1 <- c('-','*','.')
    seq.d1 <- sapply(strsplit(hla1.seq, ''), function(x) (!(x %in% a1)))[,1]
    seq.d2 <- sapply(strsplit(hla2.seq, ''), function(x) (!(x %in% a1)))[,1]

    # combine these and systematically check
    all.seq.loc <- which(seq.d1 | seq.d2)

    get.loc.score <- function(seq1, seq2, loc){
        ref.amino <- substr(ref.seq,loc,loc)
        seq1.amino <- substr(seq1,loc,loc)
        seq2.amino <- substr(seq2,loc,loc)
        if(seq1.amino == '-') seq1.amino <- ref.amino
        if(seq2.amino == '-') seq2.amino <- ref.amino
        if(seq1.amino %in% c('.','*')) return(0)
        if(seq2.amino %in% c('.','*')) return(0)
        return(score.fun(seq1.amino,seq2.amino,dist.mat))
    }

    g.scores <- sapply(all.seq.loc,
                       FUN = function(x) get.loc.score(hla1.seq, hla2.seq, x))

    # Find length of alignment
    a2 <- c('*','.')
    seq.d3 <- sapply(strsplit(hla1.seq, ''), function(x) (!(x %in% a2)))[,1]
    seq.d4 <- sapply(strsplit(hla2.seq, ''), function(x) (!(x %in% a2)))[,1]
    alg.len <- sum(seq.d3 & seq.d4)

    if(length(g.scores) > 0 & alg.len != 0){
        return(sum(g.scores)/alg.len)
    }else{
        return(0)
    }
}
