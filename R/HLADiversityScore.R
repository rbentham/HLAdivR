#' Method to calculate distance between HLA alleles using different techniques
#'
#' These functions calculate the HLA diversity between alleles using several alternative methods.
#'
#'
#' @param hla.1 HLA allele in format A*01:01
#' @param hla.2 HLA allele in format A*01:01
#' @param diversity.measure One of 'grantham', 'sandberg', 'pdist', 'PAM' or 'JTT'.
#' @param exons23 Logical of whether to only calculate the diversity of exons 2 and 3 the protein binding region of the HLA alleles.
#' @param protdist.path Path to protdist app from phylip e.g. '~/Documents/Programs/phylip-3.695/exe/protdist.app/Contents/MacOS'
#' @return HLA diversity metric
#' @name HLADiversityScore

NULL

#' @export
#' @rdname HLADiversityScore
HLADiversityScore <- function(hla.1, hla.2,
                                diversity.measure = 'grantham',
                                exons23 = TRUE,
                                protdist.path = NULL){

    dist.list <- list()
    dist.list[['grantham']] <- grantham
    dist.list[['sandberg']] <- sandberg_dist

    if(!diversity.measure %in% c('grantham','sandberg','pdist','PAM','JTT')){
        stop('diversity measure must be either "grantham", "sandberg" or "pdist"')
    }

    hla.1 <- FixAllele(hla.1)
    if(is.na(hla.1)) return(NA)
    hla.2 <- FixAllele(hla.2)
    if(is.na(hla.2)) return(NA)

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
        if(is.null(protdist.path)) stop('Path to protdist must be supplied to function')
        score <- protdistScore(hla.1,hla.2,hla.class1,exons23,method = diversity.measure, protdist.path)
        return(score)
    }

    # Get sequences for ref and two alleles
    ref.allele <- as.character(aligned_HLA_seq[[hla.class1]][1,1])
    ref.seq <- HLA_AlignedSeqGet(ref.allele)

    hla1.seq <- HLA_AlignedSeqGet(hla.1)
    hla2.seq <- HLA_AlignedSeqGet(hla.2)

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


# Function to get score from a given protein-protein distance matrix
getScore <- function(a,b,dist.mat){
    # Check amino acid is actually allowed
    if(!(a %in% c(dist.mat$FIRST,'W','X'))){
        stop(paste('Amino acid',a, 'is not recognised'))
    }
    if(!(b %in% c(dist.mat$FIRST,'W','X'))){
        stop(paste('Amino acid', b, 'is not recognised'))
    }
    if(a == b){
        return(0)
    }
    if(a == 'X' | b == 'X'){
        return(0)
    }
    a.loc <- which(dist.mat$FIRST == a)
    b.loc <- which(dist.mat$FIRST == b)
    if(length(a.loc) == 0) a.loc <- 20
    if(length(b.loc) == 0) b.loc <- 20
    if(a.loc < b.loc){
        return(as.numeric(dist.mat[a.loc, b.loc]))
    } else if(a.loc > b.loc){
        return(as.numeric(dist.mat[b.loc, a.loc]))
    }
}

# Function to calculate the binary pdistance
pdistanceGetScore <- function(a,b,dist.mat=NULL){
    # Check amino acid is actually allowed
    if(!(a %in% c(grantham$FIRST,'W','X'))){
        stop(paste('Amino acid',a, 'is not recognised'))
    }
    if(!(b %in% c(grantham$FIRST,'W','X'))){
        stop(paste('Amino acid', b, 'is not recognised'))
    }
    if(a == b){
        return(0)
    }else{
        return(1)
    }
}


# Function to select exons 2 and 3
ifExon23 <- function(seq,hla.class){
    seq.new <- seq
    seq.len <- nchar(seq.new)

    seq.end.df <- data.frame(start = c(26,26,26), end = c(220, 219, 218),
                             row.names = c('A','B','C'))

    end.loc <- min(seq.end.df[hla.class,'end'], seq.len)

    seq.new <- substr(seq, seq.end.df[hla.class,'start'], end.loc)
    return(seq.new)
}

# Function to calculate the PAM or JTT scores
protdistScore <- function(hla.1, hla.2, hla.class1, exons23, method='PAM',
                          protdist.path = '~/Documents/Programs/phylip-3.695/exe/protdist.app/Contents/MacOS'){
    # protdist.path <- '~/Documents/Programs/phylip-3.695/exe/protdist.app/Contents/MacOS'
    hla1.seq <- HLA_AlignedSeqGet(hla.1, full = TRUE)
    hla2.seq <- HLA_AlignedSeqGet(hla.2, full = TRUE)

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
    hla.phy <- phangorn::phyDat(hla.mat, type = 'AA')
    hla.ps <- Rphylip::as.proseq(hla.phy)
    hla.score <-as.numeric(Rphylip::Rprotdist(path = protdist.path,
                                              hla.ps,
                                              model = method,quiet = TRUE))
    return(hla.score)
}

# Function to pad HLA sequences to same length if requierd
HLASeqPad <- function(hla1.seq,hla2.seq){
    hla.len1 <- nchar(hla1.seq)
    hla.len2 <- nchar(hla2.seq)
    if(hla.len1 != hla.len2){
        max.len <- max(hla.len1, hla.len2)
        long.hla <- c(hla1.seq, hla2.seq)[which.max(c(hla.len1, hla.len2))]
        min.len <- min(hla.len1, hla.len2)
        short.hla <- c(hla1.seq, hla2.seq)[which.min(c(hla.len1, hla.len2))]

        short.hla <- paste(short.hla, paste(rep('*',max.len-min.len),
                                            collapse = '',sep = ''),
                           collapse = '', sep = '')
        if(hla.len1 == min.len){
            hla1.seq <- short.hla
        }else{
            hla2.seq <- short.hla
        }
    }
    return(list(hla1.seq,hla2.seq))
}

