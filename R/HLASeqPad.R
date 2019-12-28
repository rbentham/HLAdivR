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
