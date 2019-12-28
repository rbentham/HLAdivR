ifExon23 <- function(seq,hla.class){
    seq.new <- seq
    seq.len <- nchar(seq.new)

    seq.end.df <- data.frame(start = c(26,26,26), end = c(220, 219, 218),
                             row.names = c('A','B','C'))

    end.loc <- min(seq.end.df[hla.class,'end'], seq.len)

    seq.new <- substr(seq, seq.end.df[hla.class,'start'], end.loc)
    return(seq.new)
}
