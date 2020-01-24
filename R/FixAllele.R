#' Internal function to change HLA allele name to correct format
#'
#' Removes leading HLA and addition digits and checks the formatting.
#' Works for some common alternative HLA formats.
#'
#' @param allele HLA allele for genes HLA-A, HLA-B or HLA-C
#' @return HLA allele in correct format for HLADivR to function
#' @name GetDistScore

FixAllele <- function(allele){
    # function to get allele in right format
    # 1. remove begining "HLA-"
    allele.fix <- gsub('HLA-','',allele)

    # 1a. check that a '*' separates HLA class and rest
    star.loc <- grep(pattern = '[0-9]',strsplit(allele.fix,'*')[[1]])[1]-1
    if(substr(allele.fix,star.loc,star.loc)!='*'){
        allele.fix <- paste0(substr(allele.fix,1,star.loc),'*',
                             substr(allele.fix,star.loc+1,nchar(allele.fix)))
    }

    # 1b remove any remaining '-'
    allele.fix <- gsub('-','',allele.fix)

    HLA.names <- names(aligned_HLA_seq)

    # 2. check that first character is now a known HLA
    hla.class <- strsplit(allele.fix,'\\*')[[1]][1]
    if(!hla.class %in% HLA.names) {
        warning('HLA class not recognised')
        return(NA)}


    # 4. If no match, might need to specify more information - extend
    all.allele <- unique(aligned_HLA_seq[[hla.class]]$X1)
    loc1 <- which(all.allele == allele.fix)
    if(length(loc1) == 0){
        hla.split <- strsplit(allele.fix,':')[[1]]
        if(length(hla.split) >= 2){
            allele.base <- paste(hla.split[1], hla.split[2],sep = ':')
            # Choose first sequence that matches
            hla.grep <- grep(gsub('\\*','\\\\*',allele.base), all.allele)
            if(length(hla.grep) > 0){
                allele.fix <-all.allele[hla.grep[1]]
            }else{
                warning('Allele given was not matched')
                return(NA)
            }
        }else{
            warning(paste('Allele given was not matched (only A,B or C given)'))
            return(NA)
        }
    }
    return(allele.fix)
}
