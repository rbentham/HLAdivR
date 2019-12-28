# HLA.names needs to be defined in function of from internal data
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

    # 2. check that first character is now a known HLA
    hla.class <- strsplit(allele.fix,'\\*')[[1]][1]
    if(!hla.class %in% HLA.names) stop('HLA class not recognised')


    # 4. If no match, might need to specify more information - extend
    all.allele <- unique(prot.list[[hla.class]]$X1)
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
                stop('Allele given was not matched')
            }
        }else{
            stop(paste('Allele given was not matched (only A,B or C given)'))
        }
    }
    return(allele.fix)
}
