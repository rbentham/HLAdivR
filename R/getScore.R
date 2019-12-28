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
