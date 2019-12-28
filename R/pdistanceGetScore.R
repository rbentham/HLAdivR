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
