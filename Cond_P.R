Cond_P <- function(Kc=NULL, K=NULL, lambda=NULL){

    V <- rep(0,K-1)
    for(k in 1:(K-1)){
        V[k] <- rbeta(1, 1+length(which(Kc==k)), lambda+sum(length(which(Kc>k))))
    }

    Pc=NULL
    Pc_temp=NULL
    Pc_temp[1]=1
    Pc[1]=V[1]
    for(k in 2:(K-1)){
        Pc_temp[k]=Pc_temp[k-1]*(1-V[k-1])
        Pc[k]=V[k]*Pc_temp[k]
    }
    Pc[K]=Pc_temp[K-1]*(1-V[K-1])

    ## update lambda (when lambda follows Gamma(1, 1))
    lambda <- rgamma(1, 1+K-1, 1-sum(log(1-V)))

    return(list(Pc=Pc, lambda=lambda))
}

