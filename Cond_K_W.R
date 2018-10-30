Cond_K_W <- function(Wc=NULL, Mp=NULL,Wp=NULL, K=NULL, Pp=NULL, alpha=NULL, sigma=NULL){
  

  
  Pc.pre <- sapply(1:K, function(k) Pp[k]*dnorm(Wc, cbind(1,Mp,Wp)%*%(alpha[k,]) , sqrt(1/sigma[k])))
  
  zero_index <- which(rowSums(Pc.pre)==0 | is.na(rowSums(Pc.pre)))
  Pc <- Pc.pre
  Pc[zero_index,] <- 1
  
  Kc <- apply(Pc,1, function(x) sample(1:20,1, prob=x))
  #  Kc.freq <- tabulate(Kc)
  #  Kc.col <- which(Kc.freq >= 8 )
  #  Kc.ind <- which(Kc %in% which(Kc.freq < 8 ))
  
  
  #  if(length(Kc.col)!=10 & length(Kc.ind)!=0){
  #      Pc.post <- Pc[Kc.ind,Kc.col]
  #      if(length(Kc.ind)==1){if(length(Kc.col)!=1){Kc[Kc.ind] <- sample(Kc.col, 1, prob=Pc.post)}else{Kc[Kc.ind] <- Kc.col} }else{
  #          if(length(Kc.col)!=1){
  #              Kc[Kc.ind] <- apply(Pc.post,1, function(x) sample(Kc.col,1, prob=x))}else{Kc[Kc.ind] <- Kc.col }}}
  
  return(list(Kc=Kc))
}

