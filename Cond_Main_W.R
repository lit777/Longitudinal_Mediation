Cond_Main_W <- function(Wc=NULL,  Mp=NULL,WmisIndex=NULL, Wp=NULL, K=NULL,Kp=Kp, a=NULL, a_prior=NULL,tau=NULL,tau_prior=NULL, alpha=NULL,   sigma=NULL,lambda=NULL,Pp=NULL,dim.cov=NULL,COV.w=NULL, COV.a=NULL){
  
  #---- updating assignment vectors
  update0 <- Cond_K_W(Wc,Mp,Wp,K,Pp,alpha,sigma) #<-- Kc
  
  #---- updating Alpha
  update1 <- Cond_Alpha_W(Wc,Mp,Wp,update0$Kc,K,a,tau,alpha,sigma,COV.a) #<-- alpha
  
  #---- updating weights and lambda
  update2 <- Cond_P(update0$Kc,K,lambda) #<-- Pc and Lambda
  
  
  #---- Sampling Base Measure
  # Check SD of the propsal distribution (!Adaptive sampler!)
  a_prop <- NULL
  a_acc <- NULL
  COV.cov <- 2.38^2/((dim.cov+2))*(COV.w+diag(0.1^10,(dim.cov+2))) #<-- dim.cov = the dimension of the covariates.
  a_prop <- c(rmnorm(1, a, COV.cov))
  
  # a_prior$P, a_prior$mean, a_prior$sigma
  rat1 <- log(sum(a_prior$P*dmnorm(a_prop, a_prior$mean, a_prior$sigma))) + sum(sapply(1:K, function(x) sum(dnorm(update1$alpha[x,], a_prop, sqrt(1/tau),  log=TRUE)))) + dmnorm(a, a_prop, COV.cov, log=TRUE)
  rat2 <- log(sum(a_prior$P*dmnorm(a, a_prior$mean, a_prior$sigma))) + sum(sapply(1:K, function(x) sum(dnorm(update1$alpha[x,], a, sqrt(1/tau),  log=TRUE)))) + dmnorm(a_prop, a, COV.cov, log=TRUE)

  rat <- rat1 - rat2
  if(is.na(rat) | log(runif(1))>rat){
    a_prop <- a
  }else{
    a <- a_prop
  }
  
  
  
  tau_prop <- tau
  tau_acc <- NULL
  
  for(i in 1:(dim.cov+2)){
    tau_prop[i] <- rgamma(1,tau[i]*100,100)
    
    rat1 <-  sum(dgamma(tau_prop[i], 2, 1, log=TRUE))+ sum(sapply(1:K, function(x) sum(dnorm(update1$alpha[x,], a, sqrt(1/tau_prop),  log=TRUE)))) + sum(dgamma(tau[i], tau_prop[i]*100, 100, log=TRUE))
    rat2 <- sum(dgamma(tau[i], 2, 1, log=TRUE)) +  sum(sapply(1:K, function(x) sum(dnorm(update1$alpha[x,], a, sqrt(1/tau),  log=TRUE))))  + sum(dgamma(tau_prop[i], tau[i]*100, 100, log=TRUE))
    rat <- rat1 - rat2
    if(is.na(rat) | log(runif(1))>rat){
      tau_prop[i] <- tau[i]
      tau_acc <- 0
    }else{
      tau[i] <- tau_prop[i]
      tau_acc <- 1
    }
    
  }
  

  

 Wc[WmisIndex] <- rnorm(length(WmisIndex), diag(cbind(1, Mp,Wp)[WmisIndex,]%*%t(update1$alpha[update0$Kc[WmisIndex],])), sqrt(1/update1$sigma[update0$Kc[WmisIndex]]))
 

  
  return(list(Wc=Wc,alpha=update1$alpha,sigma=update1$sigma,  tau=tau,  a=a, Pc=update2$Pc, lambda=update2$lambda, Kc=update0$Kc))
}



