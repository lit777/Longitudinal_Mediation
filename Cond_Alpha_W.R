Cond_Alpha_W <- function(Wc=NULL, Mp=NULL, Wp=NULL, Kp=NULL,K=NULL, a=NULL,tau, alpha=NULL, sigma=NULL, COV.a=NULL){
    
    K_other <- setdiff(1:K, unique(Kp))
    
    alpha_prop <- alpha
    
    sigma_prop <- sigma
    
    
    for(ind in K_other){
        alpha_prop[ind,] <- rmnorm(1, a, diag(tau,dim.cov+2)/10)
        alpha[ind,] <- alpha_prop[ind,]
        
        sigma_prop[ind] <- rgamma(1,sigma[ind]*100,100)
        sigma[ind] <- sigma_prop[ind]
        
    }
    
#    Kc.freq <- tabulate(Kp)
#    Kc.valid <- which(Kc.freq >= 7 )
#    Kc.invalid <- which(Kc.freq < 7 )
    
#    for(ind in Kc.invalid){
#        alpha_prop[ind,] <- rmnorm(1, a, diag(tau,dim.cov+2)/10)
#        alpha[ind,] <- alpha_prop[ind,]
#    }

    
    # Check SD of proposal distribution (!Adaptive sampler!)
    for(ind in unique(Kp)){
        COV.alpha <- 2.38^2/((dim.cov+2))*(COV.a[(5*ind-4):(5*ind),1:5]+diag(0.1^10,(dim.cov+2)))
        alpha_prop[ind,] <- rmnorm(1, alpha[ind,], COV.alpha)
        rat1 <- dmnorm(alpha_prop[ind,], a, diag(1/tau,dim.cov+2),  log=TRUE)+dmnorm(alpha[ind,], alpha_prop[ind,], COV.alpha, log=TRUE)+sum(dnorm(Wc[which(Kp==ind)], cbind(1, Mp,Wp)[which(Kp==ind),]%*%(alpha_prop[ind,]), sqrt(1/sigma[ind]), log=TRUE))
        rat2 <- dmnorm(alpha[ind,], a, diag(1/tau,dim.cov+2),  log=TRUE)+dmnorm(alpha_prop[ind,], alpha[ind,], COV.alpha, log=TRUE)+sum(dnorm(Wc[which(Kp==ind)], cbind(1, Mp,Wp)[which(Kp==ind),]%*%(alpha[ind,]), sqrt(1/sigma[ind]), log=TRUE))
        rat <- rat1 - rat2
        if(is.na(rat) | log(runif(1))>rat){
            alpha_prop[ind,] <- alpha[ind,]
        }else{
            alpha[ind,] <- alpha_prop[ind,]
        }
        sigma_prop[ind] <- rgamma(1,sigma[ind]*1000,1000)
        rat1 <- sum(dnorm(Wc[which(Kp==ind)], cbind(1,Mp,Wp)[which(Kp==ind),]%*%(alpha[ind,]), sqrt(1/sigma_prop[ind]), log=TRUE))+  dgamma(sigma_prop[ind], 10, 1, log=TRUE) + dgamma(sigma[ind], sigma_prop[ind]*1000,1000, log=TRUE)
        rat2 <- sum(dnorm(Wc[which(Kp==ind)], cbind(1,Mp,Wp)[which(Kp==ind),]%*%(alpha[ind,]), sqrt(1/sigma[ind]), log=TRUE)) + dgamma(sigma[ind], 10, 1, log=TRUE) + dgamma(sigma_prop[ind], sigma[ind]*1000,1000, log=TRUE)
        rat <- rat1 - rat2
        if(is.na(rat) | log(runif(1))>rat){
          sigma_prop[ind] <- sigma[ind]
        }else{
          sigma[ind] <- sigma_prop[ind]
        }
        
    }
    
    return(list(alpha=alpha_prop, sigma=sigma))
}


# y = exp(x)
