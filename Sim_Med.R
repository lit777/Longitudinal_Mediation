## retrieve arguments passed from the command line. (must run in parallel)
process <- as.integer(as.character(commandArgs(trailingOnly = TRUE)))

##---------------------------------------------------------------
## Required libraries
##---------------------------------------------------------------
library(pscl)
library(mnormt)
library(DPpackage)
library(data.table)
library(MCMCpack)

source("Cond_P.R")

Cond_K <- function(Mc=NULL, Mp=NULL, Wp=NULL,K=NULL, Pp=NULL, alpha=NULL, sigma=NULL){
    Pc.pre <- sapply(1:K, function(k) Pp[k]*dnorm(Mc, cbind(1,Mp,Wp)%*%(alpha[k,]), sqrt(1/sigma[k])))
    zero_index <- which(rowSums(Pc.pre)==0 | is.na(rowSums(Pc.pre)))
    Pc <- Pc.pre
    Pc[zero_index,] <- 1
    Kc <- apply(Pc,1, function(x) sample(1:20,1, prob=x))
    return(list(Kc=Kc))
}

Cond_Alpha <- function(Mc=NULL, Mp=NULL,  Wp=NULL, Kp=NULL,K=NULL, a=NULL,tau=NULL, alpha=NULL, sigma=NULL, COV.a=NULL){
    
    K_other <- setdiff(1:K, unique(Kp))
    alpha_prop <- alpha
    alpha_acc <- NULL
    sigma_prop <- sigma
    for(ind in K_other){
        alpha_prop[ind,] <- rmnorm(1, a, diag(1/tau,dim.cov+2)/10)
        alpha[ind,] <- alpha_prop[ind,]
        sigma_prop[ind] <- rgamma(1,sigma[ind]*100,100)
        sigma[ind] <- sigma_prop[ind]
    }
    
    # Check SD of proposal distribution (!Adaptive sampler!)
    for(ind in unique(Kp)){
        COV.alpha <- 2.38^2/((dim.cov+2)^2)*(COV.a[(5*ind-4):(5*ind),1:5]+diag(0.1^10,(dim.cov+2)))
        alpha_prop[ind,] <- rmnorm(1, alpha[ind,], COV.alpha)
        rat1 <- dmnorm(alpha_prop[ind,], a, diag(1/tau,dim.cov+2),  log=TRUE)+dmnorm(alpha[ind,], alpha_prop[ind,], COV.alpha, log=TRUE)+sum(dnorm(Mc[which(Kp==ind)], cbind(1, Mp,Wp)[which(Kp==ind),]%*%(alpha_prop[ind,]), sqrt(1/sigma[ind]), log=TRUE))
        rat2 <- dmnorm(alpha[ind,], a, diag(1/tau,dim.cov+2),  log=TRUE)+dmnorm(alpha_prop[ind,], alpha[ind,], COV.alpha, log=TRUE)+sum(dnorm(Mc[which(Kp==ind)], cbind(1, Mp,Wp)[which(Kp==ind),]%*%(alpha[ind,]), sqrt(1/sigma[ind]), log=TRUE))
        rat <- rat1 - rat2
        if(is.na(rat) | log(runif(1))>rat){
            alpha_prop[ind,] <- alpha[ind,]
        }else{
            alpha[ind,] <- alpha_prop[ind,]
        }
        sigma_prop[ind] <- rgamma(1,sigma[ind]*1000,1000)
        rat1 <- sum(dnorm(Mc[which(Kp==ind)], cbind(1,Mp,Wp)[which(Kp==ind),]%*%(alpha[ind,]), sqrt(1/sigma_prop[ind]), log=TRUE))+  dgamma(sigma_prop[ind], 10, 1, log=TRUE) + dgamma(sigma[ind], sigma_prop[ind]*1000,1000, log=TRUE)
        rat2 <- sum(dnorm(Mc[which(Kp==ind)], cbind(1,Mp,Wp)[which(Kp==ind),]%*%(alpha[ind,]), sqrt(1/sigma[ind]), log=TRUE)) + dgamma(sigma[ind], 10, 1, log=TRUE) + dgamma(sigma_prop[ind], sigma[ind]*1000,1000, log=TRUE)
        rat <- rat1 - rat2
        if(is.na(rat) | log(runif(1))>rat){
          sigma_prop[ind] <- sigma[ind]
        }else{
          sigma[ind] <- sigma_prop[ind]
        }
    }
    return(list(alpha=alpha_prop, sigma=sigma))
}

Cond_Main <- function(Mc=NULL, MmisIndex=NULL, Mp=NULL,Yp=NULL, Wp=NULL, K=NULL,Kp=Kp, a=NULL, a_prior=NULL,tau=NULL,tau_prior=NULL, alpha=NULL,   sigma=NULL,lambda=NULL,Pp=NULL,dim.cov=NULL,COV.med=NULL, COV.a=NULL){
    
    #---- updating assignment vectors
    update0 <- Cond_K(Mc,Mp,Wp,K,Pp,alpha,sigma) #<-- Kc
    
    #---- updating Alpha
    update1 <- Cond_Alpha(Mc,Mp,Wp,update0$Kc,K,a,tau,alpha,sigma,COV.a) #<-- alpha
    
    #---- updating weights and lambda
    update2 <- Cond_P(update0$Kc,K,lambda) #<-- Pc and Lambda
    
    #---- Sampling Base Measure
    # Check SD of the propsal distribution (!Adaptive sampler!)
    a_prop <- NULL
    a_acc <- NULL
    COV.cov <- 2.38^2/((dim.cov+2)^2)*(COV.med+diag(0.1^10,(dim.cov+2))) #<-- dim.cov = the dimension of the covariates.
    a_prop <- c(rmnorm(1, a, COV.cov))
    
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
      tau_prop[i] <-rgamma(1, tau[i]*100, 100)
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
    
    Mc[MmisIndex] <- rnorm(length(MmisIndex), diag(cbind(1, Mp,Wp)[MmisIndex,]%*%t(update1$alpha[update0$Kc[MmisIndex],])), sqrt(1/update1$sigma[update0$Kc[MmisIndex]]))
    return(list(alpha=update1$alpha,sigma=update1$sigma, tau=tau, tau_acc=tau_acc, a=a,  Mc=Mc,Pc=update2$Pc, lambda=update2$lambda, Kc=update0$Kc))
}


#----- Load MCMC samples and Data
load("simulated_data.RData")
load(file=paste0("data_",process+1,".RData"))

##### data at time 1 : TRT=(0) #####

Mc <- as.numeric(M1)[Z1==0]
Mp <- as.numeric(Master$PM1)[Z1==0]
Wp <- cbind(Master$TTEMP3,  Master$PctUrban,  log(Master$Ter+1))[Z1==0,]
MmisIndex <- which(is.na(Mc))
Mc[MmisIndex] <- rnorm(length(MmisIndex),mean(Mc,na.rm=TRUE),sd(Mc,na.rm=TRUE))
K <- 20

Kp <- sample(1:K,length(Mc),replace=TRUE)

#### Priors ###############
dim.cov <- dim(Wp)[2]
a_prior <- list()
a_prior$P <- 1
a_prior$mean <- rep(0,dim.cov+2)
a_prior$sigma <- diag(1000, dim.cov+2)
a <- a_prior$mean
tau <- rep(0.5, dim.cov+2)

alpha <- matrix(coef(lm(Mc~Mp+Wp)), nrow=K,ncol=dim.cov+2, byrow=TRUE)

sigma <- rep(1,K)
lambda <- 2
Pp <- rep(1/K,K)
COV.med <- vcov(lm(Mc~Mp+Wp))/dim.cov

# Setting the number of MCMC interations
MCMC<-6000
COV.a <- rbind(COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med)/100

update3 <- list()
update3[[1]] <- list(Mc=Mc,
MmisIndex=MmisIndex,
Wp=Wp,
Mp=Mp,
K=K,
Kc=Kp,
a=a,
a_prior=a_prior,
tau=tau,
alpha=alpha,
sigma=sigma,
lambda=lambda,
Pc=Pp,
dim.cov=dim.cov,
COV.med=COV.med,
COV.a = COV.a)


a_storage <- NULL
alpha1_storage <- NULL
alpha2_storage <- NULL
alpha3_storage <- NULL
alpha4_storage <- NULL
alpha5_storage <- NULL
alpha6_storage <- NULL
alpha7_storage <- NULL
alpha8_storage <- NULL
alpha9_storage <- NULL
alpha10_storage <- NULL
alpha11_storage <- NULL
alpha12_storage <- NULL
alpha13_storage <- NULL
alpha14_storage <- NULL
alpha15_storage <- NULL
alpha16_storage <- NULL
alpha17_storage <- NULL
alpha18_storage <- NULL
alpha19_storage <- NULL
alpha20_storage <- NULL

for(c in 2:MCMC){
    a_storage <- rbind(a_storage, update3[[(c-1)]]$a)
    if(c > 500){
        COV.med <- cov(a_storage)
    }
    alpha1_storage <- rbind(alpha1_storage, update3[[(c-1)]]$alpha[1,])
    alpha2_storage <- rbind(alpha2_storage, update3[[(c-1)]]$alpha[2,])
    alpha3_storage <- rbind(alpha3_storage, update3[[(c-1)]]$alpha[3,])
    alpha4_storage <- rbind(alpha4_storage, update3[[(c-1)]]$alpha[4,])
    alpha5_storage <- rbind(alpha5_storage, update3[[(c-1)]]$alpha[5,])
    alpha6_storage <- rbind(alpha6_storage, update3[[(c-1)]]$alpha[6,])
    alpha7_storage <- rbind(alpha7_storage, update3[[(c-1)]]$alpha[7,])
    alpha8_storage <- rbind(alpha8_storage, update3[[(c-1)]]$alpha[8,])
    alpha9_storage <- rbind(alpha9_storage, update3[[(c-1)]]$alpha[9,])
    alpha10_storage <- rbind(alpha10_storage, update3[[(c-1)]]$alpha[10,])
    alpha11_storage <- rbind(alpha11_storage, update3[[(c-1)]]$alpha[11,])
    alpha12_storage <- rbind(alpha12_storage, update3[[(c-1)]]$alpha[12,])
    alpha13_storage <- rbind(alpha13_storage, update3[[(c-1)]]$alpha[13,])
    alpha14_storage <- rbind(alpha14_storage, update3[[(c-1)]]$alpha[14,])
    alpha15_storage <- rbind(alpha15_storage, update3[[(c-1)]]$alpha[15,])
    alpha16_storage <- rbind(alpha16_storage, update3[[(c-1)]]$alpha[16,])
    alpha17_storage <- rbind(alpha17_storage, update3[[(c-1)]]$alpha[17,])
    alpha18_storage <- rbind(alpha18_storage, update3[[(c-1)]]$alpha[18,])
    alpha19_storage <- rbind(alpha19_storage, update3[[(c-1)]]$alpha[19,])
    alpha20_storage <- rbind(alpha20_storage, update3[[(c-1)]]$alpha[20,])
    if(c > 500){
        for(i in 1:20){
            eval(parse(text=paste("COV.a[",(5*i-4),":",(5*i),",1:5] <- cov(alpha",i,"_storage)", sep="")))
        }
         }
    update3[[c]] <- Cond_Main(Mc=update3[[c-1]]$Mc,
    MmisIndex=MmisIndex,
    Mp=Mp,
    Wp=Wp,
    K=K,
    Kp=update3[[c-1]]$Kc,
    a=update3[[c-1]]$a,
    a_prior=a_prior,
    tau=update3[[c-1]]$tau,
    alpha=update3[[c-1]]$alpha,
    sigma=update3[[c-1]]$sigma,
    lambda=update3[[c-1]]$lambda,
    Pp=update3[[c-1]]$Pc,
    dim.cov=dim.cov,
    COV.med=COV.med,
    COV.a=COV.a)
}

med_time10 <- list()
for(ind in 1:2000){
    IND <- ind+4000
    med_time10[[ind]] <- list(a=update3[[IND]]$a,
    tau=update3[[IND]]$tau,
    alpha=update3[[IND]]$alpha,
    sigma=update3[[IND]]$sigma,
    Kc=update3[[IND]]$Kc)
}
rm(update3)


##### data at time 1 : TRT=(1) #####
Mc <- as.numeric(M1)[Z1==1]
Mp <- as.numeric(Master$PM1)[Z1==1]
Wp <- cbind(Master$TTEMP3,  Master$PctUrban,  log(Master$Ter+1))[Z1==1,]
MmisIndex <- which(is.na(Mc))
Mc[MmisIndex] <- rnorm(length(MmisIndex),mean(Mc,na.rm=TRUE),sd(Mc,na.rm=TRUE))
K <- 20

Kp <- sample(1:K,length(Mc),replace=TRUE)

#### Priors ###############
dim.cov <- dim(Wp)[2]
a_prior <- list()
a_prior$P <- 1
a_prior$mean <- rep(0,dim.cov+2)
a_prior$sigma <- diag(1000, dim.cov+2)
a <- a_prior$mean
tau <- rep(0.5, dim.cov+2)
alpha <- matrix(coef(lm(Mc~Mp+Wp)), nrow=K,ncol=dim.cov+2, byrow=TRUE)

sigma <- rep(1,K)
lambda <- 2
Pp <- rep(1/K,K)
COV.med <- vcov(lm(Mc~Mp+Wp))/dim.cov

# Setting the number of MCMC interations
MCMC<-6000

COV.a <- rbind(COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med)/100

update3 <- list()
update3[[1]] <- list(Mc=Mc,
                     MmisIndex=MmisIndex,
                     Wp=Wp,
                     Mp=Mp,
                     K=K,
                     Kc=Kp,
                     a=a,
                     a_prior=a_prior,
                     tau=tau,
                     alpha=alpha,
                     sigma=sigma,
                     lambda=lambda,
                     Pc=Pp,
                     dim.cov=dim.cov,
                     COV.med=COV.med,
                     COV.a = COV.a)

a_storage <- NULL
alpha1_storage <- NULL
alpha2_storage <- NULL
alpha3_storage <- NULL
alpha4_storage <- NULL
alpha5_storage <- NULL
alpha6_storage <- NULL
alpha7_storage <- NULL
alpha8_storage <- NULL
alpha9_storage <- NULL
alpha10_storage <- NULL
alpha11_storage <- NULL
alpha12_storage <- NULL
alpha13_storage <- NULL
alpha14_storage <- NULL
alpha15_storage <- NULL
alpha16_storage <- NULL
alpha17_storage <- NULL
alpha18_storage <- NULL
alpha19_storage <- NULL
alpha20_storage <- NULL

for(c in 2:MCMC){
  a_storage <- rbind(a_storage, update3[[(c-1)]]$a)
  if(c > 500){
    COV.med <- cov(a_storage)
  }
  alpha1_storage <- rbind(alpha1_storage, update3[[(c-1)]]$alpha[1,])
  alpha2_storage <- rbind(alpha2_storage, update3[[(c-1)]]$alpha[2,])
  alpha3_storage <- rbind(alpha3_storage, update3[[(c-1)]]$alpha[3,])
  alpha4_storage <- rbind(alpha4_storage, update3[[(c-1)]]$alpha[4,])
  alpha5_storage <- rbind(alpha5_storage, update3[[(c-1)]]$alpha[5,])
  alpha6_storage <- rbind(alpha6_storage, update3[[(c-1)]]$alpha[6,])
  alpha7_storage <- rbind(alpha7_storage, update3[[(c-1)]]$alpha[7,])
  alpha8_storage <- rbind(alpha8_storage, update3[[(c-1)]]$alpha[8,])
  alpha9_storage <- rbind(alpha9_storage, update3[[(c-1)]]$alpha[9,])
  alpha10_storage <- rbind(alpha10_storage, update3[[(c-1)]]$alpha[10,])
  alpha11_storage <- rbind(alpha11_storage, update3[[(c-1)]]$alpha[11,])
  alpha12_storage <- rbind(alpha12_storage, update3[[(c-1)]]$alpha[12,])
  alpha13_storage <- rbind(alpha13_storage, update3[[(c-1)]]$alpha[13,])
  alpha14_storage <- rbind(alpha14_storage, update3[[(c-1)]]$alpha[14,])
  alpha15_storage <- rbind(alpha15_storage, update3[[(c-1)]]$alpha[15,])
  alpha16_storage <- rbind(alpha16_storage, update3[[(c-1)]]$alpha[16,])
  alpha17_storage <- rbind(alpha17_storage, update3[[(c-1)]]$alpha[17,])
  alpha18_storage <- rbind(alpha18_storage, update3[[(c-1)]]$alpha[18,])
  alpha19_storage <- rbind(alpha19_storage, update3[[(c-1)]]$alpha[19,])
  alpha20_storage <- rbind(alpha20_storage, update3[[(c-1)]]$alpha[20,])
  if(c > 500){
    for(i in 1:20){
      eval(parse(text=paste("COV.a[",(5*i-4),":",(5*i),",1:5] <- cov(alpha",i,"_storage)", sep="")))
    }
  }

  update3[[c]] <- Cond_Main(Mc=update3[[c-1]]$Mc,
                            MmisIndex=MmisIndex,
                            Mp=Mp,
                            Wp=Wp,
                            K=K,
                            Kp=update3[[c-1]]$Kc,
                            a=update3[[c-1]]$a,
                            a_prior=a_prior,
                            tau=update3[[c-1]]$tau,
                            alpha=update3[[c-1]]$alpha,
                            sigma=update3[[c-1]]$sigma,
                            lambda=update3[[c-1]]$lambda,
                            Pp=update3[[c-1]]$Pc,
                            dim.cov=dim.cov,
                            COV.med=COV.med,
                            COV.a=COV.a)
}

med_time11 <- list()
for(ind in 1:2000){
  IND <- ind+4000
  med_time11[[ind]] <- list(a=update3[[IND]]$a,
                            tau=update3[[IND]]$tau,
                            alpha=update3[[IND]]$alpha,
                            sigma=update3[[IND]]$sigma,
                            Kc=update3[[IND]]$Kc)
}
rm(update3)

##### data at time 2 : TRT=(0,0) #####

Mc <- as.numeric(M2)[Z2==0]
Mp <- as.numeric(M1)[Z2==0]
Wp <- cbind(X2,  Master$PctUrban, log(Master$Ter+1))[Z2==0,]
MmisIndex <- which(is.na(Mc))
Mc[MmisIndex] <- rnorm(length(MmisIndex),mean(Mc,na.rm=TRUE),sd(Mc,na.rm=TRUE))

K <- 20
Kp <- sample(1:K,length(Mc),replace=TRUE)

#### Priors ###############
dim.cov <- dim(Wp)[2]
a_past.sample <- matrix(nrow=2000, ncol=5)
for(ind in 1:2000){
    IND <- ind+0
    a_past.sample[ind,]<-med_time10[[IND]]$a
}
a_new.sample <- matrix(nrow=2000, ncol=5)
for(ind in 1:2000){
    a_new.sample[ind,] <- a_past.sample[ind,]
}
a_prior$P <- 1
a_prior$mean <- apply(a_new.sample, 2, mean)
a_prior$sigma <- cov(a_new.sample)

a <- a_prior$mean
tau <- rep(0.5, dim.cov+2)
sigma <- rep(1,K)
alpha <- matrix(0, nrow=K,ncol=dim.cov+2, byrow=TRUE)

lambda <- 2
Pp <- rep(1/K,K)
COV.med <- vcov(lm(Mc~Mp+Wp))/(dim.cov)

# Setting the number of MCMC interations
MCMC<-6000

COV.a <- rbind(COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med)/100

update3 <- list()
update3[[1]] <- list(Mc=Mc,
MmisIndex=MmisIndex,
Mp=Mp,
Wp=Wp,
K=K,
Kc=Kp,
a=a,
a_prior=a_prior,
tau=tau,
alpha=alpha,
sigma=sigma,
lambda=lambda,
Pc=Pp,
dim.cov=dim.cov,
COV.med=COV.med,
COV.a=COV.a)


a_storage <- NULL
alpha1_storage <- NULL
alpha2_storage <- NULL
alpha3_storage <- NULL
alpha4_storage <- NULL
alpha5_storage <- NULL
alpha6_storage <- NULL
alpha7_storage <- NULL
alpha8_storage <- NULL
alpha9_storage <- NULL
alpha10_storage <- NULL
alpha11_storage <- NULL
alpha12_storage <- NULL
alpha13_storage <- NULL
alpha14_storage <- NULL
alpha15_storage <- NULL
alpha16_storage <- NULL
alpha17_storage <- NULL
alpha18_storage <- NULL
alpha19_storage <- NULL
alpha20_storage <- NULL

for(c in 2:MCMC){
    a_storage <- rbind(a_storage, update3[[(c-1)]]$a)
    if(c > 500){
        COV.med <- cov(a_storage)
    }
    alpha1_storage <- rbind(alpha1_storage, update3[[(c-1)]]$alpha[1,])
    alpha2_storage <- rbind(alpha2_storage, update3[[(c-1)]]$alpha[2,])
    alpha3_storage <- rbind(alpha3_storage, update3[[(c-1)]]$alpha[3,])
    alpha4_storage <- rbind(alpha4_storage, update3[[(c-1)]]$alpha[4,])
    alpha5_storage <- rbind(alpha5_storage, update3[[(c-1)]]$alpha[5,])
    alpha6_storage <- rbind(alpha6_storage, update3[[(c-1)]]$alpha[6,])
    alpha7_storage <- rbind(alpha7_storage, update3[[(c-1)]]$alpha[7,])
    alpha8_storage <- rbind(alpha8_storage, update3[[(c-1)]]$alpha[8,])
    alpha9_storage <- rbind(alpha9_storage, update3[[(c-1)]]$alpha[9,])
    alpha10_storage <- rbind(alpha10_storage, update3[[(c-1)]]$alpha[10,])
    alpha11_storage <- rbind(alpha11_storage, update3[[(c-1)]]$alpha[11,])
    alpha12_storage <- rbind(alpha12_storage, update3[[(c-1)]]$alpha[12,])
    alpha13_storage <- rbind(alpha13_storage, update3[[(c-1)]]$alpha[13,])
    alpha14_storage <- rbind(alpha14_storage, update3[[(c-1)]]$alpha[14,])
    alpha15_storage <- rbind(alpha15_storage, update3[[(c-1)]]$alpha[15,])
    alpha16_storage <- rbind(alpha16_storage, update3[[(c-1)]]$alpha[16,])
    alpha17_storage <- rbind(alpha17_storage, update3[[(c-1)]]$alpha[17,])
    alpha18_storage <- rbind(alpha18_storage, update3[[(c-1)]]$alpha[18,])
    alpha19_storage <- rbind(alpha19_storage, update3[[(c-1)]]$alpha[19,])
    alpha20_storage <- rbind(alpha20_storage, update3[[(c-1)]]$alpha[20,])
    if(c > 500){
        for(i in 1:20){
          eval(parse(text=paste("COV.a[",(5*i-4),":",(5*i),",1:5] <- cov(alpha",i,"_storage)", sep="")))
        }

    }
    
    update3[[c]] <- Cond_Main(Mc=update3[[c-1]]$Mc,
    MmisIndex=MmisIndex,
    Mp=Mp,
    Wp=Wp,
    K=K,
    Kp=update3[[c-1]]$Kc,
    a=update3[[c-1]]$a,
    a_prior=a_prior,
    tau=update3[[c-1]]$tau,
    alpha=update3[[c-1]]$alpha,
    sigma=update3[[c-1]]$sigma,
    lambda=update3[[c-1]]$lambda,
    Pp=update3[[c-1]]$Pc,
    dim.cov=dim.cov,
    COV.med=COV.med,
    COV.a=COV.a)
}

med_time20 <- list()
for(ind in 1:2000){
    IND <- ind+4000
    med_time20[[ind]] <- list(a=update3[[IND]]$a,
    tau=update3[[IND]]$tau,
    alpha=update3[[IND]]$alpha,
    sigma=update3[[IND]]$sigma,
    Kc=update3[[IND]]$Kc)
}
rm(update3)



##### data at time 2 : TRT=(0,1) #####
Mc <- as.numeric(M2)[Z2==1]
Mp <- as.numeric(M1)[Z2==1]
Wp <- cbind(X2,  Master$PctUrban, log(Master$Ter+1))[Z2==1,]
MmisIndex <- which(is.na(Mc))
Mc[MmisIndex] <- rnorm(length(MmisIndex),mean(Mc,na.rm=TRUE),sd(Mc,na.rm=TRUE))

K <- 20
Kp <- sample(1:K,length(Mc),replace=TRUE)

#### Priors ###############
dim.cov <- dim(Wp)[2]
a_past.sample <- matrix(nrow=2000, ncol=5)
for(ind in 1:2000){
  IND <- ind+0
  a_past.sample[ind,]<-med_time10[[IND]]$a
}
a_new.sample <- matrix(nrow=2000, ncol=5)
for(ind in 1:2000){
  a_new.sample[ind,] <- a_past.sample[ind,]
}
a_prior$P <- 1
a_prior$mean <- apply(a_new.sample, 2, mean)
a_prior$sigma <- cov(a_new.sample)

a <- a_prior$mean
tau <- rep(0.5, dim.cov+2)
sigma <- rep(1,K)
alpha <- matrix(0, nrow=K,ncol=dim.cov+2, byrow=TRUE)

lambda <- 2
Pp <- rep(1/K,K)
COV.med <- vcov(lm(Mc~Mp+Wp))/(dim.cov)

# Setting the number of MCMC interations
MCMC<-6000

COV.a <- rbind(COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med)/100

update3 <- list()
update3[[1]] <- list(Mc=Mc,
                     MmisIndex=MmisIndex,
                     Mp=Mp,
                     Wp=Wp,
                     K=K,
                     Kc=Kp,
                     a=a,
                     a_prior=a_prior,
                     tau=tau,
                     alpha=alpha,
                     sigma=sigma,
                     lambda=lambda,
                     Pc=Pp,
                     dim.cov=dim.cov,
                     COV.med=COV.med,
                     COV.a=COV.a)

a_storage <- NULL
alpha1_storage <- NULL
alpha2_storage <- NULL
alpha3_storage <- NULL
alpha4_storage <- NULL
alpha5_storage <- NULL
alpha6_storage <- NULL
alpha7_storage <- NULL
alpha8_storage <- NULL
alpha9_storage <- NULL
alpha10_storage <- NULL
alpha11_storage <- NULL
alpha12_storage <- NULL
alpha13_storage <- NULL
alpha14_storage <- NULL
alpha15_storage <- NULL
alpha16_storage <- NULL
alpha17_storage <- NULL
alpha18_storage <- NULL
alpha19_storage <- NULL
alpha20_storage <- NULL

for(c in 2:MCMC){
  a_storage <- rbind(a_storage, update3[[(c-1)]]$a)
  if(c > 500){
    COV.med <- cov(a_storage)
  }
  alpha1_storage <- rbind(alpha1_storage, update3[[(c-1)]]$alpha[1,])
  alpha2_storage <- rbind(alpha2_storage, update3[[(c-1)]]$alpha[2,])
  alpha3_storage <- rbind(alpha3_storage, update3[[(c-1)]]$alpha[3,])
  alpha4_storage <- rbind(alpha4_storage, update3[[(c-1)]]$alpha[4,])
  alpha5_storage <- rbind(alpha5_storage, update3[[(c-1)]]$alpha[5,])
  alpha6_storage <- rbind(alpha6_storage, update3[[(c-1)]]$alpha[6,])
  alpha7_storage <- rbind(alpha7_storage, update3[[(c-1)]]$alpha[7,])
  alpha8_storage <- rbind(alpha8_storage, update3[[(c-1)]]$alpha[8,])
  alpha9_storage <- rbind(alpha9_storage, update3[[(c-1)]]$alpha[9,])
  alpha10_storage <- rbind(alpha10_storage, update3[[(c-1)]]$alpha[10,])
  alpha11_storage <- rbind(alpha11_storage, update3[[(c-1)]]$alpha[11,])
  alpha12_storage <- rbind(alpha12_storage, update3[[(c-1)]]$alpha[12,])
  alpha13_storage <- rbind(alpha13_storage, update3[[(c-1)]]$alpha[13,])
  alpha14_storage <- rbind(alpha14_storage, update3[[(c-1)]]$alpha[14,])
  alpha15_storage <- rbind(alpha15_storage, update3[[(c-1)]]$alpha[15,])
  alpha16_storage <- rbind(alpha16_storage, update3[[(c-1)]]$alpha[16,])
  alpha17_storage <- rbind(alpha17_storage, update3[[(c-1)]]$alpha[17,])
  alpha18_storage <- rbind(alpha18_storage, update3[[(c-1)]]$alpha[18,])
  alpha19_storage <- rbind(alpha19_storage, update3[[(c-1)]]$alpha[19,])
  alpha20_storage <- rbind(alpha20_storage, update3[[(c-1)]]$alpha[20,])
  if(c > 500){
    for(i in 1:20){
      eval(parse(text=paste("COV.a[",(5*i-4),":",(5*i),",1:5] <- cov(alpha",i,"_storage)", sep="")))
    }
    
  }
  
  update3[[c]] <- Cond_Main(Mc=update3[[c-1]]$Mc,
                            MmisIndex=MmisIndex,
                            Mp=Mp,
                            Wp=Wp,
                            K=K,
                            Kp=update3[[c-1]]$Kc,
                            a=update3[[c-1]]$a,
                            a_prior=a_prior,
                            tau=update3[[c-1]]$tau,
                            alpha=update3[[c-1]]$alpha,
                            sigma=update3[[c-1]]$sigma,
                            lambda=update3[[c-1]]$lambda,
                            Pp=update3[[c-1]]$Pc,
                            dim.cov=dim.cov,
                            COV.med=COV.med,
                            COV.a=COV.a)
}

med_time21 <- list()
for(ind in 1:2000){
  IND <- ind+4000
  med_time21[[ind]] <- list(a=update3[[IND]]$a,
                            tau=update3[[IND]]$tau,
                            alpha=update3[[IND]]$alpha,
                            sigma=update3[[IND]]$sigma,
                            Kc=update3[[IND]]$Kc)
}
rm(update3)

##### data at time 3 : TRT=(0,0,0) #####

Mc <- as.numeric(M3)[Z3==0]
Mp <- as.numeric(M2)[Z3==0]
Wp <- cbind(X3,  Master$PctUrban,  log(Master$Ter+1))[Z3==0,]

MmisIndex <- which(is.na(Mc))
Mc[MmisIndex] <- rnorm(length(MmisIndex),mean(Mc,na.rm=TRUE),sd(Mc,na.rm=TRUE))

K <- 20

Kp <- sample(1:K,length(Mc),replace=TRUE)

#### Priors ###############
dim.cov <- dim(Wp)[2]
a_past.sample <- matrix(nrow=2000, ncol=5)
for(ind in 1:2000){
    IND <- ind+0
    a_past.sample[ind,]<-med_time20[[IND]]$a
}
a_new.sample <- matrix(nrow=2000, ncol=5)
for(ind in 1:2000){
    a_new.sample[ind,] <- a_past.sample[ind,]
}
a_prior$P <- 1
a_prior$mean <- apply(a_new.sample, 2, mean)
a_prior$sigma <- cov(a_new.sample)

a <- a_prior$mean
tau <- rep(0.5, dim.cov+2)
alpha <- matrix(0, nrow=K,ncol=dim.cov+2, byrow=TRUE)
sigma <-rep(1,K)
lambda <- 2
Pp <- rep(1/K,K)
COV.med <- vcov(lm(Mc~Mp+Wp))/(dim.cov)

# Setting the number of MCMC interations
MCMC<-6000

COV.a <- rbind(COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med)/100

update3 <- list()
update3[[1]] <- list(Mc=Mc,
MmisIndex=MmisIndex,
Mp=Mp,
Wp=Wp,
K=K,
Kc=Kp,
a=a,
a_prior=a_prior,
tau=tau,
alpha=alpha,
sigma=sigma,
lambda=lambda,
Pc=Pp,
dim.cov=dim.cov,
COV.med=COV.med,
COV.a=COV.a)

a_storage <- NULL
alpha1_storage <- NULL
alpha2_storage <- NULL
alpha3_storage <- NULL
alpha4_storage <- NULL
alpha5_storage <- NULL
alpha6_storage <- NULL
alpha7_storage <- NULL
alpha8_storage <- NULL
alpha9_storage <- NULL
alpha10_storage <- NULL
alpha11_storage <- NULL
alpha12_storage <- NULL
alpha13_storage <- NULL
alpha14_storage <- NULL
alpha15_storage <- NULL
alpha16_storage <- NULL
alpha17_storage <- NULL
alpha18_storage <- NULL
alpha19_storage <- NULL
alpha20_storage <- NULL

for(c in 2:MCMC){
    a_storage <- rbind(a_storage, update3[[(c-1)]]$a)
    if(c > 500){
        COV.med <- cov(a_storage)
    }
    alpha1_storage <- rbind(alpha1_storage, update3[[(c-1)]]$alpha[1,])
    alpha2_storage <- rbind(alpha2_storage, update3[[(c-1)]]$alpha[2,])
    alpha3_storage <- rbind(alpha3_storage, update3[[(c-1)]]$alpha[3,])
    alpha4_storage <- rbind(alpha4_storage, update3[[(c-1)]]$alpha[4,])
    alpha5_storage <- rbind(alpha5_storage, update3[[(c-1)]]$alpha[5,])
    alpha6_storage <- rbind(alpha6_storage, update3[[(c-1)]]$alpha[6,])
    alpha7_storage <- rbind(alpha7_storage, update3[[(c-1)]]$alpha[7,])
    alpha8_storage <- rbind(alpha8_storage, update3[[(c-1)]]$alpha[8,])
    alpha9_storage <- rbind(alpha9_storage, update3[[(c-1)]]$alpha[9,])
    alpha10_storage <- rbind(alpha10_storage, update3[[(c-1)]]$alpha[10,])
    alpha11_storage <- rbind(alpha11_storage, update3[[(c-1)]]$alpha[11,])
    alpha12_storage <- rbind(alpha12_storage, update3[[(c-1)]]$alpha[12,])
    alpha13_storage <- rbind(alpha13_storage, update3[[(c-1)]]$alpha[13,])
    alpha14_storage <- rbind(alpha14_storage, update3[[(c-1)]]$alpha[14,])
    alpha15_storage <- rbind(alpha15_storage, update3[[(c-1)]]$alpha[15,])
    alpha16_storage <- rbind(alpha16_storage, update3[[(c-1)]]$alpha[16,])
    alpha17_storage <- rbind(alpha17_storage, update3[[(c-1)]]$alpha[17,])
    alpha18_storage <- rbind(alpha18_storage, update3[[(c-1)]]$alpha[18,])
    alpha19_storage <- rbind(alpha19_storage, update3[[(c-1)]]$alpha[19,])
    alpha20_storage <- rbind(alpha20_storage, update3[[(c-1)]]$alpha[20,])
    if(c > 500){
        for(i in 1:20){
          eval(parse(text=paste("COV.a[",(5*i-4),":",(5*i),",1:5] <- cov(alpha",i,"_storage)", sep="")))
        }

    }
    
    update3[[c]] <- Cond_Main(Mc=update3[[c-1]]$Mc,
    MmisIndex=MmisIndex,
    Mp=Mp,
    Wp=Wp,
    K=K,
    Kp=update3[[c-1]]$Kc,
    a=update3[[c-1]]$a,
    a_prior=a_prior,
    tau=update3[[c-1]]$tau,
    alpha=update3[[c-1]]$alpha,
    sigma=update3[[c-1]]$sigma,
    lambda=update3[[c-1]]$lambda,
    Pp=update3[[c-1]]$Pc,
    dim.cov=dim.cov,
    COV.med=COV.med,
    COV.a=COV.a)
    print(c);print(update3[[c-1]]$sigma)
}

med_time30 <- list()
for(ind in 1:2000){
    IND <- ind+4000
    med_time30[[ind]] <- list(a=update3[[IND]]$a,
    tau=update3[[IND]]$tau,
    alpha=update3[[IND]]$alpha,
    sigma=update3[[IND]]$sigma,
    Kc=update3[[IND]]$Kc)
}
rm(update3)

Mc <- as.numeric(M3)[Z3==1]
Mp <- as.numeric(M2)[Z3==1]
Wp <- cbind(X3,  Master$PctUrban,  log(Master$Ter+1))[Z3==1,]

MmisIndex <- which(is.na(Mc))
Mc[MmisIndex] <- rnorm(length(MmisIndex),mean(Mc,na.rm=TRUE),sd(Mc,na.rm=TRUE))

K <- 20

Kp <- sample(1:K,length(Mc),replace=TRUE)

#### Priors ###############
dim.cov <- dim(Wp)[2]

a_past.sample <- matrix(nrow=2000, ncol=5)
for(ind in 1:2000){
  IND <- ind+0
  a_past.sample[ind,]<-med_time20[[IND]]$a
}
a_new.sample <- matrix(nrow=2000, ncol=5)
for(ind in 1:2000){
  a_new.sample[ind,] <- a_past.sample[ind,]
}
a_prior$P <- 1
a_prior$mean <- apply(a_new.sample, 2, mean)
a_prior$sigma <- cov(a_new.sample)
a <- a_prior$mean

tau <- rep(0.5, dim.cov+2)
alpha <- matrix(0, nrow=K,ncol=dim.cov+2, byrow=TRUE)


#### Starting PMs ############
sigma <-rep(1,K)
lambda <- 2
Pp <- rep(1/K,K)
COV.med <- vcov(lm(Mc~Mp+Wp))/(dim.cov)

# Setting the number of MCMC interations
MCMC<-6000

COV.a <- rbind(COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med)/100

update3 <- list()
update3[[1]] <- list(Mc=Mc,
                     MmisIndex=MmisIndex,
                     Mp=Mp,
                     Wp=Wp,
                     K=K,
                     Kc=Kp,
                     a=a,
                     a_prior=a_prior,
                     tau=tau,
                     alpha=alpha,
                     sigma=sigma,
                     lambda=lambda,
                     Pc=Pp,
                     dim.cov=dim.cov,
                     COV.med=COV.med,
                     COV.a=COV.a)


a_storage <- NULL
alpha1_storage <- NULL
alpha2_storage <- NULL
alpha3_storage <- NULL
alpha4_storage <- NULL
alpha5_storage <- NULL
alpha6_storage <- NULL
alpha7_storage <- NULL
alpha8_storage <- NULL
alpha9_storage <- NULL
alpha10_storage <- NULL
alpha11_storage <- NULL
alpha12_storage <- NULL
alpha13_storage <- NULL
alpha14_storage <- NULL
alpha15_storage <- NULL
alpha16_storage <- NULL
alpha17_storage <- NULL
alpha18_storage <- NULL
alpha19_storage <- NULL
alpha20_storage <- NULL

for(c in 2:MCMC){
  a_storage <- rbind(a_storage, update3[[(c-1)]]$a)
  if(c > 500){
    COV.med <- cov(a_storage)
  }
  alpha1_storage <- rbind(alpha1_storage, update3[[(c-1)]]$alpha[1,])
  alpha2_storage <- rbind(alpha2_storage, update3[[(c-1)]]$alpha[2,])
  alpha3_storage <- rbind(alpha3_storage, update3[[(c-1)]]$alpha[3,])
  alpha4_storage <- rbind(alpha4_storage, update3[[(c-1)]]$alpha[4,])
  alpha5_storage <- rbind(alpha5_storage, update3[[(c-1)]]$alpha[5,])
  alpha6_storage <- rbind(alpha6_storage, update3[[(c-1)]]$alpha[6,])
  alpha7_storage <- rbind(alpha7_storage, update3[[(c-1)]]$alpha[7,])
  alpha8_storage <- rbind(alpha8_storage, update3[[(c-1)]]$alpha[8,])
  alpha9_storage <- rbind(alpha9_storage, update3[[(c-1)]]$alpha[9,])
  alpha10_storage <- rbind(alpha10_storage, update3[[(c-1)]]$alpha[10,])
  alpha11_storage <- rbind(alpha11_storage, update3[[(c-1)]]$alpha[11,])
  alpha12_storage <- rbind(alpha12_storage, update3[[(c-1)]]$alpha[12,])
  alpha13_storage <- rbind(alpha13_storage, update3[[(c-1)]]$alpha[13,])
  alpha14_storage <- rbind(alpha14_storage, update3[[(c-1)]]$alpha[14,])
  alpha15_storage <- rbind(alpha15_storage, update3[[(c-1)]]$alpha[15,])
  alpha16_storage <- rbind(alpha16_storage, update3[[(c-1)]]$alpha[16,])
  alpha17_storage <- rbind(alpha17_storage, update3[[(c-1)]]$alpha[17,])
  alpha18_storage <- rbind(alpha18_storage, update3[[(c-1)]]$alpha[18,])
  alpha19_storage <- rbind(alpha19_storage, update3[[(c-1)]]$alpha[19,])
  alpha20_storage <- rbind(alpha20_storage, update3[[(c-1)]]$alpha[20,])
  if(c > 500){
    for(i in 1:20){
      eval(parse(text=paste("COV.a[",(5*i-4),":",(5*i),",1:5] <- cov(alpha",i,"_storage)", sep="")))
    }
    
  }
  
  update3[[c]] <- Cond_Main(Mc=update3[[c-1]]$Mc,
                            MmisIndex=MmisIndex,
                            Mp=Mp,
                            Wp=Wp,
                            K=K,
                            Kp=update3[[c-1]]$Kc,
                            a=update3[[c-1]]$a,
                            a_prior=a_prior,
                            tau=update3[[c-1]]$tau,
                            alpha=update3[[c-1]]$alpha,
                            sigma=update3[[c-1]]$sigma,
                            lambda=update3[[c-1]]$lambda,
                            Pp=update3[[c-1]]$Pc,
                            dim.cov=dim.cov,
                            COV.med=COV.med,
                            COV.a=COV.a)
}

med_time31 <- list()
for(ind in 1:2000){
  IND <- ind+4000
  med_time31[[ind]] <- list(a=update3[[IND]]$a,
                            tau=update3[[IND]]$tau,
                            alpha=update3[[IND]]$alpha,
                            sigma=update3[[IND]]$sigma,
                            Kc=update3[[IND]]$Kc)
}
rm(update3)


##### data at time 4 : TRT=(0,0,0,0) #####

Mc <- as.numeric(M4)[Z4==0]
Mp <- as.numeric(M3)[Z4==0]
Wp <- cbind(X4,  Master$PctUrban,  log(Master$Ter+1))[Z4==0,]
MmisIndex <- which(is.na(Mc))
Mc[MmisIndex] <- rnorm(length(MmisIndex),mean(Mc,na.rm=TRUE),sd(Mc,na.rm=TRUE))
K <- 20
Kp <- sample(1:K,length(Mc),replace=TRUE)

#### Priors ###############
dim.cov <- dim(Wp)[2]
a_past.sample <- matrix(nrow=2000, ncol=5)
for(ind in 1:2000){
    IND <- ind+0
    a_past.sample[ind,]<-med_time30[[IND]]$a
}
a_new.sample <- matrix(nrow=2000, ncol=5)
for(ind in 1:2000){
    a_new.sample[ind,] <- a_past.sample[ind,]
}
a_prior$P <- 1
a_prior$mean <- apply(a_new.sample, 2, mean)
a_prior$sigma <- cov(a_new.sample)

a <- a_prior$mean
tau <- rep(0.5, dim.cov+2)
alpha <- matrix(0, nrow=K,ncol=dim.cov+2, byrow=TRUE)


#### Starting PMs ############
sigma <- rep(1,K)
lambda <- 2
Pp <- rep(1/K,K)
COV.med <- vcov(lm(Mc~Mp+Wp))/(dim.cov)

# Setting the number of MCMC interations
MCMC<-6000

COV.a <- rbind(COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med)/100

update3 <- list()
update3[[1]] <- list(Mc=Mc,
MmisIndex=MmisIndex,
Mp=Mp,
Wp=Wp,
K=K,
Kc=Kp,
a=a,
a_prior=a_prior,
tau=tau,
alpha=alpha,
sigma=sigma,
lambda=lambda,
Pc=Pp,
dim.cov=dim.cov,
COV.med=COV.med,
COV.a=COV.a)

a_storage <- NULL
alpha1_storage <- NULL
alpha2_storage <- NULL
alpha3_storage <- NULL
alpha4_storage <- NULL
alpha5_storage <- NULL
alpha6_storage <- NULL
alpha7_storage <- NULL
alpha8_storage <- NULL
alpha9_storage <- NULL
alpha10_storage <- NULL
alpha11_storage <- NULL
alpha12_storage <- NULL
alpha13_storage <- NULL
alpha14_storage <- NULL
alpha15_storage <- NULL
alpha16_storage <- NULL
alpha17_storage <- NULL
alpha18_storage <- NULL
alpha19_storage <- NULL
alpha20_storage <- NULL

for(c in 2:MCMC){
    a_storage <- rbind(a_storage, update3[[(c-1)]]$a)
    if(c > 500){
        COV.med <- cov(a_storage)
    }
    alpha1_storage <- rbind(alpha1_storage, update3[[(c-1)]]$alpha[1,])
    alpha2_storage <- rbind(alpha2_storage, update3[[(c-1)]]$alpha[2,])
    alpha3_storage <- rbind(alpha3_storage, update3[[(c-1)]]$alpha[3,])
    alpha4_storage <- rbind(alpha4_storage, update3[[(c-1)]]$alpha[4,])
    alpha5_storage <- rbind(alpha5_storage, update3[[(c-1)]]$alpha[5,])
    alpha6_storage <- rbind(alpha6_storage, update3[[(c-1)]]$alpha[6,])
    alpha7_storage <- rbind(alpha7_storage, update3[[(c-1)]]$alpha[7,])
    alpha8_storage <- rbind(alpha8_storage, update3[[(c-1)]]$alpha[8,])
    alpha9_storage <- rbind(alpha9_storage, update3[[(c-1)]]$alpha[9,])
    alpha10_storage <- rbind(alpha10_storage, update3[[(c-1)]]$alpha[10,])
    alpha11_storage <- rbind(alpha11_storage, update3[[(c-1)]]$alpha[11,])
    alpha12_storage <- rbind(alpha12_storage, update3[[(c-1)]]$alpha[12,])
    alpha13_storage <- rbind(alpha13_storage, update3[[(c-1)]]$alpha[13,])
    alpha14_storage <- rbind(alpha14_storage, update3[[(c-1)]]$alpha[14,])
    alpha15_storage <- rbind(alpha15_storage, update3[[(c-1)]]$alpha[15,])
    alpha16_storage <- rbind(alpha16_storage, update3[[(c-1)]]$alpha[16,])
    alpha17_storage <- rbind(alpha17_storage, update3[[(c-1)]]$alpha[17,])
    alpha18_storage <- rbind(alpha18_storage, update3[[(c-1)]]$alpha[18,])
    alpha19_storage <- rbind(alpha19_storage, update3[[(c-1)]]$alpha[19,])
    alpha20_storage <- rbind(alpha20_storage, update3[[(c-1)]]$alpha[20,])
    if(c > 500){
        for(i in 1:20){
          eval(parse(text=paste("COV.a[",(5*i-4),":",(5*i),",1:5] <- cov(alpha",i,"_storage)", sep="")))
        }

    }
    
    update3[[c]] <- Cond_Main(Mc=update3[[c-1]]$Mc,
    MmisIndex=MmisIndex,
    Mp=Mp,
    Wp=Wp,
    K=K,
    Kp=update3[[c-1]]$Kc,
    a=update3[[c-1]]$a,
    a_prior=a_prior,
    tau=update3[[c-1]]$tau,
    alpha=update3[[c-1]]$alpha,
    sigma=update3[[c-1]]$sigma,
    lambda=update3[[c-1]]$lambda,
    Pp=update3[[c-1]]$Pc,
    dim.cov=dim.cov,
    COV.med=COV.med,
    COV.a=COV.a)
}

med_time40 <- list()
for(ind in 1:2000){
    IND <- ind+4000
    med_time40[[ind]] <- list(a=update3[[IND]]$a,
    tau=update3[[IND]]$tau,
    alpha=update3[[IND]]$alpha,
    sigma=update3[[IND]]$sigma,
    Kc=update3[[IND]]$Kc)
}
rm(update3)

##### data at time 4 : TRT=(0,0,0,1) #####
Mc <- as.numeric(M4)[Z4==1]
Mp <- as.numeric(M3)[Z4==1]
Wp <- cbind(X4,  Master$PctUrban,  log(Master$Ter+1))[Z4==1,]
MmisIndex <- which(is.na(Mc))
Mc[MmisIndex] <- rnorm(length(MmisIndex),mean(Mc,na.rm=TRUE),sd(Mc,na.rm=TRUE))
K <- 20
Kp <- sample(1:K,length(Mc),replace=TRUE)

#### Priors ###############
dim.cov <- dim(Wp)[2]
a_past.sample <- matrix(nrow=2000, ncol=5)
for(ind in 1:2000){
  IND <- ind+0
  a_past.sample[ind,]<-med_time30[[IND]]$a
}
a_new.sample <- matrix(nrow=2000, ncol=5)
for(ind in 1:2000){
  a_new.sample[ind,] <- a_past.sample[ind,]
}
a_prior$P <- 1
a_prior$mean <- apply(a_new.sample, 2, mean)
a_prior$sigma <- cov(a_new.sample)

a <- a_prior$mean
tau <- rep(0.5, dim.cov+2)
alpha <- matrix(0, nrow=K,ncol=dim.cov+2, byrow=TRUE)


#### Starting PMs ############
sigma <- rep(1,K)
lambda <- 2
Pp <- rep(1/K,K)
COV.med <- vcov(lm(Mc~Mp+Wp))/(dim.cov)

# Setting the number of MCMC interations
MCMC<-6000

COV.a <- rbind(COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med)/100

update3 <- list()
update3[[1]] <- list(Mc=Mc,
                     MmisIndex=MmisIndex,
                     Mp=Mp,
                     Wp=Wp,
                     K=K,
                     Kc=Kp,
                     a=a,
                     a_prior=a_prior,
                     tau=tau,
                     alpha=alpha,
                     sigma=sigma,
                     lambda=lambda,
                     Pc=Pp,
                     dim.cov=dim.cov,
                     COV.med=COV.med,
                     COV.a=COV.a)
a_storage <- NULL
alpha1_storage <- NULL
alpha2_storage <- NULL
alpha3_storage <- NULL
alpha4_storage <- NULL
alpha5_storage <- NULL
alpha6_storage <- NULL
alpha7_storage <- NULL
alpha8_storage <- NULL
alpha9_storage <- NULL
alpha10_storage <- NULL
alpha11_storage <- NULL
alpha12_storage <- NULL
alpha13_storage <- NULL
alpha14_storage <- NULL
alpha15_storage <- NULL
alpha16_storage <- NULL
alpha17_storage <- NULL
alpha18_storage <- NULL
alpha19_storage <- NULL
alpha20_storage <- NULL

for(c in 2:MCMC){
  a_storage <- rbind(a_storage, update3[[(c-1)]]$a)
  if(c > 500){
    COV.med <- cov(a_storage)
  }
  alpha1_storage <- rbind(alpha1_storage, update3[[(c-1)]]$alpha[1,])
  alpha2_storage <- rbind(alpha2_storage, update3[[(c-1)]]$alpha[2,])
  alpha3_storage <- rbind(alpha3_storage, update3[[(c-1)]]$alpha[3,])
  alpha4_storage <- rbind(alpha4_storage, update3[[(c-1)]]$alpha[4,])
  alpha5_storage <- rbind(alpha5_storage, update3[[(c-1)]]$alpha[5,])
  alpha6_storage <- rbind(alpha6_storage, update3[[(c-1)]]$alpha[6,])
  alpha7_storage <- rbind(alpha7_storage, update3[[(c-1)]]$alpha[7,])
  alpha8_storage <- rbind(alpha8_storage, update3[[(c-1)]]$alpha[8,])
  alpha9_storage <- rbind(alpha9_storage, update3[[(c-1)]]$alpha[9,])
  alpha10_storage <- rbind(alpha10_storage, update3[[(c-1)]]$alpha[10,])
  alpha11_storage <- rbind(alpha11_storage, update3[[(c-1)]]$alpha[11,])
  alpha12_storage <- rbind(alpha12_storage, update3[[(c-1)]]$alpha[12,])
  alpha13_storage <- rbind(alpha13_storage, update3[[(c-1)]]$alpha[13,])
  alpha14_storage <- rbind(alpha14_storage, update3[[(c-1)]]$alpha[14,])
  alpha15_storage <- rbind(alpha15_storage, update3[[(c-1)]]$alpha[15,])
  alpha16_storage <- rbind(alpha16_storage, update3[[(c-1)]]$alpha[16,])
  alpha17_storage <- rbind(alpha17_storage, update3[[(c-1)]]$alpha[17,])
  alpha18_storage <- rbind(alpha18_storage, update3[[(c-1)]]$alpha[18,])
  alpha19_storage <- rbind(alpha19_storage, update3[[(c-1)]]$alpha[19,])
  alpha20_storage <- rbind(alpha20_storage, update3[[(c-1)]]$alpha[20,])
  if(c > 500){
    for(i in 1:20){
      eval(parse(text=paste("COV.a[",(5*i-4),":",(5*i),",1:5] <- cov(alpha",i,"_storage)", sep="")))
    }
    
  }
  
  update3[[c]] <- Cond_Main(Mc=update3[[c-1]]$Mc,
                            MmisIndex=MmisIndex,
                            Mp=Mp,
                            Wp=Wp,
                            K=K,
                            Kp=update3[[c-1]]$Kc,
                            a=update3[[c-1]]$a,
                            a_prior=a_prior,
                            tau=update3[[c-1]]$tau,
                            alpha=update3[[c-1]]$alpha,
                            sigma=update3[[c-1]]$sigma,
                            lambda=update3[[c-1]]$lambda,
                            Pp=update3[[c-1]]$Pc,
                            dim.cov=dim.cov,
                            COV.med=COV.med,
                            COV.a=COV.a)
}

med_time41 <- list()
for(ind in 1:2000){
  IND <- ind+4000
  med_time41[[ind]] <- list(a=update3[[IND]]$a,
                            tau=update3[[IND]]$tau,
                            alpha=update3[[IND]]$alpha,
                            sigma=update3[[IND]]$sigma,
                            Kc=update3[[IND]]$Kc)
}
rm(update3)

save(med_time10,med_time20,med_time30,med_time40,med_time11,med_time21,med_time31,med_time41, file=paste0("Sim_Med_",process,".RData"))






