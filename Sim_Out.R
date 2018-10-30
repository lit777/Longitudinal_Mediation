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
library(sn)

source("Cond_P.R")

Cond_Alpha_O <- function(Yc=NULL, Mp=NULL,  Wp=NULL,Yp=NULL, Kp=NULL,K=NULL, a=NULL,tau=NULL, alpha=NULL,OFF=NULL, COV.a=NULL){
    
    K_other <- setdiff(1:K, unique(Kp))
    
    alpha_prop <- alpha
    alpha_acc <- NULL
    
    for(ind in K_other){
        alpha_prop[ind,] <- rmnorm(1, a, diag(1/tau,dim.cov+3)/10)
        alpha[ind,] <- alpha_prop[ind,]
    }
    
    # Check SD of proposal distribution (!Adaptive sampler!)
    for(ind in unique(Kp)){
        
        COV.alpha <- 2.38^2/((dim.cov+3))*(COV.a[(6*ind-5):(6*ind),1:6]+diag(0.1^10,(dim.cov+3)))
        alpha_prop[ind,] <- rmnorm(1, alpha[ind,], COV.alpha)
        
        rat1 <- dmnorm(alpha_prop[ind,], a, diag(1/tau,dim.cov+3),  log=TRUE)+dmnorm(alpha[ind,], alpha_prop[ind,], COV.alpha, log=TRUE)+sum(dpois(Yc[which(Kp==ind)], exp(log(OFF[which(Kp==ind)])+cbind(1,Mp,Yp,Wp)[which(Kp==ind),]%*%(alpha_prop[ind,])),log=TRUE))
        rat2 <- dmnorm(alpha[ind,], a, diag(1/tau,dim.cov+3),  log=TRUE)+dmnorm(alpha_prop[ind,], alpha[ind,], COV.alpha, log=TRUE)+sum(dpois(Yc[which(Kp==ind)], exp(log(OFF[which(Kp==ind)])+cbind(1,Mp,Yp,Wp)[which(Kp==ind),]%*%(alpha[ind,])),log=TRUE))
        
        rat <- rat1 - rat2
        if(is.na(rat) | log(runif(1))>rat){
            alpha_prop[ind,] <- alpha[ind,]
        }else{
            alpha[ind,] <- alpha_prop[ind,]
        }
    }
    
    return(list(alpha=alpha_prop, acc=alpha_acc))
}


Cond_K_O <- function(Yc=NULL, Mp=NULL, Wp=NULL,Yp=NULL,K=NULL, Pp=NULL, alpha=NULL, sigma=NULL, OFF=NULL){
    
    Pc.pre <- sapply(1:K, function(k) Pp[k]*dpois(Yc, exp(log(OFF)+cbind(1,Mp,Yp,Wp)%*%(alpha[k,]))))
    
    zero_index <- which(rowSums(Pc.pre)==0 | is.na(rowSums(Pc.pre)))
    Pc <- Pc.pre
    Pc[zero_index,] <- 1
    
    Kc <- apply(Pc,1, function(x) sample(1:20,1, prob=x))
    
    return(list(Kc=Kc))
}

Cond_Main_O <- function(Yc=NULL, Mp=NULL, Yp=NULL, Wp=NULL, K=NULL,Kp=Kp, a=NULL, a_prior=NULL,tau=NULL,tau_prior=NULL, alpha=NULL,   lambda=NULL,Pp=NULL,dim.cov=NULL,COV.med=NULL, OFF=NULL, COV.a=NULL, COV.tau=NULL){
    
    #---- updating Alpha
    update0 <- Cond_Alpha_O(Yc,Mp,Wp,Yp,Kp,K,a,tau,alpha,OFF,COV.a) #<-- alpha
    
    #---- updating assignment vectors
    update1 <- Cond_K_O(Yc,Mp,Wp,Yp,K,Pp,update0$alpha,sigma,OFF) #<-- Kc
    
    #---- updating weights and lambda
    update2 <- Cond_P(update1$Kc,K,lambda) #<-- Pc and Lambda
    
    #---- Sampling Base Measure
    # Check SD of the propsal distribution (!Adaptive sampler!)
    a_prop <- NULL
    a_acc <- NULL
    COV.cov <- 2.38^2/((dim.cov+3))*(COV.med+diag(0.1^10,(dim.cov+3))) #<-- dim.cov = the dimension of the covariates.
    a_prop <- c(rmnorm(1, a, COV.cov))
    
    rat1 <- log(sum(a_prior$P*dmnorm(a_prop, a_prior$mean, a_prior$sigma))) + sum(sapply(1:K, function(x) sum(dnorm(update0$alpha[x,], a_prop, sqrt(1/tau),  log=TRUE)))) + dmnorm(a, a_prop, COV.cov, log=TRUE)
    rat2 <- log(sum(a_prior$P*dmnorm(a, a_prior$mean, a_prior$sigma))) + sum(sapply(1:K, function(x) sum(dnorm(update0$alpha[x,], a, sqrt(1/tau),  log=TRUE)))) + dmnorm(a_prop, a, COV.cov, log=TRUE)
    
    rat <- rat1 - rat2
    if(is.na(rat) | log(runif(1))>rat){
        a_prop <- a
    }else{
        a <- a_prop
    }
    
    tau_prop <- tau
    tau_acc <- NULL
    
    for(i in 1:(dim.cov+3)){
    tau_prop[i] <- rgamma(1, tau[i]*100, 100)

    rat1 <-  sum(dgamma(tau_prop[i], 2, 1, log=TRUE))+ sum(sapply(1:K, function(x) sum(dnorm(update0$alpha[x,], a, sqrt(1/tau_prop),  log=TRUE)))) + dgamma(tau[i], tau_prop[i]*100, 100, log=TRUE)
    rat2 <- sum(dgamma(tau[i], 2, 1, log=TRUE)) + sum(sapply(1:K, function(x) sum(dnorm(update0$alpha[x,], a, sqrt(1/tau),  log=TRUE))))  + dgamma(tau_prop[i], tau[i]*100, 100, log=TRUE)
    rat <- rat1 - rat2
    if(is.na(rat) | log(runif(1))>rat){
        tau_prop[i] <- tau[i]
        tau_acc <- 0
    }else{
        tau[i] <- tau_prop[i]
        tau_acc <- 1
    }
}
    return(list(alpha=update0$alpha,  tau=tau, tau_acc=tau_acc, a=a, a_acc=a_acc, Pc=update2$Pc, lambda=update2$lambda, Kc=update1$Kc))
}

#----- Load Data
load("simulated_data.RData")
load(file=paste0("data_",process+1,".RData"))

##### data at time 1 : TRT=(0) #####

Yc <- Y1[Z1==0]
Mp <- as.numeric(M1[Z1==0])
Yp <- log(Master$TD1[Z1==0]+1.5)
Wp <- cbind(X1,  Master$PctUrban,  log(Master$Ter+1))[Z1==0,]
OFF <- Master$TB3[Z1==0]*3

K <- 20
Kp <- sample(1:K,length(Yc),replace=TRUE)

#### Priors ###############
dim.cov <- dim(Wp)[2]
a_prior <- list()
a_prior$P <- 1
a_prior$mean <- rep(0,6)
a_prior$sigma <- diag(1000, dim.cov+3)
a <- a_prior$mean

alpha <- matrix(coef(glm(Yc~offset(log(OFF))+Mp+Yp+Wp, family=poisson)), nrow=K,ncol=dim.cov+3, byrow=TRUE)

#### Starting PMs ############
lambda <- 2
Pp <- rep(1/K,K)
COV.med <- vcov(glm(Yc~offset(log(OFF))+Mp+Yp+Wp, family=poisson))/dim.cov

# Setting the number of MCMC interations
MCMC<-6000
tau <- rep(0.5, dim.cov+3)
COV.a <- rbind(COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med)/100

update3 <- list()
update3[[1]] <- list(Yc=Yc,
                     Mp=Mp,
                     Yp=Yp,
                     Wp=Wp,
                     K=K,
                     Kc=Kp,
                     a=a,
                     a_prior=a_prior,
                     tau=tau,
                     alpha=alpha,
                     lambda=lambda,
                     Pc=Pp,
                     dim.cov=dim.cov,
                     COV.med=COV.med,
                     COV.a=COV.a,
                     OFF=OFF)

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
          eval(parse(text=paste("COV.a[",(6*i-5),":",(6*i),",1:6] <- cov(alpha",i,"_storage)", sep="")))
      }
  }
  update3[[c]] <- Cond_Main_O(Yc=Yc,
                                   Mp=Mp,
                                   Yp=Yp,
                                   Wp=Wp,
                                   K=K,
                                   Kp=update3[[c-1]]$Kc,
                                   a=update3[[c-1]]$a,
                                   a_prior=a_prior,
                                   tau=update3[[c-1]]$tau,
                                   alpha=update3[[c-1]]$alpha,
                                   lambda=update3[[c-1]]$lambda,
                                   Pp=update3[[c-1]]$Pc,
                                   dim.cov=dim.cov,
                                   COV.med=COV.med,
                                   COV.a=COV.a,
                                   OFF=OFF)
}

out_time10 <- list()
for(ind in 1:2000){
    IND <- ind+4000
    out_time10[[ind]] <- list(a=update3[[IND]]$a,
    tau=update3[[IND]]$tau,
    alpha=update3[[IND]]$alpha,
    Kc=update3[[IND]]$Kc)
}
rm(update3)


##### data at time 1 : TRT=(1) #####

Yc <- Y1[Z1==1]
Mp <- as.numeric(M1[Z1==1])
Yp <- log(Master$TD1[Z1==1]+1.5)
Wp <- cbind(X1,  Master$PctUrban,  log(Master$Ter+1))[Z1==1,]
OFF <- Master$TB3[Z1==1]*3

K <- 20
Kp <- sample(1:K,length(Yc),replace=TRUE)

#### Priors ###############
dim.cov <- dim(Wp)[2]
a_prior <- list()
a_prior$P <- 1
a_prior$mean <- rep(0,6)
a_prior$sigma <- diag(1000, dim.cov+3)

a <- a_prior$mean

alpha <- matrix(coef(glm(Yc~offset(log(OFF))+Mp+Yp+Wp, family=poisson)), nrow=K,ncol=dim.cov+3, byrow=TRUE)

#### Starting PMs ############
lambda <- 2
Pp <- rep(1/K,K)
COV.med <- vcov(glm(Yc~offset(log(OFF))+Mp+Yp+Wp, family=poisson))/dim.cov

# Setting the number of MCMC interations
MCMC<-6000
tau <- rep(0.5, dim.cov+3)
COV.a <- rbind(COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med)/100

update3 <- list()
update3[[1]] <- list(Yc=Yc,
                     Mp=Mp,
                     Yp=Yp,
                     Wp=Wp,
                     K=K,
                     Kc=Kp,
                     a=a,
                     a_prior=a_prior,
                     tau=tau,
                     alpha=alpha,
                     lambda=lambda,
                     Pc=Pp,
                     dim.cov=dim.cov,
                     COV.med=COV.med,
                     COV.a=COV.a,
                     OFF=OFF)

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
      eval(parse(text=paste("COV.a[",(6*i-5),":",(6*i),",1:6] <- cov(alpha",i,"_storage)", sep="")))
    }
  }
  
  update3[[c]] <- Cond_Main_O(Yc=Yc,
                              Mp=Mp,
                              Yp=Yp,
                              Wp=Wp,
                              K=K,
                              Kp=update3[[c-1]]$Kc,
                              a=update3[[c-1]]$a,
                              a_prior=a_prior,
                              tau=update3[[c-1]]$tau,
                              alpha=update3[[c-1]]$alpha,
                              lambda=update3[[c-1]]$lambda,
                              Pp=update3[[c-1]]$Pc,
                              dim.cov=dim.cov,
                              COV.med=COV.med,
                              COV.a=COV.a,
                              OFF=OFF)
}

out_time11 <- list()
for(ind in 1:2000){
  IND <- ind+4000
  out_time11[[ind]] <- list(a=update3[[IND]]$a,
                            tau=update3[[IND]]$tau,
                            alpha=update3[[IND]]$alpha,
                            Kc=update3[[IND]]$Kc)
}
rm(update3)


##### data at time 2 : TRT=(0,0) #####

Yc <- Y2[Z2==0]
Mp <- as.numeric(M2[Z2==0])
Yp <- log(Y1[Z2==0]+1.5)
Wp <- cbind(X2,  Master$PctUrban,  log(Master$Ter+1))[Z2==0,]
OFF <- Master$TB5[Z2==0]*3

K <- 20
Kp <- sample(1:K,length(Yc),replace=TRUE)

#### Priors ###############
dim.cov <- dim(Wp)[2]
a_past.sample <- matrix(nrow=2000, ncol=6)
for(ind in 1:2000){
  IND <- ind+0
  a_past.sample[ind,]<-out_time10[[IND]]$a
}
a_new.sample <- matrix(nrow=2000, ncol=6)
for(ind in 1:2000){
  a_new.sample[ind,] <- a_past.sample[ind,]
}
a_prior$P <- 1
a_prior$mean <- apply(a_new.sample, 2, mean)
a_prior$sigma <- cov(a_new.sample)
a <- a_prior$mean
alpha <- matrix(0, nrow=K,ncol=dim.cov+3, byrow=TRUE)

lambda <- 2
Pp <- rep(1/K,K)
COV.med <- vcov(glm(Yc~offset(log(OFF))+Mp+Yp+Wp, family=poisson))/dim.cov
COV.a <- rbind(COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med)/100
tau <- rep(0.5, dim.cov+3)

# Setting the number of MCMC interations
MCMC<-6000

update3 <- list()
update3[[1]] <- list(Yc=Yc,
                    Mp=Mp,
                    Yp=Yp,
                    Wp=Wp,
                    K=K,
                    Kc=Kp,
                    a=a,
                    a_prior=a_prior,
                    tau=tau,
                    alpha=alpha,
                    lambda=lambda,
                    Pc=Pp,
                    dim.cov=dim.cov,
                    COV.med=COV.med,
                    COV.a=COV.a,
                    OFF=OFF)

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
            eval(parse(text=paste("COV.a[",(6*i-5),":",(6*i),",1:6] <- cov(alpha",i,"_storage)", sep="")))
        }
    }
    
    update3[[c]] <- Cond_Main_O(Yc=Yc,
                        Mp=Mp,
                        Yp=Yp,
                        Wp=Wp,
                        K=K,
                        Kp=update3[[c-1]]$Kc,
                        a=update3[[c-1]]$a,
                        a_prior=a_prior,
                        tau=update3[[c-1]]$tau,
                        alpha=update3[[c-1]]$alpha,
                        lambda=update3[[c-1]]$lambda,
                        Pp=update3[[c-1]]$Pc,
                        dim.cov=dim.cov,
                        COV.med=COV.med,
                        COV.a=COV.a,
                        OFF=OFF)
}

out_time20 <- list()
for(ind in 1:2000){
    IND <- ind+4000
    out_time20[[ind]] <- list(a=update3[[IND]]$a,
                                tau=update3[[IND]]$tau,
                                alpha=update3[[IND]]$alpha,
                                Kc=update3[[IND]]$Kc)
}
rm(update3)


##### data at time 2 : TRT=(0,1) #####

Yc <- Y2[Z2==1]
Mp <- as.numeric(M2[Z2==1])
Yp <- log(Y1[Z2==1]+1.5)
Wp <- cbind(X2,  Master$PctUrban,  log(Master$Ter+1))[Z2==1,]
OFF <- Master$TB5[Z2==1]*3
K <- 20
Kp <- sample(1:K,length(Yc),replace=TRUE)

#### Priors ###############
dim.cov <- dim(Wp)[2]
a_past.sample <- matrix(nrow=2000, ncol=6)
for(ind in 1:2000){
  IND <- ind+0
  a_past.sample[ind,]<-out_time10[[IND]]$a
}
a_new.sample <- matrix(nrow=2000, ncol=6)
for(ind in 1:2000){
  a_new.sample[ind,] <- a_past.sample[ind,]
}
a_prior$P <- 1
a_prior$mean <- apply(a_new.sample, 2, mean)
a_prior$sigma <- cov(a_new.sample)

a <- a_prior$mean
alpha <- matrix(0, nrow=K,ncol=dim.cov+3, byrow=TRUE)

lambda <- 2
Pp <- rep(1/K,K)
COV.med <- vcov(glm(Yc~offset(log(OFF))+Mp+Yp+Wp, family=poisson))/dim.cov

COV.a <- rbind(COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med)/100

tau <- rep(0.5, dim.cov+3)

# Setting the number of MCMC interations
MCMC<-6000


update3 <- list()
update3[[1]] <- list(Yc=Yc,
                     Mp=Mp,
                     Yp=Yp,
                     Wp=Wp,
                     K=K,
                     Kc=Kp,
                     a=a,
                     a_prior=a_prior,
                     tau=tau,
                     alpha=alpha,
                     lambda=lambda,
                     Pc=Pp,
                     dim.cov=dim.cov,
                     COV.med=COV.med,
                     COV.a=COV.a,
                     OFF=OFF)



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
      eval(parse(text=paste("COV.a[",(6*i-5),":",(6*i),",1:6] <- cov(alpha",i,"_storage)", sep="")))
    }
  }
  
  
  update3[[c]] <- Cond_Main_O(Yc=Yc,
                              Mp=Mp,
                              Yp=Yp,
                              Wp=Wp,
                              K=K,
                              Kp=update3[[c-1]]$Kc,
                              a=update3[[c-1]]$a,
                              a_prior=a_prior,
                              tau=update3[[c-1]]$tau,
                              alpha=update3[[c-1]]$alpha,
                              lambda=update3[[c-1]]$lambda,
                              Pp=update3[[c-1]]$Pc,
                              dim.cov=dim.cov,
                              COV.med=COV.med,
                              COV.a=COV.a,
                              OFF=OFF)
  print(c);print(update3[[c-1]]$alpha);print(table(update3[[c-1]]$Kc))
}


out_time21 <- list()
for(ind in 1:2000){
  IND <- ind+4000
  out_time21[[ind]] <- list(a=update3[[IND]]$a,
                            tau=update3[[IND]]$tau,
                            alpha=update3[[IND]]$alpha,
                            Kc=update3[[IND]]$Kc)
}
rm(update3)





##### data at time 3 : TRT=(0,0,0) #####
Yc <- Y3[Z3==0]
Mp <- as.numeric(M3[Z3==0])
Yp <- log(Y2[Z3==0]+1.5)
Wp <- cbind(X3,  Master$PctUrban,   log(Master$Ter+1))[Z3==0,]
OFF <- Master$TB7[Z3==0]*3
K <- 20
Kp <- sample(1:K,length(Yc),replace=TRUE)

#### Priors ###############
dim.cov <- dim(Wp)[2]
a_past.sample <- matrix(nrow=2000, ncol=6)
for(ind in 1:2000){
  IND <- ind+0
  a_past.sample[ind,]<-out_time20[[IND]]$a
}
a_new.sample <- matrix(nrow=2000, ncol=6)
for(ind in 1:2000){
  a_new.sample[ind,] <- a_past.sample[ind,]
}
a_prior$P <- 1
a_prior$mean <- apply(a_new.sample, 2, mean)
a_prior$sigma <- cov(a_new.sample)
a <- a_prior$mean

alpha <- matrix(0, nrow=K,ncol=dim.cov+3, byrow=TRUE)


lambda <- 2
Pp <- rep(1/K,K)
COV.med <- vcov(glm(Yc~offset(log(OFF))+Mp+Yp+Wp, family=poisson))/dim.cov
COV.a <- rbind(COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med)/100
tau <- rep(0.5, dim.cov+3)

# Setting the number of MCMC interations
MCMC<-6000

update3 <- list()
update3[[1]] <- list(Yc=Yc,
                    Mp=Mp,
                    Yp=Yp,
                    Wp=Wp,
                    K=K,
                    Kc=Kp,
                    a=a,
                    a_prior=a_prior,
                    tau=tau,
                    alpha=alpha,
                    lambda=lambda,
                    Pc=Pp,
                    dim.cov=dim.cov,
                    COV.med=COV.med,
                    COV.a=COV.a,
                    OFF=OFF)

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
            eval(parse(text=paste("COV.a[",(6*i-5),":",(6*i),",1:6] <- cov(alpha",i,"_storage)", sep="")))
        }
        
    }
    
    update3[[c]] <- Cond_Main_O(Yc=Yc,
    Mp=Mp,
    Yp=Yp,
    Wp=Wp,
    K=K,
    Kp=update3[[c-1]]$Kc,
    a=update3[[c-1]]$a,
    a_prior=a_prior,
    tau=update3[[c-1]]$tau,
    alpha=update3[[c-1]]$alpha,
    lambda=update3[[c-1]]$lambda,
    Pp=update3[[c-1]]$Pc,
    dim.cov=dim.cov,
    COV.med=COV.med,
    COV.a=COV.a,
    OFF=OFF)
}

out_time30 <- list()
for(ind in 1:2000){
    IND <- ind+4000
    out_time30[[ind]] <- list(a=update3[[IND]]$a,
    tau=update3[[IND]]$tau,
    alpha=update3[[IND]]$alpha,
    Kc=update3[[IND]]$Kc)
}
rm(update3)




##### data at time 3 : TRT=(0,0,1) #####

Yc <- Y3[Z3==1]
Mp <- as.numeric(M3[Z3==1])
Yp <- log(Y2[Z3==1]+1.5)
Wp <- cbind(X3,  Master$PctUrban,   log(Master$Ter+1))[Z3==1,]
OFF <- Master$TB7[Z3==1]*3

K <- 20
Kp <- sample(1:K,length(Yc),replace=TRUE)

#### Priors ###############
dim.cov <- dim(Wp)[2]
a_past.sample <- matrix(nrow=2000, ncol=6)
for(ind in 1:2000){
  IND <- ind+0
  a_past.sample[ind,]<-out_time20[[IND]]$a
}
a_new.sample <- matrix(nrow=2000, ncol=6)
for(ind in 1:2000){
  a_new.sample[ind,] <- a_past.sample[ind,]
}
a_prior$P <- 1
a_prior$mean <- apply(a_new.sample, 2, mean)
a_prior$sigma <- cov(a_new.sample)
a <- a_prior$mean

alpha <- matrix(0, nrow=K,ncol=dim.cov+3, byrow=TRUE)

lambda <- 2
Pp <- rep(1/K,K)
COV.med <- vcov(glm(Yc~offset(log(OFF))+Mp+Yp+Wp, family=poisson))/dim.cov
COV.a <- rbind(COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med)/100

tau <- rep(0.5, dim.cov+3)

# Setting the number of MCMC interations
MCMC<-6000

update3 <- list()
update3[[1]] <- list(Yc=Yc,
                     Mp=Mp,
                     Yp=Yp,
                     Wp=Wp,
                     K=K,
                     Kc=Kp,
                     a=a,
                     a_prior=a_prior,
                     tau=tau,
                     alpha=alpha,
                     lambda=lambda,
                     Pc=Pp,
                     dim.cov=dim.cov,
                     COV.med=COV.med,
                     COV.a=COV.a,
                     OFF=OFF)

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
      eval(parse(text=paste("COV.a[",(6*i-5),":",(6*i),",1:6] <- cov(alpha",i,"_storage)", sep="")))
    }
  }
  
  update3[[c]] <- Cond_Main_O(Yc=Yc,
                              Mp=Mp,
                              Yp=Yp,
                              Wp=Wp,
                              K=K,
                              Kp=update3[[c-1]]$Kc,
                              a=update3[[c-1]]$a,
                              a_prior=a_prior,
                              tau=update3[[c-1]]$tau,
                              alpha=update3[[c-1]]$alpha,
                              lambda=update3[[c-1]]$lambda,
                              Pp=update3[[c-1]]$Pc,
                              dim.cov=dim.cov,
                              COV.med=COV.med,
                              COV.a=COV.a,
                              OFF=OFF)
}

out_time31 <- list()
for(ind in 1:2000){
  IND <- ind+4000
  out_time31[[ind]] <- list(a=update3[[IND]]$a,
                            tau=update3[[IND]]$tau,
                            alpha=update3[[IND]]$alpha,
                            Kc=update3[[IND]]$Kc)
}
rm(update3)



##### data at time 4 : TRT=(0,0,0,0) #####

Yc <- Y4[Z4==0]
Mp <- as.numeric(M4[Z4==0])
Yp <- log(Y3[Z4==0]+1.5)
Wp <- cbind(X4,  Master$PctUrban,   log(Master$Ter+1))[Z4==0,]

OFF <- Master$TB9[Z4==0]*3
K <- 20
Kp <- sample(1:K,length(Yc),replace=TRUE)

#### Priors ###############
dim.cov <- dim(Wp)[2]
a_past.sample <- matrix(nrow=2000, ncol=6)
for(ind in 1:2000){
  IND <- ind+0
  a_past.sample[ind,]<-out_time30[[IND]]$a
}
a_new.sample <- matrix(nrow=2000, ncol=6)
for(ind in 1:2000){
  a_new.sample[ind,] <- a_past.sample[ind,]
}
a_prior$P <- 1
a_prior$mean <- apply(a_new.sample, 2, mean)
a_prior$sigma <- cov(a_new.sample)
a <- a_prior$mean

alpha <- matrix(0, nrow=K,ncol=dim.cov+3, byrow=TRUE)

lambda <- 2
Pp <- rep(1/K,K)
COV.med <- vcov(glm(Yc~offset(log(OFF))+Mp+Yp+Wp, family=poisson))/dim.cov
COV.a <- rbind(COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med)/100

tau <- rep(0.5, dim.cov+3)

# Setting the number of MCMC interations
MCMC<-6000

COV.tau <- diag(0.05, dim.cov+3)
update3 <- list()
update3[[1]] <- list(Yc=Yc,
                    Mp=Mp,
                    Yp=Yp,
                    Wp=Wp,
                    K=K,
                    Kc=Kp,
                    a=a,
                    a_prior=a_prior,
                    tau=tau,
                    alpha=alpha,
                    lambda=lambda,
                    Pc=Pp,
                    dim.cov=dim.cov,
                    COV.med=COV.med,
                    COV.a=COV.a,
                    COV.tau=COV.tau,
                    OFF=OFF)

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
            eval(parse(text=paste("COV.a[",(6*i-5),":",(6*i),",1:6] <- cov(alpha",i,"_storage)", sep="")))
        }
        
    }

    update3[[c]] <- Cond_Main_O(Yc=Yc,
    Mp=Mp,
    Yp=Yp,
    Wp=Wp,
    K=K,
    Kp=update3[[c-1]]$Kc,
    a=update3[[c-1]]$a,
    a_prior=a_prior,
    tau=update3[[c-1]]$tau,
    alpha=update3[[c-1]]$alpha,
    lambda=update3[[c-1]]$lambda,
    Pp=update3[[c-1]]$Pc,
    dim.cov=dim.cov,
    COV.med=COV.med,
    COV.a=COV.a,
    OFF=OFF)
}

out_time40 <- list()
for(ind in 1:2000){
    IND <- ind+4000
    out_time40[[ind]] <- list(a=update3[[IND]]$a,
    tau=update3[[IND]]$tau,
    alpha=update3[[IND]]$alpha,
    Kc=update3[[IND]]$Kc)
}
rm(update3)


##### data at time 4 : TRT=(0,0,0,1) #####
Yc <- Y4[Z4==1]
Mp <- as.numeric(M4[Z4==1])
Yp <- log(Y3[Z4==1]+1.5)
Wp <- cbind(X4,  Master$PctUrban,   log(Master$Ter+1))[Z4==1,]
OFF <- Master$TB9[Z4==1]*3

K <- 20
Kp <- sample(1:K,length(Yc),replace=TRUE)

#### Priors ###############
dim.cov <- dim(Wp)[2]
a_past.sample <- matrix(nrow=2000, ncol=6)
for(ind in 1:2000){
  IND <- ind+0
  a_past.sample[ind,]<-out_time30[[IND]]$a
}
a_new.sample <- matrix(nrow=2000, ncol=6)
for(ind in 1:2000){
  a_new.sample[ind,] <- a_past.sample[ind,]
}
a_prior$P <- 1
a_prior$mean <- apply(a_new.sample, 2, mean)
a_prior$sigma <- cov(a_new.sample)

a <- a_prior$mean
alpha <- matrix(0, nrow=K,ncol=dim.cov+3, byrow=TRUE)

lambda <- 2
Pp <- rep(1/K,K)
COV.med <- vcov(glm(Yc~offset(log(OFF))+Mp+Yp+Wp, family=poisson))/dim.cov
COV.a <- rbind(COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med,COV.med)/100

tau <- rep(0.5, dim.cov+3)

# Setting the number of MCMC interations
MCMC<-6000

COV.tau <- diag(0.05, dim.cov+3)

update3 <- list()
update3[[1]] <- list(Yc=Yc,
                     Mp=Mp,
                     Yp=Yp,
                     Wp=Wp,
                     K=K,
                     Kc=Kp,
                     a=a,
                     a_prior=a_prior,
                     tau=tau,
                     alpha=alpha,
                     lambda=lambda,
                     Pc=Pp,
                     dim.cov=dim.cov,
                     COV.med=COV.med,
                     COV.a=COV.a,
                     COV.tau=COV.tau,
                     OFF=OFF)

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
      eval(parse(text=paste("COV.a[",(6*i-5),":",(6*i),",1:6] <- cov(alpha",i,"_storage)", sep="")))
    }
  }
  
  update3[[c]] <- Cond_Main_O(Yc=Yc,
                              Mp=Mp,
                              Yp=Yp,
                              Wp=Wp,
                              K=K,
                              Kp=update3[[c-1]]$Kc,
                              a=update3[[c-1]]$a,
                              a_prior=a_prior,
                              tau=update3[[c-1]]$tau,
                              alpha=update3[[c-1]]$alpha,
                              lambda=update3[[c-1]]$lambda,
                              Pp=update3[[c-1]]$Pc,
                              dim.cov=dim.cov,
                              COV.med=COV.med,
                              COV.a=COV.a,
                              OFF=OFF)
}


out_time41 <- list()
for(ind in 1:2000){
  IND <- ind+4000
  out_time41[[ind]] <- list(a=update3[[IND]]$a,
                            tau=update3[[IND]]$tau,
                            alpha=update3[[IND]]$alpha,
                            Kc=update3[[IND]]$Kc)
}
rm(update3)

save(out_time10,out_time20,out_time30,out_time40,out_time11,out_time21,out_time31,out_time41, file=paste0("output/Sim_Out",process,".RData"))


