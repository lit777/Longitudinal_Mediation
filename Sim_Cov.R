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
source("Cond_K_W.R")
source("Cond_Alpha_W.R")
source("Cond_Main_W.R")

#----- Load MCMC samples and Data
load("simulated_data.RData")
library(data.table)
load(file=paste0("data_",process+1,".RData"))


##### data at time 2 : TRT=(0) #####
Wc <- X2[Z1==0]
Wp <- cbind(Master$TTEMP3,Master$PctUrban, log(Master$Ter+1))[Z1==0,]
Mp <- as.numeric(M1)[Z1==0]

WmisIndex <- which(is.na(Wc))
Wc[WmisIndex] <- rnorm(length(WmisIndex),mean(Wc,na.rm=TRUE),sd(Wc,na.rm=TRUE))
K <- 20
Kp <- sample(1:K,length(Wc),replace=TRUE)

#### Priors ###############
dim.cov <- dim(Wp)[2]
a_prior <- list()
a_prior$P <- 1
a_prior$mean <- rep(0, dim.cov+2)
a_prior$sigma <- diag(1000,dim.cov+2)
a <- a_prior$mean
tau <- rep(0.5,dim.cov+2)
alpha <- matrix(coef(lm(Wc~Mp+Wp)), nrow=K,ncol=dim.cov+2, byrow=TRUE)

sigma <-   rep(1,K)
lambda <- 2
Pp <- rep(1/K,K)
COV.w <- vcov(lm(Wc~Mp+Wp))/(dim.cov)

# Setting the number of MCMC inMaster$Terations
MCMC<-6000

COV.a <- rbind(COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w)/100

update3 <- list()
update3[[1]] <- list(Wc=Wc,
    WmisIndex=WmisIndex,
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
    COV.w=COV.w,
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
    if(c > 200){
        COV.w <- cov(a_storage)
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
    if(c > 200){
        for(i in 1:20){
            eval(parse(text=paste("COV.a[",(5*i-4),":",(5*i),",1:5] <- cov(alpha",i,"_storage)", sep="")))
        }

    }
    update3[[c]] <- Cond_Main_W(Wc=update3[[c-1]]$Wc,
    WmisIndex=WmisIndex,
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
    COV.w=COV.w,
    COV.a=COV.a)
}


cov_time10 <- list()
for(ind in 1:2000){
    IND <- ind+4000
    cov_time10[[ind]] <- list(a=update3[[IND]]$a,
    tau=update3[[IND]]$tau,
    alpha=update3[[IND]]$alpha,
    sigma=update3[[IND]]$sigma,
    Kc=update3[[IND]]$Kc)
}
rm(update3)


##### data at time 3 : TRT=(0,0) #####
Wc <- X3[Z2==0]
Wp <- cbind(X2,Master$PctUrban,  log(Master$Ter+1))[Z2==0,]
Mp <- as.numeric(M2)[Z2==0]
WmisIndex <- which(is.na(Wc))
Wc[WmisIndex] <- rnorm(length(WmisIndex),mean(Wc,na.rm=TRUE),sd(Wc,na.rm=TRUE))

K <- 20
Kp <- sample(1:K,length(Wc),replace=TRUE)

#### Priors ###############
dim.cov <- dim(Wp)[2]

a_past.sample <- matrix(nrow=1000, ncol=5)
for(ind in 1:1000){
    IND <- ind+0
    a_past.sample[ind,]<-cov_time10[[IND]]$a
}
a_new.sample <- matrix(nrow=1000, ncol=5)
for(ind in 1:1000){
    a_new.sample[ind,] <- a_past.sample[ind,]
}
a_prior$P <- 1
a_prior$mean <- apply(a_new.sample, 2, mean)
a_prior$sigma <- cov(a_new.sample)

a <- a_prior$mean
tau <- rep(0.5,dim.cov+2)
alpha <- matrix(0, nrow=K,ncol=dim.cov+2, byrow=TRUE)

sigma <- rep(1,K)
lambda <- 2
Pp <- rep(1/K,K)
COV.w <- vcov(lm(Wc~Mp+Wp))/(dim.cov)

# Setting the number of MCMC interations
MCMC<-6000

COV.a <- rbind(COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w)/100

update3 <- list()
update3[[1]] <- list(Wc=Wc,
WmisIndex=WmisIndex,
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
COV.w=COV.w,
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
    if(c > 200){
        COV.w <- cov(a_storage)
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
    if(c > 200){
        for(i in 1:20){
          eval(parse(text=paste("COV.a[",(5*i-4),":",(5*i),",1:5] <- cov(alpha",i,"_storage)", sep="")))
        }
    }
    
    update3[[c]] <- Cond_Main_W(Wc=update3[[c-1]]$Wc,
    WmisIndex=WmisIndex,
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
    COV.w=COV.w,
    COV.a=COV.a)
}

cov_time20 <- list()
for(ind in 1:2000){
    IND <- ind+4000
    cov_time20[[ind]] <- list(a=update3[[IND]]$a,
    tau=update3[[IND]]$tau,
    alpha=update3[[IND]]$alpha,
    sigma=update3[[IND]]$sigma,
    Kc=update3[[IND]]$Kc)
}
rm(update3)


##### data at time 4 : TRT=(0,0,0) #####

Wc <- X4[Z3==0]
Wp <- cbind(X3,Master$PctUrban,  log(Master$Ter+1))[Z3==0,]
Mp <- as.numeric(M3)[Z3==0]
WmisIndex <- which(is.na(Wc))
Wc[WmisIndex] <- rnorm(length(WmisIndex),mean(Wc,na.rm=TRUE),sd(Wc,na.rm=TRUE))

K <- 20
Kp <- sample(1:K,length(Wc),replace=TRUE)

#### Priors ###############
dim.cov <- dim(Wp)[2]
a_past.sample <- matrix(nrow=1000, ncol=5)
for(ind in 1:1000){
    IND <- ind+0
    a_past.sample[ind,]<-cov_time20[[IND]]$a
}
a_new.sample <- matrix(nrow=1000, ncol=5)
for(ind in 1:1000){
    a_new.sample[ind,] <- a_past.sample[ind,]
}
a_prior$P <- 1
a_prior$mean <- apply(a_new.sample, 2, mean)
a_prior$sigma <- cov(a_new.sample)
a <- a_prior$mean
tau <- rep(0.5,dim.cov+2)
alpha <- matrix(0, nrow=K,ncol=dim.cov+2, byrow=TRUE)

sigma <- rep(1,K)
lambda <- 2
Pp <- rep(1/K,K)
COV.w <- vcov(lm(Wc~Mp+Wp))/(dim.cov)

# Setting the number of MCMC interations
MCMC<-6000

COV.a <- rbind(COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w,COV.w)/100

update3 <- list()
update3[[1]] <- list(Wc=Wc,
WmisIndex=WmisIndex,
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
COV.w=COV.w,
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
    if(c > 200){
        COV.w <- cov(a_storage)
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
    if(c > 200){
        for(i in 1:20){
          eval(parse(text=paste("COV.a[",(5*i-4),":",(5*i),",1:5] <- cov(alpha",i,"_storage)", sep="")))
        }
    }
    update3[[c]] <- Cond_Main_W(Wc=update3[[c-1]]$Wc,
    WmisIndex=WmisIndex,
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
    COV.w=COV.w,
    COV.a=COV.a)
}

cov_time30<- list()
for(ind in 1:2000){
    IND <- ind+4000
    cov_time30[[ind]] <- list(a=update3[[IND]]$a,
    tau=update3[[IND]]$tau,
    alpha=update3[[IND]]$alpha,
    sigma=update3[[IND]]$sigma,
    Kc=update3[[IND]]$Kc)
}
rm(update3)

save(cov_time10,cov_time20,cov_time30, file=paste0("output/Sim_Cov_",process,".RData"))


