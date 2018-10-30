
##---------------------------------------------------------------
## Required libraries
##---------------------------------------------------------------

#----- Parallel computing
library(doParallel)

#----- Make clusters based on the number of CPU cores
cl<-makeCluster(4) # 64 cores
registerDoParallel(cl)
getDoParWorkers()

#----- Support parallel excution
library(foreach)
library(data.table)
library(mnormt)
library(sn)

load("simulated_data.RData")

Y10.T1 <- NULL
Y00.T1 <- NULL
Y11.T1 <- NULL

Y10.T2 <- NULL
Y00.T2 <- NULL
Y11.T2 <- NULL

Y10.T3 <- NULL
Y00.T3 <- NULL
Y11.T3 <- NULL

Y10.T4 <- NULL
Y00.T4 <- NULL
Y11.T4 <- NULL

main <- function(temp){
  
  eval(parse(text=paste0("load('Sim_Cov_",temp-1,".RData')")))
  eval(parse(text=paste0("load('Sim_Med_",temp-1,".RData')")))
  eval(parse(text=paste0("load('Sim_Out_",temp-1,".RData')")))

  W0 <- cbind(Master$TTEMP3,  Master$PctUrban,  log(Master$Ter+1))
  W00 <- cbind(Master$TTEMP1,  Master$PctUrban,  log(Master$Ter+1))
  
  off <- Master$TB3*3
  off0 <- Master$TB1*3
  off1 <- Master$TB3*3
  off2 <- Master$TB5*3
  off3 <- Master$TB7*3
  off4 <- Master$TB9*3
  
  MC.i <- dim(W0)[1]
  
  Y10.t1 <- NULL
  Y00.t1 <- NULL
  Y11.t1 <- NULL
  
  Y10.t2 <- NULL
  Y00.t2 <- NULL
  Y11.t2 <- NULL
  
  Y10.t3 <- NULL
  Y00.t3 <- NULL
  Y11.t3 <- NULL
  
  Y10.t4 <- NULL
  Y00.t4 <- NULL
  Y11.t4 <- NULL
  
  
  for(I in 1:200){
    ind <- I*10
    
    t.ind <- sort(unique(med_time10[[ind]]$Kc))
    prob <- (table(factor(med_time10[[ind]]$Kc, levels=1:20))/length(med_time10[[ind]]$Kc))[t.ind]
    if(length(prob)==1){
        t.ind <- rep(t.ind, MC.i)
    }else{
        t.ind <- sample(sort(unique(med_time10[[ind]]$Kc)),(MC.i), prob=prob, replace=TRUE)
    }

    m1.t1 <- rnorm(MC.i, rowSums(cbind(1,as.numeric(Master$PM1),W0)*(med_time10[[ind]]$alpha[t.ind,])),sqrt((1/med_time10[[ind]]$sigma[t.ind]))  )

    t.ind <- sort(unique(med_time11[[ind]]$Kc))
    prob <- (table(factor(med_time11[[ind]]$Kc, levels=1:20))/length(med_time11[[ind]]$Kc))[t.ind]
    if(length(prob)==1){
        t.ind <- rep(t.ind, MC.i)
    }else{
        t.ind <- sample(sort(unique(med_time11[[ind]]$Kc)),(MC.i), prob=prob, replace=TRUE)
    }

    m1 <- rnorm(MC.i, rowSums(cbind(1,as.numeric(Master$PM1),W0)*(med_time11[[ind]]$alpha[t.ind,])),sqrt((1/med_time11[[ind]]$sigma[t.ind]))  )
    
    t.ind <- sort(unique(cov_time10[[ind]]$Kc))
    
    prob <- (table(factor(cov_time10[[ind]]$Kc, levels=1:20))/length(cov_time10[[ind]]$Kc))[t.ind]
    if(length(prob)==1){
        t.ind <- rep(t.ind, MC.i)
    }else{
        t.ind <- sample(sort(unique(cov_time10[[ind]]$Kc)),(MC.i), prob=prob, replace=TRUE)
    }

    w1 <- rnorm(MC.i, rowSums(cbind(1,m1.t1,W0)*(cov_time10[[ind]]$alpha[t.ind,])),sqrt((1/cov_time10[[ind]]$sigma[t.ind]))  )
    
    W1 <- W0
    W1[,1] <- w1
    
    t.ind <- sort(unique(out_time11[[ind]]$Kc))
    
    
    prob <- (table(factor(out_time11[[ind]]$Kc, levels=1:20))/length(out_time11[[ind]]$Kc))[t.ind]
    if(length(prob)==1){
        t.ind <- rep(t.ind, MC.i)
    }else{
        t.ind <- sample(sort(unique(out_time11[[ind]]$Kc)),(MC.i), prob=prob, replace=TRUE)
    }
    y10.t1 <- rpois(length(m1.t1),exp(log(off1)+rowSums(cbind(1,m1.t1,log(Master$TD1+1.5),W0)*(out_time11[[ind]]$alpha[t.ind,]))))
    Y10.t1[I] <- mean(y10.t1/off1)
    y11.t1 <- rpois(length(m1.t1),exp(log(off1)+rowSums(cbind(1,m1,log(Master$TD1+1.5),W0)*(out_time11[[ind]]$alpha[t.ind,]))))
    Y11.t1[I] <- mean(y11.t1/off1)


    t.ind <- sort(unique(out_time10[[ind]]$Kc))

    prob <- (table(factor(out_time10[[ind]]$Kc, levels=1:20))/length(out_time10[[ind]]$Kc))[t.ind]
    if(length(prob)==1){
        t.ind <- rep(t.ind, MC.i)
    }else{
        t.ind <- sample(sort(unique(out_time10[[ind]]$Kc)),(MC.i), prob=prob, replace=TRUE)
    }
    y00.t1 <- rpois(length(m1.t1),exp(log(off1)+rowSums(cbind(1,m1.t1,log(Master$TD1+1.5),W0)*(out_time10[[ind]]$alpha[t.ind,]))))
    Y00.t1[I] <- mean(y00.t1/off1)



    
    
    ##### T=2
    
    t.ind <- sort(unique(med_time20[[ind]]$Kc))
    prob <- (table(factor(med_time20[[ind]]$Kc, levels=1:20))/length(med_time20[[ind]]$Kc))[t.ind]
    if(length(prob)==1){
        t.ind <- rep(t.ind, MC.i)
    }else{
        t.ind <- sample(sort(unique(med_time20[[ind]]$Kc)),(MC.i), prob=prob, replace=TRUE)
    }

    m2.t2 <- rnorm(MC.i, rowSums(cbind(1,m1.t1,W1)*(med_time20[[ind]]$alpha[t.ind,])),sqrt((1/med_time20[[ind]]$sigma[t.ind]))  )
    
    
    t.ind <- sort(unique(med_time21[[ind]]$Kc))
    prob <- (table(factor(med_time21[[ind]]$Kc, levels=1:20))/length(med_time21[[ind]]$Kc))[t.ind]
    if(length(prob)==1){
        t.ind <- rep(t.ind, MC.i)
    }else{
        t.ind <- sample(sort(unique(med_time21[[ind]]$Kc)),(MC.i), prob=prob, replace=TRUE)
    }

    m2 <- rnorm(MC.i, rowSums(cbind(1,m1.t1,W1)*(med_time21[[ind]]$alpha[t.ind,])),sqrt((1/med_time21[[ind]]$sigma[t.ind]))  )


    t.ind <- sort(unique(out_time21[[ind]]$Kc))
    
    prob <- (table(factor(out_time21[[ind]]$Kc, levels=1:20))/length(out_time21[[ind]]$Kc))[t.ind]
    if(length(prob)==1){
        t.ind <- rep(t.ind, MC.i)
    }else{
        t.ind <- sample(sort(unique(out_time21[[ind]]$Kc)),(MC.i), prob=prob, replace=TRUE)
    }
    y10.t2 <-rpois(length(m1.t1),exp(log(off2)+rowSums(cbind(1,m2.t2,log(y00.t1+1.5),W1)*(out_time21[[ind]]$alpha[t.ind,]))))
    Y10.t2[I] <-mean(y10.t2/off2)
    
    y11.t2 <- rpois(length(m1.t1),exp(log(off2)+rowSums(cbind(1,m2,log(y00.t1+1.5),W1)*(out_time21[[ind]]$alpha[t.ind,]))))
    Y11.t2[I] <- mean(y11.t2/off2)
    
    t.ind <- sort(unique(out_time20[[ind]]$Kc))
    
    prob <- (table(factor(out_time20[[ind]]$Kc, levels=1:20))/length(out_time20[[ind]]$Kc))[t.ind]
    if(length(prob)==1){
        t.ind <- rep(t.ind, MC.i)
    }else{
        t.ind <- sample(sort(unique(out_time20[[ind]]$Kc)),(MC.i), prob=prob, replace=TRUE)
    }
    y00.t2 <- rpois(length(m1.t1),exp(log(off2)+rowSums(cbind(1,m2.t2,log(y00.t1+1.5),W1)*(out_time20[[ind]]$alpha[t.ind,]))))
    Y00.t2[I] <- mean(y00.t2/off2)
    
    
    ####### t=3
    
    t.ind <- sort(unique(cov_time20[[ind]]$Kc))
    
    prob <- (table(factor(cov_time20[[ind]]$Kc, levels=1:20))/length(cov_time20[[ind]]$Kc))[t.ind]
    if(length(prob)==1){
        t.ind <- rep(t.ind, MC.i)
    }else{
        t.ind <- sample(sort(unique(cov_time20[[ind]]$Kc)),(MC.i), prob=prob, replace=TRUE)
    }

    w2 <- rnorm(MC.i, rowSums(cbind(1,m2.t2,W1)*(cov_time20[[ind]]$alpha[t.ind,])),sqrt((1/cov_time20[[ind]]$sigma[t.ind]))  )

    W2 <- W1
    W2[,1] <- w2
    
    t.ind <- sort(unique(med_time30[[ind]]$Kc))
    prob <- (table(factor(med_time30[[ind]]$Kc, levels=1:20))/length(med_time30[[ind]]$Kc))[t.ind]
    if(length(prob)==1){
        t.ind <- rep(t.ind, MC.i)
    }else{
        t.ind <- sample(sort(unique(med_time30[[ind]]$Kc)),(MC.i), prob=prob, replace=TRUE)
    }

    m3.t3 <- rnorm(MC.i, rowSums(cbind(1,m2.t2,W2)*(med_time30[[ind]]$alpha[t.ind,])),sqrt((1/med_time30[[ind]]$sigma[t.ind]))  )
    
    t.ind <- sort(unique(med_time31[[ind]]$Kc))
    prob <- (table(factor(med_time31[[ind]]$Kc, levels=1:20))/length(med_time31[[ind]]$Kc))[t.ind]
    if(length(prob)==1){
        t.ind <- rep(t.ind, MC.i)
    }else{
        t.ind <- sample(sort(unique(med_time31[[ind]]$Kc)),(MC.i), prob=prob, replace=TRUE)
    }
    m3 <- rnorm(MC.i, rowSums(cbind(1,m2.t2,W2)*(med_time31[[ind]]$alpha[t.ind,])),sqrt((1/med_time31[[ind]]$sigma[t.ind]))  )

    t.ind <- sort(unique(out_time31[[ind]]$Kc))
    
    prob <- (table(factor(out_time31[[ind]]$Kc, levels=1:20))/length(out_time31[[ind]]$Kc))[t.ind]
    
    if(length(prob)==1){
        t.ind <- rep(t.ind, MC.i)
    }else{
        t.ind <- sample(sort(unique(out_time31[[ind]]$Kc)),(MC.i), prob=prob, replace=TRUE)
    }
    y10.t3 <- rpois(length(m1.t1),exp(log(off3)+rowSums(cbind(1,m3.t3,log(y00.t2+1.5),W2)*(out_time31[[ind]]$alpha[t.ind,]))))
    Y10.t3[I] <- mean(y10.t3/off3)
    
    y11.t3 <- rpois(length(m1.t1),exp(log(off3)+rowSums(cbind(1,m3,log(y00.t2+1.5),W2)*(out_time31[[ind]]$alpha[t.ind,]))))
    Y11.t3[I] <- mean(y11.t3/off3)

    

    t.ind <- sort(unique(out_time30[[ind]]$Kc))

    prob <- (table(factor(out_time30[[ind]]$Kc, levels=1:20))/length(out_time30[[ind]]$Kc))[t.ind]

    if(length(prob)==1){
        t.ind <- rep(t.ind, MC.i)
    }else{
        t.ind <- sample(sort(unique(out_time30[[ind]]$Kc)),(MC.i), prob=prob, replace=TRUE)
    }

    y00.t3 <- rpois(length(m1.t1),exp(log(off3)+rowSums(cbind(1,m3.t3,log(y00.t2+1.5),W2)*(out_time30[[ind]]$alpha[t.ind,]))))
    Y00.t3[I] <- mean(y00.t3/off3)

    
    ####### t=4
    t.ind <- sort(unique(cov_time30[[ind]]$Kc))
    
    prob <- (table(factor(cov_time30[[ind]]$Kc, levels=1:20))/length(cov_time30[[ind]]$Kc))[t.ind]
    if(length(prob)==1){
        t.ind <- rep(t.ind, MC.i)
    }else{
        t.ind <- sample(sort(unique(cov_time30[[ind]]$Kc)),(MC.i), prob=prob, replace=TRUE)
    }

    w3 <- rnorm(MC.i, rowSums(cbind(1,m3.t3,W2)*(cov_time30[[ind]]$alpha[t.ind,])),sqrt((1/cov_time30[[ind]]$sigma[t.ind]))  )
    
    
    W3 <- W2
    W3[,1] <- w3
    
    
    t.ind <- sort(unique(med_time40[[ind]]$Kc))
    prob <- (table(factor(med_time40[[ind]]$Kc, levels=1:20))/length(med_time40[[ind]]$Kc))[t.ind]
    if(length(prob)==1){
        t.ind <- rep(t.ind, MC.i)
    }else{
    t.ind <- sample(sort(unique(med_time40[[ind]]$Kc)),(MC.i), prob=prob, replace=TRUE)
    }
    m4.t4 <- rnorm(MC.i, rowSums(cbind(1,m3.t3,W3)*(med_time40[[ind]]$alpha[t.ind,])),sqrt((1/med_time40[[ind]]$sigma[t.ind]))  )
    
    
    t.ind <- sort(unique(med_time41[[ind]]$Kc))
    prob <- (table(factor(med_time41[[ind]]$Kc, levels=1:20))/length(med_time41[[ind]]$Kc))[t.ind]
    if(length(prob)==1){
        t.ind <- rep(t.ind, MC.i)
    }else{
        t.ind <- sample(sort(unique(med_time41[[ind]]$Kc)),(MC.i), prob=prob, replace=TRUE)
    }
    m4 <- rnorm(MC.i, rowSums(cbind(1,m3.t3,W3)*(med_time41[[ind]]$alpha[t.ind,])),sqrt((1/med_time41[[ind]]$sigma[t.ind]))  )


    t.ind <- sort(unique(out_time41[[ind]]$Kc))
  
    
    prob <- (table(factor(out_time41[[ind]]$Kc, levels=1:20))/length(out_time41[[ind]]$Kc))[t.ind]
    
    if(length(prob)==1){
        t.ind <- rep(t.ind, MC.i)
    }else{
        t.ind <- sample(sort(unique(out_time41[[ind]]$Kc)),(MC.i), prob=prob, replace=TRUE)
    }
    y10.t4 <- rpois(length(m1.t1),exp(log(off4)+rowSums(cbind(1,m4.t4,log(y00.t3+1.5),W3)*(out_time41[[ind]]$alpha[t.ind,]))))
    Y10.t4[I] <- mean(y10.t4/off4)


    y11.t4 <-rpois(length(m1.t1),exp(log(off4)+rowSums(cbind(1,m4,log(y00.t3+1.5),W3)*(out_time41[[ind]]$alpha[t.ind,]))))
    Y11.t4[I] <- mean(y11.t4/off4)

    t.ind <- sort(unique(out_time40[[ind]]$Kc))


    prob <- (table(factor(out_time40[[ind]]$Kc, levels=1:20))/length(out_time40[[ind]]$Kc))[t.ind]

    if(length(prob)==1){
        t.ind <- rep(t.ind, MC.i)
    }else{
        t.ind <- sample(sort(unique(out_time40[[ind]]$Kc)),(MC.i), prob=prob, replace=TRUE)
    }

    y00.t4 <- rpois(length(m1.t1),exp(log(off4)+rowSums(cbind(1,m4.t4,log(y00.t3+1.5),W3)*(out_time40[[ind]]$alpha[t.ind,]))))
    Y00.t4[I] <- mean(y00.t4/off4)
    }
  
    Y11.T1 <- mean(Y11.t1)
    Y10.T1 <- mean(Y10.t1)
    Y00.T1 <- mean(Y00.t1)
  
    Y11.T2 <- mean(Y11.t2)
    Y10.T2 <- mean(Y10.t2)
    Y00.T2 <- mean(Y00.t2)
  
    Y11.T3 <- mean(Y11.t3)
    Y10.T3 <- mean(Y10.t3)
    Y00.T3 <- mean(Y00.t3)
  
    Y11.T4 <- mean(Y11.t4)
    Y10.T4 <- mean(Y10.t4)
    Y00.T4 <- mean(Y00.t4)
  
  
    nie1 <- Y11.T1-Y10.T1
    nie2 <- Y11.T2-Y10.T2
    nie3 <- Y11.T3-Y10.T3
    nie4 <- Y11.T4-Y10.T4
  
    nde1 <- Y10.T1-Y00.T1
    nde2 <- Y10.T2-Y00.T2
    nde3 <- Y10.T3-Y00.T3
    nde4 <- Y10.T4-Y00.T4
  
    te1 <- nie1 + nde1
    te2 <- nie2 + nde2
    te3 <- nie3 + nde3
    te4 <- nie4 + nde4
  
  return(c(nie1,nie2,nie3,nie4,nde1,nde2,nde3,nde4,te1,te2,te3,te4))
}


result<-foreach(temp = c(1:400), .combine = rbind) %dopar% main(temp)

save(result, file="Sim_Result.RData")
stopCluster(cl)

