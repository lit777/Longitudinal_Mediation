##---------------------------------------------------------------
## Required libraries
##---------------------------------------------------------------
library(pscl)
library(mnormt)
library(DPpackage)
library(data.table)
library(MCMCpack)

#----- Load MCMC samples and Data
load("/users/chanminkim/dropbox/Longitudinal\ Code_SO2/temp_data_long_2months_Res.RData")
library(data.table)
library(mnormt)
library(sn)


# ------ Define 7 regions by state names
Northeast = c("ME", "NH", "VT", "NY", "PA", "DE", "NJ", "MD", "DC", "VA", "MA", "CT", "RI")
IndustrialMidwest = c("WV", "OH", "KY", "IN", "IL", "WI", "MI")
Southeast = c("FL", "GA", "SC", "NC", "TN", "AL", "MS", "AR","LA")
UpperMidwest = c("MN", "IA", "MO", "KS", "NE", "SD", "ND")


#Master <- subset(Master, State.zip %in% c(Northeast, IndustrialMidwest, Southeast))
#------ After narrowing down to 3 EPA regions, there are 18573 zip codes left in the data set.

Master <- subset(Master, State.zip %in% c(Northeast))

cutoff <- summary(c(Master$Trt0,Master$Trt1,Master$Trt2,Master$Trt3,Master$Trt4, Master$Trt5, Master$Trt6, Master$Trt7, Master$Trt8, Master$Trt9, Master$Trt10))[3]

Master$TRT0 <- ifelse(Master$Trt0 < cutoff, 1, 0)
Master$TRT1 <- ifelse(Master$Trt1 < cutoff, 1, 0)
Master$TRT2 <- ifelse(Master$Trt2 < cutoff, 1, 0)
Master$TRT3 <- ifelse(Master$Trt3 < cutoff, 1, 0)
Master$TRT4 <- ifelse(Master$Trt4 < cutoff, 1, 0)
Master$TRT5 <- ifelse(Master$Trt5 < cutoff, 1, 0)
Master$TRT6 <- ifelse(Master$Trt6 < cutoff, 1, 0)
Master$TRT7 <- ifelse(Master$Trt7 < cutoff, 1, 0)
Master$TRT8 <- ifelse(Master$Trt8 < cutoff, 1, 0)
Master$TRT9 <- ifelse(Master$Trt9 < cutoff, 1, 0)
Master$TRT10 <- ifelse(Master$Trt10 < cutoff, 1, 0)

Master <- subset(Master, TB0 > 200)
Master <- subset(Master, TB1 > 200)
Master <- subset(Master, TB2 > 200)
Master <- subset(Master, TB3 > 200)
Master <- subset(Master, TB4 > 200)
Master <- subset(Master, TB5 > 200)
Master <- subset(Master, TB6 > 200)
Master <- subset(Master, TB7 > 200)
Master <- subset(Master, TB8 > 200)
Master <- subset(Master, TB9 > 200)
Master <- subset(Master, TB10 > 200)


Master <- subset(Master, TD0/TB0 < 0.02)
Master <- subset(Master, TD1/TB1 < 0.02)
Master <- subset(Master, TD2/TB2 < 0.02)
Master <- subset(Master, TD3/TB3 < 0.02)
Master <- subset(Master, TD4/TB4 < 0.02)
Master <- subset(Master, TD5/TB5 < 0.02)
Master <- subset(Master, TD6/TB6 < 0.02)
Master <- subset(Master, TD7/TB7 < 0.02)
Master <- subset(Master, TD8/TB8 < 0.02)
Master <- subset(Master, TD9/TB9 < 0.02)
Master <- subset(Master, TD10/TB10 < 0.02)



Master <- subset(Master, !is.na(as.numeric(PM0)))
Master <- subset(Master, !is.na(as.numeric(PM1)))
Master <- subset(Master, !is.na(as.numeric(PM2)))
Master <- subset(Master, !is.na(as.numeric(PM3)))
Master <- subset(Master, !is.na(as.numeric(PM4)))
Master <- subset(Master, !is.na(as.numeric(PM5)))
Master <- subset(Master, !is.na(as.numeric(PM6)))
Master <- subset(Master, !is.na(as.numeric(PM7)))
Master <- subset(Master, !is.na(as.numeric(PM8)))
Master <- subset(Master, !is.na(as.numeric(PM9)))
Master <- subset(Master, !is.na(as.numeric(PM10)))


Master$TTEMP1 <- Master$TTEMP1-273.15
Master$TTEMP2 <- Master$TTEMP2-273.15
Master$TTEMP3 <- Master$TTEMP3-273.15
Master$TTEMP4 <- Master$TTEMP4-273.15
Master$TTEMP5 <- Master$TTEMP5-273.15
Master$TTEMP6 <- Master$TTEMP6-273.15
Master$TTEMP7 <- Master$TTEMP7-273.15
Master$TTEMP8 <- Master$TTEMP8-273.15
Master$TTEMP9 <- Master$TTEMP9-273.15
Master$TTEMP10 <- Master$TTEMP10-273.15

Master$PM0 <- as.numeric(Master$PM0)
Master$PM1 <- as.numeric(Master$PM1)
Master$PM2 <- as.numeric(Master$PM2)
Master$PM3 <- as.numeric(Master$PM3)
Master$PM4 <- as.numeric(Master$PM4)
Master$PM5 <- as.numeric(Master$PM5)
Master$PM6 <- as.numeric(Master$PM6)
Master$PM7 <- as.numeric(Master$PM7)
Master$PM8 <- as.numeric(Master$PM8)
Master$PM9 <- as.numeric(Master$PM9)
Master$PM10 <- as.numeric(Master$PM10)

Master$TB0 <- Master$TB0/1000
Master$TB1 <- Master$TB1/1000
Master$TB2 <- Master$TB2/1000
Master$TB3 <- Master$TB3/1000
Master$TB4 <- Master$TB4/1000
Master$TB5 <- Master$TB5/1000
Master$TB6 <- Master$TB6/1000
Master$TB7 <- Master$TB7/1000
Master$TB8 <- Master$TB8/1000
Master$TB9 <- Master$TB9/1000
Master$TB10 <- Master$TB10/1000


#----- Load MCMC samples and Data
load("/users/chanminkim/dropbox/Longitudinal\ Code_SO2/simulation1_new/simulated_data.RData")
library(data.table)
library(mnormt)
library(sn)

Master$TB1 <- Master$TB1*3
Master$TB3 <- Master$TB3*3
Master$TB5 <- Master$TB5*3
Master$TB7 <- Master$TB7*3
Master$TB9 <- Master$TB9*3


step <- function(x){
 # x <- ifelse(x<0, 0.1, x)
#  x <- ifelse(x<5, 0, ifelse(x <10, log(x), x-10+log(10)))
  return(x)
}


for(i in 1:400){


fit1.x2 <- lm(TTEMP5~PM3+TRT3+TTEMP3+PctUrban+log(Ter+1), data= Master)
fit1.x3 <- lm(TTEMP7~PM5+TRT5+TTEMP5+PM3+TRT3+TTEMP3+PctUrban+log(Ter+1), data= Master)
fit1.x4 <- lm(TTEMP9~PM7+TRT7+TTEMP7+PM5+TRT5+TTEMP5+PM3+TRT3+TTEMP3+PctUrban+log(Ter+1), data= Master)


fit2.x2 <- lm(TTEMP5~PM3+TRT3+TTEMP3+PctUrban+log(Ter+1), data= Master)
fit2.x3 <- lm(TTEMP7~PM5+TRT5+TTEMP5+PctUrban+log(Ter+1), data= Master)
fit2.x4 <- lm(TTEMP9~PM7+TRT7+TTEMP7+PctUrban+log(Ter+1), data= Master)


fit1.m1 <- lm(PM3~PM1+TRT3+TTEMP3+PctUrban+log(Ter+1), data=Master)
fit1.m2 <- lm(PM5~PM3+TRT5+TTEMP5+PM1+TRT3+TTEMP3+ PctUrban+log(Ter+1), data=Master)
fit1.m3 <- lm(PM7~PM5+TRT7+TTEMP7+PM3+TRT5+TTEMP5+PM1+TRT3+TTEMP3+ PctUrban+log(Ter+1), data=Master)
fit1.m4 <- lm(PM9~PM7+TRT9+TTEMP9+PM5+TRT7+TTEMP7+PM3+TRT5+TTEMP5+PM1+TRT3+TTEMP3+PctUrban+log(Ter+1), data=Master)

Master$TDB1 <- Master$TD1/(Master$TB1)
Master$TDB3 <- Master$TD3/(Master$TB3)
Master$TDB5 <- Master$TD5/(Master$TB5)
Master$TDB7 <- Master$TD7/(Master$TB7)

int11 <- (Master$PM3)*Master$TRT3
int12 <- Master$TTEMP3*step(Master$PM3)

int21 <- (Master$PM5)*Master$TRT5
int22 <- Master$TTEMP5*step(Master$PM5)

int31 <- (Master$PM7)*Master$TRT7
int32 <- Master$TTEMP7*step(Master$PM7)

int41 <- (Master$PM9)*Master$TRT9
int42 <- Master$TTEMP9*step(Master$PM9)


fit1.y1 <- glm(TD3~offset(log(TB3))+step(PM3)+TRT3+log(TD1+1.5)+TTEMP3+PctUrban+log(Ter+1)+int11+int12, data=Master, family=poisson)
#  fit1.y2 <- glm(TD5~offset(log(TB5))+step(PM5)+TRT5+log(TD3+1.5)+TTEMP5+step(PM3)+TRT3+log(TD1+1.5)+TTEMP3+PctUrban+log(Ter+1)+int11+int2, data=Master, family=poisson)
#  fit1.y3 <- glm(TD7~offset(log(TB7))+step(PM7)+TRT7+log(TD5+1.5)+TTEMP7+step(PM5)+TRT5+log(TD3+1.5)+TTEMP5+step(PM3)+TRT3+log(TD1+1.5)+TTEMP3+PctUrban+log(Ter+1), data=Master, family=poisson)
#  fit1.y4 <- glm(TD9~offset(log(TB9))+step(PM9)+TRT9+log(TD7+1.5)+TTEMP9+step(PM7)+TRT7+log(TD5+1.5)+TTEMP7+step(PM5)+TRT5+log(TD3+1.5)+TTEMP5+step(PM3)+TRT3+log(TD1+1.5)+TTEMP3+PctUrban+log(Ter+1), data=Master, family=poisson)

coef.x2 <- coef(fit1.x2)
coef.x2 <- c(coef.x2[1:4], coef.x2[-(1:4)])
coef.x3 <- c(coef.x2[1],coef.x2[2:4]*0.85, coef.x2[2:4]/10, coef.x2[-(1:4)])
coef.x4 <- c(coef.x2[1],coef.x2[2:4]*0.85*0.85, coef.x2[2:4]/10, coef.x2[2:4]/50, coef.x2[-(1:4)])


coef.x21 <- coef(fit1.x2)
coef.x21 <- c(-coef.x21[1:4]/2, -coef.x21[-(1:4)])
coef.x31 <- c(coef.x21[1],coef.x21[2:4]*0.85, coef.x21[2:4]/10, coef.x21[-(1:4)])
coef.x41 <- c(coef.x21[1],coef.x21[2:4]*0.85*0.85, coef.x21[2:4]/10, coef.x21[2:4]/50, coef.x2[-(1:4)])


coef.m1 <- coef(fit1.m1)
coef.m2 <- c(coef.m1[1],coef.m1[2:4]*0.85, coef.m1[2:4]/10, coef.m1[-(1:4)])
coef.m3 <- c(coef.m1[1],coef.m1[2:4]*0.85*0.85, coef.m1[2:4]/10, coef.m1[2:4]/50, coef.m1[-(1:4)])
coef.m4 <- c(coef.m1[1],coef.m1[2:4]*0.85*0.85*0.85, coef.m1[2:4]/10, coef.m1[2:4]/50, coef.m1[2:4]/100, coef.m1[-(1:4)])

coef.y1 <- coef(fit1.y1)


coef.y1 <- c(coef(fit1.y1)[1]/4,coef(fit1.y1)[2]*2,-coef(fit1.y1)[3]*2.5,coef(fit1.y1)[4]/10,coef(fit1.y1)[c(5:7)]/2.5, coef(fit1.y1)[c(8:9)]*2)
coef.y2 <- c(coef.y1[1],coef.y1[2:5]*0.85, coef.y1[2:5]/20, coef.y1[-c(1:5)])
coef.y3 <- c(coef.y1[1],coef.y1[2:5]*0.85*0.85, coef.y1[2:5]/20, coef.y1[2:5]/50, coef.y1[-c(1:5)])
coef.y4 <- c(coef.y1[1],coef.y1[2:5]*0.85*0.85*0.85, coef.y1[2:5]/20, coef.y1[2:5]/50, coef.y1[2:5]/100, coef.y1[-c(1:5)])

coef.y11 <- coef.y1
coef.y21 <- coef.y2
coef.y31 <- coef.y3
coef.y41 <- coef.y4



# <------ Data Generation
expit <- function(x){
    exp(x)/(1+exp(x))
}


size <- dim(Master)[1]

X1 <- Master$TTEMP3
W <- cbind(Master$PctUrban,log(Master$Ter+1))




Z1 <- rbinom(size, 1, expit(-1+0.09*Master$TTEMP3+rowSums(0.05*W)))
sum(Z1)




MM1_1 <- rsn(size,cbind(1, step(Master$PM1), 1, Master$TTEMP3, W)%*%matrix(coef.m1, ncol=1),omega=1, alpha=3)
MM1_0 <- rsn(size,cbind(1, step(Master$PM1), 0, Master$TTEMP3, W)%*%matrix(coef.m1, ncol=1),omega=1, alpha=3)

M1_1 <- step(MM1_1)
M1_0 <- step(MM1_0)

I <- rbinom(size, 1, 0.5)

X2_1 <- I*rnorm(size,cbind(1,MM1_1,1,X1,W)%*%coef.x2,1)+(1-I)*rnorm(size,cbind(1,MM1_1,1,X1,W)%*%coef.x21,0.5)
X2_0 <- I*rnorm(size,cbind(1,MM1_0,0,X1,W)%*%coef.x2,1)+(1-I)*rnorm(size,cbind(1,MM1_0,0,X1,W)%*%coef.x21,0.5)


Y1_1 <- rpois(size, exp(log(Master$TB3)+cbind(1,M1_1,1,log(Master$TD1+1.5),X1,W, M1_1*1,X1*M1_1)%*%matrix(coef.y1, ncol=1)))

Y1_c <- rpois(size, exp(log(Master$TB3)+cbind(1,M1_0,1,log(Master$TD1+1.5),X1,W, M1_0*1,X1*M1_0)%*%matrix(coef.y1, ncol=1)))

Y1_0 <- rpois(size, exp(log(Master$TB3)+cbind(1,M1_0,0,log(Master$TD1+1.5),X1,W, M1_0*0,X1*M1_0)%*%matrix(coef.y1, ncol=1)))



M1 <- Z1*MM1_1+(1-Z1)*MM1_0
Y1 <- Z1*Y1_1+(1-Z1)*Y1_0
X2 <- Z1*X2_1+(1-Z1)*X2_0


Z2 <- rbinom(size, 1, expit(-0.5+0.1*Z1+0.09*X2+0.08*Master$TTEMP3+rowSums(0.05*W)))
sum(Z2)





MM2_11 <- rsn(size,cbind(1, MM1_1, 1, X2_1, (Master$PM1), 1, Master$TTEMP3,W)%*%matrix(coef.m2, ncol=1),omega=1, alpha=-3)
MM2_10 <- rsn(size,cbind(1, MM1_1, 0, X2_1, (Master$PM1), 1, Master$TTEMP3,W)%*%matrix(coef.m2, ncol=1),omega=1, alpha=-3)

MM2_01 <- rsn(size,cbind(1, MM1_0, 1, X2_0, (Master$PM1), 0, Master$TTEMP3, W)%*%matrix(coef.m2, ncol=1),omega=1, alpha=-3)
MM2_00 <- rsn(size,cbind(1, MM1_0, 0, X2_0, (Master$PM1), 0, Master$TTEMP3, W)%*%matrix(coef.m2, ncol=1),omega=1, alpha=-3)

M2_11 <- step(MM2_11)
M2_10 <- step(MM2_10)
M2_01 <- step(MM2_01)
M2_00 <- step(MM2_00)

I <- rbinom(size, 1, 0.5)

X3_11 <- I*rnorm(size,cbind(1,MM2_11,1,X2_1,MM1_1,1,X1,W)%*%coef.x3,0.5)+(1-I)*rnorm(size,cbind(1,MM2_11,1,X2_1,MM1_1,1,X1,W)%*%coef.x31,0.5)
X3_10 <- I*rnorm(size,cbind(1,MM2_10,0,X2_1,MM1_1,1,X1,W)%*%coef.x3,0.5)+(1-I)*rnorm(size,cbind(1,MM2_10,0,X2_1,MM1_1,1,X1,W)%*%coef.x31,0.5)

X3_01 <- I*rnorm(size,cbind(1,MM2_01,1,X2_0,MM1_0,0,X1,W)%*%coef.x3,0.5)+(1-I)*rnorm(size,cbind(1,MM2_01,1,X2_0,MM1_0,0,X1,W)%*%coef.x31,0.5)
X3_00 <- I*rnorm(size,cbind(1,MM2_00,0,X2_0,MM1_0,0,X1,W)%*%coef.x3,0.5)+(1-I)*rnorm(size,cbind(1,MM2_00,0,X2_0,MM1_0,0,X1,W)%*%coef.x31,0.5)


Y2_11 <- rpois(size, exp(log(Master$TB5)+cbind(1,M2_11,1,log(Y1_1+1.5),X2_1,M1_1,1,log(Master$TD1+1.5),X1,W, M2_11*1, X2_1*M2_11)%*%matrix(coef.y2, ncol=1)))
Y2_10 <- rpois(size, exp(log(Master$TB5)+cbind(1,M2_10,0,log(Y1_1+1.5),X2_1,M1_1,1,log(Master$TD1+1.5),X1,W, M2_10*0, X2_1*M2_10)%*%matrix(coef.y2, ncol=1)))
Y2_01 <- rpois(size, exp(log(Master$TB5)+cbind(1,M2_01,1,log(Y1_0+1.5),X2_0,M1_0,0,log(Master$TD1+1.5),X1,W, M2_01*1, X2_0*M2_01)%*%matrix(coef.y2, ncol=1)))

Y2_0c <- rpois(size, exp(log(Master$TB5)+cbind(1,M2_00,1,log(Y1_0+1.5),X2_0,M1_0,0,log(Master$TD1+1.5),X1,W, M2_00*1, X2_0*M2_00)%*%matrix(coef.y2, ncol=1)))

Y2_00 <- rpois(size, exp(log(Master$TB5)+cbind(1,M2_00,0,log(Y1_0+1.5),X2_0,M1_0,0,log(Master$TD1+1.5),X1,W, M2_00*0, X2_0*M2_00)%*%matrix(coef.y2, ncol=1)))



Y2 <- Z2*Z1*(Y2_11)+Z2*(1-Z1)*Y2_01+(1-Z2)*(1-Z1)*Y2_00+(1-Z2)*(Z1)*Y2_10
M2 <- Z2*Z1*(MM2_11)+Z2*(1-Z1)*MM2_01+(1-Z2)*(1-Z1)*MM2_00+(1-Z2)*(Z1)*MM2_10
X3 <- Z2*Z1*(X3_11)+Z2*(1-Z1)*X3_01+(1-Z2)*(1-Z1)*X3_00+(1-Z2)*(Z1)*X3_10


Z3 <- rbinom(size, 1, expit(0.2+0.1*Z2-0.05*X3+0.03*Master$TTEMP3+rowSums(0.05*W)))
sum(Z3)



MM3_111 <- rsn(size,cbind(1, MM2_11, 1, X3_11, MM1_1, 1, X2_1, (Master$PM1), 1, Master$TTEMP3,W)%*%matrix(coef.m3, ncol=1),omega=1, alpha=3)
MM3_110 <- rsn(size,cbind(1, MM2_11, 0, X3_11, MM1_1, 1, X2_1, (Master$PM1), 1, Master$TTEMP3,W)%*%matrix(coef.m3, ncol=1),omega=1, alpha=3)

MM3_101 <- rsn(size,cbind(1, MM2_10, 1, X3_10, MM1_1, 0, X2_1, (Master$PM1), 1, Master$TTEMP3, W)%*%matrix(coef.m3, ncol=1),omega=1, alpha=3)
MM3_100 <- rsn(size,cbind(1, MM2_10, 0, X3_10, MM1_1, 0, X2_1, (Master$PM1), 1, Master$TTEMP3,W)%*%matrix(coef.m3, ncol=1),omega=1, alpha=3)

MM3_011 <- rsn(size,cbind(1, MM2_01, 1, X3_01, MM1_0, 1, X2_0, (Master$PM1), 0, Master$TTEMP3, W)%*%matrix(coef.m3, ncol=1),omega=1, alpha=3)
MM3_010 <- rsn(size,cbind(1, MM2_01, 0, X3_01, MM1_0, 1, X2_0, (Master$PM1), 0, Master$TTEMP3,W)%*%matrix(coef.m3, ncol=1),omega=1, alpha=3)

MM3_001 <- rsn(size,cbind(1, MM2_00, 1, X3_00, MM1_0, 0, X2_0, (Master$PM1), 0, Master$TTEMP3,W)%*%matrix(coef.m3, ncol=1),omega=1, alpha=3)
MM3_000 <- rsn(size,cbind(1, MM2_00, 0, X3_00, MM1_0, 0, X2_0, (Master$PM1), 0, Master$TTEMP3, W)%*%matrix(coef.m3, ncol=1),omega=1, alpha=3)

M3_111 <- step(MM3_111)
M3_110 <- step(MM3_110)
M3_101 <- step(MM3_101)
M3_100 <- step(MM3_100)
M3_011 <- step(MM3_011)
M3_010 <- step(MM3_010)
M3_001 <- step(MM3_001)
M3_000 <- step(MM3_000)

I <- rbinom(size, 1, 0.5)


X4_111 <- I*rnorm(size,cbind(1,MM3_111,1,X3_11,MM2_11,1,X2_1,MM1_1,1,X1,W)%*%coef.x4,0.5)+(1-I)*rnorm(size,cbind(1,MM3_111,1,X3_11,MM2_11,1,X2_1,MM1_1,1,X1,W)%*%coef.x41,0.5)
X4_110 <- I*rnorm(size,cbind(1,MM3_110,0,X3_11,MM2_11,1,X2_1,MM1_1,1,X1,W)%*%coef.x4,0.5)+(1-I)*rnorm(size,cbind(1,MM3_110,0,X3_11,MM2_11,1,X2_1,MM1_1,1,X1,W)%*%coef.x41,0.5)

X4_101 <- I*rnorm(size,cbind(1,MM3_101,1,X3_10,MM2_10,0,X2_1,MM1_1,1,X1,W)%*%coef.x4,0.5)+(1-I)*rnorm(size,cbind(1,MM3_101,1,X3_10,MM2_10,0,X2_1,MM1_1,1,X1,W)%*%coef.x41,0.5)
X4_100 <- I*rnorm(size,cbind(1,MM3_100,0,X3_10,MM2_10,0,X2_1,MM1_1,1,X1,W)%*%coef.x4,0.5)+(1-I)*rnorm(size,cbind(1,MM3_100,0,X3_10,MM2_10,0,X2_1,MM1_1,1,X1,W)%*%coef.x41,0.5)

X4_011 <- I*rnorm(size,cbind(1,MM3_011,1,X3_01,MM2_01,1,X2_0,MM1_0,0,X1,W)%*%coef.x4,0.5)+(1-I)*rnorm(size,cbind(1,MM3_011,1,X3_01,MM2_01,1,X2_0,MM1_0,0,X1,W)%*%coef.x41,0.5)
X4_010 <- I*rnorm(size,cbind(1,MM3_010,0,X3_01,MM2_01,1,X2_0,MM1_0,0,X1,W)%*%coef.x4,0.5)+(1-I)*rnorm(size,cbind(1,MM3_010,0,X3_01,MM2_01,1,X2_0,MM1_0,0,X1,W)%*%coef.x41,0.5)

X4_001 <- I*rnorm(size,cbind(1,MM3_001,1,X3_00,MM2_00,0,X2_0,MM1_0,0,X1,W)%*%coef.x4,0.5)+(1-I)*rnorm(size,cbind(1,MM3_001,1,X3_00,MM2_00,0,X2_0,MM1_0,0,X1,W)%*%coef.x41,0.5)
X4_000 <- I*rnorm(size,cbind(1,MM3_000,0,X3_00,MM2_00,0,X2_0,MM1_0,0,X1,W)%*%coef.x4,0.5)+(1-I)*rnorm(size,cbind(1,MM3_000,0,X3_00,MM2_00,0,X2_0,MM1_0,0,X1,W)%*%coef.x41,0.5)




Y3_111 <- rpois(size, exp(log(Master$TB7)+cbind(1,M3_111,1,log(Y2_11+1.5),X3_11,M2_11,1,log(Y1_1+1.5),X2_1,M1_1,1,log(Master$TD1+1.5),X1,W, M3_111*1, X3_11*M3_111)%*%matrix(coef.y3, ncol=1)))
Y3_110 <- rpois(size, exp(log(Master$TB7)+cbind(1,M3_110,0,log(Y2_11+1.5),X3_11,M2_11,1,log(Y1_1+1.5),X2_1,M1_1,1,log(Master$TD1+1.5),X1,W, M3_110*0, X3_11*M3_110)%*%matrix(coef.y3, ncol=1)))

Y3_101 <- rpois(size, exp(log(Master$TB7)+cbind(1,M3_101,1,log(Y2_10+1.5),X3_10,M2_10,0,log(Y1_1+1.5),X2_1,M1_1,1,log(Master$TD1+1.5),X1,W, M3_101*1, X3_10*M3_101)%*%matrix(coef.y3, ncol=1)))
Y3_100 <- rpois(size, exp(log(Master$TB7)+cbind(1,M3_100,0,log(Y2_10+1.5),X3_10,M2_10,0,log(Y1_1+1.5),X2_1,M1_1,1,log(Master$TD1+1.5),X1,W, M3_100*0, X3_10*M3_100)%*%matrix(coef.y3, ncol=1)))

Y3_011 <- rpois(size, exp(log(Master$TB7)+cbind(1,M3_011,1,log(Y2_01+1.5),X3_01,M2_01,1,log(Y1_0+1.5),X2_0,M1_0,0,log(Master$TD1+1.5),X1,W, M3_011*1, X3_01*M3_011)%*%matrix(coef.y3, ncol=1)))
Y3_010 <- rpois(size, exp(log(Master$TB7)+cbind(1,M3_010,0,log(Y2_01+1.5),X3_01,M2_01,1,log(Y1_0+1.5),X2_0,M1_0,0,log(Master$TD1+1.5),X1,W, M3_010*0, X3_01*M3_010)%*%matrix(coef.y3, ncol=1)))

Y3_001 <- rpois(size, exp(log(Master$TB7)+cbind(1,M3_001,1,log(Y2_00+1.5),X3_00,M2_00,0,log(Y1_0+1.5),X2_0,M1_0,0,log(Master$TD1+1.5),X1,W, M3_001*1, X3_00*M3_001)%*%matrix(coef.y3, ncol=1)))

Y3_00c <- rpois(size, exp(log(Master$TB7)+cbind(1,M3_000,1,log(Y2_00+1.5),X3_00,M2_00,0,log(Y1_0+1.5),X2_0,M1_0,0,log(Master$TD1+1.5),X1,W, M3_000*1, X3_00*M3_000)%*%matrix(coef.y3, ncol=1)))

Y3_000 <- rpois(size, exp(log(Master$TB7)+cbind(1,M3_000,0,log(Y2_00+1.5),X3_00,M2_00,0,log(Y1_0+1.5),X2_0,M1_0,0,log(Master$TD1+1.5),X1,W, M3_000*0, X3_00*M3_000)%*%matrix(coef.y3, ncol=1)))



Y3 <- Z3*Z2*Z1*(Y3_111)+Z3*(1-Z2)*Z1*(Y3_101)+Z3*Z2*(1-Z1)*(Y3_011)+(1-Z3)*Z2*Z1*(Y3_110)+
(1-Z3)*(1-Z2)*Z1*(Y3_100)+(1-Z3)*Z2*(1-Z1)*(Y3_010)+Z3*(1-Z2)*(1-Z1)*(Y3_001)+(1-Z3)*(1-Z2)*(1-Z1)*(Y3_000)

M3 <- Z3*Z2*Z1*(MM3_111)+Z3*(1-Z2)*Z1*(MM3_101)+Z3*Z2*(1-Z1)*(MM3_011)+(1-Z3)*Z2*Z1*(MM3_110)+
(1-Z3)*(1-Z2)*Z1*(MM3_100)+(1-Z3)*Z2*(1-Z1)*(MM3_010)+Z3*(1-Z2)*(1-Z1)*(MM3_001)+(1-Z3)*(1-Z2)*(1-Z1)*(MM3_000)

X4 <- Z3*Z2*Z1*(X4_111)+Z3*(1-Z2)*Z1*(X4_101)+Z3*Z2*(1-Z1)*(X4_011)+(1-Z3)*Z2*Z1*(X4_110)+
(1-Z3)*(1-Z2)*Z1*(X4_100)+(1-Z3)*Z2*(1-Z1)*(X4_010)+Z3*(1-Z2)*(1-Z1)*(X4_001)+(1-Z3)*(1-Z2)*(1-Z1)*(X4_000)



Z4 <- rbinom(size, 1, expit(0.8+0.1*Z3+0.05*X4+0.01*Master$TTEMP3+rowSums(0.05*W)))
sum(Z4)



MM4_1111 <- rsn(size,cbind(1, MM3_111, 1, X4_111, MM2_11, 1, X3_11, MM1_1, 1, X2_1, (Master$PM1), 1, Master$TTEMP3,W)%*%matrix(coef.m4, ncol=1),omega=1, alpha=-3)
MM4_1110 <- rsn(size,cbind(1, MM3_111, 0, X4_111, MM2_11, 1, X3_11, MM1_1, 1, X2_1, (Master$PM1), 1, Master$TTEMP3,W)%*%matrix(coef.m4, ncol=1),omega=1, alpha=-3)

MM4_1101 <- rsn(size,cbind(1, MM3_110, 1, X4_110, MM2_11, 0, X3_11, MM1_1, 1, X2_1, (Master$PM1), 1, Master$TTEMP3,W)%*%matrix(coef.m4, ncol=1),omega=1, alpha=-3)
MM4_1100 <- rsn(size,cbind(1, MM3_110, 0, X4_110, MM2_11, 0, X3_11, MM1_1, 1, X2_1, (Master$PM1), 1, Master$TTEMP3,W)%*%matrix(coef.m4, ncol=1),omega=1, alpha=-3)

MM4_1011 <- rsn(size,cbind(1, MM3_101, 1, X4_101, MM2_10, 1, X3_10,MM1_1, 0, X2_1, (Master$PM1), 1, Master$TTEMP3, W)%*%matrix(coef.m4, ncol=1),omega=1, alpha=-3)
MM4_1010 <- rsn(size,cbind(1, MM3_101, 0, X4_101, MM2_10, 1, X3_10,MM1_1, 0, X2_1, (Master$PM1), 1, Master$TTEMP3, W)%*%matrix(coef.m4, ncol=1),omega=1, alpha=-3)

MM4_1001 <- rsn(size,cbind(1, MM3_100, 1, X4_100, MM2_10, 0, X3_10, MM1_1, 0, X2_1, (Master$PM1), 1, Master$TTEMP3,W)%*%matrix(coef.m4, ncol=1),omega=1, alpha=-3)
MM4_1000 <- rsn(size,cbind(1, MM3_100, 0, X4_100, MM2_10, 0, X3_10, MM1_1, 0, X2_1, (Master$PM1), 1, Master$TTEMP3,W)%*%matrix(coef.m4, ncol=1),omega=1, alpha=-3)

MM4_0111 <- rsn(size,cbind(1, MM3_011, 1, X4_011,MM2_01, 1, X3_01,MM1_0, 1, X2_0, (Master$PM1), 0, Master$TTEMP3, W)%*%matrix(coef.m4, ncol=1),omega=1, alpha=-3)
MM4_0110 <- rsn(size,cbind(1, MM3_011, 0, X4_011,MM2_01, 1, X3_01,MM1_0, 1, X2_0, (Master$PM1), 0, Master$TTEMP3, W)%*%matrix(coef.m4, ncol=1),omega=1, alpha=-3)

MM4_0101 <- rsn(size,cbind(1, MM3_010, 1, X4_010, MM2_01, 0, X3_01, MM1_0, 1, X2_0, (Master$PM1), 0, Master$TTEMP3,W)%*%matrix(coef.m4, ncol=1),omega=1, alpha=-3)
MM4_0100 <- rsn(size,cbind(1, MM3_010, 0, X4_010, MM2_01, 0, X3_01, MM1_0, 1, X2_0, (Master$PM1), 0, Master$TTEMP3,W)%*%matrix(coef.m4, ncol=1),omega=1, alpha=-3)

MM4_0011 <- rsn(size,cbind(1, MM3_001, 1, X4_001, MM2_00,1, X3_00, MM1_0, 0, X2_0, (Master$PM1), 0, X1, W)%*%matrix(coef.m4, ncol=1),omega=1, alpha=-3)

MM4_0010 <- rsn(size,cbind(1, MM3_001, 0, X4_001, MM2_00,1, X3_00, MM1_0, 0, X2_0, (Master$PM1), 0, X1, W)%*%matrix(coef.m4, ncol=1),omega=1, alpha=-3)

MM4_0001 <- rsn(size,cbind(1, MM3_000, 1, X4_000, MM2_00, 0, X3_00, MM1_0, 0, X2_0, (Master$PM1), 0, Master$TTEMP3,W)%*%matrix(coef.m4, ncol=1),omega=1, alpha=-3)
MM4_0000 <- rsn(size,cbind(1, MM3_000, 0, X4_000, MM2_00, 0, X3_00, MM1_0, 0, X2_0, (Master$PM1), 0, Master$TTEMP3,W)%*%matrix(coef.m4, ncol=1),omega=1, alpha=-3)



M4_1000 <- step(MM4_1000)
M4_0100 <- step(MM4_0100)
M4_0010 <- step(MM4_0010)
M4_0001 <- step(MM4_0001)

M4_1100 <- step(MM4_1100)
M4_0110 <- step(MM4_0110)
M4_0011 <- step(MM4_0011)
M4_1001 <- step(MM4_1001)

M4_1010 <- step(MM4_1010)
M4_0101 <- step(MM4_0101)
M4_1110 <- step(MM4_1110)
M4_0111 <- step(MM4_0111)

M4_1011 <- step(MM4_1011)
M4_1101 <- step(MM4_1101)
M4_1111 <- step(MM4_1111)
M4_0000 <- step(MM4_0000)




Y4_1111 <- rpois(size, exp(log(Master$TB9)+cbind(1,M4_1111,1,log(Y3_111+1.5),X4_111,M3_111,1,log(Y2_11+1.5),X3_11,M2_11,1,log(Y1_1+1.5),X2_1,M1_1,1,log(Master$TD1+1.5),X1,W, M4_1111*1, X4_111*M4_1111)%*%matrix(coef.y4, ncol=1)))
Y4_1110 <- rpois(size, exp(log(Master$TB9)+cbind(1,M4_1110,0,log(Y3_111+1.5),X4_111,M3_111,1,log(Y2_11+1.5),X3_11,M2_11,1,log(Y1_1+1.5),X2_1,M1_1,1,log(Master$TD1+1.5),X1,W, M4_1110*0, X4_111*M4_1110)%*%matrix(coef.y4, ncol=1)))

Y4_1101 <- rpois(size, exp(log(Master$TB9)+cbind(1,M4_1101,1,log(Y3_110+1.5),X4_110,M3_110,0,log(Y2_11+1.5),X3_11,M2_11,1,log(Y1_1+1.5),X2_1,M1_1,1,log(Master$TD1+1.5),X1,W, M4_1101*1, X4_110*M4_1101)%*%matrix(coef.y4, ncol=1)))
Y4_1100 <- rpois(size, exp(log(Master$TB9)+cbind(1,M4_1100,0,log(Y3_110+1.5),X4_110,M3_110,0,log(Y2_11+1.5),X3_11,M2_11,1,log(Y1_1+1.5),X2_1,M1_1,1,log(Master$TD1+1.5),X1,W, M4_1100*0, X4_110*M4_1100)%*%matrix(coef.y4, ncol=1)))

Y4_1011 <- rpois(size, exp(log(Master$TB9)+cbind(1,M4_1011,1,log(Y3_101+1.5),X4_101,M3_101,1,log(Y2_10+1.5),X3_10,M2_10,0,log(Y1_1+1.5),X2_1,M1_1,1,log(Master$TD1+1.5),X1,W, M4_1011*1, X4_101*M4_1011)%*%matrix(coef.y4, ncol=1)))
Y4_1010 <- rpois(size, exp(log(Master$TB9)+cbind(1,M4_1010,0,log(Y3_101+1.5),X4_101,M3_101,1,log(Y2_10+1.5),X3_10,M2_10,0,log(Y1_1+1.5),X2_1,M1_1,1,log(Master$TD1+1.5),X1,W, M4_1010*0, X4_101*M4_1010)%*%matrix(coef.y4, ncol=1)))

Y4_1001 <- rpois(size, exp(log(Master$TB9)+cbind(1,M4_1001,1,log(Y3_100+1.5),X4_100,M3_100,0,log(Y2_10+1.5),X3_10,M2_10,0,log(Y1_1+1.5),X2_1,M1_1,1,log(Master$TD1+1.5),X1,W, M4_1001*1, X4_100*M4_1001)%*%matrix(coef.y4, ncol=1)))
Y4_1000 <- rpois(size, exp(log(Master$TB9)+cbind(1,M4_1000,0,log(Y3_100+1.5),X4_100,M3_100,0,log(Y2_10+1.5),X3_10,M2_10,0,log(Y1_1+1.5),X2_1,M1_1,1,log(Master$TD1+1.5),X1,W, M4_1000*0, X4_100*M4_1000)%*%matrix(coef.y4, ncol=1)))

Y4_0111 <- rpois(size, exp(log(Master$TB9)+cbind(1,M4_0111,1,log(Y3_011+1.5),X4_011,M3_011,1,log(Y2_01+1.5),X3_01,M2_01,1,log(Y1_0+1.5),X2_0,M1_0,0,log(Master$TD1+1.5),X1,W, M4_0111*1, X4_011*M4_0111)%*%matrix(coef.y4, ncol=1)))
Y4_0110 <- rpois(size, exp(log(Master$TB9)+cbind(1,M4_0110,0,log(Y3_011+1.5),X4_011,M3_011,1,log(Y2_01+1.5),X3_01,M2_01,1,log(Y1_0+1.5),X2_0,M1_0,0,log(Master$TD1+1.5),X1,W, M4_0110*0, X4_011*M4_0110)%*%matrix(coef.y4, ncol=1)))

Y4_0101 <- rpois(size, exp(log(Master$TB9)+cbind(1,M4_0101,1,log(Y3_010+1.5),X4_010,M3_010,0,log(Y2_01+1.5),X3_01,M2_01,1,log(Y1_0+1.5),X2_0,M1_0,0,log(Master$TD1+1.5),X1,W, M4_0101*1, X4_010*M4_0101)%*%matrix(coef.y4, ncol=1)))
Y4_0100 <- rpois(size, exp(log(Master$TB9)+cbind(1,M4_0100,0,log(Y3_010+1.5),X4_010,M3_010,0,log(Y2_01+1.5),X3_01,M2_01,1,log(Y1_0+1.5),X2_0,M1_0,0,log(Master$TD1+1.5),X1,W, M4_0100*0, X4_010*M4_0100)%*%matrix(coef.y4, ncol=1)))

Y4_0011 <- rpois(size, exp(log(Master$TB9)+cbind(1,M4_0011,1,log(Y3_001+1.5),X4_001,M3_001,1,log(Y2_00+1.5),X3_00,M2_00,0,log(Y1_0+1.5),X2_0,M1_0,0,log(Master$TD1+1.5),X1,W, M4_0011*1, X4_001*M4_0011)%*%matrix(coef.y4, ncol=1)))
Y4_0010 <- rpois(size, exp(log(Master$TB9)+cbind(1,M4_0010,0,log(Y3_001+1.5),X4_001,M3_001,1,log(Y2_00+1.5),X3_00,M2_00,0,log(Y1_0+1.5),X2_0,M1_0,0,log(Master$TD1+1.5),X1,W, M4_0010*0, X4_001*M4_0010)%*%matrix(coef.y4, ncol=1)))

Y4_0001 <- rpois(size, exp(log(Master$TB9)+cbind(1,M4_0001,1,log(Y3_000+1.5),X4_000,M3_000,0,log(Y2_00+1.5),X3_00,M2_00,0,log(Y1_0+1.5),X2_0,M1_0,0,log(Master$TD1+1.5),X1,W, M4_0001*1, X4_000*M4_0001)%*%matrix(coef.y4, ncol=1)))
Y4_000c <- rpois(size, exp(log(Master$TB9)+cbind(1,M4_0000,1,log(Y3_000+1.5),X4_000,M3_000,0,log(Y2_00+1.5),X3_00,M2_00,0,log(Y1_0+1.5),X2_0,M1_0,0,log(Master$TD1+1.5),X1,W, M4_0000*1, X4_000*M4_0000)%*%matrix(coef.y4, ncol=1)))
Y4_0000 <- rpois(size, exp(log(Master$TB9)+cbind(1,M4_0000,0,log(Y3_000+1.5),X4_000,M3_000,0,log(Y2_00+1.5),X3_00,M2_00,0,log(Y1_0+1.5),X2_0,M1_0,0,log(Master$TD1+1.5),X1,W, M4_0000*0, X4_000*M4_0000)%*%matrix(coef.y4, ncol=1)))


Y4 <- Z1*Z2*Z3*Z4*Y4_1111+
(1-Z1)*Z2*Z3*Z4*Y4_0111+Z1*(1-Z2)*Z3*Z4*Y4_1011+Z1*Z2*(1-Z3)*Z4*Y4_1101+Z1*Z2*Z3*(1-Z4)*Y4_1110+
(1-Z1)*(1-Z2)*Z3*Z4*Y4_0011+(1-Z1)*Z2*(1-Z3)*Z4*Y4_0101+(1-Z1)*Z2*Z3*(1-Z4)*Y4_0110+Z1*(1-Z2)*(1-Z3)*Z4*Y4_1001+Z1*(1-Z2)*Z3*(1-Z4)*Y4_1010+Z1*Z2*(1-Z3)*(1-Z4)*Y4_1100+
Z1*(1-Z2)*(1-Z3)*(1-Z4)*Y4_1000+(1-Z1)*Z2*(1-Z3)*(1-Z4)*Y4_0100+(1-Z1)*(1-Z2)*Z3*(1-Z4)*Y4_0010+(1-Z1)*(1-Z2)*(1-Z3)*Z4*Y4_0001+(1-Z1)*(1-Z2)*(1-Z3)*(1-Z4)*Y4_0000

M4 <- Z1*Z2*Z3*Z4*MM4_1111+
(1-Z1)*Z2*Z3*Z4*MM4_0111+Z1*(1-Z2)*Z3*Z4*MM4_1011+Z1*Z2*(1-Z3)*Z4*MM4_1101+Z1*Z2*Z3*(1-Z4)*MM4_1110+
(1-Z1)*(1-Z2)*Z3*Z4*MM4_0011+(1-Z1)*Z2*(1-Z3)*Z4*MM4_0101+(1-Z1)*Z2*Z3*(1-Z4)*MM4_0110+Z1*(1-Z2)*(1-Z3)*Z4*MM4_1001+Z1*(1-Z2)*Z3*(1-Z4)*MM4_1010+Z1*Z2*(1-Z3)*(1-Z4)*MM4_1100+
Z1*(1-Z2)*(1-Z3)*(1-Z4)*MM4_1000+(1-Z1)*Z2*(1-Z3)*(1-Z4)*MM4_0100+(1-Z1)*(1-Z2)*Z3*(1-Z4)*MM4_0010+(1-Z1)*(1-Z2)*(1-Z3)*Z4*MM4_0001+(1-Z1)*(1-Z2)*(1-Z3)*(1-Z4)*MM4_0000




summary(Y1)
summary(Y2)
summary(Y3)
summary(Y4)

sum(Z1)
sum(Z2)
sum(Z3)
sum(Z4)

mean(Y1_1/Master$TB3-Y1_c/Master$TB3)
mean(Y2_01/Master$TB5-Y2_0c/Master$TB5)
mean(Y3_001/Master$TB7-Y3_00c/Master$TB7)
mean(Y4_0001/Master$TB9-Y4_000c/Master$TB9)

mean(Y1_c/Master$TB3-Y1_0/Master$TB3)
mean(Y2_0c/Master$TB5-Y2_00/Master$TB5)
mean(Y3_00c/Master$TB7-Y3_000/Master$TB7)
mean(Y4_000c/Master$TB9-Y4_0000/Master$TB9)


# Per 10,000

save(Y1,Y2,Y3,Y4,M1,M2,M3,M4,X1,X2,X3,X4, Z1,Z2,Z3,Z4, file=paste0("~/dropbox/longitudinal code_so2/simulation1_new/simulated_data/data_",i,".RData"))

}



