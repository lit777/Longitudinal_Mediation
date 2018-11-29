#----- Load MCMC samples and Data
load("simulated_data.RData")
library(data.table)
library(mnormt)
library(sn)

Master$TB1 <- Master$TB1*3
Master$TB3 <- Master$TB3*3
Master$TB5 <- Master$TB5*3
Master$TB7 <- Master$TB7*3
Master$TB9 <- Master$TB9*3


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

# Time 1
X1 <- Master$TTEMP3
W <- cbind(Master$PctUrban,log(Master$Ter+1))
Z1 <- rbinom(size, 1, expit(-1+0.09*Master$TTEMP3+rowSums(0.05*W)))


MM1_1 <- rsn(size,cbind(1, step(Master$PM1), 1, Master$TTEMP3, W)%*%matrix(coef.m1, ncol=1),omega=1, alpha=3)
MM1_0 <- rsn(size,cbind(1, step(Master$PM1), 0, Master$TTEMP3, W)%*%matrix(coef.m1, ncol=1),omega=1, alpha=3)

M1_1 <- (MM1_1)
M1_0 <- (MM1_0)

I <- rbinom(size, 1, 0.5)

X2_1 <- I*rnorm(size,cbind(1,MM1_1,1,X1,W)%*%coef.x2,1)+(1-I)*rnorm(size,cbind(1,MM1_1,1,X1,W)%*%coef.x21,0.5)
X2_0 <- I*rnorm(size,cbind(1,MM1_0,0,X1,W)%*%coef.x2,1)+(1-I)*rnorm(size,cbind(1,MM1_0,0,X1,W)%*%coef.x21,0.5)

Y1_1 <- rpois(size, exp(log(Master$TB3)+cbind(1,M1_1,1,log(Master$TD1+1.5),X1,W, M1_1*1,X1*M1_1)%*%matrix(coef.y1, ncol=1)))

Y1_c <- rpois(size, exp(log(Master$TB3)+cbind(1,M1_0,1,log(Master$TD1+1.5),X1,W, M1_0*1,X1*M1_0)%*%matrix(coef.y1, ncol=1)))

Y1_0 <- rpois(size, exp(log(Master$TB3)+cbind(1,M1_0,0,log(Master$TD1+1.5),X1,W, M1_0*0,X1*M1_0)%*%matrix(coef.y1, ncol=1)))

M1 <- Z1*MM1_1+(1-Z1)*MM1_0
Y1 <- Z1*Y1_1+(1-Z1)*Y1_0
X2 <- Z1*X2_1+(1-Z1)*X2_0

# Time 2

Z2 <- rbinom(size, 1, expit(-0.5+0.1*Z1+0.09*X2+0.08*Master$TTEMP3+rowSums(0.05*W)))

MM2_11 <- rsn(size,cbind(1, MM1_1, 1, X2_1, (Master$PM1), 1, Master$TTEMP3,W)%*%matrix(coef.m2, ncol=1),omega=1, alpha=-3)
MM2_10 <- rsn(size,cbind(1, MM1_1, 0, X2_1, (Master$PM1), 1, Master$TTEMP3,W)%*%matrix(coef.m2, ncol=1),omega=1, alpha=-3)

MM2_01 <- rsn(size,cbind(1, MM1_0, 1, X2_0, (Master$PM1), 0, Master$TTEMP3, W)%*%matrix(coef.m2, ncol=1),omega=1, alpha=-3)
MM2_00 <- rsn(size,cbind(1, MM1_0, 0, X2_0, (Master$PM1), 0, Master$TTEMP3, W)%*%matrix(coef.m2, ncol=1),omega=1, alpha=-3)

M2_11 <- (MM2_11)
M2_10 <- (MM2_10)
M2_01 <- (MM2_01)
M2_00 <- (MM2_00)

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


# Time 3

Z3 <- rbinom(size, 1, expit(0.2+0.1*Z2-0.05*X3+0.03*Master$TTEMP3+rowSums(0.05*W)))

MM3_111 <- rsn(size,cbind(1, MM2_11, 1, X3_11, MM1_1, 1, X2_1, (Master$PM1), 1, Master$TTEMP3,W)%*%matrix(coef.m3, ncol=1),omega=1, alpha=3)
MM3_110 <- rsn(size,cbind(1, MM2_11, 0, X3_11, MM1_1, 1, X2_1, (Master$PM1), 1, Master$TTEMP3,W)%*%matrix(coef.m3, ncol=1),omega=1, alpha=3)

MM3_101 <- rsn(size,cbind(1, MM2_10, 1, X3_10, MM1_1, 0, X2_1, (Master$PM1), 1, Master$TTEMP3, W)%*%matrix(coef.m3, ncol=1),omega=1, alpha=3)
MM3_100 <- rsn(size,cbind(1, MM2_10, 0, X3_10, MM1_1, 0, X2_1, (Master$PM1), 1, Master$TTEMP3,W)%*%matrix(coef.m3, ncol=1),omega=1, alpha=3)

MM3_011 <- rsn(size,cbind(1, MM2_01, 1, X3_01, MM1_0, 1, X2_0, (Master$PM1), 0, Master$TTEMP3, W)%*%matrix(coef.m3, ncol=1),omega=1, alpha=3)
MM3_010 <- rsn(size,cbind(1, MM2_01, 0, X3_01, MM1_0, 1, X2_0, (Master$PM1), 0, Master$TTEMP3,W)%*%matrix(coef.m3, ncol=1),omega=1, alpha=3)

MM3_001 <- rsn(size,cbind(1, MM2_00, 1, X3_00, MM1_0, 0, X2_0, (Master$PM1), 0, Master$TTEMP3,W)%*%matrix(coef.m3, ncol=1),omega=1, alpha=3)
MM3_000 <- rsn(size,cbind(1, MM2_00, 0, X3_00, MM1_0, 0, X2_0, (Master$PM1), 0, Master$TTEMP3, W)%*%matrix(coef.m3, ncol=1),omega=1, alpha=3)

M3_111 <- (MM3_111)
M3_110 <- (MM3_110)
M3_101 <- (MM3_101)
M3_100 <- (MM3_100)
M3_011 <- (MM3_011)
M3_010 <- (MM3_010)
M3_001 <- (MM3_001)
M3_000 <- (MM3_000)

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


# Time 4

Z4 <- rbinom(size, 1, expit(0.8+0.1*Z3+0.05*X4+0.01*Master$TTEMP3+rowSums(0.05*W)))

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


M4_1000 <- (MM4_1000)
M4_0100 <- (MM4_0100)
M4_0010 <- (MM4_0010)
M4_0001 <- (MM4_0001)

M4_1100 <- (MM4_1100)
M4_0110 <- (MM4_0110)
M4_0011 <- (MM4_0011)
M4_1001 <- (MM4_1001)

M4_1010 <- (MM4_1010)
M4_0101 <- (MM4_0101)
M4_1110 <- (MM4_1110)
M4_0111 <- (MM4_0111)

M4_1011 <- (MM4_1011)
M4_1101 <- (MM4_1101)
M4_1111 <- (MM4_1111)
M4_0000 <- (MM4_0000)

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



save(Y1,Y2,Y3,Y4,M1,M2,M3,M4,X1,X2,X3,X4, Z1,Z2,Z3,Z4, file=paste0("simulated_data/data_",i,".RData"))

}



