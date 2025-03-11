rm(list = ls())
library(MASS)
library(mnormt)
library(nnet)
#install.packages('LaplacesDemon')
library(LaplacesDemon)
#install.packages('Matching')
library(Matching)
library(dplyr)
library(cluster)
library(e1071)
#install.packages('ramsvm')
library(ramsvm)
library(caret)
library(survival)
library(parallel)
#install.packages('snowfall')
library(snowfall)
library(randomForestSRC)



set.seed(2023)
outcome_setting <- "coxph" 
set.seed(2023)

if (!(outcome_setting %in% c("lognormal", "coxph", "coxph2"))) {
  print("No such setting!")
}
cen_coxph <- 2.7  # length of study
cen_coxph2 <- 9.1
cen_lognormal <- 11.7

theta_coxph <- 0.2  # for censoring rate
theta_coxph2 <- 0.1
theta_lognormal <- 0.08

if (outcome_setting == "coxph") {
  cen <- cen_coxph
  theta <- theta_coxph
} else if (outcome_setting == "coxph2") {
  cen <- cen_coxph2
  theta <- theta_coxph2
} else {
  cen <- cen_lognormal
  theta <- theta_lognormal
}

####################################################################################
# Toy dataset for multiple sensitive features (S1, S2)
####################################################################################
n <- 1000
S1 <- rbinom(n, 1, 0.5)
S2 <- rbinom(n, 1, 0.5)
X11 <- rnorm(sum(S1), 1, 1/1.5)
X10 <- rnorm(n - sum(S1), -1, 1/1.5)
X21 <- rnorm(sum(S2), 1, 1/1.5)
X20 <- rnorm(n - sum(S2), -1, 1/1.5)
X3 <- rnorm(n, 0, 1)
X <- matrix(NA, nrow = n, ncol = 3)
X[S1 == 1, 1] <- X11
X[S1 == 0, 1] <- X10
X[S2 == 1, 2] <- X21
X[S2 == 0, 2] <- X20
X[, 3] <- X3

expit <- function(x) { exp(x) / (1 + exp(x)) }
pa <- expit(X3)
A <- rbinom(n, 1, pa)
e <- rnorm(n, 0, 0.1)
tau <- X[,1]^3/4 + X[,2]/2

# Define R1 and R0
R1 <- X[,1]^3 / 4 + X[,2] / 2 + X[,3] + e
R0 <- X[,3] + e

# Generate R based on A
R_continuous <- A * R1 + (1 - A) * R0

# Define outcome_setting
outcome_setting <- "coxph"  # survival outcome 

R <- rep(0, n)

for (i in 1:n) {
  if (outcome_setting == "coxph") {
    # Scale parameter influenced by R_continuous
    scale_param <- exp(-(R_continuous[i]))^(-1/2)
    R[i] <- min(rweibull(1, shape = 2, scale = scale_param), cen)
  } else if (outcome_setting == "coxph2") {
    R[i] <- min(coxph2(A = A[i], R_continuous = R_continuous[i]), cen)
  } else {
    R[i] <- min(exp(R_continuous[i] + rnorm(1)), cen)
  }
}

censoringtime <- rexp(n, theta)
event <- as.numeric(R <= censoringtime)
observeR <- R * event + censoringtime * (1 - event)

train = as.data.frame(cbind(X,S1,S2,A,censoringtime,event,observeR))
colnames(train) <- c("X1", "X2", "X3", "S1", "S2", "A", "C","censor","observeR")

condition_survival <- function(S.hat, Y.grid) {
  Y.diff <- diff(c(0, Y.grid))
  Q.hat <- matrix(NA, nrow(S.hat), ncol(S.hat))
  dot.products <- sweep(S.hat[, 1:(ncol(S.hat) - 1)], 2, Y.diff[2:ncol(S.hat)], "*")
  Q.hat[, 1] <- rowSums(dot.products)
  for (i in 2:(ncol(Q.hat) - 1)) {
    Q.hat[, i] <- Q.hat[, i - 1] - dot.products[, i - 1]
  }
  Q.hat <- Q.hat/S.hat
  Q.hat <- sweep(Q.hat, 2, Y.grid, "+") 
  Q.hat[, ncol(Q.hat)] <- cen
  
  Q.hat0 <- matrix(expected_survival(S.hat, Y.grid), nrow(S.hat), 1)
  Q.hat <- cbind(Q.hat0, Q.hat)
  return(Q.hat)
}

#imputation E[T|A,X]
expected_survival <- function(S.hat, Y.grid) {
  grid.diff <- diff(c(0, Y.grid, max(Y.grid)))
  return(c(cbind(1, S.hat) %*% grid.diff))
}

traindata <- train
rsf_fit <- rfsrc(Surv(observeR, censor) ~ ., data = traindata, block.size = 1)

#calculate imputation by Cui et al., 2017, but not used in the paper
Q.hat <- condition_survival(S.hat = rsf_fit$survival.oob[which(traindata$censor == 0), ], Y.grid = rsf_fit$time.interest)
impute_all <- rep(NA, nrow(Q.hat))
Y <- traindata[which(traindata$censor == 0), ]$observeR 
for (i in 1:nrow(Q.hat)) {
  Y.index <- findInterval(Y[i], rsf_fit$time.interest, rightmost.closed = T) + 1
  Q.Y.hat <- Q.hat[i, Y.index]
  impute_all[i] <- Q.Y.hat
}

#imputation by RSF
impute_R1 <- expected_survival(S.hat = rsf_fit$survival.oob, Y.grid = rsf_fit$time.interest)
simudata <- cbind(traindata, imputeR1 = impute_R1, imputeR2 = traindata$observeR, treat =traindata$A )
simudata$imputeR2[which(traindata$censor == 0)] <- impute_all  

#we use the imputation by RSF in the paper, and only replace the censored observations, i.e. R*=ΔT+(1-Δ)E(T|A,X)
simudata <- cbind(simudata, impute = simudata$imputeR1) 
simudata$impute[which(traindata$censor == 1)] <- simudata$observeR[which(traindata$censor == 1)]  

dat= simudata
colnames(dat) <- c("X1", "X2", "X3", "S1", "S2", "A","C", "censor", "observeR", "imputeR1", "imputeR2", "treat", "Y")
dat <- dat[, -c(7:12)]

table(dat$S1)
#install.packages('latex2exp')
library(latex2exp)
#install.packages('matrixcalc')
library(matrixcalc)
#install.packages('mbend')
library(mbend)
library(Matrix)
#install.packages('Rmosek')
#library(Rmosek)
library(ggplot2)
#install.packages('ranger')
library(ranger)
library(quadprog)
source("/Users/.../utils.R")



####################################################################################
# fair and robust CATE estimation
####################################################################################

nm.s.features <- c("S1")
#install.packages('Rmosek')
#Rmosek::mosek_attachbuilder("/Users/..../mosek/10.1/tools/platform/osxaarch64/bin")
#install.rmosek()
#install.packages('ranger')
#library(ranger)
#install.packages('Matrix')
#Sys.getenv("MOSEKLM_LICENSE_FILE")
#Sys.setenv(MOSEKLM_LICENSE_FILE ="..../Downloads/mosek/mosek.lic")
#library(Matrix)
#source("/Users/..../utils.R")



# Estimate CATE with fairness constraints
fr.cate0.1 <-fr_cate(dat, nm.s.features = c("S1"),
                     nm.l.factors = c(NA,NA),
                     solver="quadprog",
                     nuisance.est='rf',
                     delta=0.005, sqr=TRUE, interactions=TRUE)
b.mat0.1<- fr.cate0.1$b.mat
beta.hat0.1 <- fr.cate0.1$beta.hat
tau.hat0.005 <- b.mat0.1 %*% beta.hat0.1
tau0.005 = tau.hat0.005 




fr.cate10 <-fr_cate(dat, nm.s.features = c("S1"),
                    nm.l.factors = c(NA,NA),
                    solver="quadprog",
                    nuisance.est='rf',
                    delta=4, sqr=TRUE, interactions=TRUE)
b.mat10<- fr.cate10$b.mat
beta.hat10 <- fr.cate10$beta.hat
tau.hat10 <- b.mat10 %*% beta.hat10
tau10 = tau.hat10



###画出对比图
## CATE density and S1 
pal1 <- c("#FF9933", "#9966FF")
#pal2 <- c("#99CC99", "#CCCCCC")
d1.s0.inf <- density(tau10[S1==0])
d1.s1.inf <- density(tau10[S1==1])
#d2.s0.inf <- density(tau10[S2==0])
#d2.s1.inf <- density(tau10[S2==1])
y.max <- max(c(d1.s0.inf$y,d1.s1.inf$y))
y.min <- min(c(d1.s0.inf$y,d1.s1.inf$y))
x.max <- max(c(d1.s0.inf$x,d1.s1.inf$x))
x.min <- min(c(d1.s0.inf$x,d1.s1.inf$x))
den <- density(tau10)
par(mfrow=c(2,1), mar=c(5, 4, 4, 2) + 0.1)
# CATE density and S1 
plot(den, xlab = TeX(r'(${\hat{\widetilde{\tau}}}$)'), ylab="", main="", cex.main=1.5, cex.lab=1.5)
tau1 <- tau10[S1==1]
tau0 <- tau10[S1==0]
points(tau1, rep(-0.002,length(tau1)), col=alpha(pal1[1],0.5), pch=19)
points(tau0, rep(0.084,length(tau0)), col=alpha(pal1[2],0.5), pch=19)
legend("topleft", legend= TeX(r'($\delta= \infty$)'))
legend("topright", legend=c(TeX(r'($S=1$)'), TeX(r'($S=0$)')),
       col=pal1, cex=1, pch=19)
abline(v=weighted.mean(d1.s0.inf$x, d1.s0.inf$y), lty=2, lwd=3, col=pal1[2], xlim=c(x.min,x.max))
abline(v=weighted.mean(d1.s1.inf$x, d1.s1.inf$y), lty=3, lwd=3, col=pal1[1], xlim=c(x.min,x.max))
#fair
d1.s0.inf <- density(tau0.005[S1==0])
d1.s1.inf <- density(tau0.005[S1==1])
#d2.s0.inf <- density(tau0.005[S2==0])
#d2.s1.inf <- density(tau0.005[S2==1])
y.max <- max(c(d1.s0.inf$y,d1.s1.inf$y))
y.min <- min(c(d1.s0.inf$y,d1.s1.inf$y))
x.max <- max(c(d1.s0.inf$x,d1.s1.inf$x))
x.min <- min(c(d1.s0.inf$x,d1.s1.inf$x))
den <- density(tau0.005)
# CATE density and S1 
plot(den, xlab = TeX(r'(${\hat{\widetilde{\tau}}}$)'), ylab="", main="", cex.main=1.5, cex.lab=1.5)
tau1 <- tau0.005[S1==1]
tau0 <- tau0.005[S1==0]
points(tau1, rep(-0.002,length(tau1)), col=alpha(pal1[1],0.5), pch=19)
points(tau0, rep(0.092,length(tau0)), col=alpha(pal1[2],0.5), pch=19)
legend("topleft", legend= TeX(r'($\delta= 0$)'))
legend("topright", legend=c(TeX(r'($S=1$)'), TeX(r'($S=0$)')),
       col=pal1, cex=1, pch=19)
abline(v=weighted.mean(d1.s0.inf$x, d1.s0.inf$y), lty=2, lwd=3, col=pal1[2], xlim=c(x.min,x.max))
abline(v=weighted.mean(d1.s1.inf$x, d1.s1.inf$y), lty=3, lwd=3, col=pal1[1], xlim=c(x.min,x.max))
####


# bar plot: OTR vs sensitive features
par(mfrow=c(2,1))
par(mar = c(3.5, 4.5, 3, 7))
# S1
ts1 = sum(S1==1 & tau10 > 0)
ts0 = sum(S1==0 & tau10 > 0)
uts1 = sum(S1==1 & tau10 < 0)
uts0 = sum(S1==0 & tau10 < 0)
data <- matrix(c(ts1, ts0, uts1, uts0) , nrow=2)
data_percentage <- apply(data, 2, function(x){x*100/sum(x,na.rm=T)})
barplot(data_percentage, col=pal1 , border="white", ylab="% of each subgroup", cex.names=1.5, cex.lab=1.5, cex.main=2,
        names=c(TeX(r'($\hat{g}=1(\hat{\widetilde{\tau}}>0)$)'),TeX(r'($\hat{g}=0(\hat{\widetilde{\tau}}<0)$)')), main =TeX(r'($\delta = \infty$)'))
legend( x = "right", bty = "n",
        inset = c(-0.25, 0),
        legend=c(TeX(r'($S=1$)'), TeX(r'($S=0$)')),
        col=pal1, 
        pch=c(19, 19), cex = 1.5, xpd = TRUE)

# S1
ts1 = sum(S1==1 & tau0.005 > 0)
ts0 = sum(S1==0 & tau0.005 > 0)
uts1 = sum(S1==1 & tau0.005 < 0)
uts0 = sum(S1==0 & tau0.005 < 0)
data <- matrix(c(ts1, ts0, uts1, uts0) , nrow=2)
data_percentage <- apply(data, 2, function(x){x*100/sum(x,na.rm=T)})
barplot(data_percentage, col=pal1 , border="white", ylab="% of each subgroup", cex.names=1.5, cex.lab=1.5, cex.main=2,
        names=c(TeX(r'($\hat{g}=1(\hat{\widetilde{\tau}}>0)$)'),TeX(r'($\hat{g}=0(\hat{\widetilde{\tau}}<0)$)')), main = TeX(r'($\delta = 0$)'))
legend( x = "right", bty = "n",
        inset = c(-0.25, 0),
        legend=c(TeX(r'($S=1$)'), TeX(r'($S=0$)')),
        col=pal1, 
        pch=c(19, 19), cex = 1.5, xpd = TRUE)














