rm(list = ls())
library(MASS)
library(mnormt)
library(nnet)
library(LaplacesDemon)
library(Matching)
library(dplyr)
library(cluster)
library(e1071)
library(ramsvm)
library(caret)
library(survival)
library(parallel)
library(snowfall)
library(randomForestSRC)
library(latex2exp)
#install.packages('matrixcalc')
library(matrixcalc)
#install.packages('mbend')
library(mbend)
library(Matrix)
#install.packages('Rmosek')
library(Rmosek)
library(ggplot2)
#install.packages('ranger')
library(ranger)
library(quadprog)
library(ranger)
#install.packages('Matrix')
library(Matrix)


source("/Users/.../utils.R")
source("/Users.../ITR-covgpsmatch-funcs.r")
hcc_original <- read.csv("/Users/../hcc-simplified.csv")
set.seed(8884)
hcc <- hcc_original

################### Basic parameters #############################
K <- 3  #treatment 
fold <- 5  #5-fold cross validation
calipernum <- NULL
lambda_param <- c(1e-06, 1e-05, 1e-04, 0.001, 0.01, 0.1, 1, 5, 10, 20, 50, 100, 200) #penalty tuning parameter
kernel_param <- c(1)  #Gaussian kernel parameter
############### Preprocess the raw hcc data #######################
hcc <- preprocessing(hcc) #deal with survival time and censoring indicator
#categorical covariates: dummy variable
dataX_f <- data.frame(Gender = factor(hcc$Gender), hypertension = factor(hcc$hypertension), diabetes = factor(hcc$diabetes.type2), Group = factor(hcc$Group),
                      MELD = factor(hcc$MELD.grading), Child = factor(hcc$Child.Pugh.grading), AFP = factor(hcc$AFP3))
dmy <- dummyVars(~., data = dataX_f, fullRank = T) 
dataX_f <- data.frame(predict(dmy, newdata = dataX_f))
#continuous covariates: standardized
dataX_c <- data.frame(Age = hcc$Age, BMI = hcc$BMI, ALBI = hcc$ALBI.score, APRI = log(hcc$APRI.score + 10^-5))
dataX_c <- data.frame(scale(dataX_c, center = T, scale = T))
#cbind all the covariates
dataX <- cbind(dataX_f, dataX_c)
p <- dim(dataX)[2]  #dimension

#right censored data, here we focus on overall survival time
tao <- 2000  #set the observation window
hcc$day_truncate <- hcc$hss.day
hcc$day_truncate[which(hcc$hss.day > tao)] <- tao  #truncation
hcc$censor_truncate <- hcc$hss.censor 
hcc$censor_truncate[which(hcc$hss.day > tao)] <- 0 #re-code the censoring indicator, set censor=0 for survival time > tao
dataset_truncate <- data.frame(dataX, treat = hcc$operation, censor = hcc$censor_truncate, truncateR = hcc$day_truncate)

############### Random survival forest imputation####################
rsf_fit <- rfsrc(formula(Surv(truncateR, censor) ~ .), data = dataset_truncate, block.size = 1) #fit the RSF
Y.grid <- rsf_fit$time.interest #event time
R1 <- expected_survival(rsf_fit$survival.oob, Y.grid) #expected survival time
hss_impute <- dataset_truncate
imputeR1 <- dataset_truncate$truncateR
imputeR1[which(hss_impute$censor == 0)] <- R1[which(hss_impute$censor == 0)] #replace the censored observations with imputation
hss_impute$impute <- imputeR1 #the final imputation
hss_impute <- hss_impute[ , -which(names(hss_impute) == "truncateR")]


#write.csv(hss_impute,"/Users/hongniwang/dissertation/fair-HTE/hss_impute.csv")

table(hss_impute$treat)
table(hss_impute$Gender.1)

set.seed(2024)
data_no3 <- hss_impute[hss_impute$treat != 3, ]
table(data_no3$treat)
table(data_no3$Gender.1)

sum(data_no3$treat==1 & data_no3$Gender.1==1)
sum(data_no3$treat==1 & data_no3$Gender.1==0)
sum(data_no3$treat==2 & data_no3$Gender.1==1)
sum(data_no3$treat==2 & data_no3$Gender.1==0)

indices_to_delete_pool <- which(data_no3$treat == 2 & data_no3$Gender.1 == 1)
print(length(indices_to_delete_pool)) 

set.seed(123) 
indices_to_delete <- sample(indices_to_delete_pool, 145)

data_no3 <- data_no3[-indices_to_delete, ]
remaining_count <- sum(data_no3$treat ==2 & data_no3$Gender.1 == 1)
print(remaining_count) 

sum(data_no3$treat==2& data_no3$Gender.1==1)
#######
indices_to_delete_pool2 <- which(data_no3$treat == 1 & data_no3$Gender.1 == 1)
#print(length(indices_to_delete_pool2)) 
set.seed(123) 
indices_to_delete2 <- sample(indices_to_delete_pool2, 90)
data_no3 <- data_no3[-indices_to_delete2, ]
remaining_count <- sum(data_no3$treat ==1 & data_no3$Gender.1 == 1)
sum(data_no3$treat==1& data_no3$Gender.1==1)
######
data_no3new <- data_no3 %>%
  rename(
    S1 = Gender.1,
    S2 = Age,
    Y = impute,
    A = treat
  )

# 获取原有列名
original_colnames <- colnames(data_no3new)

# 除了指定的列之外，其他的都改为 X1, X2, ...
new_colnames <- original_colnames
other_cols <- setdiff(original_colnames, c("S1", "S2", "Y", "A"))
new_colnames[which(original_colnames %in% other_cols)] <- paste0("X", seq_along(other_cols))
colnames(data_no3new) <- new_colnames

pal2 <- c("turquoise2", "yellow")

##动作，之前1还是1对应mu1；之前2赋值为0对应mu0.这样就直接对应mu1-mu0为最优g*
data_no3new <- data_no3new %>%
  mutate(A = ifelse(A == 1, 1, ifelse(A ==2, 0, A)))



####################################################################################################


# Estimate CATE without fairness constraints
##fair
fr.cate0.005 <-fr_cate(data_no3new, nm.s.features = c("S1"),
                     nm.l.factors = c(NA,NA),
                     solver="quadprog",
                     nuisance.est='rf',
                     delta=0.00005, sqr=TRUE, interactions=TRUE)
b.mat0.005<- fr.cate0.005$b.mat
beta.hat0.005 <- fr.cate0.005$beta.hat
tau.hat0.005 <- b.mat0.005 %*% beta.hat0.005
tau0.005 = tau.hat0.005 



##unfair
fr.cate10 <-fr_cate(data_no3new, nm.s.features = c("S1"),
                    nm.l.factors = c(NA,NA),
                    solver="quadprog",
                    nuisance.est='rf',
                    delta=90, sqr=TRUE, interactions=TRUE)
b.mat10<- fr.cate10$b.mat
beta.hat10 <- fr.cate10$beta.hat
tau.hat10 <- b.mat10 %*% beta.hat10
tau10 = tau.hat10 



###画出对比图
## CATE density and S1 
#pal1 <- c("#FF9933", "#9966FF")
#pal2 <- c("#99CC99", "#CCCCCC")
#pal1 <- c("pink", "darkblue")
#pal1 <- c("tomato4", "cyan3")
#pal2 <- c("orange","cyan3" )
x.range=c(-800,800)
d1.s0.inf <- density(tau10[data_no3new$S1==0])
d1.s1.inf <- density(tau10[data_no3new$S1==1])
#d2.s0.inf <- density(tau10[S2==0])
#d2.s1.inf <- density(tau10[S2==1])
y.max <- max(c(d1.s0.inf$y,d1.s1.inf$y))
y.min <- min(c(d1.s0.inf$y,d1.s1.inf$y))
x.max <- max(c(d1.s0.inf$x,d1.s1.inf$x))
x.min <- min(c(d1.s0.inf$x,d1.s1.inf$x))
den <- density(tau10)
par(mfrow=c(2,1), mar=c(5, 4, 4, 2) + 0.1)
# CATE density and S1 
plot(den, xlab = TeX(r'(${\hat{\widetilde{\tau}}}$)'), ylab="", main='', cex.main=1.5, cex.lab=1.5,xlim=x.range)
tau1 <- tau10[data_no3new$S1==1]
tau0 <- tau10[data_no3new$S1==0]
points(tau1, rep(-0.00001,length(tau1)), col=alpha(pal2[1],0.5), pch=19)
points(tau0, rep(0.000068,length(tau0)), col=alpha(pal2[2],0.5), pch=19)
legend("topleft", legend= TeX(r'($\delta= \infty$)'))
legend("topright", legend=c(TeX(r'($S=1$)'), TeX(r'($S=0$)')),
       col=pal2, cex=1.2, pch=19)
abline(v=weighted.mean(d1.s0.inf$x, d1.s0.inf$y), lty=2, lwd=3, col=pal2[2], xlim=c(x.min,x.max))
abline(v=weighted.mean(d1.s1.inf$x, d1.s1.inf$y), lty=3, lwd=3, col=pal2[1], xlim=c(x.min,x.max))
#fair
d1.s0.inf <- density(tau0.005[data_no3new$S1==0])
d1.s1.inf <- density(tau0.005[data_no3new$S1==1])
#d2.s0.inf <- density(tau0.005[data_no3new$S2==0])
#d2.s1.inf <- density(tau0.005[data_no3new$S2==1])
y.max <- max(c(d1.s0.inf$y,d1.s1.inf$y))
y.min <- min(c(d1.s0.inf$y,d1.s1.inf$y))
x.max <- max(c(d1.s0.inf$x,d1.s1.inf$x))
x.min <- min(c(d1.s0.inf$x,d1.s1.inf$x))
den <- density(tau0.005)
# CATE density and S1 
plot(den, xlab = TeX(r'(${\hat{\widetilde{\tau}}}$)'), ylab="", main='', cex.main=1.5, cex.lab=1.5,xlim=c(-1000,1000))
tau1 <- tau0.005[data_no3new$S1==1]
tau0 <- tau0.005[data_no3new$S1==0]
points(tau1, rep(-0.00001,length(tau1)), col=alpha(pal2[1],0.5), pch=19)
points(tau0, rep(0.000068,length(tau0)), col=alpha(pal2[2],0.5), pch=19)
legend("topleft", legend= TeX(r'($\delta= 0$)'))
legend("topright", legend=c(TeX(r'($S=1$)'), TeX(r'($S=0$)')),
       col=pal2, cex=1.2, pch=19)
abline(v=weighted.mean(d1.s0.inf$x, d1.s0.inf$y), lty=2, lwd=3, col=pal2[2], xlim=c(x.min,x.max))
abline(v=weighted.mean(d1.s1.inf$x, d1.s1.inf$y), lty=3, lwd=3, col=pal2[1], xlim=c(x.min,x.max))
####


# bar plot: OTR vs sensitive features
par(mfrow=c(2,1))
par(mar = c(3.5, 4.5, 3, 7))
# S1
ts1 = sum(data_no3new$S1==1 & tau10 > 0)
ts0 = sum(data_no3new$S1==0 & tau10 > 0)
uts1 = sum(data_no3new$S1==1 & tau10 < 0)
uts0 = sum(data_no3new$S1==0 & tau10 < 0)
data <- matrix(c(ts1, ts0, uts1, uts0) , nrow=2)
data_percentage <- apply(data, 2, function(x){x*100/sum(x,na.rm=T)})
barplot(data_percentage, col=pal2 , border="white", ylab="% of each subgroup", cex.names=1.5, cex.lab=1.5, cex.main=2,
        names=c(TeX(r'($\hat{g}=1(\hat{\widetilde{\tau}}>0)$)'),TeX(r'($\hat{g}=0(\hat{\widetilde{\tau}}<0)$)')), main = TeX(r'($\delta = \infty$)'))
legend( x = "right", bty = "n",
        inset = c(-0.25, 0),
        legend=c(TeX(r'($S=1$)'), TeX(r'($S=0$)')),
        col=pal2, 
        pch=c(19, 19), cex = 1.5, xpd = TRUE)

# S1
ts1 = sum(data_no3new$S1==1 & tau0.005 > 0)
ts0 = sum(data_no3new$S1==0 & tau0.005 > 0)
uts1 = sum(data_no3new$S1==1 & tau0.005 < 0)
uts0 = sum(data_no3new$S1==0 & tau0.005 < 0)
data <- matrix(c(ts1, ts0, uts1, uts0) , nrow=2)
data_percentage <- apply(data, 2, function(x){x*100/sum(x,na.rm=T)})
barplot(data_percentage, col=pal2 , border="white", ylab="% of each subgroup", cex.names=1.5, cex.lab=1.5, cex.main=2,
        names=c(TeX(r'($\hat{g}=1(\hat{\widetilde{\tau}}>0)$)'),TeX(r'($\hat{g}=0(\hat{\widetilde{\tau}}<0)$)')), main = TeX(r'($\delta = 0$)'))
legend( x = "right", bty = "n",
        inset = c(-0.25, 0),
        legend=c(TeX(r'($S=1$)'), TeX(r'($S=0$)')),
        col=pal2, 
        pch=c(19, 19), cex = 1.5, xpd = TRUE)

