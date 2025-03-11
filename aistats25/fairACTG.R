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
actg_original <- read.csv("/Users/.../ACTG_175_combined.csv")
actg_original
actg <- actg_original


###suvival time
actg_ori_C1 <- subset(actg_original,actg_original$cid==1)
actg_ori_C0 <- subset(actg_original, actg_original$cid==0)
t0s1_ori_C1  <- actg_ori_C1[actg_ori_C1$treat == 0 & actg_ori_C1$gender == 1, ]
t0s0_ori_C1  <- actg_ori_C1[actg_ori_C1$treat == 0 & actg_ori_C1$gender == 0, ]
t1s1_ori_C1  <- actg_ori_C1[actg_ori_C1$treat == 1 & actg_ori_C1$gender == 1, ]
t1s0_ori_C1  <- actg_ori_C1[actg_ori_C1$treat == 1 & actg_ori_C1$gender == 0, ]
mean(t0s1_ori_C1$time)
mean(t0s0_ori_C1$time)
mean(t1s1_ori_C1$time)
mean(t1s0_ori_C1$time)

#######select data
actg <- actg_original
table(actg$gender)
table(actg$race)
table(actg$treat)
sum(actg$treat==0 & actg$gender==1)
sum(actg$treat==0 & actg$gender==0)
sum(actg$treat==1 & actg$gender==1)
sum(actg$treat==1 & actg$gender==0)

indices_to_delete_pool <- which(actg$treat == 1 & actg$gender == 1)
print(length(indices_to_delete_pool)) 

#set.seed(123) 
indices_to_delete <- sample(indices_to_delete_pool, 1300)

actg <- actg[-indices_to_delete, ]
remaining_count <- sum(actg$treat ==1 & actg$gender == 1)
print(remaining_count) 

sum(actg$treat==1& actg$gender==1)


t0s1 <- actg[actg$treat == 0 & actg$gender == 1, ]
t0s0 <- actg[actg$treat == 0 & actg$gender == 0, ]
t1s1 <- actg[actg$treat == 1 & actg$gender == 1, ]
t1s0 <- actg[actg$treat == 1 & actg$gender == 0, ]
start_time <- Sys.time()  # 记录开始时间



############## Random survival forest imputation###################
# 创建生存对象并执行随机森林生存分析
table(actg$cid)
# 拟合随机生存森林模型
rsf_model <- rfsrc(Surv(time, cid) ~ ., data = actg, block.size = 1)
# 提取时间网格
time_grid <- rsf_model$time.interest

# 获取生存概率矩阵
surv_prob <- rsf_model$survival.oob  # 或者 rsf_model$survival

# 定义计算预期生存时间的函数
expected_survival_time <- function(surv_probs, times) {
  # 计算生存概率的差值
  surv_diff <- -diff(c(1, surv_probs))
  # 计算预期生存时间
  exp_time <- sum(surv_diff * times)
  return(exp_time)
}

# 对每个观测计算预期生存时间
expected_times <- apply(surv_prob, 1, expected_survival_time, times = time_grid)

# 创建插补后的时间向量
imputed_time <- actg$time

# 对于删失观测，替换生存时间为预期生存时间
imputed_time[actg$cid == 0] <- expected_times[actg$cid == 0]

# 将插补后的时间添加到数据集中
actg$imputed_time <- imputed_time

# 
actg <- actg[ , -which(names(actg) == "time")]
actg <- actg[ , -which(names(actg) == "trt")]

# 重命名指定的变量
names(actg)[names(actg) == 'imputed_time'] <- 'Y'
names(actg)[names(actg) == 'gender'] <- 'S1'
names(actg)[names(actg) == 'race'] <- 'S2'
names(actg)[names(actg) == 'treat'] <- 'A'

# 获取需要重命名的其他变量名称
variables_to_rename <- setdiff(names(actg), c('Y', 'S1', 'S2', 'A'))

# 将其他变量重命名为 X1, X2, X3, ..., XP
names(actg)[names(actg) %in% variables_to_rename] <- paste0('X', seq_along(variables_to_rename))

set.seed(123)  # 您可以选择任意整数
n_total <- nrow(actg)
indices <- sample(n_total)
n_train <- floor(0.75* n_total)
train_indices <- indices[1:n_train]
test_indices <- indices[(n_train + 1):n_total]
train_data <- actg[train_indices, ]
test_data <- actg[test_indices, ]

cat("训练集样本数：", nrow(train_data), "\n")
cat("测试集样本数：", nrow(test_data), "\n")
actg<-train_data




sum(actg$A==0 & actg$S1==0)
sum(actg$A==0 & actg$S1==1)
# 定义颜色
#pal2 <- c("darkgreen", "#444")
pal2 <- c("pink", "#CCCCCC")
# 设置绘图区域
par(mfrow=c(1,1))
par(mar = c(3.5, 4.5, 3, 7))

# S1
ts1 = sum(actg$S1==1 & actg$A==0 )
ts0 = sum(actg$S1==0 & actg$A==0 )
uts1 = sum(actg$S1==1 & actg$A==1 )
uts0 = sum(actg$S1==0 & actg$A==1 )
data <- matrix(c(uts1, uts0, ts1, ts0), nrow = 2)
data_percentage <- apply(data, 2, function(x){ x * 100 / sum(x, na.rm = TRUE) })

barplot(
  data_percentage,
  col = pal2,
  border = "white",
  ylab = "% of each subgroup",
  cex.names = 1.5,
  cex.lab = 1.5,
  cex.main = 2,
  names = c(
    TeX(r'(${g^*}=1($\tau$>0$))'),
    TeX(r'(${g^*}=0($\tau$<0$))')
  )
)

legend(
  x = "right",
  bty = "n",
  inset = c(-0.25, 0),
  legend = c(TeX(r'($S=1$)'), TeX(r'($S=0$)')),
  col = pal2,
  pch = c(19, 19),
  cex = 1.5,
  xpd = TRUE
)

####################################################################################################

# Estimate CATE with fairness constraints
fr.cate0.005 <-fr_cate(actg, nm.s.features = c("S1"),
                     nm.l.factors = c(NA,NA),
                     solver="quadprog",
                     nuisance.est='rf',
                     delta=0.005, sqr=TRUE, interactions=TRUE)
b.mat0.005<- fr.cate0.005$b.mat
beta.hat0.005 <- fr.cate0.005$beta.hat
tau.hat0.005 <- b.mat0.005 %*% beta.hat0.005
tau0.005 = tau.hat0.005 


fr.cate10 <-fr_cate(actg, nm.s.features = c("S1"),
                    nm.l.factors = c(NA,NA),
                    solver="quadprog",
                    nuisance.est='rf',
                    delta=70, sqr=TRUE, interactions=TRUE)
b.mat10<- fr.cate10$b.mat
beta.hat10 <- fr.cate10$beta.hat
tau.hat10 <- b.mat10 %*% beta.hat10
tau10 = tau.hat10


# 设置绘图区域
par(mfrow=c(1,1))
par(mar = c(3.5, 4.5, 3, 7))

# S1
ts1 = sum(actg$S1==1 & tau0.005  > 0)
ts0 = sum(actg$S1==0 & tau0.005 >0)
uts1 = sum(actg$S1==1 & tau0.005  < 0)
uts0 = sum(actg$S1==0 & tau0.005  < 0)
data <- matrix(c(ts1, ts0, uts1, uts0) , nrow=2)
data_percentage <- apply(data, 2, function(x){x*100/sum(x,na.rm=T)})
barplot(data_percentage, col=pal2 , border="white", ylab="% of each subgroup", cex.names=1.5, cex.lab=1.5, cex.main=2,
        names=c(TeX(r'($\hat{g}=1(\hat{\widetilde{\tau}}>0)$)'),TeX(r'($\hat{g}=0(\hat{\widetilde{\tau}}<0)$)')))
legend( x = "right", bty = "n",
        inset = c(-0.25, 0),
        legend=c(TeX(r'($S=1$)'), TeX(r'($S=0$)')),
        col=pal2, 
        pch=c(19, 19), cex = 1.5, xpd = TRUE)




###compare density
x.range=c(-350,350)
d1.s0.inf <- density(tau10[actg$S1==0])
d1.s1.inf <- density(tau10[actg$S1==1])
y.max <- max(c(d1.s0.inf$y,d1.s1.inf$y))
y.min <- min(c(d1.s0.inf$y,d1.s1.inf$y))
x.max <- max(c(d1.s0.inf$x,d1.s1.inf$x))
x.min <- min(c(d1.s0.inf$x,d1.s1.inf$x))
den <- density(tau10)
par(mfrow=c(2,1), mar=c(5, 4, 4, 2)-0.3 )
# CATE density and S1 
plot(den, xlab = TeX(r'(${\hat{\widetilde{\tau}}}$)'), ylab="", main="", cex.main=1.5, cex.lab=1.5,xlim=x.range)
tau1 <- tau10[actg$S1==1]
tau0 <- tau10[actg$S1==0]
points(tau1, rep(-0.000002,length(tau1)), col=alpha(pal2[1],0.5), pch=19)
points(tau0, rep(0.00016,length(tau0)), col=alpha(pal2[2],0.5), pch=19)
legend("topleft", legend= TeX(r'($\delta= \infty$)'))
legend("topright", legend=c(TeX(r'($S=1$)'), TeX(r'($S=0$)')),
       col=pal2, cex=1, pch=19)
abline(v=weighted.mean(d1.s0.inf$x, d1.s0.inf$y), lty=2, lwd=3, col=pal2[2], xlim=c(x.min,x.max))
abline(v=weighted.mean(d1.s1.inf$x, d1.s1.inf$y), lty=3, lwd=3, col=pal2[1], xlim=c(x.min,x.max))
#fair
d1.s0.inf <- density(tau0.005[actg$S1==0])
d1.s1.inf <- density(tau0.005[actg$S1==1])
y.max <- max(c(d1.s0.inf$y,d1.s1.inf$y))
y.min <- min(c(d1.s0.inf$y,d1.s1.inf$y))
x.max <- max(c(d1.s0.inf$x,d1.s1.inf$x))
x.min <- min(c(d1.s0.inf$x,d1.s1.inf$x))
den <- density(tau0.005)
# CATE density and S1 
plot(den, xlab = TeX(r'(${\hat{\widetilde{\tau}}}$)'), ylab="",main="",cex.main=1.5, cex.lab=1.5,xlim=x.range)
tau1 <- tau0.005[actg$S1==1]
tau0 <- tau0.005[actg$S1==0]
points(tau1, rep(-0.000002,length(tau1)), col=alpha(pal2[1],0.5), pch=19)
points(tau0, rep(0.00016,length(tau0)), col=alpha(pal2[2],0.5), pch=19)
legend("topleft", legend= TeX(r'($\delta= 0$)'))
legend("topright", legend=c(TeX(r'($S=1$)'), TeX(r'($S=0$)')),
       col=pal2, cex=1.1, pch=19)
abline(v=weighted.mean(d1.s0.inf$x, d1.s0.inf$y), lty=2, lwd=3, col=pal2[2], xlim=c(x.min,x.max))
abline(v=weighted.mean(d1.s1.inf$x, d1.s1.inf$y), lty=3, lwd=3, col=pal2[1], xlim=c(x.min,x.max))
####

















