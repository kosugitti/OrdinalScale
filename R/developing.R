library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())

# Purpose -----------------------------------------------------------------

## 1,双対尺度法の順序尺度
## 2.多次元展開法
## 3.GRMの確率モデル
## 4.GPCMの確率モデル



# Dual Scaling ------------------------------------------------------------
Wisky <- matrix(c(1, 3, 6, 3, 5, 3, 6, 2, 0), nrow = 3)
DualScalingMA(Wisky)
DualScaling(Wisky)


# Dual Scaling for Ordinal Scale ------------------------------------------

ord.data <- matrix(c(
  3, 2, 3, 1, 2,
  6, 7, 5, 8, 7,
  8, 5, 7, 6, 4,
  1, 3, 2, 3, 5,
  4, 1, 1, 2, 1,
  5, 8, 6, 7, 8,
  2, 4, 4, 5, 3,
  7, 6, 8, 4, 6
), nrow = 5)
DualScalingOrdinal(ord.data)

source("http://aoki2.si.gunma-u.ac.jp/R/src/ro.dual.R", encoding = "euc-jp")
ro.dual(ord.data)


# Dual Scaling for Multiple choice ----------------------------------------


# Dual Scaling for paired comparison --------------------------------------


pair.data <- matrix(c(
  1, 0, 2, 1, 1, 2, 1, 1,
  2, 2, 2, 2, 2, 1, 2, 2,
  1, 2, 2, 0, 2, 1, 1, 2,
  2, 2, 2, 2, 2, 1, 1, 2,
  2, 1, 2, 1, 0, 2, 1, 2,
  1, 1, 1, 1, 1, 1, 1, 2
), ncol = 6)


# Multidimensional Unfolding model ---------------------------------------------------------

## simulation data
D <- 2
M <- 9
N <- 50
Fsc <- matrix(nrow=N,ncol=D)
Apos <- matrix(nrow=M,ncol=D)
for(n in 1:N){
  for(d in 1:D){
    Fsc[n,d] = rnorm(1,0,1);
  }
}
Apos <- matrix(c(
                  0,-1,
                  1,-1,
                 -1, 0,
                  0, 0,
                  1, 0,
                  0, 1,
                 -1, 1,
                 -1,-1,
                  1, 1),ncol=2,byrow=T)
beta0 <- 10
beta1 <- 5
sig <- 3

DataMat <- matrix(nrow=N,ncol=M)
for(n in 1:N){
  for(m in 1:M){
    tmp = 0
    for(d in 1:D){
      tmp = tmp + (Fsc[n,d]-Apos[m,d])^2
    }
    tmp = sqrt(tmp)
    DataMat[n,m] = beta0 - beta1*tmp + rnorm(1,mean=0,sd=sig)
  }
}

DataMat %>% as.data.frame %>% mutate(ID=row_number()) %>%
  tidyr::gather(key,val,-ID,factor_key=TRUE) %>%
  mutate(key=as.numeric(key)) -> dat.long

datalist <- list(
  L = NROW(dat.long), M = max(dat.long$key), N = max(dat.long$ID),
  Pid = dat.long$ID, Jid = dat.long$key, val = dat.long$val, D= 2
)

UNFLDconst <- stan_model("R/UNFLD_const.stan")
UNFLD <- stan_model("R/UNFLD.stan")
fitL <- sampling(UNFLD,datalist)
fitLc <- sampling(UNFLDconst,datalist)
print(fitL,pars=c('beta0','beta1','sig'))
print(fitLc,pars=c('beta0','beta1','sig'))

traceplot(fitL,pars='a')
traceplot(fitLc,pars='a')

print(fitL,pars='a')
print(fitLc,pars='a')
### 多次元展開法は座標をuniformにすべきか？
### MDSのように象限を固定してやるべきか？？

# GRM ---------------------------------------------------------------------
library(ltm)
GRMmodel <- stan_model("R/GRM.stan")



## Sample Data
dat <- Science[c(1, 3, 4, 7)]
result.grm <- grm(dat,IRT.param = TRUE)
dat %>%
  mutate(ID = row_number()) %>%
  tidyr::gather(key, val, -ID) %>%
  mutate(key = ifelse(key=="Comfort",1,
                      ifelse(key=="Work",2,
                             ifelse(key=="Future",3,4)))) %>%
  mutate(val = ifelse(val=="strongly disagree",1,
                      ifelse(val=="disagree",2,
                             ifelse(val=="agree",3,4))))  %>% arrange(ID)-> dat.long
datalist <- list(
  L = NROW(dat.long), C = max(dat.long$val), M = max(dat.long$key), N = max(dat.long$ID),
  Pid = dat.long$ID, Jid = dat.long$key, val = dat.long$val, D= 1
)

fit <- vb(GRMmodel,datalist)
result.grm
print(fit,pars=c("a","b","loc"))
fit %>% rstan::extract() %>% as.data.frame() %>%
  dplyr::select(starts_with("theta")) %>% tidyr::gather(key,val,factor_key=TRUE) %>%
  group_by(key) %>% summarise(EAP=mean(val)) %>%
  mutate(SC=factor.scores.grm(result.grm,resp.patterns = dat)$score.dat$z1) %>%
  ggplot(aes(x=EAP,y=SC))+geom_point()

## Sample Data2
datYG <- read_table("R/YG外向.txt") %>% rename(val=ABCDEFGHIJ) %>%
  mutate(A = str_sub(.$val,start=1, end=1) %>% as.numeric(),
         B = str_sub(.$val,start=2, end=2)%>% as.numeric(),
         C = str_sub(.$val,start=3, end=3)%>% as.numeric(),
         D = str_sub(.$val,start=4, end=4)%>% as.numeric(),
         E = str_sub(.$val,start=5, end=5)%>% as.numeric(),
         F = str_sub(.$val,start=6, end=6)%>% as.numeric(),
         G = str_sub(.$val,start=7, end=7)%>% as.numeric(),
         H = str_sub(.$val,start=8, end=8)%>% as.numeric(),
         I = str_sub(.$val,start=9, end=9)%>% as.numeric(),
         J = str_sub(.$val,start=10, end=10)%>% as.numeric()) %>%
  dplyr::select(-"val") %>%
  mutate(ID=row_number()) %>%
  tidyr::gather(key,val,-ID) %>%
  mutate(key=as.numeric(as.factor(key))) %>%
  mutate(val=val+1)

datalist <- list(
  L = NROW(datYG), C = max(datYG$val), M = max(datYG$key), N = max(datYG$ID),
  Pid = datYG$ID, Jid = datYG$key, val = datYG$val
)
fit <- sampling(GRMmodel,datalist)
print(fit,pars=c("a","loc"))

# GPCM --------------------------------------------------------------------
GPCM <- stan_model("R/GPCM.stan")
fit <- sampling(GPCM,datalist)
print(fit,pars=c("a","b","loc"))


# 多次元GRM ------------------------------------------------------------------
psych::bfi %>% dplyr::select(starts_with("C"),starts_with("E"),-education) %>%
  mutate(ID=row_number()) %>%
  # 逆転項目の処理
  mutate(C4 = 7-C4, C5 = 7-C5, E1 = 7-E1, E2=7-E2) %>%
  # 具合の良さそうな項目だけ残す
  dplyr::select(C1,C2,C3,E1,E2,ID) %>%
  # 数を減らす
  dplyr::slice(1:100) %>%
  tidyr::gather(key,val,-ID,factor_key=TRUE) %>% na.omit %>%
  mutate(key=as.numeric(key)) -> bfi.dat


#read_csv("R/CE.csv") %>% mutate(ID = row_number()) %>%
#  tidyr::gather(key,val,-ID,factor_key=TRUE) %>%
#  mutate(key = as.numeric(key),val=val+1) -> bfi.dat

datalist <- list(
  L = NROW(bfi.dat), C = max(bfi.dat$val), M = max(bfi.dat$key), N = max(bfi.dat$ID),
  Pid = bfi.dat$ID, Jid = bfi.dat$key, val = bfi.dat$val,D=2
)

init_list <- list(a=matrix(c(1.5,1.7,1.3,0.1,0.3,0.3,0.3,0.1,1.5,2.3),ncol=2))

GRMD <- stan_model("R/GRMD.stan")
fit_D <- sampling(GRMD,datalist,init=list(init_list,init_list,init_list,init_list))
print(fit_D,pars=c('a','b'))
# 綺麗なラベルスイッチング
traceplot(fit_D,pars=c("a[1,1]","a[1,2]"))
traceplot(fit_D,pars=c("a[2,1]","a[2,2]"))
traceplot(fit_D,pars=c("a[3,1]","a[3,2]"))
traceplot(fit_D,pars=c("a[4,1]","a[4,2]"))
traceplot(fit_D,pars=c("a[5,1]","a[5,2]"))

# 多次元GPCM -----------------------------------------------------------------
GPCMD <- stan_model("R/GPCMD.stan")
fit_pd <- sampling(GPCMD,datalist)
print(fit_pd,pars=c('a','b'))
