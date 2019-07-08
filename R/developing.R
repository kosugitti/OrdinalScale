library(tidyverse)
library(rstan)
options(mc.cores = parallel::detectCores())

# Purpose -----------------------------------------------------------------

## 1,双対尺度法の順序尺度
## 2.展開法の確率モデル
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


# Unfolding model ---------------------------------------------------------



# GRM ---------------------------------------------------------------------
GRMmodel <- stan_model("GRM.stan")
# GPCM --------------------------------------------------------------------


# 多次元GRM ------------------------------------------------------------------


# 多次元GPCM -----------------------------------------------------------------





