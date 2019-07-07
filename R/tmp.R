Wisky <- matrix(c(1, 3, 6, 3, 5, 3, 6, 2, 0), nrow = 3)
DualScalingMA(Wisky)
DualScaling(Wisky)



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



ord.data2 <- matrix(c(
  2, 1, 1, 2, 1, 1,
  1, 2, 2, 1, 3, 2,
  3, 4, 3, 4, 2, 3,
  4, 3, 4, 3, 4, 4
), ncol = 4)

Height.Blood <- matrix(c(2, 1, 2, 1, 4, 0, 1, 1, 3), nrow = 3)
colnames(Height.Blood) <- c("High", "Mid", "Short")
rownames(Height.Blood) <- c("Hi-pre", "Normal", "Low-pre")

pair.data <- matrix(c(
  1, 0, 2, 1, 1, 2, 1, 1,
  2, 2, 2, 2, 2, 1, 2, 2,
  1, 2, 2, 0, 2, 1, 1, 2,
  2, 2, 2, 2, 2, 1, 1, 2,
  2, 1, 2, 1, 0, 2, 1, 2,
  1, 1, 1, 1, 1, 1, 1, 2
), ncol = 6)



source("R/make_dominance.R")
tmp <- make.dominance(pair.data)
d.data <- make.dominance(ord.data)
dat.tmp <- d.data
