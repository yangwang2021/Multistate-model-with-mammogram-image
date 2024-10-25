options(warn=-1)
rm(list=ls())
library(MFPCA);library(funData);library(MASS); 
library(msm); library(readxl); library(msm); library(readr); library(tidyverse); 
library(msm.stacked); library(ggplot2);library(survival); library(SDMTools);
library(RColorBrewer);library(gplots)

source("Auxiliary.R")

K=10
n1=n2=40
nsub=500
R=matrix(0, nrow = 10, ncol = 10)
diag(R) <- c(0.5, 0.4, 0.3, 0.2, 0.1, 0.08, 0.05, 0.03, 0.02, 0.01)
mu=rep(0,10)
seed=1

Scoresv.image1 = Scoresv.image2 = Scoresv.image3 = NULL;

for (i in (seed*100 + 1):(seed*100 + 500)) {
  set.seed(i)
  print(paste('iter', i))
  
  Imag <- Image_data(K, nsub, mu, R)
  
  tscorev <- Imag$tscore
  
  num.e1 <- Imag$choose1
  num.e2 <- Imag$choose2
  num.e3 <- Imag$choose3
  eigf1 <- Imag$eigf[c(1:num.e1), , ]
  eigf2 <- Imag$eigf[c(1:num.e2), , ]
  eigf3 <- Imag$eigf[c(1:num.e3), , ]
  
  meigf1 <- t(apply(eigf1, 1, c))
  meigf2 <- t(apply(eigf2, 1, c))
  meigf3 <- t(apply(eigf3, 1, c))
  
  mscorev1 <- Imag$mscore[ ,c(1:num.e1)]
  mscorev2 <- Imag$mscore[ ,c(1:num.e2)]
  mscorev3 <- Imag$mscore[ ,c(1:num.e3)]
  
 
  Q.sim <- matrix(c(-0.32, 0.32, 0,
                    0, -0.53, 0.53,
                    0, 0, 0), ncol = 3, nrow = 3, byrow = T)
  
  coef.true0 <- list(age = c(0.1, 0.15),
                     image.1=c(0,0),image.2=c(0,0),image.3=c(0,0),
                     image.4=c(0,0),image.5=c(0,0),image.6=c(0,0),
                     image.7=c(0,0),image.8=c(0,0),image.9=c(0,0),
                     image.10=c(0,0))
 
  #set.seed(300)
  sim.data1 <- gene_multistat(Q.0 = Q.sim, coef.value = coef.true0, tscore = tscorev, num.score = num.e1, mscore = mscorev1)
  sim.data2 <- gene_multistat(Q.0 = Q.sim, coef.value = coef.true0, tscore = tscorev, num.score = num.e2, mscore = mscorev2)
  sim.data3 <- gene_multistat(Q.0 = Q.sim, coef.value = coef.true0, tscore = tscorev, num.score = num.e3, mscore = mscorev3)
  
  q = 0.05
  Q.init = rbind(c(-q, q, 0),
                 c(0,-q, q),
                 c(0, 0, 0))
  
  testsvalue.image1 <- score_image_fun(data.type = sim.data1, Q.1=Q.init, num.score = num.e1)
  testsvalue.image2 <- score_image_fun(data.type = sim.data2, Q.1=Q.init, num.score = num.e2)
  testsvalue.image3 <- score_image_fun(data.type = sim.data3, Q.1=Q.init, num.score = num.e3)
  
  ## save results
  Scoresv.image1 = rbind(Scoresv.image1, testsvalue.image1)
  Scoresv.image2 = rbind(Scoresv.image2, testsvalue.image2)
  Scoresv.image3 = rbind(Scoresv.image3, testsvalue.image3)

  save(Scoresv.image1, Scoresv.image2, Scoresv.image3, file=paste0('Resultimage/','seedmqTtl',seed,'.RData'))
}

