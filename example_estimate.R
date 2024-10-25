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

Endsta1 = Endsta2 = Endsta3 = NULL
Est1 = Est2 = Est3 = NULL


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
  
  coef.true1 <- list(age = c(0.1, 0.15),
                     image.1=c(0.055, 0.05),image.2=c(0.05, 0.04),image.3=c(0.04, 0.03),
                     image.4=c(0.03, 0.025),image.5=c(0.025, 0.02),image.6=c(0.02, 0.015),
                     image.7=c(0.015, 0.01),image.8=c(0.01, 0.008),image.9=c(0.008, 0.005),
                     image.10=c(0.005, 0.003)
  )
  
  
  sim.data1 <- gene_multistat(Q.0 = Q.sim, coef.value = coef.true1, tscore = tscorev, num.score = num.e1, mscore = mscorev1)
  sim.data2 <- gene_multistat(Q.0 = Q.sim, coef.value = coef.true1, tscore = tscorev, num.score = num.e2, mscore = mscorev2)
  sim.data3 <- gene_multistat(Q.0 = Q.sim, coef.value = coef.true1, tscore = tscorev, num.score = num.e3, mscore = mscorev3)
  
  endsta1 <- endsta.fun(data=sim.data1)
  endsta2 <- endsta.fun(data=sim.data2)
  endsta3 <- endsta.fun(data=sim.data3)
  
  q = 0.05
  Q.init = rbind(c(-q, q, 0),
                 c(0,-q, q),
                 c(0, 0, 0))

  estimate1 <- estimate.fun(data.type = sim.data1, Q.1=Q.init, num.score = num.e1, meigf = meigf1)
  estimate2 <- estimate.fun(data.type = sim.data2, Q.1=Q.init, num.score = num.e2, meigf = meigf2)
  estimate3 <- estimate.fun(data.type = sim.data3, Q.1=Q.init, num.score = num.e3, meigf = meigf3)
  
  ## save results
  Endsta1 = rbind(Endsta1, endsta1)
  Endsta2 = rbind(Endsta2, endsta2)
  Endsta3 = rbind(Endsta3, endsta3)
  
  Est1<-c(Est1, estimate1)
  Est2<-c(Est2, estimate2)
  Est3<-c(Est3, estimate3)
  
  save(Endsta1, Endsta2, Endsta3, Est1, Est2, Est3,
       file=paste0('Resultimage/','seedmehTtl',seed,'.RData'))
}

