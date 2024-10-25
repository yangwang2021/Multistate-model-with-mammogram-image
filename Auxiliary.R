
### obtain two fourier eigenfuntion for generating image data
eigen_fun<-function(K){
  Twobasis <-  simMultiFunData(type = "weighted", 
                               argvals = list(list(seq(0,1,length.out=n1), seq(0,1,length.out=n2))), 
                               M = list(c(K, K)), eFunType = list(c("Fourier", "Fourier")), 
                               eValType = "exponential", N=1)
  #Trueigen <- Twobasis$trueFuns[[1]]@X[c(2,3,6,7,8,9,11,12,17,19), , ]
  Trueigen <- Twobasis$trueFuns[[1]]@X[c(2,3,5,11,13,15,21,23,24,25), , ]
  #Trueigen <- Twobasis$trueFuns[[1]]@X[c(2,3,5,6,11,13,15,21,23,25), , ]
  #Trueigen <- Twobasis$trueFuns[[1]]@X[c(3,5,11,12,13,14,15,21,22,23), , ]
  
  return(Trueigen)
}


### simulate image data for two-dimensional domains
Image_data<- function(K, nsub, mu, R){
  ### generate score from multiple normal
  tscore=MASS::mvrnorm(nsub, mu = mu, Sigma = R)
  
  ### get true image eigenfunction
  trueigen <- eigen_fun(K/2)
  
  ### get simulated image data
  Meigf <- t(apply(trueigen, 1, c))
  Gimag <- tscore %*% Meigf
  
  ### trans the generated image data to the format which can suite the MFPCA
  ### scale image data in order to the mean equal to zero
  ### scale data based on patients
  SGimag <- scale(Gimag, center=TRUE, scale=FALSE)
  
  AGimag <- array(SGimag, dim = c(nsub, n1, n2))
  Simag <- funData(argvals = list(seq(0,1,length.out=n1), seq(0,1,length.out=n2)), AGimag)
  SMimag <- as.multiFunData(Simag)
  SMimag1 <- SMimag - meanFunction(SMimag)
  ### add error for observed image data
  Obsimag <- addError(SMimag1, sd=0.2) 
  mfpca <- MFPCA(Obsimag, M = 12, uniExpansions = list(list(type = "splines2D", k=c(12,12))))
  Cpr <- round(summary(mfpca)[3,],6)
  
  #choos.ind <- which.min(abs(Cpr - 0.991))
  #choos.ind <- length(which((Cpr - 0.99) < 0))+1 ## more than 99%
  choos.ind1 <- length(which((Cpr - 0.75) < 0))+1
  choos.ind2 <- length(which((Cpr - 0.85) < 0))+1
  choos.ind3 <- length(which((Cpr - 0.95) < 0))+1
  
  # eigfun <- mfpca$functions[[1]]@X[c(1:choos.ind),,]
  # scoreV <- mfpca$scores[ ,c(1:choos.ind)]
  
  eigfun <- mfpca$functions[[1]]@X
  scoreV <- mfpca$scores
  
  return(list(eigf=eigfun, mscore=scoreV, tscore = tscore, 
              choose1=choos.ind1, choose2=choos.ind2, choose3=choos.ind3))
}


### generate multi-state data using true score 
gene_multistat <- function(Q.0, coef.value, tscore, num.score, mscore){
  
  # generate multi-state data using true score 
  funcov <- sapply(1:K, function(x){rep(tscore[,x], rep(6,nsub))})
  ### time should random from poisson process
  ### generate arrival times of a Poisson process with 
  ### lambda = rate * T = 0.5 * mean(T) = 0.5 * 5 * mean(X) = 0.5*5*(1/0.5)=5,
  ### where X is the interval time ~ exp(rate=0.5).
  
  timeobs <- rep()
  for(i in 1:nsub){
    tim <- c(0,round(cumsum(rexp(5, rate = 0.5)),3))
    while (max(tim) >= 15) {
      tim <- c(0,round(cumsum(rexp(5, rate = 0.5)),3))
    }
    timeobs <- c(timeobs,tim)
  }
  
  sim.df <- data.frame(subject = rep(1:nsub, rep(6,nsub)), 
                       time = timeobs, 
                       age = rep(rnorm(nsub, mean=0, sd=1), rep(6,nsub)), 
                       image = funcov)
  
  simudata <- simmulti.msm(sim.df, qmatrix = Q.0, start = 1, 
                           covariates = coef.value)
  # statetable.msm(state, subject, data = simudata)
  # to
  # from 1   2   3
  # 1   87 153 347
  # 2    0  76 150
  simudata.add <- simudata %>%
    group_by(subject) %>%
    mutate(as.data.frame(sapply(1:num.score, function(x){rep(mscore[subject[1],x],length(subject))}))) %>% 
    dplyr::select(-c(keep, image.1, image.2, image.3, image.4, image.5, image.6, image.7, image.8, image.9, image.10))

  return(simudata.add)
}



# obtain the estimator for inference
estimate.fun <- function(data.type, Q.1, num.score, meigf){
  
  #### get the name of covariates
  cov.name <- "age"
  for (k in 1:num.score) {
    score.name <- paste0("V",k)
    cov.name <- c(cov.name, score.name)
  }
  
  fitd = msm(state ~ time,
             covariates = as.formula(paste("~", paste(cov.name, collapse = "+"))),
             subject = subject, data = data.type, qmatrix = Q.1,
             control=list(trace=0, fnscale=5000, reltol=1e-8, maxit=5000),
             center = T)
  #print(fitd)
  
  params <- fitd$estimates
  params.cil <- c(log(fitd$ci[c(1,2),1]), fitd$ci[-c(1,2),1])
  params.cir <- c(log(fitd$ci[c(1,2),2]), fitd$ci[-c(1,2),2])
  
  params.se <- NULL
  for(i in 1:(2 + num.score)){
    se<- c(fitd$QmatricesSE[[i]][1,2],fitd$QmatricesSE[[i]][2,3])
    params.se <- c(params.se, se)
  }
  
  coef <- fitd$estimates[-c(1:4)]
  cfcs1 <- coef[seq(1,length(coef), by=2)]
  cfcs2 <- coef[seq(2,length(coef), by=2)]
  
  ecs1 <- cfcs1 %*% meigf
  ecs2 <- cfcs2 %*% meigf
  
  secs1 <- params.se[seq(5, length(params.se), by=2)] %*% meigf
  secs2 <- params.se[seq(6, length(params.se), by=2)] %*% meigf
  
  return(list(par=params, par.se=params.se, par.cil=params.cil, 
              par.cir=params.cir, ecs1=ecs1, ecs2=ecs2, secs1=secs1, secs2=secs2))
}


# obtain score value for image
score_image_fun <- function(data.type, Q.1, num.score){
  fitd0 = msm(state ~ time,
              covariates = ~ age,
              subject = subject, data = data.type, qmatrix = Q.1, 
              control=list(fnscale=5000, reltol=1e-8, maxit=5000), 
              center = T)
  ######## center=TRUE
  #fitd0$estimates.t
  coe.age <- unname(fitd0$estimates.t[c(3,4)])
  
  covini <- list(age=coe.age)
  for (k in 1:num.score) {
    score.k <- list(c(0,0))
    names(score.k) <- paste0("V",k)
    covini <- c(covini, score.k)
  }
  
  q1=fitd0$estimates.t[1]; q2=fitd0$estimates.t[2]; 
  Q.2 = rbind(c(-q1, q1, 0),
              c(0,-q2, q2),
              c(0, 0, 0))
  
  fitd1 = msm(state ~ time,
              covariates = as.formula(paste("~", paste(names(covini), collapse = "+"))),
              covinits = covini,
              subject = subject, data = data.type, qmatrix = Q.2, 
              control=list(trace=1, maxit=0), center = T)
  
  der <- fitd1$deriv/2
  cova <- fitd1$covmat
  testscore0 <- t(der) %*% cova %*% der
  testscore <- t(der[-c(1:4)]) %*% cova[-c(1:4), -c(1:4)] %*% der[-c(1:4)]
  
  pseq1 <- seq(5,length(der),by=2)
  pseq2 <- seq(6,length(der),by=2)
  dfm1 <- length(pseq1)
  dfm2 <- length(pseq2)
  testscore1 <- t(der[pseq1]) %*% cova[pseq1, pseq1] %*% der[pseq1]
  testscore2 <- t(der[pseq2]) %*% cova[pseq2, pseq2] %*% der[pseq2]
  p_test1<- pchisq(testscore1, df=dfm1, lower.tail = TRUE, log.p = FALSE)
  p_test2<- pchisq(testscore2, df=dfm2, lower.tail = TRUE, log.p = FALSE)
  
  return(as.numeric(c(testscore0, testscore, testscore1, testscore2, p_test1, p_test2, dfm1, dfm2)))
}







