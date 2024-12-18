# Sim no ties

source("Final_Consensus_SEARCH_updated.r")

#Load the computed results
# load("Final_results.rdata")

#...or compute them by yourself!


iterazioni <- 100

Final_results <- matrix(NA, nrow= 9, ncol= 11)
colnames(Final_results) <- c("Eltime_us","Eltime_naive","%us","%naive","AVGdist_us","AVGdist_naive","MEDIANdist_us","MEDIANdist_naive","D_us<D_naive","D_us=D_naive","D_us>_Dnaive")
rownames(Final_results)<- c("n=10,m=10",
                            "n=20,m=10",
                            "n=30,m=10",
                            "n=20,m=20",
                            "n=40,m=20",
                            "n=60,m=20",
                            "n=50,m=50",
                            "n=100,m=50",
                            "n=150,m=50")

###############################################################################################################
#######################################à



# 10 items
# 10 giudici

nitems <- 10
njudges <- 10


Sigma=matrix(0,nrow=nitems,ncol=nitems)
diag(Sigma)=1

results_100_100 <- matrix(NA, nrow = iterazioni, ncol=7)
TF_100_100 <- rep(NA,iterazioni)


for (i in 1:iterazioni) {
  
  set.seed(i)
  mu <- rnorm(nitems, 0, 1)
  true_conf <- MASS::mvrnorm(n=njudges,mu,Sigma=Sigma)
  
  rank <- t(apply(true_conf,1,rank,ties.method="min"))
  approv <-round(runif(njudges,min=0,max=nitems),0)
  
  approv_m <- matrix(nrow=nrow(rank),ncol=ncol(rank))
  
  for (j in 1:nrow(rank)) {
    approv_m[j,] <-  ifelse(rank[j,]<=approv[j],1,0)
  }
  
  
  X = data.frame(cbind(rank,approv_m))
  
  
  
  colmeans<-apply(X,2,mean)
  colmeans
  
  naive_r <- rank(colmeans[1:(length(colmeans)/2)])
  naive_a <- ifelse(colmeans[(length(colmeans)/2+1):length(colmeans)]>0.5,1,0)
  
  
  naive_consensus <- matrix(c(naive_r,naive_a),nrow=1)
  
  #è un preference-approval valido?
  TF_100_100[i] <- controlla(upperTriangle(score_r(naive_r)), upperTriangle(score_a(naive_a)))
  
  FC <- Final_Consensus(X,algorithm = "quick")
  
  results_100_100[i,1] <- FC$D_lambda
  results_100_100[i,2]  <- FC$Elapsed
  results_100_100[i,3]  <- mean(Pref_dist2(as.matrix(X),as.matrix(naive_consensus)))

  
  print(i)
}

summary(results_100_100[,1])
summary(results_100_100[,2])
summary(results_100_100[,3])
summary(results_100_100[,4])
summary(results_100_100[,5])
sum(TF_100_100)


Final_results[1,] <- c(mean(results_100_100[,2]),0,0,sum(TF_100_100),mean(results_100_100[which(TF_100_100),1]),mean(results_100_100[which(TF_100_100),3]),
                        median(results_100_100[which(TF_100_100),1]),median(results_100_100[which(TF_100_100),3]),
                        100- round(sum(results_100_100[which(TF_100_100),1] == results_100_100[which(TF_100_100),3])),
                        round(sum(results_100_100[which(TF_100_100),1] == results_100_100[which(TF_100_100),3])),
                        0)

###############################################################################################################
#######################################à

# 10 items
# 20 giudici

nitems <- 10
njudges <- 20


Sigma=matrix(0,nrow=nitems,ncol=nitems)
diag(Sigma)=1

results_100_200 <- matrix(NA, nrow = iterazioni, ncol=7)
TF_100_200 <- rep(NA,iterazioni)


for (i in 1:iterazioni) {
  
  set.seed(1000+i)
  mu <- rnorm(nitems, 0, 1)
  true_conf <- MASS::mvrnorm(n=njudges,mu,Sigma=Sigma)
  
  rank <- t(apply(true_conf,1,rank,ties.method="min"))
  approv <-round(runif(njudges,min=0,max=nitems),0)
  
  approv_m <- matrix(nrow=nrow(rank),ncol=ncol(rank))
  
  for (j in 1:nrow(rank)) {
    approv_m[j,] <-  ifelse(rank[j,]<=approv[j],1,0)
  }
  
  
  X = data.frame(cbind(rank,approv_m))
  
  
  
  colmeans<-apply(X,2,mean)
  colmeans
  
  naive_r <- rank(colmeans[1:(length(colmeans)/2)])
  naive_a <- ifelse(colmeans[(length(colmeans)/2+1):length(colmeans)]>0.5,1,0)
  
  
  naive_consensus <- matrix(c(naive_r,naive_a),nrow=1)
  
  #è un preference-approval valido?
  TF_100_200[i] <- controlla(upperTriangle(score_r(naive_r)), upperTriangle(score_a(naive_a)))
  
  FC <- Final_Consensus(X,algorithm = "quick")

  results_100_200[i,1] <- FC$D_lambda
  results_100_200[i,2]  <- FC$Elapsed
  results_100_200[i,3]  <- mean(Pref_dist2(as.matrix(X),as.matrix(naive_consensus)))
  
  
  
  print(i)
}

summary(results_100_200[,1])
summary(results_100_200[,2])
summary(results_100_200[,3])
summary(results_100_200[,4])
summary(results_100_200[,5])
sum(TF_100_200)


Final_results[2,] <- c(mean(results_100_200[,2]),0,0,sum(TF_100_200),mean(results_100_200[which(TF_100_200),1]),mean(results_100_200[which(TF_100_200),3]),
                        median(results_100_200[which(TF_100_200),1]),median(results_100_200[which(TF_100_200),3]),
                        100- round(sum(results_100_200[which(TF_100_200),1] == results_100_200[which(TF_100_200),3])),
                        round(sum(results_100_200[which(TF_100_200),1] == results_100_200[which(TF_100_200),3])),
                        0)



###############################################################################################################
#######################################à

# 10 items
# 30 giudici

nitems <- 10
njudges <- 30


Sigma=matrix(0,nrow=nitems,ncol=nitems)
diag(Sigma)=1

results_100_300 <- matrix(NA, nrow = iterazioni, ncol=7)
TF_100_300 <- rep(NA,iterazioni)


for (i in 1:iterazioni) {
  
  set.seed(i)
  mu <- rnorm(nitems, 0, 1)
  true_conf <- MASS::mvrnorm(n=njudges,mu,Sigma=Sigma)
  
  rank <- t(apply(true_conf,1,rank,ties.method="min"))
  approv <-round(runif(njudges,min=0,max=nitems),0)
  
  approv_m <- matrix(nrow=nrow(rank),ncol=ncol(rank))
  
  for (j in 1:nrow(rank)) {
    approv_m[j,] <-  ifelse(rank[j,]<=approv[j],1,0)
  }
  
  
  X = data.frame(cbind(rank,approv_m))
  
  
  
  colmeans<-apply(X,2,mean)
  colmeans
  
  naive_r <- rank(colmeans[1:(length(colmeans)/2)])
  naive_a <- ifelse(colmeans[(length(colmeans)/2+1):length(colmeans)]>0.5,1,0)
  
  
  naive_consensus <- matrix(c(naive_r,naive_a),nrow=1)
  
  #è un preference-approval valido?
  TF_100_300[i] <- controlla(upperTriangle(score_r(naive_r)), upperTriangle(score_a(naive_a)))
  
  FC <- Final_Consensus(X,algorithm = "quick")
  
  results_100_300[i,1] <- FC$D_lambda
  results_100_300[i,2]  <- FC$Elapsed
  results_100_300[i,3]  <- mean(Pref_dist2(as.matrix(X),as.matrix(naive_consensus)))

  
  
  print(i)
}

summary(results_100_300[,1])
summary(results_100_300[,2])
summary(results_100_300[,3])
summary(results_100_300[,4])
summary(results_100_300[,5])
sum(TF_100_300)


Final_results[3,] <- c(mean(results_100_300[,2]),0,0,sum(TF_100_300),mean(results_100_300[which(TF_100_300),1]),mean(results_100_300[which(TF_100_300),3]),
                        median(results_100_300[which(TF_100_300),1]),median(results_100_300[which(TF_100_300),3]),
                        100- round(sum(results_100_300[which(TF_100_300),1] == results_100_300[which(TF_100_300),3])),
                        round(sum(results_100_300[which(TF_100_300),1] == results_100_300[which(TF_100_300),3])),
                        0)



#######################################à

# 20 items
# 20 giudici

nitems <- 20
njudges <- 20


Sigma=matrix(0,nrow=nitems,ncol=nitems)
diag(Sigma)=1

results_20_20 <- matrix(NA, nrow = iterazioni, ncol=7)
TF_20_20 <- rep(NA,iterazioni)

for (i in 1:iterazioni) {
  
  set.seed(i)
  mu <- rnorm(nitems, 0, 1)
  true_conf <- MASS::mvrnorm(n=njudges,mu,Sigma=Sigma)
  
  rank <- t(apply(true_conf,1,rank,ties.method="min"))
  approv <-round(runif(njudges,min=0,max=nitems),0)
  approv_m <- matrix(nrow=nrow(rank),ncol=ncol(rank))
  
  for (j in 1:nrow(rank)) {
    approv_m[j,] <-  ifelse(rank[j,]<=approv[j],1,0)
  }
  
  
  X = data.frame(cbind(rank,approv_m))
  
  
  
  colmeans<-apply(X,2,mean)
  colmeans
  
  naive_r <- rank(colmeans[1:(length(colmeans)/2)])
  naive_a <- ifelse(colmeans[(length(colmeans)/2+1):length(colmeans)]>0.5,1,0)
  
  
  naive_consensus <- matrix(c(naive_r,naive_a),nrow=1)
  
  #è un preference-approval valido?
  TF_20_20[i] <- controlla(upperTriangle(score_r(naive_r)), upperTriangle(score_a(naive_a)))
  
  FC <-   Final_Consensus(X,algorithm = "quick")
  
  rank_mu <- rank(mu,ties.method="min")
  approv_mu <-  ifelse(rank_mu<=round(sqrt(nitems),0),1,0)
  dgp <- data.frame(matrix(c(rank_mu,approv_mu),nrow=1))
  
  results_20_20[i,1] <- FC$D_lambda
  results_20_20[i,2]  <- FC$Elapsed
  results_20_20[i,3]  <- mean(Pref_dist2(as.matrix(X),as.matrix(naive_consensus)))
  results_20_20[i,4] <-  min(Pref_dist2(dgp,FC$Consensus))
  results_20_20[i,5]  <- Pref_dist2(dgp,naive_consensus)
  
  
  print(i)
}


Final_results[4,] <- c(mean(results_20_20[,2]),0,0,sum(TF_20_20),mean(results_20_20[which(TF_20_20),1]),mean(results_20_20[which(TF_20_20),3]),
                       median(results_20_20[which(TF_20_20),1]),median(results_20_20[which(TF_20_20),3]),
                       100- round(sum(results_20_20[which(TF_20_20),1] == results_20_20[which(TF_20_20),3])),
                       round(sum(results_20_20[which(TF_20_20),1] == results_20_20[which(TF_20_20),3])),
                       0)
###############################################################################################################
#######################################à

# 20 items
# 40 giudici

nitems <- 20
njudges <- 40


Sigma=matrix(0,nrow=nitems,ncol=nitems)
diag(Sigma)=1

results_20_40 <- matrix(NA, nrow = iterazioni, ncol=7)
TF_20_40 <- rep(NA,iterazioni)

for (i in 1:iterazioni) {
  
  set.seed(i)
  mu <- rnorm(nitems, 0, 1)
  true_conf <- MASS::mvrnorm(n=njudges,mu,Sigma=Sigma)
  
  rank <- t(apply(true_conf,1,rank,ties.method="min"))
  approv <-round(runif(njudges,min=0,max=nitems),0)
  approv_m <- matrix(nrow=nrow(rank),ncol=ncol(rank))
  
  for (j in 1:nrow(rank)) {
    approv_m[j,] <-  ifelse(rank[j,]<=approv[j],1,0)
  }
  
  rank[1,]
  approv_m[1,]
  
  X = data.frame(cbind(rank,approv_m))
  
  
  
  colmeans<-apply(X,2,mean)
  colmeans
  
  naive_r <- rank(colmeans[1:(length(colmeans)/2)])
  naive_a <- ifelse(colmeans[(length(colmeans)/2+1):length(colmeans)]>0.5,1,0)
  
  
  naive_consensus <- matrix(c(naive_r,naive_a),nrow=1)
  
  #è un preference-approval valido?
  TF_20_40[i] <- controlla(upperTriangle(score_r(naive_r)), upperTriangle(score_a(naive_a)))
  
  FC <- Final_Consensus(X,algorithm = "quick")

  
  rank_mu <- rank(mu,ties.method="min")
  approv_mu <-  ifelse(rank_mu<=round(sqrt(nitems),0),1,0)
  dgp <- data.frame(matrix(c(rank_mu,approv_mu),nrow=1))
  
  results_20_40[i,1] <- FC$D_lambda
  results_20_40[i,2]  <- FC$Elapsed
  results_20_40[i,3]  <- mean(Pref_dist2(as.matrix(X),as.matrix(naive_consensus)))
  results_20_40[i,4] <-  min(Pref_dist2(dgp,FC$Consensus))
  results_20_40[i,5]  <- Pref_dist2(dgp,naive_consensus)
  
  
  print(i)
}

summary(results_20_40[,1])
summary(results_20_40[,2])
summary(results_20_40[,3])
summary(results_20_40[,4])
summary(results_20_40[,5])
sum(TF_20_40)


Final_results[5,] <- c(mean(results_20_40[,2]),0,0,sum(TF_20_40),mean(results_20_40[which(TF_20_40),1]),mean(results_20_40[which(TF_20_40),3]),
                       median(results_20_40[which(TF_20_40),1]),median(results_20_40[which(TF_20_40),3]),
                       100- round(sum(results_20_40[which(TF_20_40),1] == results_20_40[which(TF_20_40),3])),
                       round(sum(results_20_40[which(TF_20_40),1] == results_20_40[which(TF_20_40),3])),
                       0)
###############################################################################################################
#######################################à

# 20 items
# 60 giudici

nitems <- 20
njudges <- 60


Sigma=matrix(0,nrow=nitems,ncol=nitems)
diag(Sigma)=1

results_20_60 <- matrix(NA, nrow = iterazioni, ncol=7)
TF_20_60 <- rep(NA,iterazioni)

for (i in 1:iterazioni) {
  
  set.seed(i)
  mu <- rnorm(nitems, 0, 1)
  true_conf <- MASS::mvrnorm(n=njudges,mu,Sigma=Sigma)
  
  rank <- t(apply(true_conf,1,rank,ties.method="min"))
  approv <-round(runif(njudges,min=0,max=nitems),0)
  approv_m <- matrix(nrow=nrow(rank),ncol=ncol(rank))
  
  for (j in 1:nrow(rank)) {
    approv_m[j,] <-  ifelse(rank[j,]<=approv[j],1,0)
  }
  
  rank[1,]
  approv_m[1,]
  
  X = data.frame(cbind(rank,approv_m))
  
  
  
  colmeans<-apply(X,2,mean)
  colmeans
  
  naive_r <- rank(colmeans[1:(length(colmeans)/2)])
  naive_a <- ifelse(colmeans[(length(colmeans)/2+1):length(colmeans)]>0.5,1,0)
  
  
  naive_consensus <- matrix(c(naive_r,naive_a),nrow=1)
  
  #è un preference-approval valido?
  TF_20_60[i] <- controlla(upperTriangle(score_r(naive_r)), upperTriangle(score_a(naive_a)))
  
  FC <- Final_Consensus(X,algorithm = "quick")

  
  rank_mu <- rank(mu,ties.method="min")
  approv_mu <-  ifelse(rank_mu<=round(sqrt(nitems),0),1,0)
  dgp <- data.frame(matrix(c(rank_mu,approv_mu),nrow=1))
  
  results_20_60[i,1] <- FC$D_lambda
  results_20_60[i,2]  <- FC$Elapsed
  results_20_60[i,3]  <- mean(Pref_dist2(as.matrix(X),as.matrix(naive_consensus)))
  results_20_60[i,4] <-  min(Pref_dist2(dgp,FC$Consensus))
  results_20_60[i,5]  <- Pref_dist2(dgp,naive_consensus)
  

  
  print(i)
}

summary(results_20_60[,1])
summary(results_20_60[,2])
summary(results_20_60[,3])
summary(results_20_60[,4])
summary(results_20_60[,5])
sum(TF_20_60)


Final_results[6,] <- c(mean(results_20_60[,2]),0,0,sum(TF_20_60),mean(results_20_60[which(TF_20_60),1]),mean(results_20_60[which(TF_20_60),3]),
                       median(results_20_60[which(TF_20_60),1]),median(results_20_60[which(TF_20_60),3]),
                       100- round(sum(results_20_60[which(TF_20_60),1] == results_20_60[which(TF_20_60),3])),
                       round(sum(results_20_60[which(TF_20_60),1] == results_20_60[which(TF_20_60),3])),
                       0)
###############################################################################################################
#######################################à

# 50 items
# 50 giudici

nitems <- 50
njudges <- 50


Sigma=matrix(0,nrow=nitems,ncol=nitems)
diag(Sigma)=1

results_50_50 <- matrix(NA, nrow = iterazioni, ncol=7)
TF_50_50 <- rep(NA,iterazioni)

for (i in 1:iterazioni) {
  
  set.seed(i)
  mu <- rnorm(nitems, 0, 1)
  true_conf <- MASS::mvrnorm(n=njudges,mu,Sigma=Sigma)
  
  rank <- t(apply(true_conf,1,rank,ties.method="min"))
  approv <-round(runif(njudges,min=0,max=nitems),0)
  approv_m <- matrix(nrow=nrow(rank),ncol=ncol(rank))
  
  for (j in 1:nrow(rank)) {
    approv_m[j,] <-  ifelse(rank[j,]<=approv[j],1,0)
  }
  
  
  X = data.frame(cbind(rank,approv_m))
  
  
  
  colmeans<-apply(X,2,mean)
  colmeans
  
  naive_r <- rank(colmeans[1:(length(colmeans)/2)])
  naive_a <- ifelse(colmeans[(length(colmeans)/2+1):length(colmeans)]>0.5,1,0)
  
  
  naive_consensus <- matrix(c(naive_r,naive_a),nrow=1)
  
  #è un preference-approval valido?
  TF_50_50[i] <- controlla(upperTriangle(score_r(naive_r)), upperTriangle(score_a(naive_a)))
  
  FC <- Final_Consensus(X,algorithm = "quick")
  
  rank_mu <- rank(mu,ties.method="min")
  approv_mu <-  ifelse(rank_mu<=round(sqrt(nitems),0),1,0)
  dgp <- data.frame(matrix(c(rank_mu,approv_mu),nrow=1))
  
  results_50_50[i,1] <- FC$D_lambda
  results_50_50[i,2]  <- FC$Elapsed
  results_50_50[i,3]  <- mean(Pref_dist2(as.matrix(X),as.matrix(naive_consensus)))
  results_50_50[i,4] <-  min(Pref_dist2(dgp,FC$Consensus))
  results_50_50[i,5]  <- Pref_dist2(dgp,naive_consensus)
  
  
  print(i)
}

summary(results_50_50[,1])
summary(results_50_50[,2])
summary(results_50_50[,3])
summary(results_50_50[,4])
summary(results_50_50[,5])
sum(TF_50_50)


Final_results[7,] <- c(median(results_50_50[,2]),0,0,sum(TF_50_50),mean(results_50_50[which(TF_50_50),1]),mean(results_50_50[which(TF_50_50),3]),
                       median(results_50_50[which(TF_50_50),1]),median(results_50_50[which(TF_50_50),3]),
                       100- round(sum(results_50_50[which(TF_50_50),1] == results_50_50[which(TF_50_50),3])),
                       round(sum(results_50_50[which(TF_50_50),1] == results_50_50[which(TF_50_50),3])),
                       0)


###############################################################################################################
#######################################à

# 50 items
# 100 giudici

nitems <- 50
njudges <- 100


Sigma=matrix(0,nrow=nitems,ncol=nitems)
diag(Sigma)=1

results_50_100 <- matrix(NA, nrow = iterazioni, ncol=7)
TF_50_100 <- rep(NA,iterazioni)

for (i in 1:iterazioni) {
  
  set.seed(i)
  mu <- rnorm(nitems, 0, 1)
  true_conf <- MASS::mvrnorm(n=njudges,mu,Sigma=Sigma)
  
  rank <- t(apply(true_conf,1,rank,ties.method="min"))
  approv <-round(runif(njudges,min=0,max=nitems),0)
  approv_m <- matrix(nrow=nrow(rank),ncol=ncol(rank))
  
  for (j in 1:nrow(rank)) {
    approv_m[j,] <-  ifelse(rank[j,]<=approv[j],1,0)
  }
  
  
  X = data.frame(cbind(rank,approv_m))
  
  
  
  colmeans<-apply(X,2,mean)
  colmeans
  
  naive_r <- rank(colmeans[1:(length(colmeans)/2)])
  naive_a <- ifelse(colmeans[(length(colmeans)/2+1):length(colmeans)]>0.5,1,0)
  
  
  naive_consensus <- matrix(c(naive_r,naive_a),nrow=1)
  
  #è un preference-approval valido?
  TF_50_100[i] <- controlla(upperTriangle(score_r(naive_r)), upperTriangle(score_a(naive_a)))
  
  FC <- Final_Consensus(X,algorithm = "quick")

  
  results_50_100[i,1] <- FC$D_lambda
  results_50_100[i,2]  <- FC$Elapsed
  results_50_100[i,3]  <- mean(Pref_dist2(as.matrix(X),as.matrix(naive_consensus)))

  
  
  print(i)
}

summary(results_50_100[,1])
summary(results_50_100[,2])
summary(results_50_100[,3])
summary(results_50_100[,4])
summary(results_50_100[,5])
sum(TF_50_100)


Final_results[8,] <- c(median(results_50_100[,2]),0,0,sum(TF_50_100),mean(results_50_100[which(TF_50_100),1]),mean(results_50_100[which(TF_50_100),3]),
                       median(results_50_100[which(TF_50_100),1]),median(results_50_100[which(TF_50_100),3]),
                       100- round(sum(results_50_100[which(TF_50_100),1] == results_50_100[which(TF_50_100),3])),
                       round(sum(results_50_100[which(TF_50_100),1] == results_50_100[which(TF_50_100),3])),
                       0)

###############################################################################################################
#######################################à

# 50 items
# 150 giudici

nitems <- 50
njudges <- 150


Sigma=matrix(0,nrow=nitems,ncol=nitems)
diag(Sigma)=1

results_50_150 <- matrix(NA, nrow = iterazioni, ncol=7)
TF_50_150 <- rep(NA,iterazioni)

for (i in 1:iterazioni) {
  
  set.seed(i)
  mu <- rnorm(nitems, 0, 1)
  true_conf <- MASS::mvrnorm(n=njudges,mu,Sigma=Sigma)
  
  rank <- t(apply(true_conf,1,rank,ties.method="min"))
  approv <-round(runif(njudges,min=0,max=nitems),0)
  approv_m <- matrix(nrow=nrow(rank),ncol=ncol(rank))
  
  for (j in 1:nrow(rank)) {
    approv_m[j,] <-  ifelse(rank[j,]<=approv[j],1,0)
  }
  
  
  X = data.frame(cbind(rank,approv_m))
  
  
  
  colmeans<-apply(X,2,mean)
  colmeans
  
  naive_r <- rank(colmeans[1:(length(colmeans)/2)])
  naive_a <- ifelse(colmeans[(length(colmeans)/2+1):length(colmeans)]>0.5,1,0)
  
  
  naive_consensus <- matrix(c(naive_r,naive_a),nrow=1)
  
  #è un preference-approval valido?
  TF_50_150[i] <- controlla(upperTriangle(score_r(naive_r)), upperTriangle(score_a(naive_a)))
  
  FC <- Final_Consensus(X,algorithm = "quick")
  
  rank_mu <- rank(mu,ties.method="min")
  approv_mu <-  ifelse(rank_mu<=round(sqrt(nitems),0),1,0)
  dgp <- data.frame(matrix(c(rank_mu,approv_mu),nrow=1))
  
  results_50_150[i,1] <- FC$D_lambda
  results_50_150[i,2]  <- FC$Elapsed
  results_50_150[i,3]  <- mean(Pref_dist2(as.matrix(X),as.matrix(naive_consensus)))
  results_50_150[i,4] <-  min(Pref_dist2(dgp,FC$Consensus))
  results_50_150[i,5]  <- Pref_dist2(dgp,naive_consensus)
  
  
  
  print(i)
}

summary(results_50_150[,1])
summary(results_50_150[,2])
summary(results_50_150[,3])
summary(results_50_150[,4])
summary(results_50_150[,5])
sum(TF_50_150)


Final_results[9,] <- c(median(results_50_150[,2]),0,0,sum(TF_50_150),mean(results_50_150[which(TF_50_150),1]),mean(results_50_150[which(TF_50_150),3]),
                       median(results_50_150[which(TF_50_150),1]),median(results_50_150[which(TF_50_150),3]),
                       100- round(sum(results_50_150[which(TF_50_150),1] == results_50_150[which(TF_50_150),3])),
                       round(sum(results_50_150[which(TF_50_150),1] == results_50_150[which(TF_50_150),3])),
                       0)








