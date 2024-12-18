# MALLOW SIMULATIONS

source("Final_Consensus_SEARCH_updated.R")


#Load the computed results
#load("mallow2.rdata")

#...or compute them by yourself!




fix_preference_approval = function(x){
  xa = x[((length(x)/2)+1):length(x)]
  xr = x[1:(length(x)/2)]
  w1 <- which(xa==1)
  w0 <- which(xa==0)
  
  xr[w1] <- floor(runif(length(w1),1,10))
  xr[w0] <- floor(runif(length(w0),11,20))
  return(c(rank(xr,ties.method = "min"),xa))
}



fix_new <- function(x,base){
  
  xa = x[((length(x)/2)+1):length(x)]
  xr = x[1:(length(x)/2)]
  w1 <- which(xa==1)
  w0 <- which(xa==0)
  
  xr[w1] <- floor(runif(length(w1),1,10))
  xr[w0] <- floor(runif(length(w0),11,20))
  return(c(rank(xr,ties.method = "min"),xa))
  
}

#####################################################
# 10 items

##########################################################################################

nitems <- 10
njudges <- 10


###################################################
# SPAMM for generating preference approvals
###################################################

# Low noise
xr <- 1:10
xa <- rep(c(1,0),c(5,5))

njudges <- 2000
Sigma=matrix(0,nrow=nitems,ncol=nitems)
diag(Sigma)=1
mu <- 1:10/3

true_conf <- MASS::mvrnorm(n=njudges,mu,Sigma=Sigma)
ties <- rpois(nrow(true_conf),nitems*0.1)
rank <- t(apply(true_conf,1,rank,ties.method="min"))

approv <-rpois(nrow(true_conf),nitems)
approv_m <- matrix(nrow=nrow(rank),ncol=ncol(rank))

for (j in 1:nrow(rank)) {
  approv_m[j,] <-  ifelse(rank[j,]<=approv[j],1,0)
  rank[j,] <- ifelse(rank[j,]<=ties[j],min(rank[j,rank[j,] <= ties[j]]),rank[j,])
}

X1_10 = data.frame(cbind(rank,approv_m))
di_10 <- Pref_dist2(as.matrix(X1_10),c(xr,xa))
# summary(di_10)
# head(sort(di_10),20)
# plot(ecdf(di_10))


# High noise
nitems <- 10
njudges <- 2000

Sigma=matrix(0,nrow=nitems,ncol=nitems)
diag(Sigma)=1
set.seed(i)

mu <- 10:1/3

true_conf <- MASS::mvrnorm(n=njudges,mu,Sigma=Sigma)
ties <- rpois(nrow(true_conf),nitems*0.1)
rank <- t(apply(true_conf,1,rank,ties.method="min"))

approv <-rpois(nrow(true_conf),nitems)
approv_m <- matrix(nrow=nrow(rank),ncol=ncol(rank))

for (j in 1:nrow(rank)) {
  approv_m[j,] <-  ifelse(rank[j,]<=approv[j],1,0)
  rank[j,] <- ifelse(rank[j,]<=ties[j],min(rank[j,rank[j,] <= ties[j]]),rank[j,])
}

X2_10 = data.frame(cbind(rank,approv_m))
di2_10 <- Pref_dist2(as.matrix(X2_10),c(xr,xa))
# summary(di2_10)
# head(sort(di2_10),20)
# plot(ecdf(di2_10))


######### Medium noise
nitems <- 10
njudges <- 2000

Sigma=matrix(0,nrow=nitems,ncol=nitems)
diag(Sigma)=1

order <- c(10,2,3,7,6,5,4,8,9,1)


mu <-  1:10/3
mu <- mu[order(order)]


true_conf <- MASS::mvrnorm(n=njudges,mu,Sigma=Sigma)
ties <- rpois(nrow(true_conf),nitems*0.1)
rank <- t(apply(true_conf,1,rank,ties.method="min"))

approv <-rpois(nrow(true_conf),nitems)
approv_m <- matrix(nrow=nrow(rank),ncol=ncol(rank))

for (j in 1:nrow(rank)) {
  approv_m[j,] <-  ifelse(rank[j,]<=approv[j],1,0)
  rank[j,] <- ifelse(rank[j,]<=ties[j],min(rank[j,rank[j,] <= ties[j]]),rank[j,])
}

X3_10 = data.frame(cbind(rank,approv_m))
di3_10 <- Pref_dist2(as.matrix(X3_10),c(xr,xa))
#summary(di3_10)


# Building reduced universe of preference approvals

# Choosing central permutation

xr_10 <- 1:10
xa_10 <- rep(c(1,0),each=5)

dists_10 <- c(di_10,di2_10,di3_10)

universe_10 <- rbind(X1_10,X2_10,X3_10)
indici_10 <- seq(0,1, by=0.01)

probs0_10 <- exp(-0.1*indici_10)/ sum(exp(-0.1*indici_10))
probs5_10 <- exp(-5*indici_10)/ sum(exp(-5*indici_10))
probs10_10 <- exp(-10*indici_10)/ sum(exp(-1*indici_10))


Final_10_items_new <- matrix(NA,nrow=9,ncol=3)
df_10 = expand.grid(c(0.01,0.5,1),c(10,20,30))
rownames(Final_10_items_new) <- paste("lambda=",df_10[,1], "n=", df_10[,2]) 

r10 <- list()
for (k in 1:9) {
  results_10 = matrix(NA,nrow=100,ncol=3)
  if(df_10$Var1[k]==0.01) probs <- probs0_10
  if(df_10$Var1[k]==0.5) probs <- probs5_10
  if(df_10$Var1[k]==1) probs <- probs10_10
  
  
  for (j in 1:100) {
    
    set.seed(j)
    gen0_10 <- indici_10[sample(1:length(indici_10), 
                                df_10$Var2[k] ,
                                replace = TRUE,
                                prob = probs)]
    
    generati <- matrix(NA, length(gen0_10), 20)
    
    
    for (i in 1:length(gen0_10)) {
      chi <- which.min(abs(dists_10-gen0_10[i]))
      generati[i,] <- as.numeric(universe_10[chi,])
    }
    
    invisible(capture.output(FC <- Final_Consensus(generati,algo="quick")))
    
    Pc <- mean(Pref_dist2(as.matrix(generati),c(xr_10,xa_10)))
    
    results_10[j,1] <- round(Pc - FC$D_lambda,4)
    results_10[j,2] <- round((Pc - FC$D_lambda)/Pc,4)
    results_10[j,3] <- round((Pc - FC$D_lambda)/FC$D_lambda,4)
    
    
    cat("Iterazione\n",j)
  }
  
  
  Final_10_items_new[k,1] <- mean(results_10[,1])
  Final_10_items_new[k,2] <- mean(results_10[,2])
  Final_10_items_new[k,3] <- mean(results_10[,3])
  
  cat("macro iterazione", k) 
  r10[[k]] <- results_10
}



###############################
# 20 items

##################

###################################################
# SPAMM for generating preference approvals
###################################################

xr <- 1:20
xa <- rep(c(1,0),c(10,10))

# Low noise
nitems <- 20
njudges <- 2000

Sigma=matrix(0,nrow=nitems,ncol=nitems)
diag(Sigma)=1
set.seed(i)
mu <- seq(0,1,by=0.05)

mu <- 1:20/2

true_conf <- MASS::mvrnorm(n=njudges,mu,Sigma=Sigma)
ties <- rpois(nrow(true_conf),nitems*0.1)
rank <- t(apply(true_conf,1,rank,ties.method="min"))

approv <-rpois(nrow(true_conf),nitems)
approv_m <- matrix(nrow=nrow(rank),ncol=ncol(rank))

for (j in 1:nrow(rank)) {
  approv_m[j,] <-  ifelse(rank[j,]<=approv[j],1,0)
  rank[j,] <- ifelse(rank[j,]<=ties[j],min(rank[j,rank[j,] <= ties[j]]),rank[j,])
}

X1 = data.frame(cbind(rank,approv_m))
di <- Pref_dist2(as.matrix(X1),c(xr,xa))


# High noise
nitems <- 20
njudges <- 2000

Sigma=matrix(0,nrow=nitems,ncol=nitems)
diag(Sigma)=1
set.seed(i)
mu <- seq(0,1,by=0.05)

mu <- 20:1/2

true_conf <- MASS::mvrnorm(n=njudges,mu,Sigma=Sigma)
ties <- rpois(nrow(true_conf),nitems*0.1)
rank <- t(apply(true_conf,1,rank,ties.method="min"))

approv <-rpois(nrow(true_conf),nitems)
approv_m <- matrix(nrow=nrow(rank),ncol=ncol(rank))

for (j in 1:nrow(rank)) {
  approv_m[j,] <-  ifelse(rank[j,]<=approv[j],1,0)
  rank[j,] <- ifelse(rank[j,]<=ties[j],min(rank[j,rank[j,] <= ties[j]]),rank[j,])
}

X2 = data.frame(cbind(rank,approv_m))
di2 <- Pref_dist2(as.matrix(X2),c(xr,xa))


######### Medium noise
nitems <- 20
njudges <- 2000

Sigma=matrix(0,nrow=nitems,ncol=nitems)
diag(Sigma)=1
set.seed(i)

order <-  c(20,19,3,4,16,6,7,8,9,10,
            11,12,13,14,15,5,17,18,2,1)


mu <-  1:20/2
mu <- mu[order(order)]



true_conf <- MASS::mvrnorm(n=njudges,mu,Sigma=Sigma)
ties <- rpois(nrow(true_conf),nitems*0.1)
rank <- t(apply(true_conf,1,rank,ties.method="min"))

approv <-rpois(nrow(true_conf),nitems)
approv_m <- matrix(nrow=nrow(rank),ncol=ncol(rank))

for (j in 1:nrow(rank)) {
  approv_m[j,] <-  ifelse(rank[j,]<=approv[j],1,0)
  rank[j,] <- ifelse(rank[j,]<=ties[j],min(rank[j,rank[j,] <= ties[j]]),rank[j,])
}

X3 = data.frame(cbind(rank,approv_m))
di3 <- Pref_dist2(as.matrix(X3),c(xr,xa))


summary(c(di,di2,di3))



# Creatint the reduced universe

# Choosing the central permutation 

xr_20 <- 1:20
xa_20 <- rep(c(1,0),each=10)

dists_20 <- c(di,di2,di3)

universe_20 <- rbind(X1,X2,X3)


indici_20 <- seq(0,1, by=0.01)

probs0_20 <- exp(-0.1*indici_20)/ sum(exp(-0.1*indici_20))
probs5_20 <- exp(-5*indici_20)/ sum(exp(-5*indici_20))
probs10_20 <- exp(-10*indici_20)/ sum(exp(-1*indici_20))


Final_20_items_new <- matrix(NA,nrow=9,ncol=3)
df_20 = expand.grid(c(0.01,0.5,1),c(20,40,60))
rownames(Final_20_items_new) <- paste("lambda=",df_20[,1], "n=", df_20[,2]) 

r20 = list()
for (k in 1:9) {
  results_20 = matrix(NA,nrow=100,ncol=3)
  if(df_20$Var1[k]==0.01) probs <- probs0_20
  if(df_20$Var1[k]==0.5) probs <- probs5_20
  if(df_20$Var1[k]==1) probs <- probs10_20
  
  
  for (j in 1:100) {
    
    set.seed(j)
    gen0_50 <- indici_20[sample(1:length(indici_20), df_20$Var2[k] , replace = TRUE, prob = probs)]
    generati <- matrix(NA, length(gen0_50), 40)
    
    
    for (i in 1:length(gen0_50)) {
      chi <- which.min(abs(dists_20-gen0_50[i]))
      generati[i,] <- as.numeric(universe_20[chi,])
    }
    
    invisible(capture.output(FC <- Final_Consensus(generati,algo="quick")))
    
    Pc <- mean(Pref_dist2(as.matrix(generati),c(xr_20,xa_20)))
    
    results_20[j,1] <- round(Pc - FC$D_lambda,4)
    results_20[j,2] <- round((Pc - FC$D_lambda)/Pc,4)
    results_20[j,3] <- round((Pc - FC$D_lambda)/FC$D_lambda,4)
    
    
    cat("Iterazione\n",j)
    
  }
  
  
  Final_20_items_new[k,1] <- mean(results_20[,1])
  Final_20_items_new[k,2] <- mean(results_20[,2])
  Final_20_items_new[k,3] <- mean(results_20[,3])
  
  cat("macro iterazione", k) 
  r20[[k]] <- results_20
}


#################################

# 50 items 

##########################################################################################


###################################################
# SPAMM for generating preference approvals
###################################################

xr <- 1:50
xa <- rep(c(1,0),c(25,25))

# Low noise
nitems <- 50
njudges <- 2000

Sigma=matrix(0,nrow=nitems,ncol=nitems)
diag(Sigma)=1
set.seed(i)


mu <- 1:50

true_conf <- MASS::mvrnorm(n=njudges,mu,Sigma=Sigma)
ties <- rpois(nrow(true_conf),nitems*0.1)
rank <- t(apply(true_conf,1,rank,ties.method="min"))

approv <-rpois(nrow(true_conf),nitems)
approv_m <- matrix(nrow=nrow(rank),ncol=ncol(rank))

for (j in 1:nrow(rank)) {
  approv_m[j,] <-  ifelse(rank[j,]<=approv[j],1,0)
  rank[j,] <- ifelse(rank[j,]<=ties[j],min(rank[j,rank[j,] <= ties[j]]),rank[j,])
}

X1_50 = data.frame(cbind(rank,approv_m))
di_50 <- Pref_dist2(as.matrix(X1_50),c(xr,xa))


# High noise
nitems <- 50
njudges <- 2000

Sigma=matrix(0,nrow=nitems,ncol=nitems)
diag(Sigma)=1

mu <- 50:1/2

true_conf <- MASS::mvrnorm(n=njudges,mu,Sigma=Sigma)
ties <- rpois(nrow(true_conf),nitems*0.1)
rank <- t(apply(true_conf,1,rank,ties.method="min"))

approv <-rpois(nrow(true_conf),nitems)
approv_m <- matrix(nrow=nrow(rank),ncol=ncol(rank))

for (j in 1:nrow(rank)) {
  approv_m[j,] <-  ifelse(rank[j,]<=approv[j],1,0)
  rank[j,] <- ifelse(rank[j,]<=ties[j],min(rank[j,rank[j,] <= ties[j]]),rank[j,])
}

X2_50 = data.frame(cbind(rank,approv_m))
di2_50 <- Pref_dist2(as.matrix(X2_50),c(xr,xa))


######### Medium noise 1
nitems <- 50
njudges <- 2000

Sigma=matrix(0,nrow=nitems,ncol=nitems)
diag(Sigma)=1


order <- c(50,49,48,47,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,
           26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,4,3,2,1)


mu <-  1:50/20
mu <- mu[order(order)]


true_conf <- MASS::mvrnorm(n=njudges,mu,Sigma=Sigma)
ties <- rpois(nrow(true_conf),nitems*0.1)
rank <- t(apply(true_conf,1,rank,ties.method="min"))

approv <-rpois(nrow(true_conf),nitems)
approv_m <- matrix(nrow=nrow(rank),ncol=ncol(rank))

for (j in 1:nrow(rank)) {
  approv_m[j,] <-  ifelse(rank[j,]<=approv[j],1,0)
  rank[j,] <- ifelse(rank[j,]<=ties[j],min(rank[j,rank[j,] <= ties[j]]),rank[j,])
}

X3_50 = data.frame(cbind(rank,approv_m))
di3_50 <- Pref_dist2(as.matrix(X3_50),c(xr,xa))



######### Medium noise 2
nitems <- 50
njudges <- 2000

Sigma=matrix(0,nrow=nitems,ncol=nitems)
diag(Sigma)=1


order <- c(50,49,48,47,46,45,44,43,42,41,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,
           26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,10,9,8,7,6,5,4,3,2,1)



mu <-  1:50/10
mu <- mu[order(order)]


true_conf <- MASS::mvrnorm(n=njudges,mu,Sigma=Sigma)
ties <- rpois(nrow(true_conf),nitems*0.1)
rank <- t(apply(true_conf,1,rank,ties.method="min"))

approv <-rpois(nrow(true_conf),nitems)
approv_m <- matrix(nrow=nrow(rank),ncol=ncol(rank))

for (j in 1:nrow(rank)) {
  approv_m[j,] <-  ifelse(rank[j,]<=approv[j],1,0)
  rank[j,] <- ifelse(rank[j,]<=ties[j],min(rank[j,rank[j,] <= ties[j]]),rank[j,])
}

X4_50 = data.frame(cbind(rank,approv_m))
di4_50 <- Pref_dist2(as.matrix(X4_50),c(xr,xa))


summary(c(di_50,di2_50,di3_50,di4_50))



# Creating reduced universe

# Choosing central permutation

xr_50 <- 1:50
xa_50 <- rep(c(1,0),each=25)

dists_50 <- c(di_50,di2_50,di3_50,di4_50)

universe_50 <- rbind(X1_50,X2_50,X3_50,X4_50)


indici_50 <- seq(0,1, by=0.01)

probs0_50 <- exp(-0.1*indici_50)/ sum(exp(-0.1*indici_50))
probs5_50 <- exp(-5*indici_50)/ sum(exp(-5*indici_50))
probs10_50 <- exp(-10*indici_50)/ sum(exp(-1*indici_50))


Final_50_items_new <- matrix(NA,nrow=9,ncol=3)
df_50 = expand.grid(c(0.01,0.5,1),c(50,100,150))
rownames(Final_50_items_new) <- paste("lambda=",df_50[,1], "n=", df_50[,2]) 


r50 <- list()
for (k in 1:9) {
  results_50 = matrix(NA,nrow=100,ncol=3)
  if(df_50$Var1[k]==0.01) probs <- probs0_50
  if(df_50$Var1[k]==0.5) probs <- probs5_50
  if(df_50$Var1[k]==1) probs <- probs10_50
  
  
  for (j in 1:100) {
    
    set.seed(j)
    gen0_50 <- indici_50[sample(1:length(indici_50), df_50$Var2[k] , replace = TRUE, prob = probs)]
    generati <- matrix(NA, length(gen0_50), 100)
    
    
    for (i in 1:length(gen0_50)) {
      chi <- which.min(abs(dists_50-gen0_50[i]))
      generati[i,] <- as.numeric(universe_50[chi,])
    }
    
    invisible(capture.output(FC <- Final_Consensus(generati,algo="quick")))
    
    Pc <- mean(Pref_dist2(generati,c(xr_50,xa_50)))
    
    results_50[j,1] <- round(Pc - FC$D_lambda,4)
    results_50[j,2] <- round((Pc - FC$D_lambda)/Pc,4)
    results_50[j,3] <- round((Pc - FC$D_lambda)/FC$D_lambda,4)
    
    
    cat("Iterazione\n",j)
    r50[[k]] <- results_50
  }
  
  
  Final_50_items_new[k,1] <- mean(results_50[,1])
  Final_50_items_new[k,2] <- mean(results_50[,2])
  Final_50_items_new[k,3] <- mean(results_50[,3])
  
  cat("macro iterazione", k) 
}





###################
# Summary of the results
dt <- round(cbind(Final_10_items_new[,1],
                  Final_10_items_new[,2]*100,
                  Final_20_items_new[,1],
                  Final_20_items_new[,2]*100, 
                  Final_50_items_new[,1],
                  Final_50_items_new[,2]*100),3)

colnames(dt) <- c(10,10,20,20,50,50)


