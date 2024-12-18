# FULL consensus

library(ConsRank)
library(gtools)
library(gdata)
library(tidyr)
library(dplyr)
library(hms)

#############
# UTILITIES #
#############
binary_dist <- function (x, y = NULL) 
{
  if (is(y, "NULL")) {
    N <- nrow(x)
    indice <- combinations(N, 2, repeats.allowed = T)
    kx <- mat.or.vec(N, N)
    for (j in 1:nrow(indice)) {
      kx[indice[j, 1], indice[j, 2]] = kx[indice[j, 2], 
                                          indice[j, 1]] <- sum(abs(x[indice[j, 1], ] - x[indice[j, 2], ]))
    }
    kx
  }
  else {
    if (is(x, "numeric") & !is(x, "matrix")) {
      x <- matrix(x, ncol = length(x))
    }
    if (is(y, "numeric") & !is(y, "matrix")) {
      y <- matrix(y, ncol = length(y))
    }
    nx <- nrow(x)
    ny <- nrow(y)
    indice <- data.frame(crossing(1:nx, 1:ny))
    kx <- matrix(NA, nx, ny)
    for (j in 1:nrow(indice)) {
      kx[indice[j, 1], indice[j, 2]] <- sum(abs(x[indice[j,1], ] - y[indice[j, 2], ]))
    }
    kx
  }
}

Pref_dist2 <- function(X,Y=NULL,lambda=0.5){
  if (is(X, "numeric") & !is(X, "matrix")) {
    X <- matrix(X, ncol = length(X))
  }
  x_r <- X[,1:(ncol(X)/2),drop=F]
  x_a <- X[,(ncol(X)/2+1):ncol(X),drop=F]
  if(is.null(Y)){
    y_r <- NULL
    y_a <- NULL
  }
  else{
    if (is(Y, "numeric") & !is(Y, "matrix")) {
      Y <- matrix(Y, ncol = length(Y))
    }
    y_r <- Y[,1:(ncol(Y)/2),drop=F]
    y_a <- Y[,(ncol(Y)/2+1):ncol(Y),drop=F] 
  }
  
  fkx <-
    lambda * as.matrix(kemenyd(X = x_r, Y = y_r)) / (ncol(x_r) * (ncol(x_r) -
                                                                    1)) + (1 - lambda) * binary_dist(x = x_a, y = y_a) / ncol(x_r)
  
  return(fkx)
  
}

score_r=function (X) 
{
  itemnames <- names(X)
  if (is(X, "numeric") & !is(X, "matrix")) {
    X <- matrix(X, ncol = length(X))
  }
  c <- ncol(X)
  sm <- matrix(0, c, c)
  colnames(sm) <- itemnames
  row.names(sm) <- itemnames
  X <- as.numeric(X)
  for (j in 1:c) {
    diffs <- sign(X[j] - X[setdiff(1:c, j)])
    ind <- setdiff(1:c, j)
    sm[j, ind] <- diffs
  }
  idn <- is.na(sm)
  sm <- -sm
  sm[idn] <- 0
  sm[sm==0]=0.5
  diag(sm)=0
  sm
}

score_a=function (X) 
{
  itemnames <- names(X)
  if (is(X, "numeric") & !is(X, "matrix")) {
    X <- matrix(X, ncol = length(X))
  }
  c <- ncol(X)
  sm <- matrix(0, c, c)
  colnames(sm) <- itemnames
  row.names(sm) <- itemnames
  X <- as.numeric(X)
  for (j in 1:c) {
    diffs <- sign(X[j] - X[setdiff(1:c, j)])
    ind <- setdiff(1:c, j)
    sm[j, ind] <- diffs
  }
  idn <- is.na(sm)
  sm[idn] <- 0
  sm[sm==0]=0.5
  diag(sm)=0
  sm
}


controlla = function(X,Y){
  logic=rep(NA,length(X))
  for (i in 1:length(X)) {
    x=X[i]
    y=Y[i]
    
    if(x==1){
      logic[i]<-(y>=0.5)
    }
    if(x==0.5){
      logic[i]<-(y==0.5)
    }
    if(x==-1){
      logic[i]<-(y<=0.5)
    }
  }
  all(logic)
}


FindApproval <- function(x){
  x <- as.numeric(x)
  buckets <- length(table(x))
  approval <- matrix(0,nrow=(buckets+1), ncol=length(x))
  rank <- matrix(x,nrow=(buckets+1), ncol=length(x),byrow=T)
  
  for (i in 1:buckets) {
    wb <- sort(unique(x))
    which(x <= wb[i])
    approval[(i+1),which(x <= wb[i])] <- 1
  }
  cbind.data.frame(rank,approval)
}

genera_preference = function(m) {
  start=Sys.time()
  pa <- do.call(rbind.data.frame,apply(univranks(m)$Cuniv$R,1,FindApproval))
  cat(paste("Generated", nrow(pa), "preferences in", as_hms(Sys.time()-start)))
  
  colnames(pa) = paste("V",1:(2*m),sep="")
  return(
    pa
  )

}

############################################
# Core function used by the DIVA Algorithm #
############################################

pa_cons<- function(X, x_a, algorithm = "decor"){
  if(is.null(colnames(X)))colnames(X) = paste("V",1:ncol(X),sep="")
  #Simple case
  #If there is a PA of only 1 (All approved) or 0 (Nothing approved) we need just to find a simple consensus
  if(length(names(table(x_a))) != 2){ 
    # Finding the consensus ranking 
    cr_1 <- as.data.frame((consrank( X, algorithm = algorithm)$Consensus))
    cr_0 <- NULL
  }
  else{
    #Selecting all the "approved" or "not approved" items
    x_a1 = which(x_a == 1)
    x_a0 = which(x_a == 0)
    
    #Simple case
    #If there is just one "approved" item, this will be the top in the consensus approval
    if(length(x_a1) == 1){
      cr_1 = as.data.frame(1)
      colnames(cr_1) = colnames(X)[x_a1]
    }
    else{
      #Creating a data matrix that store the preferences but remove the not considered items
      x1 = reordering(as.matrix(X[,x_a1]))
      #Compute the consensus ranking for "approved" items
      if(length(x_a1) == 2){
        
        
        cr_1 <- as.data.frame(matrix(rank(colSums(x1),ties.method = "min"),nrow=1))
        colnames(cr_1) = colnames(X)[x_a1]  
      }else{
        cr_1 <- as.data.frame((consrank(x1, algorithm = algorithm)$Consensus)) }
    }
    
    #Simple case
    #If there is just one "non approved" item, this will be the last in the consensus approval
    if(length(x_a0) == 1){
      cr_0 = as.data.frame(1)
      colnames(cr_0) = colnames(X)[x_a0]
    }
    else{
      #Creating a data matrix that store the preferences but remove the not considered items
      x0 = reordering(as.matrix(X[,x_a0]))
      #Compute the consensus ranking for "non approved" items
      if(length(x_a0) == 2){
        
        cr_0 <- as.data.frame(matrix(rank(colSums(x0),ties.method = "min"),nrow=1))
        colnames(cr_0) = colnames(X)[x_a0] 
      }else{
        cr_0 <- as.data.frame((consrank( x0, algorithm = algorithm)$Consensus))}
    }
  }
  
  if(!is.null(cr_0)){
    #This return all the possible combinations of obtained consensus from the "approved" and "non approved" component
    Consensus = data.frame(matrix(NA, nrow = nrow(cr_1) * nrow(cr_0), ncol = ncol(X)+length(x_a)))
    colnames(Consensus) = paste("V", 1:length(Consensus), sep="")
    count = 1
    
    
    for(i in 1:nrow(cr_1)){
      for(j in 1:nrow(cr_0)){
        # The consensus ranking is the merge of cr_1 and cr_0 
        cr = unlist(c(cr_1[i,], cr_0[j,]+length(x_a1)))
        names(cr) = c(colnames(cr_1), colnames(cr_0))
        # Reordering the items to make them follow the same initial order
        cr = reordering(matrix(cr[colnames(X)], ncol = ncol(X)))
        cr = unlist(as.data.frame(cr))
        # The consensus approval is the previously computed consensus with the approval component concateneted in tail
        Consensus[count,] = c(cr, x_a)
        count = count+1
      }
    }
  }
  
  
  else{
    Consensus = data.frame(matrix(NA, nrow = nrow(cr_1), ncol = ncol(X)+length(x_a)))
    names(Consensus) = paste("V", 1:length(Consensus), sep="")
    Consensus[,1:ncol(cr_1)] = cr_1
    Consensus[,(ncol(cr_1)+1):ncol(Consensus)] = x_a
  }
  
  return(list(Consensus=Consensus))
}


#####################################
# Main function: the DIVA Algorithm #
#####################################
Final_Consensus <- function(X, search=FALSE, algorithm="decor",lambda=0.5){
  
  start = Sys.time()
  # First check to see if there is a simple solution with a valid preference approval
  
  m <- ncol(X)/2
  #Selecting the rankings component
  x_r <- X[,1:(ncol(X)/2)]
  #Selecting the approvals component
  x_a <- X[,(ncol(X)/2+1):ncol(X)]
  # Computing consensus ranking 
  cr <- as.data.frame((consrank( x_r, algorithm = algorithm)$Consensus))
  # Finding the  consensus approval 
  amean <- apply(x_a,2,mean)
  num <- sum(amean == 0.5)
  
  
  
  if(search & num>0){
    who <- which(amean == 0.5)
    dispo = permutations(2, num, repeats.allowed = T) - 1
    cr2 <-  matrix(rep(amean,2^num), nrow= 2^num, ncol= m, byrow=T)
    cr2[,who] <- dispo
    cr2 <- ifelse(cr2>0.5 ,1 ,0)} 
  else{
    cr2 <- data.frame(t(ifelse(apply(x_a,2,mean)>0.5 ,1 ,0)))
  }
  
  
  
  #If they are compatible: STOP!
  
  check <- matrix(NA, nrow = nrow(cr), ncol= nrow(cr2))
  
  for (j in 1:nrow(cr2)) {
    for (i in 1:nrow(cr)) {
      check[i,j] <-
        controlla(upperTriangle(score_r(cr[i, ])), upperTriangle(score_a(cr2[j,])))
      
    }
  }
  
  
  if (any(check)) {
    print("Solution 0")
    
    
    select <- which(check == T, arr.ind = TRUE)
    
    if(nrow(select)==1){
      candidati =  data.frame(matrix(c(as.numeric(cr[(select[,1]),]),as.numeric(cr2[(select[,2]),])),ncol=(m*2)))
    }
    else{
      candidati =  data.frame(cr[(select[,1]),],cr2[(select[,2]),])}
    
    
    mat_dist <- Pref_dist2(as.matrix(X),as.matrix(candidati),lambda = lambda)
    
    Consensus= candidati[near(apply(mat_dist,2,mean) , min(apply(mat_dist,2,mean))),]
    D_lambda = min(apply(mat_dist,2,mean))
    
  }
  else{
    # First Branch
    # FROM RANK TO APPROVAL
    print("Solution 1")
    approval = apply(cr, 1, FindApproval)
    candidati <- do.call(rbind.data.frame, approval)
    colnames(candidati)<- paste("V",1:(2*m),sep="")
    
    
    # Computing the distance of each item from the judges set
    mat_dist <- Pref_dist2(as.matrix(X),as.matrix(candidati),lambda = lambda)
    
    Consensus_1= candidati[near(apply(mat_dist,2,mean) , min(apply(mat_dist,2,mean))),]
    D_lambda_1 = min(apply(mat_dist,2,mean))
    
    # Second Branch
    # FROM APPROVAL TO RANK
    print("Solution 2")
    
    c2 = list()
    for (i in 1:nrow(cr2)) {
      c2[[i]] <- pa_cons(X = x_r, x_a = as.numeric(cr2[i,]),algorithm = algorithm)$Consensus
    }
    
    candidati2 <- do.call(rbind.data.frame, c2)
    colnames(candidati2)<- paste("V",1:(2*m),sep="")
    
    mat_dist2 <- Pref_dist2(as.matrix(X),as.matrix(candidati2),lambda = lambda)
    Consensus_2 <- candidati2[near(apply(mat_dist2,2,mean) , min(apply(mat_dist2,2,mean))),]
    D_lambda_2 <- min(apply(mat_dist2,2,mean))
    
    
    D_lambda= min(D_lambda_1,D_lambda_2)
    
    if (near(D_lambda_1, D_lambda_2)) {
      Consensus = rbind(Consensus_1, Consensus_2)
    }
    else{
      if (D_lambda_1 > D_lambda_2) {
        Consensus = Consensus_2
      }
      else{
        Consensus = Consensus_1
      }
    } 
     
  }
  
  
  El <- as_hms(Sys.time()-start)
  
  return(list(D_lambda=D_lambda, Consensus=Consensus, Elapsed = El))
  
}

#Simple example of DIVA Algorithm

# universe=genera_preference(5)
# 
# X=universe[sample(1:nrow(universe),20,replace=T),]
# Final_Consensus(X, search = TRUE)
# Final_Consensus(X)
# Final_Consensus(X,algorithm = "fast")
# Final_Consensus(X,algorithm = "quick")



######################################################
# DIVA Algorithm version for the sensibility analysis#
######################################################
Consensus_sensibility <- function(X, search=FALSE, algorithm="quick",lambda=0.5){
  
# First check to see if there is a simple solution with a valid preference approval
  
  m <- ncol(X)/2
  #Selecting the rankings component
  x_r <- X[,1:(ncol(X)/2)]
  #Selecting the approvals component
  x_a <- X[,(ncol(X)/2+1):ncol(X)]
  # Computing consensus ranking 
  cr <- as.data.frame((consrank( x_r, algorithm = algorithm)$Consensus))
  # Finding the  consensus approval  
  amean <- apply(x_a,2,mean)
  num <- sum(amean == 0.5)
  
  if(search & num>0){
    who <- which(amean == 0.5)
    dispo = permutations(2, num, repeats.allowed = T) - 1
    cr2 <-  matrix(rep(amean,2^num), nrow= 2^num, ncol= m, byrow=T)
    cr2[,who] <- dispo
    cr2 <- ifelse(cr2>0.5 ,1 ,0)} else{
    cr2 <- data.frame(t(ifelse(apply(x_a,2,mean)>0.5 ,1 ,0)))
  }
  
  #If they are compatible: STOP!
  
  check <- matrix(NA, nrow = nrow(cr), ncol= nrow(cr2))
  
  for (j in 1:nrow(cr2)) {
    for (i in 1:nrow(cr)) {
      check[i,j] <-
        controlla(upperTriangle(score_r(cr[i, ])), upperTriangle(score_a(cr2[j,])))
      
    }
  }
  
  
  if (any(check)) {
    print("Solution 0")
    
    
    select <- which(check == T, arr.ind = TRUE)
    
    if(nrow(select)==1){
      candidati =  data.frame(matrix(c(as.numeric(cr[(select[,1]),]),as.numeric(cr2[(select[,2]),])),ncol=(m*2)))
    }
    else{
      candidati =  data.frame(cr[(select[,1]),],cr2[(select[,2]),])
      candidati = candidati[1,]}
    
    
    mat_dist_0 <- Pref_dist2(as.matrix(X),as.matrix(candidati),lambda = 0)
    mat_dist_1 <- Pref_dist2(as.matrix(X),as.matrix(candidati),lambda = 0.1)
    mat_dist_2 <- Pref_dist2(as.matrix(X),as.matrix(candidati),lambda = 0.2)
    mat_dist_3 <- Pref_dist2(as.matrix(X),as.matrix(candidati),lambda = 0.3)
    mat_dist_4 <- Pref_dist2(as.matrix(X),as.matrix(candidati),lambda = 0.4)
    mat_dist_5 <- Pref_dist2(as.matrix(X),as.matrix(candidati),lambda = 0.5)
    mat_dist_6 <- Pref_dist2(as.matrix(X),as.matrix(candidati),lambda = 0.6)
    mat_dist_7 <- Pref_dist2(as.matrix(X),as.matrix(candidati),lambda = 0.7)
    mat_dist_8 <- Pref_dist2(as.matrix(X),as.matrix(candidati),lambda = 0.8)
    mat_dist_9 <- Pref_dist2(as.matrix(X),as.matrix(candidati),lambda = 0.9)
    mat_dist_10 <- Pref_dist2(as.matrix(X),as.matrix(candidati),lambda = 1)
    
    D_lambda_final = c(min(apply(mat_dist_0,2,mean)),
                 min(apply(mat_dist_1,2,mean)),
                 min(apply(mat_dist_2,2,mean)),
                 min(apply(mat_dist_3,2,mean)),
                 min(apply(mat_dist_4,2,mean)),
                 min(apply(mat_dist_5,2,mean)),
                 min(apply(mat_dist_6,2,mean)),
                 min(apply(mat_dist_7,2,mean)),
                 min(apply(mat_dist_8,2,mean)),
                 min(apply(mat_dist_9,2,mean)),
                 min(apply(mat_dist_10,2,mean)))
    
  }
  else{
    # First Branch
    # FROM RANK TO APPROVAL
    print("Solution 1")
    approval = apply(cr, 1, FindApproval)
    candidati <- do.call(rbind.data.frame, approval)
    colnames(candidati)<- paste("V",1:(2*m),sep="")
    

    mat_dist_0 <- Pref_dist2(as.matrix(X),as.matrix(candidati),lambda = 0)
    mat_dist_1 <- Pref_dist2(as.matrix(X),as.matrix(candidati),lambda = 0.1)
    mat_dist_2 <- Pref_dist2(as.matrix(X),as.matrix(candidati),lambda = 0.2)
    mat_dist_3 <- Pref_dist2(as.matrix(X),as.matrix(candidati),lambda = 0.3)
    mat_dist_4 <- Pref_dist2(as.matrix(X),as.matrix(candidati),lambda = 0.4)
    mat_dist_5 <- Pref_dist2(as.matrix(X),as.matrix(candidati),lambda = 0.5)
    mat_dist_6 <- Pref_dist2(as.matrix(X),as.matrix(candidati),lambda = 0.6)
    mat_dist_7 <- Pref_dist2(as.matrix(X),as.matrix(candidati),lambda = 0.7)
    mat_dist_8 <- Pref_dist2(as.matrix(X),as.matrix(candidati),lambda = 0.8)
    mat_dist_9 <- Pref_dist2(as.matrix(X),as.matrix(candidati),lambda = 0.9)
    mat_dist_10 <- Pref_dist2(as.matrix(X),as.matrix(candidati),lambda = 1)
    
    D_lambda_1 = c(min(apply(mat_dist_0,2,mean)),
                 min(apply(mat_dist_1,2,mean)),
                 min(apply(mat_dist_2,2,mean)),
                 min(apply(mat_dist_3,2,mean)),
                 min(apply(mat_dist_4,2,mean)),
                 min(apply(mat_dist_5,2,mean)),
                 min(apply(mat_dist_6,2,mean)),
                 min(apply(mat_dist_7,2,mean)),
                 min(apply(mat_dist_8,2,mean)),
                 min(apply(mat_dist_9,2,mean)),
                 min(apply(mat_dist_10,2,mean)))

    # Second Branch
    # FROM APPROVAL  TO RANK
    print("Solution 2")
    
    c2 = list()
    for (i in 1:nrow(cr2)) {
      c2[[i]] <- pa_cons(X = x_r, x_a = as.numeric(cr2[i,]),algorithm = algorithm)$Consensus
    }
    
  
    
    candidati2 <- do.call(rbind.data.frame, c2)
    colnames(candidati2)<- paste("V",1:(2*m),sep="")
  
    mat_dist_0_2 <- Pref_dist2(as.matrix(X),as.matrix(candidati2),lambda = 0)
    mat_dist_1_2 <- Pref_dist2(as.matrix(X),as.matrix(candidati2),lambda = 0.1)
    mat_dist_2_2 <- Pref_dist2(as.matrix(X),as.matrix(candidati2),lambda = 0.2)
    mat_dist_3_2 <- Pref_dist2(as.matrix(X),as.matrix(candidati2),lambda = 0.3)
    mat_dist_4_2 <- Pref_dist2(as.matrix(X),as.matrix(candidati2),lambda = 0.4)
    mat_dist_5_2 <- Pref_dist2(as.matrix(X),as.matrix(candidati2),lambda = 0.5)
    mat_dist_6_2 <- Pref_dist2(as.matrix(X),as.matrix(candidati2),lambda = 0.6)
    mat_dist_7_2 <- Pref_dist2(as.matrix(X),as.matrix(candidati2),lambda = 0.7)
    mat_dist_8_2 <- Pref_dist2(as.matrix(X),as.matrix(candidati2),lambda = 0.8)
    mat_dist_9_2 <- Pref_dist2(as.matrix(X),as.matrix(candidati2),lambda = 0.9)
    mat_dist_10_2 <- Pref_dist2(as.matrix(X),as.matrix(candidati2),lambda = 1)
    
    D_lambda_2 = c(min(apply(mat_dist_0_2,2,mean)),
                   min(apply(mat_dist_1_2,2,mean)),
                   min(apply(mat_dist_2_2,2,mean)),
                   min(apply(mat_dist_3_2,2,mean)),
                   min(apply(mat_dist_4_2,2,mean)),
                   min(apply(mat_dist_5_2,2,mean)),
                   min(apply(mat_dist_6_2,2,mean)),
                   min(apply(mat_dist_7_2,2,mean)),
                   min(apply(mat_dist_8_2,2,mean)),
                   min(apply(mat_dist_9_2,2,mean)),
                   min(apply(mat_dist_10_2,2,mean)))
    
    
    D_lambda_final <- apply(cbind(D_lambda_1,D_lambda_2),1,min)
   
    
    
  }
  
  return(list(D_lambda=D_lambda_final))
}


