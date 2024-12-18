################################################################################
source("Final_Consensus_SEARCH_updated.r")
#Iterazioni 


library(SimDesign)
library(ggplot2)
library(tidyr)
library(dplyr)

###############################################
# SETTING 1
# LAMBDA IMPACT ON THE CONSENSUS

################################################################

# 10 items 10 judges

# Iterations number
n_iter <- 500

nitems <- 10
njudges <- 10

# Matrices to store the results
SCR_results_10_30 <- matrix(NA, nrow = n_iter, ncol = 9)
SMA_results_10_30 <- matrix(NA, nrow = n_iter, ncol = 9)



Sigma=matrix(0,nrow=nitems,ncol=nitems)
diag(Sigma)=9.8


for (i in 1:n_iter) {
  # Data simulation
  set.seed(600+i)
  mu <- rnorm(nitems, 10, 1)
  true_conf <- MASS::mvrnorm(n = njudges, mu, Sigma = Sigma)
  rank <- t(apply(true_conf, 1, rank, ties.method = "min"))
  
  
  approv <- rpois(njudges,5)
  ties <- sapply(approv, function(x) sample(0:x, 1))
  approv_m <- matrix(nrow = nrow(rank), ncol = ncol(rank))
  
  
  for (j in 1:nrow(rank)) {
      rank[j, ] <- ifelse(rank[j, ] <= ties[j], min(rank[j, rank[j, ] <= ties[j]]), rank[j, ])
      approv_m[j, ] <- ifelse(rank[j, ] <= approv[j], 1, 0)
   
  }
  
  X <- data.frame(cbind(rank, approv_m))
  
  M <- list()
  
  for (lambda_index in 1:9) {
    lambda <- lambda_index / 10
    M[[lambda_index]] <- quiet(Final_Consensus(X,
                                               algorithm = "quick", 
                                               lambda = lambda))$Consensus
  }
  
  CR <- quiet(consrank(X[, 1:nitems], algorithm = "quick")$Consensus)
  CA <- ifelse(apply(X[, (nitems + 1):(2 * nitems)], 2, mean) > 0.5, 1, 0)
  
  res = res2 =  rep(NA,9)
  for (j in 1:9) {
    res[j] <- min(kemenyd(as.matrix(CR), as.matrix(M[[j]][, 1:nitems])))/ (nitems * (nitems - 1))
    res2[j] <- min(binary_dist(CA, as.matrix(M[[j]][, (nitems + 1):(2 * nitems)])))/ nitems
  }
  
  
  SCR_results_10_30[i, ] <- res
  SMA_results_10_30[i, ] <- res2
  
  
  if(!all(res == sort(res,decreasing = T))){
    print("STOP!")
    break
  }
  
  print(i)
}

# Computing the average of the obtained results 
SCR_mean <- colMeans(SCR_results_10_30, na.rm = TRUE)
SMA_mean <- colMeans(SMA_results_10_30, na.rm = TRUE)


plot(1:9/10,SCR_mean, type="o")



data_10_30 <- data.frame(
  mean = c(
    apply(as.data.frame(SCR_results_10_30), 2, mean, na.rm = TRUE),
    apply(as.data.frame(SMA_results_10_30), 2, mean, na.rm = TRUE)
  ),
  
  lower = c(
    apply(as.data.frame(SCR_results_10_30), 2, function(x) quantile(x, probs = 0.025, na.rm = TRUE)),
    apply(as.data.frame(SMA_results_10_30), 2, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
  ),
  
  upper = c(
    apply(as.data.frame(SCR_results_10_30), 2, function(x) quantile(x, probs = 0.975, na.rm = TRUE)),
    apply(as.data.frame(SMA_results_10_30), 2, function(x) quantile(x, probs = 0.975, na.rm = TRUE))
  ),
  
  Type = rep(c("SCR", "SMA"), each = 9),
  lambda = rep(seq(0.1, 0.9, by = 0.1), 2)
)

# 20 items 20 giudici

# Numero di iterazioni
n_iter <- 500

nitems <- 20
njudges <- 20

# Matrici per raccogliere i risultati
SCR_results_20_30 <- matrix(NA, nrow = n_iter, ncol = 9)
SMA_results_20_30 <- matrix(NA, nrow = n_iter, ncol = 9)


Sigma=matrix(0,nrow=nitems,ncol=nitems)
diag(Sigma)=10


for (i in 1:n_iter) {
  # Simulazione dei dati
  set.seed(600+i)
  mu <- rnorm(nitems, 10, 1)
  true_conf <- MASS::mvrnorm(n = njudges, mu, Sigma = Sigma)
  rank <- t(apply(true_conf, 1, rank, ties.method = "min"))
  
  approv <- rpois(njudges,10)
  ties <- sapply(approv, function(x) sample(0:x, 1))
  approv_m <- matrix(nrow = nrow(rank), ncol = ncol(rank))
  
  
  for (j in 1:nrow(rank)) {
    rank[j, ] <- ifelse(rank[j, ] <= ties[j], min(rank[j, rank[j, ] <= ties[j]]), rank[j, ])
    approv_m[j, ] <- ifelse(rank[j, ] <= approv[j], 1, 0)

  }
  
  X <- data.frame(cbind(rank, approv_m))
  
  M <- list()
  
  for (lambda_index in 1:9) {
    lambda <- lambda_index / 10
    M[[lambda_index]] <- quiet(Final_Consensus(X, algorithm = "quick", lambda = lambda))$Consensus
  }
  
  CR <- quiet(consrank(X[, 1:nitems], algorithm = "quick")$Consensus)
  CA <- ifelse(apply(X[, (nitems + 1):(2 * nitems)], 2, mean) > 0.5, 1, 0)
  
  res = res2 =  rep(NA,9)
  for (j in 1:9) {
    res[j] <- min(kemenyd(as.matrix(CR), as.matrix(M[[j]][, 1:nitems])))/ (nitems * (nitems - 1))
    res2[j] <- min(binary_dist(CA, as.matrix(M[[j]][, (nitems + 1):(2 * nitems)])))/ nitems
  }
  
  
  SCR_results_20_30[i, ] <- res
  SMA_results_20_30[i, ] <- res2
  
  
  if(!all(res == sort(res,decreasing = T))){
    print("stoppami")
    break
  }
  
  print(i)
}

# Calcolo della media dei risultati
SCR_mean <- colMeans(SCR_results_20_30, na.rm = TRUE)
SMA_mean <- colMeans(SMA_results_20_30, na.rm = TRUE)

plot(1:9/10,SCR_mean, type="o")


data_20_30 <- data.frame(
  mean = c(
    apply(as.data.frame(SCR_results_20_30), 2, mean, na.rm = TRUE),
    apply(as.data.frame(SMA_results_20_30), 2, mean, na.rm = TRUE)
  ),
  
  lower = c(
    apply(as.data.frame(SCR_results_20_30), 2, function(x) quantile(x, probs = 0.025, na.rm = TRUE)),
    apply(as.data.frame(SMA_results_20_30), 2, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
  ),
  
  upper = c(
    apply(as.data.frame(SCR_results_20_30), 2, function(x) quantile(x, probs = 0.975, na.rm = TRUE)),
    apply(as.data.frame(SMA_results_20_30), 2, function(x) quantile(x, probs = 0.975, na.rm = TRUE))
  ),
  
  Type = rep(c("SCR", "SMA"), each = 9),
  lambda = rep(seq(0.1, 0.9, by = 0.1), 2)
)

################################################################

# 50 items 50 giudici

# Numero di iterazioni
n_iter <- 500

nitems <- 50
njudges <- 50

# Matrici per raccogliere i risultati
SCR_results_50_50 <- matrix(NA, nrow = n_iter, ncol = 9)
SMA_results_50_50 <- matrix(NA, nrow = n_iter, ncol = 9)


Sigma=matrix(0,nrow=nitems,ncol=nitems)
diag(Sigma)=10


for (i in 1:n_iter) {
  # Simulazione dei dati
  set.seed(600+i)
  mu <- rnorm(nitems, 10, 1)
  true_conf <- MASS::mvrnorm(n = njudges, mu, Sigma = Sigma)
  rank <- t(apply(true_conf, 1, rank, ties.method = "min"))
  
  approv <- rpois(njudges,25)
  ties <- sapply(approv, function(x) sample(0:x, 1))
  approv_m <- matrix(nrow = nrow(rank), ncol = ncol(rank))
  
  
  for (j in 1:nrow(rank)) {
    rank[j, ] <- ifelse(rank[j, ] <= ties[j], min(rank[j, rank[j, ] <= ties[j]]), rank[j, ])
    approv_m[j, ] <- ifelse(rank[j, ] <= approv[j], 1, 0)
    
  }
  
  X <- data.frame(cbind(rank, approv_m))
  
  
  M <- list()
  
  for (lambda_index in 1:9) {
    lambda <- lambda_index / 10
    M[[lambda_index]] <- quiet(Final_Consensus(X, algorithm = "quick", lambda = lambda))$Consensus
  }
  
  CR <- quiet(consrank(X[, 1:nitems], algorithm = "quick")$Consensus)
  CA <- ifelse(apply(X[, (nitems + 1):(2 * nitems)], 2, mean) > 0.5, 1, 0)
  
  res = res2 =  rep(NA,9)
  for (j in 1:9) {
    res[j] <- min(kemenyd(as.matrix(CR), as.matrix(M[[j]][, 1:nitems])))/ (nitems * (nitems - 1))
    res2[j] <- min(binary_dist(CA, as.matrix(M[[j]][, (nitems + 1):(2 * nitems)])))/ nitems
  }
  
  
  SCR_results_50_50[i, ] <- res
  SMA_results_50_50[i, ] <- res2
  
  
  if(!all(res == sort(res,decreasing = T))){
    print("stoppami")
    break
  }
  
  print(i)
}

# Calcolo della media dei risultati
SCR_mean_50_50 <- colMeans(SCR_results_50_50, na.rm = TRUE)
SMA_mean_50_50 <- colMeans(SMA_results_50_50, na.rm = TRUE)


plot(1:9/10,SCR_mean_50_50, type="o")


data_50_50 <- data.frame(
  mean = c(
    apply(as.data.frame(SCR_results_50_50), 2, mean, na.rm = TRUE),
    apply(as.data.frame(SMA_results_50_50), 2, mean, na.rm = TRUE)
  ),
  
  lower = c(
    apply(as.data.frame(SCR_results_50_50), 2, function(x) quantile(x, probs = 0.025, na.rm = TRUE)),
    apply(as.data.frame(SMA_results_50_50), 2, function(x) quantile(x, probs = 0.025, na.rm = TRUE))
  ),
  
  upper = c(
    apply(as.data.frame(SCR_results_50_50), 2, function(x) quantile(x, probs = 0.975, na.rm = TRUE)),
    apply(as.data.frame(SMA_results_50_50), 2, function(x) quantile(x, probs = 0.975, na.rm = TRUE))
  ),
  
  Type = rep(c("SCR", "SMA"), each = 9),
  lambda = rep(seq(0.1, 0.9, by = 0.1), 2)
)

#load("sensitivity1_500ite.Rdata")

ALL <-  bind_rows(
data_10_30 |> 
  mutate(setting = "10"),
 data_20_30 |> 
   mutate(setting = "20"),
data_50_50 |> 
  mutate(setting="50")

)

library(ggplot2)
library(tidyverse)

# Save plot to PDF
#pdf("Setting1_ties_500.pdf", width = 7, height = 4.5)
ALL %>%
  ggplot(aes(x = lambda, y = mean, shape = setting, color = setting)) +
  geom_point(size = 2) +  # Larger points with slight transparency
  geom_line(size = 1, linetype = "solid", alpha = 0.8) +
  geom_vline(xintercept = 0.5, col = "orange", linetype = "dashed") +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  labs(
    x = expression(lambda),
    y = "Average",
    shape = expression(n),
    color = expression(n)
  ) +
  scale_x_continuous(
    breaks = c(0.1, 0.3, 0.5, 0.7, 0.9),
    labels = c(0.1, 0.3, 0.5, 0.7, 0.9)
  ) +
  facet_wrap(~Type) +
  ylim(0, 0.2) +
  theme_minimal(base_size = 12) +
  theme(
    legend.text = element_text(size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    panel.spacing = unit(2, "lines"),
    strip.text = element_text(face = "bold", size = 12)
  )
#dev.off()




##################################################
# IMPACT OF LAMBDA ON THE GOODNESS OF THE SOLUTION
##################################################
  
  
########################################################################  
  # Setting 1
  # Ranking: low noise
  # approval: low noise
######################################

 # 1.1 item 10
  
  
  nitems= 10
  njudges = 10
  Sigma=matrix(0,nrow=nitems,ncol=nitems)
  diag(Sigma)=1
  n_iter <- 100
  
  consensus_result_1_10_iter <- matrix(NA, nrow=n_iter, ncol= 11)
  
  for (i in 1:n_iter) {
    set.seed(i)
    # Simulazione dei dati
    mu <- rnorm(nitems, 10, 5.5)
    true_conf <- MASS::mvrnorm(n = njudges, mu, Sigma = Sigma)
    rank <- t(apply(true_conf, 1, rank, ties.method = "min"))


    approv <- rpois(njudges,0.6)
    approv_m <- matrix(nrow = nrow(rank), ncol = ncol(rank))
    
    for (j in 1:nrow(rank)) {
      approv_m[j, ] <- ifelse(rank[j, ] <= approv[j], 1, 0)
    }
    
    X <- data.frame(cbind(rank, approv_m))
    
    
    
    
    rd <- mean(kemenyd(rank)) / 90
    ad <- mean(as.dist(binary_dist(approv_m))) / 10
    condition <- between(rd,0.07,0.09) * between(ad,0.07,0.09)
    

    
    if (condition ==0){
      print("NO SOLUTION")
    }else{
      consensus_result_1_10_iter[i,] <- quiet(Consensus_sensibility(X)$D_lambda)
    }
    
    print(i)
  }
  
  
  ###############################################################  

  # 1.2 item 20
  nitems= 20
  njudges = 20
  
  Sigma=matrix(0,nrow=nitems,ncol=nitems)
  diag(Sigma)=1
  
  
  n_iter <- 100
  
  consensus_result_1_20_iter <- matrix(NA, nrow=n_iter, ncol= 11)
  for (i in 1:n_iter) {
  set.seed(i)
  # Simulazione dei dati
  mu <- rnorm(nitems, 10, 5.5)
  true_conf <- MASS::mvrnorm(n = njudges, mu, Sigma = Sigma)
  rank <- t(apply(true_conf, 1, rank, ties.method = "min"))
  
  
  
  approv <- rpois(njudges,1.8)
  approv_m <- matrix(nrow = nrow(rank), ncol = ncol(rank))
  
  for (j in 1:nrow(rank)) {
    approv_m[j, ] <- ifelse(rank[j, ] <= approv[j], 1, 0)
  }
  
  X <- data.frame(cbind(rank, approv_m))
  
  
  rd <- mean(kemenyd(rank)) / 380
  ad <- mean(as.dist(binary_dist(approv_m))) / 20
  condition <- between(rd,0.07,0.09) * between(ad,0.07,0.09)
  
  if (condition ==0){
    print("NO SOLUTION")
  }else{
    consensus_result_1_20_iter[i,] <- quiet(Consensus_sensibility(X)$D_lambda)
  }

print(i)
}

###############################################################  

# 1.3 item 50
nitems = 50
njudges = 50
  

n_iter <- 100
Sigma=matrix(0,nrow=nitems,ncol=nitems)
diag(Sigma)=1

consensus_result_1_50_iter <- matrix(NA, nrow=n_iter, ncol= 11)

for (i in 1:n_iter) {
  set.seed(i)
  # Simulazione dei dati
  mu <- rnorm(nitems, 10, 5.5)
  true_conf <- MASS::mvrnorm(n = njudges, mu, Sigma = Sigma)
  rank <- t(apply(true_conf, 1, rank, ties.method = "min"))

  
  approv <- rpois(njudges,8.5)
  approv_m <- matrix(nrow = nrow(rank), ncol = ncol(rank))
  
  for (j in 1:nrow(rank)) {
    approv_m[j, ] <- ifelse(rank[j, ] <= approv[j], 1, 0)
  }
  
  X <- data.frame(cbind(rank, approv_m))
  
  (rd <- mean(kemenyd(rank)) / 2450)
  
  (ad <- mean(as.dist(binary_dist(approv_m))) / 50)
  
  condition <- between(rd,0.07,0.09) * between(ad,0.07,0.09)
  
  if (condition ==0){
    print("NO SOLUTION")
  }else{
    consensus_result_1_50_iter[i,] <- quiet(Consensus_sensibility(X)$D_lambda)
  }
  

print(i)

}


###############################################################    
########################################################################  
# Setting 4
# Ranking: moderate noise
# approval: low noise
######################################
  
  
# 4.1 10 items

nitems = 10
njudges = 10 

n_iter <- 300

consensus_result_2_10_iter <- matrix(NA, nrow=n_iter, ncol= 11)


Sigma=matrix(0,nrow=nitems,ncol=nitems)
diag(Sigma)=1

for (i in 1:n_iter) {
set.seed(i)
# Simulazione dei dati
mu <- rnorm(nitems, 10, 1)
true_conf <- MASS::mvrnorm(n = njudges, mu, Sigma = Sigma)
rank <- t(apply(true_conf, 1, rank, ties.method = "min"))
  
  
  
approv <- rbinom(njudges,nitems,0.975)
approv_m <- matrix(nrow = nrow(rank), ncol = ncol(rank))
  
  for (j in 1:nrow(rank)) {
    approv_m[j, ] <- ifelse(rank[j, ] <= approv[j], 1, 0)
  }
  
X <- data.frame(cbind(rank, approv_m))

(rd <- mean(kemenyd(rank)) / 90)

(ad <- mean(as.dist(binary_dist(approv_m))) / 10)

condition <- between(rd,0.325,0.345) * between(ad,0.025,0.055)

if (condition ==0){
  print("NO SOLUTION")
}else{
  consensus_result_2_10_iter[i,] <- quiet(Consensus_sensibility(X)$D_lambda)
}
  
print(i)
}

###############################################################    

# 4.2 20 items

nitems=20
njudges = 20 

n_iter <- 200

consensus_result_2_20_iter <- matrix(NA, nrow=n_iter, ncol= 11)

Sigma=matrix(0,nrow=nitems,ncol=nitems)
diag(Sigma)=1

for (i  in 1:n_iter) {
# Simulazione dei dati
  set.seed(i)
mu <- rnorm(nitems, 10, 1)
true_conf <- MASS::mvrnorm(n = njudges, mu, Sigma = Sigma)
rank <- t(apply(true_conf, 1, rank, ties.method = "min"))




approv <- rbinom(njudges,nitems,0.975)
approv_m <- matrix(nrow = nrow(rank), ncol = ncol(rank))

for (j in 1:nrow(rank)) {
  approv_m[j, ] <- ifelse(rank[j, ] <= approv[j], 1, 0)
}

X <- data.frame(cbind(rank, approv_m))

(rd <- mean(kemenyd(rank)) / 380)

(ad <- mean(as.dist(binary_dist(approv_m))) / 20)


condition <- between(rd,0.325,0.345) * between(ad,0.025,0.055)

if (condition ==0){
  print("NO SOLUTION")
}else{
  consensus_result_2_20_iter[i,] <- quiet(Consensus_sensibility(X)$D_lambda)
}
print(i)
}


############################################################### 

# 4.3 50 items

nitems = 50
njudges = 50
n_iter = 100 

Sigma=matrix(0,nrow=nitems,ncol=nitems)
diag(Sigma)=1
# Simulazione dei dati
consensus_result_2_50_iter <- matrix(NA, nrow=n_iter, ncol= 11)

for (i in 1:n_iter) {

mu <- rnorm(nitems, 10, 1)
true_conf <- MASS::mvrnorm(n = njudges, mu, Sigma = Sigma)
rank <- t(apply(true_conf, 1, rank, ties.method = "min"))


approv <- rbinom(50,nitems,0.975)
approv_m <- matrix(nrow = nrow(rank), ncol = ncol(rank))

for (j in 1:nrow(rank)) {
  approv_m[j, ] <- ifelse(rank[j, ] <= approv[j], 1, 0)
}

X <- data.frame(cbind(rank, approv_m))

(rd <- mean(kemenyd(rank)) / 2450)

(ad <- mean(as.dist(binary_dist(approv_m))) / 50)

condition <- between(rd,0.325,0.345) * between(ad,0.03,0.05)

if (condition ==0){
  print("NO SOLUTION")
}else{
  consensus_result_2_50_iter[i,] <- quiet(Consensus_sensibility(X)$D_lambda)
}
print(i)
}

############################################################### 

########################################################################  
# Setting 3
# Ranking: low noise
# approval: moderate noise
######################################  

# 3.1 10 item

nitems= 10
njudges = 10
n_iter = 300

Sigma=matrix(0,nrow=nitems,ncol=nitems)
diag(Sigma)=1
  
  
  consensus_result_3_10_iter <- matrix(NA, nrow=n_iter, ncol= 11)
  for (i in 1:n_iter) {
set.seed(i)
      # Simulazione dei dati
  mu <- rnorm(nitems, 10, 10.5)
  true_conf <- MASS::mvrnorm(n = njudges, mu, Sigma = Sigma)
  rank <- t(apply(true_conf, 1, rank, ties.method = "min"))
  
  
  
  approv <- round(runif(njudges, min = 0, max = nitems), 0)
  approv_m <- matrix(nrow = nrow(rank), ncol = ncol(rank))
  
  for (j in 1:nrow(rank)) {
    approv_m[j, ] <- ifelse(rank[j, ] <= approv[j], 1, 0)
  }
  
  X <- data.frame(cbind(rank, approv_m))
  
  (rd <- mean(kemenyd(rank)) / 90)
  
  (ad <- mean(as.dist(binary_dist(approv_m))) / 10)
  
  condition <- between(ad,0.325,0.345) * between(rd,0.03,0.05)
  
  if (condition ==0){
    print("NO SOLUTION")
  }else{
    consensus_result_3_10_iter[i,] <- quiet(Consensus_sensibility(X)$D_lambda)
  }
  
  print(i)
  }
  
  ###############################################################   
  
  # 3.2 20 item
  
  nitems= 20
  njudges = 20
  n_iter = 100
  
  Sigma=matrix(0,nrow=nitems,ncol=nitems)
  diag(Sigma)=1
  
  
  consensus_result_3_20_iter <- matrix(NA, nrow=n_iter, ncol= 11)
 
   for (i in 1:n_iter) {
  set.seed(i)
  Sigma=matrix(0,nrow=nitems,ncol=nitems)
  diag(Sigma)=1
  # Simulazione dei dati
  mu <- rnorm(nitems, 10, 10.5)
  true_conf <- MASS::mvrnorm(n = njudges, mu, Sigma = Sigma)
  rank <- t(apply(true_conf, 1, rank, ties.method = "min"))
  
  
  approv <- round(runif(njudges, min = 0, max = nitems), 0)
  approv_m <- matrix(nrow = nrow(rank), ncol = ncol(rank))
  
  for (j in 1:nrow(rank)) {
    approv_m[j, ] <- ifelse(rank[j, ] <= approv[j], 1, 0)
  }
  
  X <- data.frame(cbind(rank, approv_m))
  
  (rd <- mean(kemenyd(rank)) / 380)
  
  (ad <- mean(as.dist(binary_dist(approv_m))) / 20)
  
  condition <- between(ad,0.325,0.345) * between(rd,0.03,0.05)
  

  if (condition ==0){
    print("NO SOLUTION")
  }else{
    consensus_result_3_20_iter[i,] <- quiet(Consensus_sensibility(X)$D_lambda)
  }
  print(i)
   }
  
################################################################
  
  # 3.3 50 item
  
  nitems= 50
  njudges = 50
  n_iter= 100
  Sigma=matrix(0,nrow=nitems,ncol=nitems)
  diag(Sigma)=1
  # Simulazione dei dati
  consensus_result_3_50_iter <- matrix(NA, nrow=n_iter, ncol= 11)
  for (i in 1:n_iter) {
  mu <- rnorm(nitems, 10, 12)
  true_conf <- MASS::mvrnorm(n = njudges, mu, Sigma = Sigma)
  rank <- t(apply(true_conf, 1, rank, ties.method = "min"))
  
  
  approv <- round(runif(njudges, min = 0, max = nitems), 0)
  approv_m <- matrix(nrow = nrow(rank), ncol = ncol(rank))
  
  
  
  for (j in 1:nrow(rank)) {
    approv_m[j, ] <- ifelse(rank[j, ] <= approv[j], 1, 0)
  }

    X <- data.frame(cbind(rank, approv_m))

    (rd <- mean(kemenyd(rank)) / 2450)
    
    (ad <- mean(as.dist(binary_dist(approv_m))) / 50)
    
    condition <- between(ad,0.325,0.345) * between(rd,0.03,0.05)

    if (condition ==0){
      print("NO SOLUTION")
    }else{
      consensus_result_3_50_iter[i,] <- quiet(Consensus_sensibility(X)$D_lambda)
    }

 print(i)
  
  }

  ################################################################

  
  ########################################################################  
  # Setting 2
  # Ranking: moderate noise
  # approval: moderate noise
  ######################################  
  
  
  # 2.1 10 item
  nitems = 10
  njudges = 10
  n_iter = 300
  Sigma=matrix(0,nrow=nitems,ncol=nitems)
  diag(Sigma)=1
  
  consensus_result_4_10_iter <- matrix(NA, nrow=n_iter, ncol= 11) 
  for (i in 1:n_iter) {
    set.seed(1995+i)
  # Simulazione dei dati
  mu <- rnorm(nitems, 10, 1)
  true_conf <- MASS::mvrnorm(n = njudges, mu, Sigma = Sigma)
  rank <- t(apply(true_conf, 1, rank, ties.method = "min"))
  
  
  
  approv <- rpois(njudges,3.1)
  approv_m <- matrix(nrow = nrow(rank), ncol = ncol(rank))
  
  for (j in 1:nrow(rank)) {
    approv_m[j, ] <- ifelse(rank[j, ] <= approv[j], 1, 0)
  }
  
  X <- data.frame(cbind(rank, approv_m))
  (rd <- mean(kemenyd(rank)) / 90)
  
  (ad <- mean(as.dist(binary_dist(approv_m))) / 10)
  
  (condition <- between(ad,0.325,0.345) * between(rd,0.324,0.345))

  
  if (condition ==0){
    print("NO SOLUTION")
  }else{
  
    consensus_result_4_10_iter[i,] <- quiet(Consensus_sensibility(X)$D_lambda)
  }
  
  
  print(i)
  }

  
  ################################################################
  

  
  # 2.2 20 item
  nitems=20
  njudges = 20
  n_iter = 200
  
  Sigma=matrix(0,nrow=nitems,ncol=nitems)
  diag(Sigma)=1
  # Simulazione dei dati
  
  
  consensus_result_4_20_iter <- matrix(NA, nrow=n_iter, ncol= 11) 
  
  for (i in 1:n_iter) {
  set.seed(i)
  mu <- rnorm(nitems, 10, 1)
  true_conf <- MASS::mvrnorm(n = njudges, mu, Sigma = Sigma)
  rank <- t(apply(true_conf, 1, rank, ties.method = "min"))
  
  
  
  approv <- rpois(njudges,7)
  approv_m <- matrix(nrow = nrow(rank), ncol = ncol(rank))
  
  for (j in 1:nrow(rank)) {
    approv_m[j, ] <- ifelse(rank[j, ] <= approv[j], 1, 0)
  }
  
  X <- data.frame(cbind(rank, approv_m))
  
  (rd <- mean(kemenyd(rank)) / 380)
  
  (ad <- mean(as.dist(binary_dist(approv_m))) / 20)
  
  (condition <- between(ad,0.325,0.345) * between(rd,0.324,0.345))
  
 
  
  
  if (condition ==0){
    print("NO SOLUTION")
  }else{
    consensus_result_4_20_iter[i,] <- quiet(Consensus_sensibility(X)$D_lambda)
  }
  print(i)
  }
  
  
  ################################################################

  # 2.3 50 item
  nitems=50
  njudges = 50
  n_iter = 100
  
  Sigma=matrix(0,nrow=nitems,ncol=nitems)
  diag(Sigma)=1
  # Simulazione dei dati
  
  consensus_result_4_50_iter <- matrix(NA, nrow=n_iter, ncol= 11) 
  
  for (i in 1:n_iter) {
  set.seed(i)
  mu <- rnorm(nitems, 10, 1)
  true_conf <- MASS::mvrnorm(n = njudges, mu, Sigma = Sigma)
  rank <- t(apply(true_conf, 1, rank, ties.method = "min"))
  
  
  
  approv <- rpois(njudges,18)
  approv_m <- matrix(nrow = nrow(rank), ncol = ncol(rank))
  
  for (j in 1:nrow(rank)) {
    approv_m[j, ] <- ifelse(rank[j, ] <= approv[j], 1, 0)
  }
  
  X <- data.frame(cbind(rank, approv_m))
  
  (rd <- mean(kemenyd(rank)) / 2450)
  
  (ad <- mean(as.dist(binary_dist(approv_m))) / 50)
  
  (condition <- between(ad,0.325,0.345) * between(rd,0.324,0.345))
  
  
  if (condition ==0){
    print("NO SOLUTION")
  }else{
    consensus_result_4_50_iter[i,] <- quiet(Consensus_sensibility(X)$D_lambda)
  }

  print(i)
  
  }
  
  ################################################################ 
  
  #load(file="setting2_sensitivity_iter.rdata")
  
  #############################
  # Plotting iter 
  library(tibble)
  library(dplyr)
  library(ggplot2)
  
  combined_data <- bind_rows(
    # Setting S1
    tibble(Badness_of_fit = colMeans(consensus_result_1_10_iter,na.rm = T), setting = "S1",
           nitems = 10, 
           lambda = 0:10 / 10),
    tibble(Badness_of_fit = colMeans(consensus_result_1_20_iter,na.rm = T), 
           setting = "S1",
           nitems = 20,
           lambda = 0:10 / 10),
    tibble(Badness_of_fit = colMeans(consensus_result_1_50_iter,na.rm = T),
           setting = "S1",
           nitems = 50,
           lambda = 0:10 / 10),#,
    
    # Setting S2
     tibble(Badness_of_fit = colMeans(consensus_result_4_10_iter,na.rm = T), setting = "S2", nitems = 10, lambda = 0:10 / 10),
     tibble(Badness_of_fit = colMeans(consensus_result_4_20_iter,na.rm = T), setting = "S2", nitems = 20, lambda = 0:10 / 10),
     tibble(Badness_of_fit = colMeans(consensus_result_4_50_iter,na.rm = T), setting = "S2", nitems = 50, lambda = 0:10 / 10),
    # 
    # # Setting S3
     tibble(Badness_of_fit = colMeans(consensus_result_3_10_iter,na.rm = T), setting = "S3", nitems = 10, lambda = 0:10 / 10),
     tibble(Badness_of_fit = colMeans(consensus_result_3_20_iter,na.rm = T), setting = "S3", nitems = 20, lambda = 0:10 / 10),
     tibble(Badness_of_fit = colMeans(consensus_result_3_50_iter,na.rm = T), setting = "S3", nitems = 50, lambda = 0:10 / 10),
    # 
    # # Setting S4
     tibble(Badness_of_fit = colMeans(consensus_result_2_10_iter,na.rm = T), setting = "S4", nitems = 10, lambda = 0:10 / 10),
     tibble(Badness_of_fit = colMeans(consensus_result_2_20_iter,na.rm = T), setting = "S4", nitems = 20, lambda = 0:10 / 10),
     tibble(Badness_of_fit = colMeans(consensus_result_2_50_iter,na.rm = T), setting = "S4", nitems = 50, lambda = 0:10 / 10)
  )
  
  y_min <- min(combined_data$Badness_of_fit) * 0.9
  y_max <- max(combined_data$Badness_of_fit) * 1.1
  
  library(tibble)
  library(dplyr)
  library(ggplot2)
  
  setting_labels <- c("S1" = "Setting 1", 
                      "S2" = "Setting 2", 
                      "S3" = "Setting 3", 
                      "S4" = "Setting 4")
  
  
  
  #pdf("Setting2_new.pdf",width = 6, height=4.5)

  ggplot(combined_data, aes(x = lambda, y = Badness_of_fit, color = factor(nitems), shape = factor(nitems), group = nitems)) +
    geom_point(size = 2) +  
    geom_line(size = 1, linetype = "solid", alpha = 0.8) +
    scale_y_continuous(limits = c(y_min, y_max)) +
    scale_x_continuous(
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c(0, 0.25, 0.5, 0.75, 1)
    )  +
    scale_color_brewer(palette = "Set1") +
    scale_fill_brewer(palette = "Set1") +
    labs(
      x = expression(lambda),  
      y = expression(bar(d)[lambda]),
      color = expression(n),
      shape = expression(n)
    ) +
    facet_wrap(~setting, ncol = 2, labeller = labeller(setting = setting_labels)) + 
    theme_minimal() +
    theme(
      legend.text = element_text(size = 10),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.y = element_text(margin = margin(r = 10)),
      panel.grid.major = element_line(color = "grey80"),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      panel.spacing = unit(0.5, "lines"),
      strip.text = element_text(face = "bold", size = 12)
    )
  #dev.off()
  
  
  

  
  
  

  
  

  
  