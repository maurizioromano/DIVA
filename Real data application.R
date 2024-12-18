# Real data
###
# DATA PREPROCESSING
source("Final_Consensus_SEARCH_updated.r")
df <- read.table("00029-00000001.dat",sep=",")
names <- c("Approve of Francois Bayrou",
           "Rating of Francois Bayrou",
           "Approve of Olivier Besancenot",
           "Rating of Olivier Besancenot",
           "Approve of Christine Boutin",
           "Rating of Christine Boutin",
           "Approve of Jacques Cheminade",
           "Rating of Jacques Cheminade",
           "Approve of Jean-Pierre Chevenement",
           "Rating of Jean-Pierre Chevenement",
           "Approve of Jacques Chirac",
           "Rating of Jacques Chirac",
           "Approve of Robert Hue",
           "Rating of Robert Hue",
           "Approve of Lionel Jospin",
           "Rating of Lionel Jospin",
           "Approve of Arlette Laguiller",
           "Rating of Arlette Laguiller",
           "Approve of Brice Lalonde",
           "Rating of Brice Lalonde",
           "Approve of Corine Lepage",
           "Rating of Corine Lepage",
           "Approve of Jean-Marie Le Pen",
           "Rating of Jean-Marie Le Pen",
           "Approve of Alain Madelin",
           "Rating of Alain Madelin",
           "Approve of Noel Mamere",
           "Rating of Noel Mamere",
           "Approve of Bruno Maigret",
           "Rating of Bruno Maigret")

length(names)
dim(df)
colnames(df) <-  names

df_rank <- df[,1:15*2]
df_approv <- df[,seq(1,29,by=2)]

X <- matrix(NA, nrow=nrow(df_rank),ncol=30)
TF <- rep(NA, nrow(df_rank))

for (i in 1:nrow(df_rank)) {
  
  r <- df_rank[i,]
  a <- df_approv[i,]
  
  r[r < 0 & a==1] <-  max(-1,max(r[a==0]))+0.01
  
  X[i,] <- c(as.numeric(rank(-r,ties.method="min")),as.numeric(a))
  TF[i] <- controlla(upperTriangle(score_r(rank(-r,ties.method="min"))), upperTriangle(score_a(a)))
  
  print(i)
}

table(TF)
which(!TF)



df_rank[372,]
df_approv[372,]


colMeans(df)
colMeans(df_rank)
colMeans(X)



X_1 <- X[which(TF),]



##########################
# FINE PREPARAZIONE
#########################
# CONSENSUS COMPARISONS

Fc <- Final_Consensus(X_1,algorithm="quick")

colmeans<-apply(X_1,2,mean)
round(colmeans,2)

naive_r <- rank(colmeans[1:(length(colmeans)/2)])
naive_a <- ifelse(colmeans[(length(colmeans)/2+1):length(colmeans)]>0.5,1,0)

naive_consensus <- matrix(c(naive_r,naive_a),nrow=1)


# DIVA D_LAMBDA

Fc$D_lambda

# Naive (RPA) 
mean(Pref_dist2(naive_consensus,X_1))

#########################################
# Visualizing the pref-appro of candidates


n <- gsub("Approve of ","", names[grep("Approve of",names)])
n2 <- paste0("x_",1:16)

colnames(X) <- rep(n,2)

# Candidates name
n[which(Fc$Consensus==1)[which(Fc$Consensus==1)<16]]
# [1] "Lionel Jospin"
n[which(Fc$Consensus==2)[which(Fc$Consensus==2)<16]]
# [1] "Francois Bayrou"         "Jean-Pierre Chevenement"
# [3] "Jacques Chirac"          "Noel Mamere"            
n[which(Fc$Consensus==3)[which(Fc$Consensus==3)<16]]
 # [1] "Olivier Besancenot" "Christine Boutin"   "Jacques Cheminade" 
 # [4] "Robert Hue"         "Arlette Laguiller"  "Brice Lalonde"     
 # [7] "Corine Lepage"      "Jean-Marie Le Pen"  "Alain Madelin"     
 # [10] "Bruno Maigret"


# ID

n2[which(Fc$Consensus==1)[which(Fc$Consensus==1)<16]]
# [1] "Lionel Jospin"
n2[which(Fc$Consensus==2)[which(Fc$Consensus==2)<16]]
# [1] "Francois Bayrou"         "Jean-Pierre Chevenement"
# [3] "Jacques Chirac"          "Noel Mamere"            
n2[which(Fc$Consensus==3)[which(Fc$Consensus==3)<16]]
# [1] "Olivier Besancenot" "Christine Boutin"   "Jacques Cheminade" 
# [4] "Robert Hue"         "Arlette Laguiller"  "Brice Lalonde"     
# [7] "Corine Lepage"      "Jean-Marie Le Pen"  "Alain Madelin"     
# [10] "Bruno Maigret"


# Consensus ID

n2[order(naive_r)]
