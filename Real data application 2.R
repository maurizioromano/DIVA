# Real data
###
#----------------------------------------
# Funzione per leggere i dati da preflib, sia soc che soi. E' necessario che siano strutturati due file .csv/testo cos? strutturati
# description esempio: 
# 1,Bielefeld University
# 2,Alexandria University
# 3,Dalhousie University
# ...
# raWdata esempio:
# 418,243,425,124,531, ...
# 518,395,16,490,101,480, ...
# il parametro soi, indica se i dati sono incompleti (default=FALSE (dunque SOC))
# i dati soi vengono trattati inserendo gli oggetti non valutati da ciascun giudice, come ex aequo all'ultimo posto
read_preflib_soc = function(file_name_description, file_name_rawData, soi=FALSE){
  require(ConsRank)
  data_raw <- read.csv(file_name_rawData, header=F)
  #la riga seguente non funziona come dovrebbe
  data_clean <- do.call(cbind, data_raw) #cosi' leva gli spazi
  #importa codifica da file excel
  varcode <- read.csv(file_name_description, header=F)
  colnames(varcode) = c("varCode","itemCode")
  idvarcode <- as.numeric(varcode$varCode) #coding originale
  itemscode <- as.matrix(varcode$itemCode) #etichette associate al coding
  #first, trasformiamo gli ordinamenti originali in rankings
  for(i in 1:nrow(data_clean)){
    data_clean[i,] = as.character(data_clean[i,])
  }
  
  if(soi){
    
    if(ncol(data_clean) < length(varcode$varCode)){
      data_clean = cbind(data_clean, matrix(NA, nrow=nrow(data_clean), ncol=length(varcode$varCode)-ncol(data_clean)))
      colnames(data_clean) = paste("V", 1:ncol(data_clean), sep="")
    }
    
    for(i in 1:nrow(data_clean)){
      to_add = as.character(varcode$varCode[!(as.numeric(varcode$varCode) %in% as.numeric(data_clean[i,]))])
      if(length(to_add) > 0){
        if(length(to_add) > 1){
          to_add[1] = paste("{",to_add[1],sep="")
          to_add[length(to_add)] = paste(to_add[length(to_add)],"}",sep="")
        }
        data_clean[i,is.na(data_clean[i,])] = to_add
      }
      
    }
  }
  
  data_rank = order2rank(data_clean)
  labsranks <- colnames(data_rank) #ora ogni colonna è un item e la matrice contene i ranks
  idlabs <- as.numeric(labsranks) #trasforma le etichette in numeri
  objects <- itemscode[idlabs] #well, order the real sushi names according to the order of
  #                            the labels in coding format (scusa, è scappato l'inglese)
  colnames(data_rank) <- objects #assegna alle colonne i nomi dei sushi
  
  
  
  return(data_rank)
}
####
# PREPARAZIONE DATI
source("Final_Consensus_SEARCH_updated.r")

#Università 2012 soc
#df_rank = read_preflib_soc(file_name_description = "00046-00000001_description.soc", file_name_rawData = "00046-00000001_rawData.soc")

#Università 2013 soi (n mooooolto più grande)
#df_rank = read_preflib_soc(file_name_description = "00046-00000002_description.soi", file_name_rawData = "00046-00000002_rawData.soi", soi=T)

#breakfast dataset (n piccolo)
#df_rank = read_preflib_soc(file_name_description = "00035-00000002_description.soc", file_name_rawData = "00035-00000002_rawData.soc")

#F1 italia 2020 dataset (n piccolo)
#df_rank = read_preflib_soc(file_name_description = "00053-00000447_description.soi", file_name_rawData = "00053-00000447_rawData.soi", soi=T)

#Formula 1 Seasons 1950
df_rank = read_preflib_soc(file_name_description = "00052-00000001_description.soi", file_name_rawData = "00052-00000001_rawData.soi", soi=T)



dim(df_rank)

df_approv = matrix(0, nrow = nrow(df_rank), ncol = ncol(df_rank))

#prova con la top10
for(i in 1:nrow(df_approv)){
  approved = rep(0, ncol(df_rank))
  approved[order(df_rank[i,])[1:5]] = 1
  df_approv[i,] = approved
}

# for(i in 1:nrow(df_approv)){
#   approved = FindApproval(df_rank[i,])
#   set.seed(i)
#   df_approv[i,] = as.numeric(approved[sample(1:nrow(approved), 1), ((ncol(approved)/2)+1):ncol(approved)])
# }

colnames(df_approv) = paste("Approve of", colnames(df_rank)) 
  
df = cbind(df_rank, df_approv)

colMeans(df)
colMeans(df_rank)

#########################################################################################
# FINE PREPARAZIONE
#########################
# cONFRONTO CONSENSUS

Fc <- Final_Consensus(df,algorithm="quick")
Fc
#Final_Consensus(X_1,algorithm="fast")
#Approved
approved = which(Fc$Consensus[(ncol(df_rank)+1):ncol(df)]==1)
approved[unlist(Fc$Consensus[approved])]

approved[order(unlist(Fc$Consensus[approved]))]
colnames(df)[approved[ order(unlist(Fc$Consensus[approved]))]]

#Non approved
non_approved = which(Fc$Consensus[(ncol(df_rank)+1):ncol(df)]==0)
table(unlist(Fc$Consensus[non_approved]))

non_approved[order(unlist(Fc$Consensus[non_approved]))]


colmeans<-apply(df,2,mean)
round(colmeans,2)

naive_r <- rank(colmeans[1:(length(colmeans)/2)])
naive_a <- ifelse(colmeans[(length(colmeans)/2+1):length(colmeans)]>0.5,1,0)

naive_consensus <- matrix(c(naive_r,naive_a),nrow=1)

#Approved
approved = which(naive_consensus[(ncol(df_rank)+1):ncol(df)]==1)
approved[order(naive_consensus[approved])]
colnames(df)[approved[ order(naive_consensus[approved])]]

#Non approved
non_approved = which(naive_consensus[(ncol(df_rank)+1):ncol(df)]==0)
table(naive_consensus[non_approved])

non_approved[order(naive_consensus[non_approved])]


# IL NOSTRO 

Fc$D_lambda

# Naive 
mean(Pref_dist2(naive_consensus,df))

controlla(upperTriangle(score_r(naive_r)), upperTriangle(score_a(naive_a)))

#########################################
# visualizzare il pref-appro dei candidati

#Final_Consensus(X_1,algorithm="quick")

#which(!TF)

# tibble(Politician=n,
# Rating = as.numeric(df_rank[10,]),
# Approval = as.numeric(df_approv[10,])
# )


# nome candidati
colnames(df)[which(Fc$Consensus==1)[which(Fc$Consensus==1)<=(ncol(df)/2)]]
# [1] "Lionel Jospin"
colnames(df)[which(Fc$Consensus==2)[which(Fc$Consensus==2)<=(ncol(df)/2)]]
# [1] "Francois Bayrou"         "Jean-Pierre Chevenement"
# [3] "Jacques Chirac"          "Noel Mamere"            
colnames(df)[which(Fc$Consensus==3)[which(Fc$Consensus==3)<=(ncol(df)/2)]]
 # [1] "Olivier Besancenot" "Christine Boutin"   "Jacques Cheminade" 
 # [4] "Robert Hue"         "Arlette Laguiller"  "Brice Lalonde"     
 # [7] "Corine Lepage"      "Jean-Marie Le Pen"  "Alain Madelin"     
 # [10] "Bruno Maigret"


# identificativo

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


# identificativo consensus

n2[order(naive_r)]

library(ConsRank)
citation("ConsRank")
