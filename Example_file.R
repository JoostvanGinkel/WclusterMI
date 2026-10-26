setwd("C:/directory")

#Packages
library(MGMM)
library(MASS)
library(cluster)
source("WclusterMI.r")

J=4                          # number of variables
N=2000                       # sample size
numberofclusters=4           #number of clusters
numberofgroups = 3	     #number of groups for generating MAR
numberofstarts = 10	     #number of start sets for K-means clustering
maxiterations = 1000         #maximum number of iterations


NumberOfImp=5                                            #number of imputed datasets

percent=0.3    # proportions of missing


#Creating parameters for simulation model
se=2.4
mu = list()
sigma = list() 
for (j in 1:numberofclusters) {
  sigma[[j]] = matrix(0.03,J,J)
  mu[[j]] = rep(0, times = J)
}
for (j in 1:numberofclusters) {
  mu[[j]] = 2*rep(j*se, times=J)/numberofclusters
  sigma[[j]] = 3*se*j*sigma[[j]]/numberofclusters
  diag(sigma[[j]]) = 1
}
 
###########
# FULL data
###########
ID = (numberofgroups-1)*100 + 1
Seed = ID + 329 
set.seed(Seed) 
prob <- c(1:numberofclusters)/sum(c(1:numberofclusters))
cluster = sample(1:numberofclusters,N,prob=prob,replace=TRUE)
datacom=matrix(,N,J)
for(i in 1:numberofclusters){
  datacom[cluster==i,] = mvrnorm(n=nrow(datacom[cluster==i,]),mu=rep(0,times=J),Sigma=diag(J))
  eS <- eigen(sigma[[i]], symmetric = TRUE)
  ev <- eS$values
  datacom[cluster==i,] = mu[[i]] + eS$vectors%*%diag(sqrt(pmax(ev, 0)), J)%*%t(datacom[cluster==i,]) 
}
  
#Gap statistic for full data
numberofstarts2 = numberofstarts

gap_stat = (clusGap(datacom, FUNcluster = kmeans, nstart=numberofstarts2, K.max = numberofclusters, iter.max=maxiterations, B=50, d.power = 2, verbose=FALSE))$Tab[,3]

  
#K-means clustering with optimal K      
bestK = which(gap_stat == max(gap_stat))
clustering = kmeans(datacom,centers=bestK, iter.max=maxiterations, nstart=numberofstarts)
solution = clustering$cluster   
          
         
###########
# With Missing
###########
Seed = ID + 3*212
set.seed(Seed) 

randommatrix = runif(nrow(datacom)*ncol(datacom)) 
randommatrix = matrix(randommatrix,nrow(datacom),ncol(datacom))
missingmatrix = matrix(1,nrow(datacom),ncol(datacom))
Numberofmissing = round(percent*(nrow(datacom)*ncol(datacom)))
Candidates = sample(1:numberofgroups,nrow(datacom),prob=rep(1/numberofgroups,times=numberofgroups),replace=TRUE)

Numberofmissingvec = rep(0, times=numberofgroups)  
for (v in 1:length(Numberofmissingvec)) Numberofmissingvec[v] = round(percent*(sum(Candidates==v)*ncol(datacom)))
difference = Numberofmissing - sum(Numberofmissingvec)
if (difference < 0) {
  for (v in 1:length(Numberofmissingvec)) {
    if (Numberofmissingvec[v] == max(Numberofmissingvec)) {
      Numberofmissingvec[v] = Numberofmissingvec[v] + difference
      break
    }
  }
}
if (difference > 0) {
  for (v in 1:length(Numberofmissingvec)) {
    if (Numberofmissingvec[v] == min(Numberofmissingvec)) {
      Numberofmissingvec[v] = Numberofmissingvec[v] + difference
      break
    }
  }
}

weightmatrix = matrix(, nrow(datacom), ncol(datacom))
weightmatrix[Candidates==1,] <- 1
weightmatrix[Candidates==2,] <- datacom[Candidates==2,ncol(datacom)]-min(datacom[Candidates==2,ncol(datacom)])+0.01
weightmatrix[Candidates==2,ncol(datacom)] <- 0
weightmatrix[Candidates==3,] <- datacom[Candidates==3,1]-min(datacom[Candidates==3,1])+0.01 
weightmatrix[Candidates==3,1] <- 0      


randommatrix <- randommatrix*weightmatrix
for (v in 1:numberofgroups) {
  stoploop = FALSE
  for (k in 1:Numberofmissingvec[v]){
    for (ii in 1:nrow(randommatrix[Candidates==v,])){
      for (jj in 1:ncol(randommatrix[Candidates==v,])){
        if (randommatrix[Candidates==v,][ii,jj] == max(randommatrix[Candidates==v,]) & sum(1*(randommatrix[Candidates==v,][ii,-jj] != 0))>1) {
          missingmatrix[Candidates==v,][ii,jj] = 0
          randommatrix[Candidates==v,][ii,jj] = 0 
          stoploop = TRUE
          break             
        }
      }
      if (stoploop == TRUE) {
        stoploop = FALSE
        break
      }
    }
  }          
}
datamis <- datacom
datamis[missingmatrix==0] = NA

dimnames(datamis)[2] <- list(c(paste0("x", 1:J,sep="")))


#Imputing the data
Seed = ID + 7
set.seed(Seed)

lambda = rep(1/numberofclusters, times = numberofclusters)
  
warnmessage = FALSE
impwarns = 0

fit = FitGMM(datamis, init_covs=sigma, init_means=mu, init_props=lambda, k=numberofclusters, report=FALSE)
dataimp = list()
for (m in 1:NumberOfImp) {
  dataimp[[m]] = GenImputation(fit)
}

Seed = ID*2 + 15

#Gap statistic for imputed datasets using Gap - Wbar
set.seed(Seed)        
gap_solution1 = (clusGapMI(dataimp, K.max = numberofclusters, NumberOfImp=NumberOfImp, poolmethod="Wbar", n.start=numberofstarts, n=N, B=50, d.power = 2))
bestK1 = gap_solution1$bestK

#Cluster analysis for imputed datasets with optimal K using Wbar
solution1 = Wbar(data=dataimp,numberofclusters=bestK1,numberofstarts=numberofstarts,N=N,J=J,NumberOfImp=NumberOfImp) 

#Gap statistic for imputed datasets using Gap - W*
set.seed(Seed)  
gap_solution2 = (clusGapMI(dataimp, K.max = numberofclusters, NumberOfImp=NumberOfImp, poolmethod="Wstar", n.start=numberofstarts, n=N, B=50, d.power = 2))
bestK2 = gap_solution2$bestK

#Cluster analysis for imputed datasets with optimal K using Wbar
solution2 = Wstar(data=dataimp,numberofclusters=bestK2,numberofstarts=numberofstarts,N=N,J=J,NumberOfImp=NumberOfImp)
