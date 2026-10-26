# WclusterMI
Calculates a pooled cluster solution from a cluster analysis on multiply imputed datasets

Wbar = function(data,
                numberofclusters,
                N,
                J,
                NumberOfImp,
                maxiterations = 100,
                numberofstarts = 20,
                Prob = rep(1/numberofclusters, times = numberofclusters)))

Wstar = function(data,
                 numberofclusters,
                 N,
                 J,
                 NumberOfImp,
                 maxiterations = 100,
                 numberofstarts = 20,
                 Prob =  rep(1/numberofclusters, times = numberofclusters)))

Input:
- data: List object of size NumberOfImp where each element contains an imputed dataset of size N x J.
- numberofclusters: Number of presumed clusters, K.
- N: Sample size of each imputed dataset.
- J: Number of variables in each imputed dataset.
- NumberOfImp: Number of imputed datasets.
- maxiterations: Maximum number of iterations of the K-means clustering algorithm. Default: 100
- numberofstarts: Number of starting values for the K-means clustering algorithm.
- Prob: vector of cluster proportions. Default: equal proportions for all clusters.

clusGapMI <- function (x, 
                       K.max, 
                       NumberOfImp, 
                       poolmethod, 
                       n.start=20, 
                       B = 100, 
                       d.power, 
                       spaceH0 = "scaledPCA", 
                       n)
Input:
- x: List object of size NumberOfImp where each element contains an imputed dataset of size N x J.
- K.max: maximum number of clusters.
- NumberOfImp: Number of imputed datasets.
- poolmethod: method for obtaining the pooled cluster solution. Can either be "Wbar" or "Wstar".
- n.start: Number of starting values for the K-means clustering algorithm.
- B: Number of bootstrap samples for the simulating the reference distribution for the Gap statistic. Default: 100.
- d.power: a positive integer specifying the power which is applied to the euclidean distances (dist) before they are summed up to give W(k). Default: "d.power = 1", 
- spaceH0: Method to generate reference data ("scaledPCA" or "original").
- n: Sample size of each imputed dataset.
