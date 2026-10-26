#setwd("C:/Users/ginkeljrvan/OneDrive - Universiteit Leiden/#COPY_P_DRIVE_9-DEC-2024/ginkeljrvan/VanGinkelLovik3")
Wbar = function(data,
                numberofclusters,
                N,
                J,
                NumberOfImp,
                maxiterations = 100,
                numberofstarts = 20,
                Prob = rep(1/numberofclusters, times = numberofclusters)) {
          huge = 1E10
          tiny = 1E-7
          centroids = list()
          W=c()
          sclusters=c()
          dyn.load("wcalc.dll")
          count = rep(0, times = numberofstarts)
          for (s in 1:numberofstarts){
            minclustersize = 1
            while (minclustersize <= 1) {
              scluster = sample(1:numberofclusters,N,prob=Prob,replace=TRUE)
              W2 = huge
              for (it in 1:maxiterations) {
                for (c in 1:numberofclusters) {
                  centroids[[c]] = matrix(0,NumberOfImp,J)     
                  for (m in 1:NumberOfImp) {
                   centroids[[c]][m,] = colMeans(data[[m]][scluster==c,])
                  }
                }

                data_array <- array(0, dim = c(J, N, NumberOfImp))
                
                for (m in 1:NumberOfImp) {
                  data_array[,,m] <- as.matrix(data[[m]])
                }

                centroid_array <- array(0, dim = c(numberofclusters, J, NumberOfImp))
                for (c in 1:numberofclusters) centroid_array[c,,] <- centroids[[c]]

                result <- .Fortran("w_bar",
                   data = as.double(data_array),
                   centroids = as.double(centroid_array),
                   n = as.integer(N),
                   p = as.integer(J),
                   m = as.integer(NumberOfImp),
                   c = as.integer(numberofclusters),
                   cluster = integer(N))
                
                scluster2 <- result$cluster

                scluster = scluster2
                if (length(unique(scluster)) < numberofclusters) {
                  minclustersize <- 0
                } else {
                  minclustersize <- min(table(scluster2))
                }

                if (minclustersize > 1) {
                  result <- .Fortran("w1calc",
                   data = data_array,
                   scluster = as.integer(scluster),
                   n = as.integer(N),
                   p = as.integer(J),
                   nimp = as.integer(NumberOfImp),
                   ncluster = as.integer(numberofclusters),
                   w1 = double(1))

                  W1 <- result$w1
                  if (abs(W2-W1) > tiny) {
                    W2 = W1 
                  } else {
                    break
                  }
                } else {
                  break
                }  
              }
            }
            sclusters=cbind(sclusters,scluster)  
            W=cbind(W,W2)
          }
          indexbest = which(W==min(W))
          if (length(indexbest > 1)) indexbest = indexbest[1]
          bestcluster = sclusters[,indexbest]
          resultslist <- list()
          resultslist$centroids_m <- centroid_array
          resultslist$centroids <- apply(centroid_array, c(1, 2), mean)
          resultslist$solution <- bestcluster
          resultslist
}

Wstar = function(data,
                 numberofclusters,
                 N,
                 J,
                 NumberOfImp,
                 maxiterations = 100,
                 numberofstarts = 20,
                 Prob = rep(1/numberofclusters, times = numberofclusters)) {
          huge = 1E10
          tiny = 1E-7
          centroids = list()
          centroids_m <- list()

          W=c()
          sclusters=c()
          dyn.load("wcalc.dll") 
          for (s in 1:numberofstarts){
            minclustersize = 1
            while (minclustersize <= 1) {
              scluster = sample(1:numberofclusters,N,prob=Prob,replace=TRUE)
              W2 = huge
              for (it in 1:maxiterations) {
                for (c in 1:numberofclusters) {
                  centroids[[c]] = rep(0, times = J)
                  centroids_m[[c]] = matrix(0,NumberOfImp,J)  
                  for (m in 1:NumberOfImp) {
                   centroids[[c]] = centroids[[c]] + colMeans(data[[m]][scluster==c,])
                   centroids_m[[c]][m,] = colMeans(data[[m]][scluster==c,])
                  }
                  centroids[[c]] = centroids[[c]]/NumberOfImp
                }

                data_array <- array(0, dim = c(J, N, NumberOfImp))
                
                for (m in 1:NumberOfImp) {
                  data_array[,,m] <- as.matrix(data[[m]])
                }

                centroid_array <- array(0, dim = c(numberofclusters, J))
                for (c in 1:numberofclusters) centroid_array[c,] <- centroids[[c]]

                centroid_array_m <- array(0, dim = c(numberofclusters, J, NumberOfImp))
                for (c in 1:numberofclusters) centroid_array_m[c,,] <- centroids_m[[c]]

                result <- .Fortran("w_star",
                   data = as.double(data_array),
                   centroids = as.double(centroid_array),
                   n = as.integer(N),
                   p = as.integer(J),
                   m = as.integer(NumberOfImp),
                   c = as.integer(numberofclusters),
                   cluster = integer(N))

                scluster2 <- result$cluster

                scluster = scluster2
                if (length(unique(scluster)) < numberofclusters) {
                  minclustersize <- 0
                } else {
                  minclustersize <- min(table(scluster2))
                }

                if (minclustersize > 1) {
                  result <- .Fortran("w1calc",
                   data = data_array,
                   scluster = as.integer(scluster),
                   n = as.integer(N),
                   p = as.integer(J),
                   nimp = as.integer(NumberOfImp),
                   ncluster = as.integer(numberofclusters),
                   w1 = double(1))

                  W1 <- result$w1
                  if (abs(W2-W1) > tiny) {
                    W2 = W1 
                  } else {
                    break
                  }
                } else {
                  break
                }  
              }
            }
            sclusters=cbind(sclusters,scluster)  
            W=cbind(W,W2)
          }
          indexbest = which(W==min(W))
          if (length(indexbest > 1)) indexbest = indexbest[1]
          bestcluster = sclusters[,indexbest]
          resultslist <- list()
          resultslist$solution <- bestcluster
          resultslist$centroids_m <- centroid_array_m
          resultslist$centroids <- centroid_array
          resultslist
}

clusGapMI <- function (x, 
                       K.max, 
                       NumberOfImp, 
                       poolmethod, 
                       n.start=20, 
                       B = 100, 
                       d.power, 
                       spaceH0 = "scaledPCA", 
                       n)
{
    if (B != (B. <- as.integer(B)) || (B <- B.) <= 0) 
        stop("'B' has to be a positive integer")
    cl. <- match.call()
    if (is.data.frame(x)) 
        x <- as.matrix(x)
    ii <- seq_len(n)
    W.k <- function(X, kk, M, poolmethod, Prob) {
        clus <- if (kk > 1) {
            if (poolmethod == "Wbar") 
              Wbar(data=X,numberofclusters=kk,N=nrow(X[[1]]),J=ncol(X[[1]]),numberofstarts=n.start,NumberOfImp=NumberOfImp,Prob=Prob)$solution
            else if (poolmethod == "Wstar") 
              Wstar(data=X,numberofclusters=kk,N=nrow(X[[1]]),J=ncol(X[[1]]),numberofstarts=n.start,NumberOfImp=NumberOfImp,Prob=Prob)$solution             
        } else { 
          rep.int(1L, nrow(X[[1]]))
        }
        xs <- 0
        for (m in 1:M) {
          xs <- xs + 0.5 * sum(vapply(split(ii, clus), function(I) {
              xsm <- X[[m]][I, , drop = FALSE]
              sum(dist(xsm)^d.power/nrow(xsm))
          }, 0))
        }
        xs/M
    }
    logW <- E.logW <- SE.sim <- testje <- numeric(K.max)
    for (k in 1:K.max) {
      Prob = rep(1/k, times = k)
      logW[k] <- log(W.k(x, k, NumberOfImp, poolmethod, Prob))
    }

    spaceH0 <- "scaledPCA"
    
    xs <- lapply(x, function(C) {scale(C, center = TRUE, scale = FALSE)})

    m.x <- lapply(xs, function(C) {rep(attr(C, "scaled:center"), each = n)})
    switch(spaceH0, scaledPCA = {
        for (m in 1:NumberOfImp) {
          V.sx <- svd(xs[[m]], nu = 0)$v 
          xs[[m]] <- xs[[m]] %*% V.sx
        }
    }, original = {
    }, stop("invalid 'spaceH0':", spaceH0))

    rng.xl <- list()
    for (m in 1:NumberOfImp) {
      rng.xlm <- apply(xs[[m]], 2L, range)
      rng.xl[[m]] <- rng.xlm
    }
    logWks <- matrix(0, B, K.max)

    z <- list()
    for (b in 1:B) {
      zcore <- runif(n*ncol(x[[1]]))
      zcore <- matrix(zcore,n,ncol(x[[1]]))
      for (m in 1:NumberOfImp) {
        minimumz <- matrix(rep(1, times = n))%*%rng.xl[[m]][1,]
        maximumz <- matrix(rep(1, times = n))%*%rng.xl[[m]][2,]

        z1 <- minimumz + zcore*(maximumz-minimumz)

        zm <- switch(spaceH0, scaledPCA = tcrossprod(z1, V.sx), 
            original = z1) + m.x[[m]]
        z[[m]] <- zm
      }
      for (k in 1:K.max) {
          Prob = rep(1/k, times = k)
          logWks[b, k] <-  log(W.k(z, k, NumberOfImp, poolmethod, Prob))
      }
    }
    gap <- E.logW - logW
    resultslist <- list()
    resultslist$E.logW <- colMeans(logWks)
    resultslist$logW <- logW
    resultslist$SE.sim <- sqrt((1 + 1/B) * apply(logWks, 2, var))
    resultslist$gap <- gap
    resultslist$bestK <- which(gap == max(gap))
    resultslist
}
