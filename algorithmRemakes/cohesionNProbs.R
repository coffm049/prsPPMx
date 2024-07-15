library(tidyverse)

lcohesion <- function(Mprecision, clusterSize) {
  return(log(Mprecision * gamma(clusterSize)))
}


assignIcluster <- function(i, data, Mprecision) {
  # remove the ith row from data
  # might just pop head off and add to tail
  tempData <- data[-i, ]

  clusterSizes <- table(tempData[, "cluster"])
  clusterMeans <- tapply(tempData$x, tempData$cluster, mean)
  # initialize "empty" vector to store probs
  clusterllikes <- clusterSizes
  for (j in names(clusterSizes)) {
    if (clusterSizes[j] == 0) {
      clusterllikes[j] <- lcohesion(Mprecision, 1)
    } else {
      clusterllikes[j] <- lcohesion(Mprecision, clusterSizes[j] + 1) / lcohesion(Mprecision, clusterSizes[j])
      clusterllikes[j] <- clusterllikes[j] + dnorm(data[i, "x"], clusterMeans[j], log = TRUE)
    }
  }
  # replace missing vluaes inclusterllikes with 0
  is.na(clusterllikes) <- 0
  data[i, "cluster"] <- sample(names(clusterSizes), size = 1, prob = exp(clusterllikes))
  return(data)
}

# loop assignIcluster over all i subjects
assignClusters <- function(data, Mprecision) {
  for (i in nrow(data)) {
    data <- assignIcluster(i, data, Mprecision)
  }
  return(data)
}

# Make data to test
data <- data.frame("clusterTrue" = rep(c(1, 15, 30), each = 10))
data["x"] <- rnorm(n = nrow(data), data$cluster)
data$cluster <- rep(1, each = nrow(data))
Mprecision <- 1
mNeal <- 1
data$cluster <- factor(data$cluster, levels = 1:6)
data2 <- assignIcluster(1, data, 1)

for (i in 1:1000) {
  data <- assignClusters(data, Mprecision)
}
