setwd('C:/Users/RYU/Documents/R/TRSD_R')
data <- read.csv("data/checkerboard.csv", header=F)
colnames(data) <- c('y', 'X1', 'X2')

for (i in 7:9) {
  data <- read.csv("data/mnist.csv", header=F, sep=",")
  colnames(data)[1] <- 'y'
  
  data[which(data$y != i), "y"] <- -1
  data[which(data$y == i), "y"] <- 1
  write.table(data, paste("data/mnist_", i, ".csv", sep=""), sep = ",", row.names = F, col.names = F)
}

data_0 <- data[which(data$y == 0), ]
data_1 <- data[which(data$y != 1), ]
data_2 <- data[which(data$y != 2), ]

rownames(data_2) <- 1:nrow(data_2)

write.csv(data_2, "data/waveform-2.csv", row.names = F)
data(iris3)

setosa <- iris[which(iris$Species != "setosa"), ]
versicolor <- iris[which(iris$Species != "versicolor"), ]
virginica <- iris[which(iris$Species != "virginica"), ]

write.csv(virginica, "data/iris-virginica.csv", row.names = F)
