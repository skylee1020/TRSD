library(tree);library(rpart);library(e1071);library(MASS);library(RANN);
suppressMessages(library(Rfast));library(dict);library(xlsx);
suppressMessages(library(caret));library(optparse);library(progress);
setwd('C:/Users/RYU/Documents/R/TRSD_R')
source('algorithm.R')

opt_list <- list(
  make_option(c("-d", "--data"), type="character", default="checkerboard"),
  make_option(c("-b", "--beta"), type="double", default=0.1),
  make_option(c("--seed"), type="double", default=42),
  make_option(c("--nfold"), type="double", default=5),
  make_option(c("-m", "--method"), type="character", default="TRSD",
              help="[TRSD, DTFLD, CVS, SVDISC, RVS, FULL]"),
  make_option(c("-p", "--pca"), type="logical", action="store_true", default=F),
  make_option(c("--pdim"), type="double", default=3),
  make_option(c("-s", "--scale"), type="logical", action="store_true", default=F),
  make_option(c("--shuffle"), type="logical", action="store_true", default=FALSE)
)

opt_parser <- OptionParser(option_list = opt_list)
opt <- parse_args(opt_parser)

file_path <- paste('data/', opt$data, '.csv', sep='')
data <- read.csv(file_path)

gammas <- unique(c(10^(-3:0), 1/(ncol(data) - 1)))
Cs <- 10^(-1:2)
beta <- opt$beta

mincut <- 1 / beta; mindev <- 1e-2
data <- remove_constant_features(data)
if (opt$pca) {
  data <- cbind(data[, 1], prcomp(data[, -1], center=TRUE, scale.=TRUE)$x[, 1:min(3, ncol(data) - 1)])
}

colnames(data) <- c('y', 1:(ncol(data)-1))
data <- data.frame(data)
if(opt$shuffle) {
  shuffle_index <- sample(nrow(data), nrow(data))
  data <- data[shuffle_index, ]
}
if (opt$scale) {data[, 2:ncol(data)] <- apply(data[, 2:ncol(data)], 2, scale)}

set.seed(opt$seed)
folds <- createFolds(factor(data$y), k = opt$nfold, list = TRUE)

if(opt$method != "FULL") {algorithm <- get(opt$method)}

best_accuracy <- 0; best_gamma <- NULL; best_C <- NULL;
pb <- progress_bar$new(format = sprintf("[%s | %s] [:bar] (:current/:total)", opt$method, opt$data), 
                       total = opt$nfold * length(gammas) * length(Cs), show_after = 0)
capture.output(pb$tick(0), file="NUL")
for (gamma in gammas) {
  for (C in Cs) {
    accuracies <- c(); times <- c(); percents <- c()
    for (i in 1:length(folds)) {
      test_index <- folds[[i]]
      train <- data[-test_index,]
      test <- data[test_index,]
      t1 <- Sys.time()
      if (opt$method == "FULL") {
        candidate_index <- 1:nrow(train)
      }
      else {
        suppressWarnings(candidate_index <- algorithm(train, mincut, mindev, beta))
      }
      candidate_train <- train[candidate_index, ]
      candidate_train <- remove_constant_features(candidate_train)
      flag <- TRUE
      tryCatch(svm.m <- svm(as.factor(y)~., data=candidate_train, gamma=gamma, cost=C, scale=FALSE),
               error = function(e){flag <<- FALSE})
      pb$tick()
      if(!flag) {next}
      t2 <- Sys.time()
      pred <- predict(svm.m, newdata=test, type="class")
      acc <- 100 * mean(pred == test[, 1])
      accuracies <- c(accuracies, acc)
      times <- c(times, as.numeric(difftime(t2, t1, units="secs")))
      percents <- c(percents, length(candidate_index)/nrow(train))
    }
    if (best_accuracy < mean(accuracies)) {
      best_accuracy <- mean(accuracies)
      best_gamma <- gamma
      best_C <- C
    }
  }
}
if(length(accuracies) == 0) {cat(sprintf("This method doesn't work on [%s]\n", data_name)); stop()}
cat(sprintf("[%s | %s] acc: %.3f%%, std: %.3f, time: %.3f, percent: %.3f, best_gamma: %.3f, best_C: %.3f\n",
            opt$method, opt$data, mean(accuracies), sd(accuracies), mean(times), mean(percents), best_gamma, best_C))