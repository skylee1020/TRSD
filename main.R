library(tree);library(rpart);library(e1071);library(MASS);library(RANN);
suppressMessages(library(Rfast));library(dict);library(xlsx);
suppressMessages(library(caret));library(optparse);library(progress);
suppressMessages(library(randomForest));
setwd('C:/Users/RYU/Documents/R/TRSD_R')
source('algorithm.R')

opt_list <- list(
  make_option(c("-d", "--data"), type="character", default="fourclass"),
  make_option(c("-b", "--beta"), type="double", default=0.1),
  make_option(c("-g", "--gamma"), type="double"),
  make_option(c("--C"), type="double", default=1.0),
  make_option(c("--seed"), type="double", default=42),
  make_option(c("--nfold"), type="double", default=5),
  make_option(c("-m", "--method"), type="character", default="pRTRSD",
              help="[TRSD, DTFLD, CVS, SVDISC, RVS, FULL]"),
  make_option(c("-p", "--pca"), type="logical", action="store_true", default=F),
  make_option(c("--pdim"), type="double", default=3),
  make_option(c("-s", "--scale"), type="logical", action="store_true", default=F),
  make_option(c("--shuffle"), type="logical", action="store_true", default=FALSE)
)

opt_parser <- OptionParser(option_list = opt_list)
opt <- parse_args(opt_parser)
file_path <- paste('data/', opt$data, '.csv', sep='')
data <- read.csv(file_path, header=FALSE)
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
if (opt$scale) {data[, 2:ncol(data)] <- apply(data[, 2:ncol(data)], 2, normalize)}

data$y <- as.factor(data$y)
set.seed(opt$seed)
folds <- createFolds(factor(data$y), k = opt$nfold, list = TRUE)

if(opt$method != "FULL") {algorithm <- get(opt$method)}
if (is.null(opt$gamma)) {gamma <- 1/(ncol(data) - 1)} else {gamma <- opt$gamma}
beta <- opt$beta; C <- opt$C
mincut <- 1 / beta; mindev <- 1e-2

cat(sprintf("start running [%s] algorithm on [%s] dataset with [beta=%.3f, gamma=%.3f, C=%.3f].\n", 
            opt$method, opt$data, opt$beta, gamma, opt$C))
accuracies <- c(); percents <- c(); times <- c();
precisions <- c(); recalls <- c(); f1scores <- c();
pb <- progress_bar$new(format = sprintf("[%s | %s] [:bar] (:current/:total)", opt$method, opt$data), 
                       total = opt$nfold, show_after = 0)
capture.output(pb$tick(0), file="NUL")
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
  t2 <- Sys.time()
  candidate_train <- train[candidate_index, ]
  candidate_train <- remove_constant_features(candidate_train)
  tryCatch(svm.m <- svm(y~., data=candidate_train, gamma=gamma, cost=C),
           error = function(e){print(e); stop()})
  t3 <- Sys.time()
  pred <- predict(svm.m, newdata=test, type="class")
  xtab <- table(pred, test$y)
  result <- confusionMatrix(xtab)
  accuracies <- c(accuracies, 100 * result$overall["Accuracy"])
  precisions <- c(precisions, 100 * result$byClass["Precision"])
  recalls <- c(recalls, 100 * result$byClass["Recall"])
  f1scores <- c(f1scores, 100 * result$byClass["F1"])
  index_times <- c(times, as.numeric(difftime(t2, t1, units="secs")))
  svm_times <- c(times, as.numeric(difftime(t3, t2, units="secs")))
  total_times <- c(times, as.numeric(difftime(t3, t1, units="secs")))
  percents <- c(percents, length(candidate_index)/nrow(train))
  pb$tick()
}
cat(sprintf(paste("[%s | %s] acc: %.3f%%, std: %.3f, precision: %.3f, recall: %.3f, f1_score: %.3f\n", 
            "index_time: %.3f, svm_time: %.3f, total_time: %.3f, percent: %.3f\n", sep=""),
            opt$method, opt$data, mean(accuracies), sd(accuracies), mean(precisions), mean(recalls), mean(f1scores),
            mean(index_times), mean(svm_times), mean(total_times), mean(percents)))