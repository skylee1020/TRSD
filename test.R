library(tree);library(rpart);library(e1071);library(MASS);library(RANN);
suppressMessages(library(Rfast));library(dict);library(xlsx);
suppressMessages(library(caret));library(optparse);library(progress);
setwd("C:/Users/RYU/Documents/R/TRSD_R")
source("algorithm.R")

opt_list <- list(
  make_option(c("--gamma"), type="double"),
  make_option(c("--C"), type="double", default=1.0),
  make_option(c("-m", "--method"), type="character", default="TRSD",
              help="[TRSD, DTFLD, CVS, SVDISC, RVS, FULL]"),
  make_option(c("--shuffle"), type="logical", action="store_true", default=FALSE)
)

seeds <- c(102, 7519, 92241)
betas <- c(0.1, 0.2, 0.3)
opt_parser <- OptionParser(option_list = opt_list)
opt <- parse_args(opt_parser)
if(opt$method != "FULL") {method <- get(opt$method)}
file_names <- dir("data")

for (beta in betas) {
  if(method == "FULL" && beta %in% c(0.1, 0.2)) {next}
  C <- opt$C; mincut <- 1 / beta; mindev <- 1e-2;
  for (file_name in file_names) {
    data_name <- gsub(".csv", "", file_name)
    if (data_name %in% c("rna", "covtype")) {next}
    file_path <- paste("data/", file_name, sep="")
    data <- read.csv(file_path, header=FALSE)
    colnames(data) <- c("y", 1:(ncol(data)-1))
    data <- data.frame(data)
    if(opt$shuffle) {
      shuffle_index <- sample(nrow(data), nrow(data))
      data <- data[shuffle_index, ]
    }
    data <- remove_constant_features(data)
    if (is.null(opt$gamma)) {gamma <- 1/(ncol(data) - 1)} else {gamma <- opt$gamma}
    cat(sprintf("start running [%s] algorithm on [%s] dataset with [beta=%.2f, gamma=%.2f, C=%.2f].\n", 
                opt$method, data_name, beta, gamma, opt$C))
    accuracies <- c(); percents <- c(); 
    index_times <- c(); svm_times <- c(); total_times <- c()
    for (seed in seeds) {
      set.seed(seed)
      folds <- createFolds(factor(data$y), k = 5, list = TRUE)
      for (k in 1:length(folds)) {
        test_index <- folds[[k]]
        train <- data[-test_index,]
        test <- data[test_index,]
        t1 <- Sys.time()
        if (opt$method == "FULL") {
          candidate_index <- 1:nrow(train)
        }
        else {
          candidate_index <- method(train, mincut, mindev, beta)
        }
        t2 <- Sys.time()
        candidate_train <- train[candidate_index, ]
        candidate_train <- remove_constant_features(candidate_train)
        flag <- TRUE
        tryCatch(svm.m <- svm(as.factor(y)~., data=candidate_train, gamma=gamma, cost=C),
                 error = function(e){flag <<- FALSE})
        if(!flag) {next}
        t3 <- Sys.time()
        pred <- predict(svm.m, newdata=test, type="class")
        acc <- 100 * mean(pred == test$y)
        accuracies <- c(accuracies, acc)
        index_times <- c(times, as.numeric(difftime(t2, t1, units="secs")))
        svm_times <- c(times, as.numeric(difftime(t3, t2, units="secs")))
        total_times <- c(times, as.numeric(difftime(t3, t1, units="secs")))
        percents <- c(percents, length(candidate_index)/nrow(train))
      }
    }
    if(length(accuracies) == 0) {cat(sprintf("This method doesn't work on [%s]\n", data_name)); next}
    
    cat(sprintf("acc: %.3f%%, std: %.3f, index_time: %.3f, svm_time: %.3f, total_time: %.3f, percent: %.3f\n",
                mean(accuracies), sd(accuracies), mean(index_times), 
                mean(svm_times), mean(total_times), mean(percents)))
    accuracies <- c(accuracies, mean(accuracies))
    accuracies <- c(accuracies, sd(accuracies))
    times <- c(times, mean(times))
    times <- c(times, sd(times))
    percents <- c(percents, mean(percents))
    percents <- c(percents, sd(percents))
    if (dir.exists('experiments')) {dir.create('experiments')}
    df <- data.frame(accuracy = accuracies, time = times, percent = percents)
    write.xlsx(df, file=paste("experiments/", opt$method, "-", beta, ".xlsx", sep=""),
               sheetName=data_name, append=TRUE,  row.names=F)
  }
}
