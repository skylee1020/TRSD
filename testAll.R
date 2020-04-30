library(tree);library(rpart);library(e1071);library(MASS);library(RANN);
suppressMessages(library(Rfast));library(dict);library(xlsx);
suppressMessages(library(caret));library(optparse);library(progress);
suppressMessages(library(randomForest));
setwd("C:/Users/RYU/Documents/R/TRSD_R")
source("algorithm.R")

opt_list <- list(
  make_option(c("--gamma"), type="double"),
  make_option(c("--C"), type="double", default=1.0),
  make_option(c("--nfold"), type="double", default=3),
  make_option(c("--shuffle"), type="logical", action="store_true", default=FALSE),
  make_option(c("-s", "--scale"), type="logical", action="store_true", default=F)
)

methods <- c("SVIS")
seeds <- c(502124, 7519, 92241)
betas <- c(0.1)
opt_parser <- OptionParser(option_list = opt_list)
opt <- parse_args(opt_parser)
file_names <- dir("data")

for (method in methods) {
  if (opt$scale) {
    params <- read.xlsx('experiments_scaled/grid_search.xlsx', sheetName = method)
  } else {
    params <- read.xlsx('experiments/grid_search.xlsx', sheetName = method)
  }
  if(method != "FULL") {algorithm <- get(method)}
  for (beta in betas) {
    append <- TRUE
    if(method == "FULL" && beta %in% c(0.1, 0.2)) {next}
    mincut <- 1 / beta; mindev <- 1e-2; #C <- opt$C; 
    for (i in 1:length(file_names)) {
      data_name <- gsub(".csv", "", file_names[i])
      #if (!data_name %in% c("a9a", "banana", "w8a")) {next}
      gamma <- params[which(params[, "data"] == data_name), 3]
      C <- params[which(params[, "data"] == data_name), 4]
      file_path <- paste("data/", file_names[i], sep="")
      data <- read.csv(file_path, header=FALSE)
      colnames(data) <- c("y", 1:(ncol(data)-1))
      data <- data.frame(data)
      if(opt$shuffle) {
        shuffle_index <- sample(nrow(data), nrow(data))
        data <- data[shuffle_index, ]
      }
      data <- remove_constant_features(data)
      if (opt$scale) {data[, 2:ncol(data)] <- apply(data[, 2:ncol(data)], 2, scale)}
      data$y <- as.factor(data$y)
      cat(sprintf("start running [%s] algorithm on [%s] dataset with [beta=%.3f, gamma=%.3f, C=%.3f].\n", 
                  method, data_name, beta, gamma, C))
      
      pb <- progress_bar$new(format = sprintf("[%s | %s] [:bar] (:current/:total)", method, data_name), 
                             total = opt$nfold * length(seeds), show_after = 0)
      pb$tick(0)
      
      accuracies <- c(); percents <- c(); 
      index_times <- c(); svm_times <- c(); total_times <- c()
      for (seed in seeds) {
        set.seed(seed)
        folds <- createFolds(factor(data$y), k = opt$nfold, list = TRUE)
        for (k in 1:length(folds)) {
          test_index <- folds[[k]]
          train <- data[-test_index,]
          test <- data[test_index,]
          t1 <- Sys.time()
          if (method == "FULL") {
            candidate_index <- 1:nrow(train)
          }
          else {
            flag <- TRUE
            tryCatch(candidate_index <- algorithm(train, mincut, mindev, beta),
                     error = function(e){flag <<- FALSE})
            if(!flag) {next}
          }
          t2 <- Sys.time()
          candidate_train <- train[candidate_index, ]
          candidate_train <- remove_constant_features(candidate_train)
          flag <- TRUE
          tryCatch(svm.m <- svm(y~., data=candidate_train, gamma=gamma, cost=C),
                   error = function(e){flag <<- FALSE})
          if(!flag) {next}
          t3 <- Sys.time()
          pred <- predict(svm.m, newdata=test, type="class")
          acc <- 100 * mean(pred == test$y)
          accuracies <- c(accuracies, acc)
          index_times <- c(index_times, as.numeric(difftime(t2, t1, units="secs")))
          svm_times <- c(svm_times, as.numeric(difftime(t3, t2, units="secs")))
          total_times <- c(total_times, as.numeric(difftime(t3, t1, units="secs")))
          percents <- c(percents, length(candidate_index)/nrow(train))
          pb$tick()
        }
      }
      if(length(accuracies) == 0) {cat(sprintf("This method doesn't work on [%s]\n", data_name)); next}
      
      cat(sprintf("acc: %.3f%%, std: %.3f, index_time: %.3f, svm_time: %.3f, total_time: %.3f, percent: %.3f\n",
                  mean(accuracies), sd(accuracies), mean(index_times), 
                  mean(svm_times), mean(total_times), mean(percents)))
      accuracies <- c(accuracies, mean(accuracies))
      accuracies <- c(accuracies, sd(accuracies))
      index_times <- c(index_times, mean(index_times))
      index_times <- c(index_times, sd(index_times))
      svm_times <- c(svm_times, mean(svm_times))
      svm_times <- c(svm_times, sd(svm_times))
      total_times <- c(total_times, mean(total_times))
      total_times <- c(total_times, sd(total_times))
      percents <- c(percents, mean(percents))
      percents <- c(percents, sd(percents))
      df <- data.frame(accuracy = accuracies, index_times=index_times, 
                       svm_times=svm_times, total_times = total_times, percent = percents)
      if (opt$scale) {exp_folder <- "experiments_scaled"} else {exp_folder <- "experiments"}
      if (!dir.exists(exp_folder)) {dir.create(exp_folder)}
      if (!file.exists(paste(exp_folder, "/", method, "-", beta, ".xlsx", sep=""))) {
        file.create(paste(exp_folder, "/", method, "-", beta, ".xlsx", sep=""))
      }
      if (i > 1) {append <- TRUE}
      write.xlsx(df, file=paste(exp_folder, "/", method, "-", beta, ".xlsx", sep=""),
                 sheetName=data_name, append=append,  row.names=F)
    }
  }
}

