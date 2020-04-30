library(tree);library(rpart);library(e1071);library(MASS);library(RANN);
suppressMessages(library(Rfast));library(dict);library(xlsx);
suppressMessages(library(caret));library(optparse);library(progress);
suppressMessages(library(randomForest));
setwd("C:/Users/RYU/Documents/R/TRSD_R")
source("algorithm.R")

opt_list <- list(
  make_option(c("--beta"), type="double", default=0.1),
  make_option(c("--seed"), type="double", default=42),
  make_option(c("--nfold"), type="double", default=2),
  make_option(c("-s", "--scale"), type="logical", action="store_true", default=F)
)

opt_parser <- OptionParser(option_list = opt_list)
opt <- parse_args(opt_parser)

file_names <- dir("data")

methods <- c("SVIS")
beta <- opt$beta
mincut <- 1 / beta; mindev <- 1e-2

append <- TRUE
for (i in 1:length(methods)) {
  method <- methods[i]
  bestAccs <- c(); dataNames <- c(); bestGams <- c(); bestCs <-c();
  for (file_name in file_names) {
    data_name <- gsub(".csv", "", file_name)
    #if (!data_name %in% c("a9a", "banana", "w8a")) {next}
    file_path <- paste("data/", file_name, sep="")
    data <- read.csv(file_path, header=FALSE)
    colnames(data) <- c("y", 1:(ncol(data)-1))
    data <- data.frame(data)
    data <- remove_constant_features(data)
    if (opt$scale) {data[, 2:ncol(data)] <- apply(data[, 2:ncol(data)], 2, scale)}
    data$y <- as.factor(data$y)
    set.seed(opt$seed)
    folds <- createFolds(factor(data$y), k = opt$nfold, list = TRUE)
    
    best_accuracy <- 0; best_gamma <- NULL; best_C <- NULL;
    algorithm <- get(method)
    
    gammas <- unique(c(10^(-3:0), 1/(ncol(data) - 1)))
    Cs <- 10^(-1:2)
    
    pb <- progress_bar$new(format = sprintf("[%s | %s] [:bar] (:current/:total)", method, data_name), 
                           total = opt$nfold * length(gammas) * length(Cs), show_after = 0)
    pb$tick(0)
    for (gamma in gammas) {
      for (C in Cs) {
        accuracies <- c(); times <- c(); percents <- c()
        for (j in 1:length(folds)) {
          test_index <- folds[[j]]
          train <- data[-test_index,]
          test <- data[test_index,]
          t1 <- Sys.time()
          suppressWarnings(candidate_index <- algorithm(train, mincut, mindev, beta))
          candidate_train <- train[candidate_index, ]
          candidate_train <- remove_constant_features(candidate_train)
          flag <- TRUE
          tryCatch(svm.m <- svm(y~., data=candidate_train, gamma=gamma, cost=C),
                   error = function(e){flag <<- FALSE})
          if(!flag) {next}
          t2 <- Sys.time()
          pred <- predict(svm.m, newdata=test, type="class")
          acc <- 100 * mean(pred == test[, 1])
          accuracies <- c(accuracies, acc)
          times <- c(times, as.numeric(difftime(t2, t1, units="secs")))
          percents <- c(percents, length(candidate_index)/nrow(train))
          pb$tick()
        }
        if(length(accuracies) == 0) {next}
        if (best_accuracy < mean(accuracies)) {
          best_accuracy <- mean(accuracies)
          best_time <- mean(times)
          best_percent <- mean(percents)
          best_gamma <- gamma
          best_C <- C
        }
      }
    }
    if(is.null(best_C)) {cat(sprintf("[%s] doesn't work on [%s]\n", method, data_name)); next}
    pb$terminate()
    cat(sprintf("[%s | %s] acc: %.3f%%, time: %.3f, percent: %.3f, best_gamma: %.3f, best_C: %.3f\n",
                method, data_name, best_accuracy, best_time, best_percent, best_gamma, best_C))
    dataNames <- c(dataNames, data_name)
    bestAccs <- c(bestAccs, best_accuracy)
    bestGams <- c(bestGams, best_gamma)
    bestCs <- c(bestCs, best_C)
  }
  df <- data.frame(data = dataNames, best_Acc = bestAccs, best_Gamma = bestGams, best_C = bestCs)
  if (opt$scale) {
    save_name <- "grid_search_scaled_add.xlsx"
  } else {
    save_name <- "grid_search.xlsx"
  }
  if (!dir.exists('experiments')) {dir.create('experiments')}
  if (!file.exists(paste("experiments/", save_name, sep = ""))) {
    file.create(paste("experiments/", save_name, sep = ""))
  }
  if (i > 1) {append <- TRUE}
  write.xlsx(df, file=paste("experiments/", save_name, sep = ""),
             sheetName=method, append=append,  row.names=F)
}
