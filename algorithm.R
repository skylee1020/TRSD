get_neighbors <- function(train, tree_matrix) {
  leaf_nodes <- which(tree_matrix$var=='leaf')
  n_leaves <- length(leaf_nodes)
  n_nodes <- n_leaves * 2 - 1
  
  d <- dict()
  for(i in 1:n_nodes) {
    d[[tree_matrix$node[i]]] <- i
  }
  
  B <- c(as.numeric(apply(train[,-1], 2, min)), as.numeric(apply(train[, -1], 2, max)))
  M <- matrix(0,nrow=n_nodes, ncol=2*ncol(train[, -1]))
  M[1,] <- B
  
  n <- length(B) / 2
  for (i in 2:n_nodes) {
    p <- d[[tree_matrix$node[i] %/% 2]] # parent node to matrix row mapping
    if (tree_matrix$node[i] %% 2 == 0) {
      bound <- M[d[[tree_matrix$node[p]]],]
      bound[n + as.numeric(tree_matrix$var[p])] <- tree_matrix$threshold[p]
      M[d[[tree_matrix$node[i]]],] <- bound
    } else {
      bound <- M[d[[tree_matrix$node[p]]],]
      bound[as.numeric(tree_matrix$var[p])] <- tree_matrix$threshold[p]
      M[d[[tree_matrix$node[i]]],] <- bound
    }
  }
  
  neighbors <- dict()
  for (i in 1:n_leaves) {
    neighbor <- c()
    for (j in (1:n_leaves)[-i]) {
      if (is_neighbor(M[leaf_nodes[i],], M[leaf_nodes[j],])
          && tree_matrix$label[leaf_nodes[i]] != tree_matrix$label[leaf_nodes[j]]) {
        neighbor <- c(neighbor, leaf_nodes[j])
      }
    }
    if (length(neighbor) == 0) {
      next
    }
    neighbors[[leaf_nodes[i]]] <- neighbor
  }
  return(neighbors)
}

get_pairwise_neighbors <-function(train, tree_matrix) {
  
  leaf_nodes <- which(tree_matrix$var=='leaf')
  n_leaves <- length(leaf_nodes)
  n_nodes <- n_leaves * 2 - 1
  
  d <- dict()
  for(i in 1:n_nodes) {
    d[[tree_matrix$node[i]]] <- i
  }
  
  B <- c(as.numeric(apply(train[,-1], 2, min)), as.numeric(apply(train[, -1], 2, max)))
  M <- matrix(0,nrow=n_nodes, ncol=2*ncol(train[, -1]))
  M[1,] <- B
  
  n <- length(B) / 2
  for (i in 2:n_nodes) {
    p <- d[[tree_matrix$node[i] %/% 2]] # parent node to matrix row mapping
    if (tree_matrix$node[i] %% 2 == 0) {
      bound <- M[d[[tree_matrix$node[p]]],]
      bound[n + as.numeric(tree_matrix$var[p])] <- tree_matrix$threshold[p]
      M[d[[tree_matrix$node[i]]],] <- bound
    } else {
      bound <- M[d[[tree_matrix$node[p]]],]
      bound[as.numeric(tree_matrix$var[p])] <- tree_matrix$threshold[p]
      M[d[[tree_matrix$node[i]]],] <- bound
    }
  }
  
  neighbor <- c()
  for (i in 1:(n_leaves-1)) {
    for (j in (i+1):n_leaves) {
      if (is_neighbor(M[leaf_nodes[i],], M[leaf_nodes[j],]) 
          && tree_matrix$label[leaf_nodes[i]] != tree_matrix$label[leaf_nodes[j]]) {
        neighbor <- c(neighbor, c(leaf_nodes[i], leaf_nodes[j]))
      }
    }
  }
  neighbor <- matrix(neighbor, ncol=2, byrow=T)
  return(neighbor)
}


is_neighbor <- function(a, b) {
  n <- length(a) %/% 2
  cond = FALSE
  for (i in 1:n) {
    if (a[n + i] == b[i] || a[i] == b[n + i]) {
      cond = TRUE
      k = i
      break
    }
  }
  if (!cond) {
    return(FALSE)
  }
  
  cond = FALSE
  for (i in (1:n)[-k]) {
    if ((a[i] <= b[i] && b[i] <= a[n + i]) || 
        (a[i] <= b[n + i] && b[n + i] <= a[n + i])) {
      return(TRUE)
    }
  }
  return(FALSE)
}

get_cp <- function(rtree) {
  first <- TRUE
  for (i in 1:nrow(rtree$cptable)) {
    if (rtree$cptable[i, "nsplit"] != 0 ) {
      if (first) {
        min_index <- i
        min_xerror <- rtree$cptable[i, "xerror"]
        first <- FALSE
      } else {
        if (min_xerror > rtree$cptable[i, "xerror"]) {
          min_index <- i
          min_xerror <- rtree$cptable[i, "xerror"]
        }
      }
    }
  }
  cp <- rtree$cptable[min_index, "CP"]
  return (cp)
}

get_splits <- function(ptree) {
  split_values <- c()
  ncompetes <- ptree$frame[which(ptree$frame$var != "<leaf>"), "ncompete"] + 1
  nonzero_splits <- ptree$splits[which(ptree$splits[, "count"] != 0), ]
  split_size <- nrow(nonzero_splits)
  index <- 1
  for (i in 1:length(ncompetes)) {
    split_values <- c(split_values, nonzero_splits[index, 4])
    index <- index + ncompetes[i]
    if (index > split_size) {
      break
    }
  }
  
  k <- 0
  threshold <- c()
  for (i in 1:nrow(ptree$frame)) {
    if (ptree$frame$var[i] != "<leaf>") {
      k <- k + 1
      threshold <- c(threshold, as.character(split_values[k]))
    } else {
      threshold <- c(threshold, "")
    }
  }
  return(threshold)
}

get_height <- function(tree_matrix) {
  height <- 0
  node <- max(tree_matrix$node)
  while (node > 1) {
    height <- height + 1
    node <- node %/% 2
  }
  return(height)
}

get_candidates_v1 <- function(train, dtree, tree_matrix, neighbor, beta) {
  candidate_index <- c()
  
  for (i in 1:nrow(neighbor)) {
    first_node <- neighbor[i, 1]
    second_node <- neighbor[i, 2]
    first_majority_class <- tree_matrix$label[first_node]
    second_majority_class <- tree_matrix$label[second_node]
    first_index <- intersect(which(dtree$where == first_node), which(train$y == first_majority_class))
    second_index <- intersect(which(dtree$where == second_node), which(train$y == second_majority_class))
    cat_index <- c(first_index, second_index)
    x1 <- train[first_index, ]
    x2 <- train[second_index, ]
    x <- rbind(x1, x2)
    m1 <- as.numeric(apply(x1[, -1], 2, mean))
    m2 <- as.numeric(apply(x2[, -1], 2, mean))
    W <- m1 - m2
    M <- (m1 + m2) / 2
    r1 <- as.numeric(sqrt(apply((x1[, -1] - m1) ^ 2, 1, sum)))
    r2 <- as.numeric(sqrt(apply((x2[, -1] - m2) ^ 2, 1, sum)))
    r <- c(r1, r2)
    h <- abs(as.matrix(x[, -1] - M) %*% as.matrix(W)) / sqrt(sum(W ^ 2))
    r <- normalize(r)
    h <- normalize(h)
    t <- sigmoid(r) / (h + 1e-8)
    count <- round(length(t) * beta)
    index <- order(t, decreasing = T)[1:count]
    candidate_index <- c(candidate_index, cat_index[index])
  }
  candidate_index <- unique(candidate_index)
  return(candidate_index)
}


get_candidates_v2 <- function(train, dtree, tree_matrix, neighbors, beta) {
  leaf_nodes <- which(tree_matrix$var=='leaf')
  n_leaves <- length(leaf_nodes)
  
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  candidate_index <- c()
  for (i in 1:n_leaves) {
    main_node <- leaf_nodes[i]
    if (!main_node %in% neighbors$keys()) {
      next
    }
    
    main_train <- train[which(dtree$where == main_node), ]
    main_majority_class <- getmode(main_train$y)
    
    neighbor_nodes <- neighbors[[leaf_nodes[i]]]
    neighbor_train <- train[which(dtree$where %in% neighbor_nodes), ]
    
    neighbor_majority_class <- getmode(neighbor_train$y)
    main_index <- intersect(which(dtree$where == main_node), 
                            which(train$y == main_majority_class))
    neighbor_index <- intersect(which(dtree$where %in% neighbor_nodes), 
                                which(train$y == neighbor_majority_class))
    
    x1 <- train[main_index, ]
    x2 <- train[neighbor_index, ]
    
    m1 <- as.numeric(apply(x1[, -1], 2, mean))
    m2 <- as.numeric(apply(x2[, -1], 2, mean))
    W <- m1 - m2
    M <- (m1 + m2) / 2
    
    main_index <- which(dtree$where == main_node)
    x1 <- train[main_index, ]
    
    r <- as.numeric(sqrt(apply((x1[, -1] - m1) ^ 2, 1, sum)))
    h <- abs(as.matrix(x1[, -1] - M) %*% as.matrix(W)) / sqrt(sum(W ^ 2))
    
    r <- normalize(r)
    h <- normalize(h)
    t <- sigmoid(r) / (h + 1e-8)
    
    count <- round(length(t) * beta)
    
    index <- order(t, decreasing = T)[1:count]
    candidate_index <- c(candidate_index, main_index[index])
  }
  candidate_index <- unique(candidate_index)
  return(candidate_index)
}


get_candidates_v3 <- function(train, dtree, tree_matrix, neighbors, beta) {
  leaf_nodes <- which(tree_matrix$var=='leaf')
  n_leaves <- length(leaf_nodes)
  
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  candidate_index <- c()
  for (i in 1:n_leaves) {
    
    main_node <- leaf_nodes[i]
    if (!main_node %in% neighbors$keys()) {
      next
    }
    
    can_index <- c()
    main_train <- train[which(dtree$where == main_node), ]
    main_majority_class <- getmode(main_train$y)
    
    neighbor_nodes <- neighbors[[leaf_nodes[i]]]
    neighbor_train <- train[which(dtree$where %in% neighbor_nodes), ]
    neighbor_majority_class <- getmode(neighbor_train$y)
    
    main_index <- intersect(which(dtree$where == main_node),
                            which(train$y == main_majority_class))
    neighbor_index <- intersect(which(dtree$where %in% neighbor_nodes),
                                which(train$y == neighbor_majority_class))
    
    x1 <- train[main_index, ]
    o1 <- train[neighbor_index, ]
    
    m1 <- as.numeric(apply(x1[, -1], 2, mean))
    
    dist_center1 <- as.numeric(sqrt(apply((x1[, -1] - m1) ^ 2, 1, sum)))
    dist_opposite1 <- nn2(o1[,-1], x1[,-1], k=1, treetype = "kd")$nn.dist
    
    dist_center1 <- scale(dist_center1)
    dist_opposite1 <- scale(dist_opposite1)
    
    t1 <- sigmoid(dist_center1 / dist_opposite1)
    count <- round(length(t1) * beta)
    index1 <- order(t1, decreasing = T)[1:count]
    
    can_index <- c(can_index, main_index[index1])
    
    
    main_opposite_index <- intersect(which(dtree$where == main_node),
                                     which(train$y == neighbor_majority_class))
    neighbor_opposite_index <- intersect(which(dtree$where %in% neighbor_nodes),
                                         which(train$y == main_majority_class))
    
    if (length(main_opposite_index) * length(neighbor_opposite_index) != 0) {
      x2 <- train[main_opposite_index, ]
      o2 <- train[neighbor_opposite_index, ]
      
      m2 <- as.numeric(apply(x2[, -1], 2, mean))
      
      dist_center2 <- as.numeric(sqrt(apply((x2[, -1] - m2) ^ 2, 1, sum)))
      dist_opposite2 <- nn2(o2[,-1], x2[,-1], k=1, treetype = "kd")$nn.dist
      
      dist_center2 <- scale(dist_center2)
      dist_opposite2 <- scale(dist_opposite2)
      t2 <- sigmoid(dist_center2 / dist_opposite2)
      count2 <- round(length(t2) * beta)
      index2 <- order(t2, decreasing = T)[1:count2]
      can_index <- c(can_index, main_opposite_index[index2])
    }
    
    candidate_index <- c(candidate_index, can_index)
  }
  candidate_index <- unique(candidate_index)
  return(candidate_index)
}


get_candidates_v4 <- function(train, ptree, tree_matrix, neighbors, beta) {
  leaf_nodes <- which(tree_matrix$var=='leaf')
  n_leaves <- length(leaf_nodes)
  
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  candidate_index <- c()
  for (i in 1:n_leaves) {
    main_node <- leaf_nodes[i]
    if (!main_node %in% neighbors$keys()) {
      next
    }
    neighbor_nodes <- neighbors[[main_node]]
    main_train <- train[which(ptree$where == main_node), ]
    main_majority_class <- getmode(main_train$y)
    main_index <- intersect(which(ptree$where == main_node), 
                            which(train$y == main_majority_class))
    x1 <- train[main_index, ]
    m1 <- as.numeric(apply(x1[, -1], 2, mean))
    main_index <- which(ptree$where == main_node)
    x1 <- train[main_index, ]
    r <- as.numeric(sqrt(apply((x1[, -1] - m1) ^ 2, 1, sum)))
    total_r <- c(); total_h <- c(); total_index <- c()
    for (j in 1:length(neighbor_nodes)) {
      neighbor_train <- train[which(ptree$where == neighbor_nodes[j]), ]
      neighbor_majority_class <- getmode(neighbor_train$y)
      neighbor_index <- intersect(which(ptree$where == neighbor_nodes[j]), 
                                  which(train$y == neighbor_majority_class))
      
      x2 <- train[neighbor_index, ]
      m2 <- as.numeric(apply(x2[, -1], 2, mean))
      W <- m1 - m2
      M <- (m1 + m2) / 2
      h <- abs(as.matrix(x1[, -1] - M) %*% as.matrix(W)) / sqrt(sum(W ^ 2))
      
      total_r <- c(total_r, r)
      total_h <- c(total_h, h)
      total_index <- c(total_index, main_index)
    }
    t <- sigmoid(total_r) / (normalize(total_h) + 1e-8)
    max_count <- round(length(main_index) * beta)
    sorted <- order(t, decreasing = T)
    index <- sorted[1:max_count]
    count <- max_count
    can_index <- unique(total_index[index])
    while (length(can_index) < max_count) {
      count <- count + 1
      new_index <- total_index[sorted[count]]
      if (!new_index %in% can_index) {
        can_index <- c(can_index, new_index)
      }
    }
    candidate_index <- c(candidate_index, can_index)
  }
  candidate_index <- unique(candidate_index)
  return(candidate_index)
}


get_candidates_with_random <- function(train, ptree, tree_matrix, neighbors, beta) {
  sample_index <- sample(nrow(train), round(nrow(train) * beta / 2))
  
  leaf_nodes <- which(tree_matrix$var=='leaf')
  n_leaves <- length(leaf_nodes)
  
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  candidate_index <- sample_index
  for (i in 1:n_leaves) {
    main_node <- leaf_nodes[i]
    if (!main_node %in% neighbors$keys()) {
      next
    }
    neighbor_nodes <- neighbors[[main_node]]
    main_train <- train[which(ptree$where == main_node), ]
    main_majority_class <- getmode(main_train$y)
    main_index <- intersect(which(ptree$where == main_node), 
                            which(train$y == main_majority_class))
    x1 <- train[main_index, ]
    m1 <- as.numeric(apply(x1[, -1], 2, mean))
    main_index <- which(ptree$where == main_node)
    x1 <- train[main_index, ]
    r <- as.numeric(sqrt(apply((x1[, -1] - m1) ^ 2, 1, sum)))
    total_r <- c(); total_h <- c(); total_index <- c()
    for (j in 1:length(neighbor_nodes)) {
      neighbor_train <- train[which(ptree$where == neighbor_nodes[j]), ]
      neighbor_majority_class <- getmode(neighbor_train$y)
      neighbor_index <- intersect(which(ptree$where == neighbor_nodes[j]), 
                                  which(train$y == neighbor_majority_class))
      
      x2 <- train[neighbor_index, ]
      m2 <- as.numeric(apply(x2[, -1], 2, mean))
      W <- m1 - m2
      M <- (m1 + m2) / 2
      h <- abs(as.matrix(x1[, -1] - M) %*% as.matrix(W)) / sqrt(sum(W ^ 2))
      
      total_r <- c(total_r, r)
      total_h <- c(total_h, h)
      total_index <- c(total_index, main_index)
    }
    t <- sigmoid(total_r) / (normalize(total_h) + 1e-8)
    max_count <- round(length(main_index) * beta / 2)
    sorted <- order(t, decreasing = T)
    index <- sorted[1:max_count]
    count <- max_count
    can_index <- unique(total_index[index])
    while (length(can_index) < max_count) {
      count <- count + 1
      new_index <- total_index[sorted[count]]
      if (!new_index %in% can_index && !new_index %in% candidate_index) {
        can_index <- c(can_index, new_index)
      }
    }
    candidate_index <- c(candidate_index, can_index)
  }
  candidate_index <- unique(candidate_index)
  return(candidate_index)
}


get_relative_distance <- function(train, where, tree_matrix, neighbors, beta) {
  leaf_nodes <- which(tree_matrix$var=='leaf')
  n_leaves <- length(leaf_nodes)
  
  relative_distance <- rep(0, nrow(train))
  for (i in 1:n_leaves) {
    main_node <- leaf_nodes[i]
    main_index <- which(where == main_node)
    
    if (!main_node %in% neighbors$keys() || length(main_index) == 0) {next}
    
    neighbor_nodes <- neighbors[[main_node]]
    main_train <- train[main_index, ]
    neighbor_train <- train[which(where %in% neighbor_nodes), ]
    distance <- as.numeric(sqrt(apply(scale(main_train[, -1])^2, 1, sum)))
    relative_distance[main_index] <- distance
  }
  return(relative_distance)
}

remove_redundant_features <- function(data) {
  y_value <- unique(data$y)
  redundant_index <- c()
  for (i in 2:ncol(data)) {
    for (j in 1:length(y_value)) {
      if (length(unique(data[data$y==y_value[j], i])) == 1) {
        redundant_index <- c(redundant_index, i)
        break
      }
    }
  }
  if(length(redundant_index) != 0) {
    data <- data[, -redundant_index]
  }
  return(data)
}

# FLDSVM
FLD <- function(train, dtree, tree_matrix, neighbor, beta) {
  
  candidate_index <- c()
  for (i in 1:nrow(neighbor)) {
    first_node <- neighbor[i, 1]
    second_node <- neighbor[i, 2]
    first_index <- which(dtree$where == first_node)
    second_index <- which(dtree$where == second_node)
    cat_index <- c(first_index, second_index)
    x1 <- train[first_index, ]
    x2 <- train[second_index, ]
    x <- rbind(x1, x2)
    lda.m <- lda(y~.,data=remove_redundant_features(x))
    lda_pred <- predict(lda.m, x)$x
    count <- round(length(cat_index) * beta)
    index <- order(abs(lda_pred))[1:count]
    candidate_index <- c(candidate_index, cat_index[index])
  }
  candidate_index <- unique(candidate_index)
  return(candidate_index)
}

# Revised FLDSVM
RFLD <- function(train, dtree, tree_matrix, neighbors, beta) {
  leaf_nodes <- which(tree_matrix$var=='leaf')
  n_leaves <- length(leaf_nodes)
  
  candidate_index <- c()
  for (i in 1:n_leaves) {
    main_node <- leaf_nodes[i]
    main_index <- which(dtree$where == main_node)
    main_train <- train[main_index, ]
    if (!main_node %in% neighbors$keys()) {
      next
    }
    neighbor_nodes <- neighbors[[leaf_nodes[i]]]
    total_dist <- c(); total_index <- c();
    for (j in 1:length(neighbor_nodes)) {
      neighbor_train <- train[which(dtree$where == neighbor_nodes[j]), ]
      x <- rbind(main_train, neighbor_train)
      lda.m <- lda(y~.,data=remove_redundant_features(x))
      lda.pred <- predict(lda.m, x)$x[1:length(main_index)]
      total_dist <- c(total_dist, lda.pred)
      total_index <- c(total_index, main_index)
    }
    max_count <- round(length(main_index) * beta)
    sorted <- order(abs(total_dist))
    index <- sorted[1:max_count]
    count <- max_count
    can_index <- unique(total_index[index])
    while (length(can_index) < max_count) {
      count <- count + 1
      new_index <- total_index[sorted[count]]
      if (!new_index %in% can_index) {
        can_index <- c(can_index, new_index)
      }
    }
    candidate_index <- c(candidate_index, can_index)
  }
  candidate_index <- unique(candidate_index)
  return(candidate_index)
}

# Improved FLDSVM
IFLD_v1 <- function(train, dtree, tree_matrix, neighbors, beta) {
  leaf_nodes <- which(tree_matrix$var=='leaf')
  n_leaves <- length(leaf_nodes)
  
  candidate_index <- c()
  for (i in 1:n_leaves) {
    main_node <- leaf_nodes[i]
    if (!main_node %in% neighbors$keys()) {
      next
    }
    
    main_index <- which(dtree$where == main_node)
    neighbor_nodes <- neighbors[[leaf_nodes[i]]]
    neighbor_index <- which(dtree$where %in% neighbor_nodes)
    
    x1 <- train[main_index, ]
    x2 <- train[neighbor_index, ]
    
    x <- rbind(x1, x2)
    lda.m <- lda(y~.,data=remove_redundant_features(x))
    lda_pred <- predict(lda.m, x)$x[1:length(main_index)]
    count <- round(length(main_index) * beta)
    index <- order(abs(lda_pred))[1:count]
    candidate_index <- c(candidate_index, main_index[index])
  }
  candidate_index <- unique(candidate_index)
  return(candidate_index)
}


IFLD_v2 <- function(train, dtree, tree_matrix, neighbors, beta) {
  leaf_nodes <- which(tree_matrix$var=='leaf')
  n_leaves <- length(leaf_nodes)
  
  candidate_index <- c()
  for (i in 1:n_leaves) {
    main_node <- leaf_nodes[i]
    if (!main_node %in% neighbors$keys()) {
      next
    }
    neighbor_nodes <- neighbors[[leaf_nodes[i]]]
    main_index <- which(dtree$where == main_node)
    x1 <- train[main_index, ]
    total_pred <- c(); total_index <- c()
    for (j in 1:length(neighbor_nodes)) {
      neighbor_index <- which(dtree$where == neighbor_nodes[j])
      x2 <- train[neighbor_index, ]
      
      x <- rbind(x1, x2)
      lda_cls <- lda(y~.,data=remove_redundant_features(x))
      pred <- predict(lda_cls, x)$x[1:length(main_index)]
      
      total_pred <- c(total_pred, pred)
      total_index <- c(total_index, main_index)
    }
    
    max_count <- round(length(main_index) * beta)
    sorted <- order(total_pred, decreasing = T)
    index <- sorted[1:max_count]
    count <- max_count
    can_index <- unique(total_index[index])
    while (length(can_index) < max_count) {
      count <- count + 1
      new_index <- total_index[sorted[count]]
      if (!new_index %in% can_index) {
        can_index <- c(can_index, new_index)
      }
    }
    candidate_index <- c(candidate_index, can_index)
  }
  candidate_index <- unique(candidate_index)
  return(candidate_index)
}


# CVS implementation
get_cvs <- function(train, beta) {
  y_value <- unique(train$y)
  first_index <- which(train$y==y_value[1])
  second_index <- which(train$y==y_value[2])
  
  x1 <- train[first_index, ]
  x2 <- train[second_index, ]
  
  dist1 <- nn2(x2[,-1], x1[,-1], k=1, treetype = "kd") 
  dist2 <- nn2(x1[,-1], x2[,-1], k=1, treetype = "kd") 
  
  r1 <- (dist1[[2]] + 1)/(dist2[[2]][dist1[[1]], 1] + 1)
  r2 <- (dist2[[2]] + 1)/(dist1[[2]][dist2[[1]], 1] + 1)
  
  index1 <- order(r1)[1:round(length(r1)*beta)]
  index2 <- order(r2)[1:round(length(r2)*beta)]
  
  candidate_index <- c(first_index[index1], second_index[index2])
  return(candidate_index)
}

get_dtcvs <- function(train, dtree, tree_matrix, neighbors, beta) {
  leaf_nodes <- which(tree_matrix$var=='leaf')
  n_leaves <- length(leaf_nodes)
  
  candidate_index <- c()
  for (i in 1:n_leaves) {
    main_node <- leaf_nodes[i]
    main_index <- which(dtree$where == main_node)
    main_train <- train[main_index, ]
    if (!main_node %in% neighbors$keys()) {
      next
    }
    neighbor_nodes <- neighbors[[leaf_nodes[i]]]
    neighbor_train <- train[which(dtree$where %in% neighbor_nodes), ]
    dist1 <- nn2(neighbor_train[,-1], main_train[,-1], k=1, treetype = "kd") 
    dist2 <- nn2(main_train[,-1], neighbor_train[,-1], k=1, treetype = "kd")
    
    r1 <- (dist1[[2]] + 1)/(dist2[[2]][dist1[[1]], 1] + 1)
    r2 <- (dist2[[2]] + 1)/(dist1[[2]][dist2[[1]], 1] + 1)
    
    index <- order(r1)[1:round(length(r1)*beta)]
    candidate_index <- c(candidate_index, main_index[index])
  }
  candidate_index <- unique(candidate_index)
  return(candidate_index)
}

normalize <- function(x) {
  x <- (x - min(x)) / (max(x) - min(x) + 1e-8)
  return(x)
}

get_initial_selection <- function(train, p_0=0.1, tau_l=0.1, tau_u=0.25) {
  n <- nrow(train)
  y_label <- unique(train$y)
  positive_train <- train[train$y == y_label[1], ]
  negative_train <- train[train$y == y_label[2], ]
  n_pos <- nrow(positive_train)
  n_neg <- nrow(negative_train)
  p <- min(n_pos, n_neg) / n
  
  if (n_neg < n_pos) {
    if (p < tau_l) {
      positive_train <- positive_train[sample(n_pos, n_neg), ]
    }
    else if (p < tau_u) {
      positive_train <- positive_train[sample(n_pos, round(p * n_pos)), ]
      negative_train <- negative_train[sample(n_neg, round((1 - p) * n_neg)), ]
    } else {
      positive_train <- positive_train[sample(n_pos, round(p_0 * n_pos)), ]
      negative_train <- negative_train[sample(n_neg, round(p_0 * n_neg)), ]
    }
  } else {
    if (p < tau_l) {
      negative_train <- negative_train[sample(n_neg, n_pos), ]
    }
    else if (p < tau_u) {
      positive_train <- positive_train[sample(n_pos, round((1 - p) * n_pos)), ]
      negative_train <- negative_train[sample(n_neg, round(p * n_neg)), ]
    } else {
      positive_train <- positive_train[sample(n_pos, round(p_0 * n_pos)), ]
      negative_train <- negative_train[sample(n_neg, round(p_0 * n_neg)), ]
    }
  }
  sub_train <- rbind(positive_train, negative_train)
  return(sub_train)
}


TRSD_v1 <- function(train, mincut=5, mindev=1e-5, beta=0.1) {
  dtree <- tree(y~., data=train, mincut=mincut, mindev=mindev)
  
  tree_matrix <- data.frame(node=as.numeric(rownames(dtree$frame)),
                            var=paste(gsub("[Xx<+>]","",dtree$frame$var)),
                            threshold=as.numeric(gsub("[<+>]","",dtree$frame$splits[, 1])),
                            label=paste(dtree$frame$yval))
  
  neighbor <- get_pairwise_neighbors(train, tree_matrix)
  candidate_index <- get_candidates_v1(train, dtree, tree_matrix, neighbor, beta)
  return(candidate_index)
}

TRSD_v2 <- function(train, mincut=5, mindev=1e-5, beta=0.1) {
  dtree <- tree(y~., data=train, mincut=mincut, mindev=mindev)
  
  tree_matrix <- data.frame(node=as.numeric(rownames(dtree$frame)),
                            var=paste(gsub("[X<+>]","",dtree$frame$var)),
                            threshold=as.numeric(gsub("[<+>]","",dtree$frame$splits[, 1])),
                            label=paste(dtree$frame$yval))
  
  neighbors <- get_neighbors(train, tree_matrix)
  candidate_index <- get_candidates_v2(train, dtree, tree_matrix, neighbors, beta)
  return(candidate_index)
}

TRSD_v3 <- function(train, mincut=5, mindev=1e-5, beta=0.1) {
  dtree <- tree(y~., data=train, mincut=mincut, mindev=mindev)
  
  tree_matrix <- data.frame(node=as.numeric(rownames(dtree$frame)),
                            var=paste(gsub("[X<+>]","",dtree$frame$var)),
                            threshold=as.numeric(gsub("[<+>]","",dtree$frame$splits[, 1])),
                            label=paste(dtree$frame$yval))
  
  neighbors <- get_neighbors(train, tree_matrix)
  candidate_index <- get_candidates_v3(train, dtree, tree_matrix, neighbors, beta)
  return(candidate_index)
}

TRSD <- function(train, mincut=5, mindev=1e-5, beta=0.1) {
  dtree <- tree(y~., data=train, mincut=mincut, mindev=mindev)
  
  tree_matrix <- data.frame(node=as.numeric(rownames(dtree$frame)),
                            var=paste(gsub("[X<+>]","",dtree$frame$var)),
                            threshold=as.numeric(gsub("[<+>]","",dtree$frame$splits[, 1])),
                            label=paste(dtree$frame$yval))
  
  neighbors <- get_neighbors(train, tree_matrix)
  candidate_index <- get_candidates_v4(train, dtree, tree_matrix, neighbors, beta)
  return(candidate_index)
}

pTRSD <- function(train, mincut=5, mindev=1e-5, beta) {
  rtree <- rpart(y~., data=train, minsplit=mincut, method="class")
  ptree <- prune(rtree, cp=get_cp(rtree))
  ptree$frame$splits <- get_splits(ptree)
  tree_matrix <- data.frame(node=as.numeric(rownames(ptree$frame)),
                            var=paste(gsub("[X<+>]","",ptree$frame$var)),
                            threshold=as.numeric(ptree$frame$splits),
                            label=paste(ptree$frame$yval))

  neighbors <- get_neighbors(train, tree_matrix)
  candidate_index <- get_candidates_v4(train, ptree, tree_matrix, neighbors, beta)
  return(candidate_index)
}

pRTRSD <- function(train, mincut=5, mindev=1e-5, beta) {
  rtree <- rpart(y~., data=train, minsplit=mincut, method="class")
  ptree <- prune(rtree, cp=get_cp(rtree))
  ptree$frame$splits <- get_splits(ptree)
  tree_matrix <- data.frame(node=as.numeric(rownames(ptree$frame)),
                            var=paste(gsub("[X<+>]","",ptree$frame$var)),
                            threshold=as.numeric(ptree$frame$splits),
                            label=paste(ptree$frame$yval))
  
  neighbors <- get_neighbors(train, tree_matrix)
  candidate_index <- get_candidates_with_random(train, ptree, tree_matrix, neighbors, beta)
  return(candidate_index)
}

pTRSD_v2 <- function(train, mincut=5, mindev=1e-5, beta) {
  rtree <- rpart(y~., data=train, minsplit=mincut, method="class")
  ptree <- prune(rtree, cp=get_cp(rtree))
  ptree$frame$splits <- get_splits(ptree)
  tree_matrix <- data.frame(node=as.numeric(rownames(ptree$frame)),
                            var=paste(gsub("[X<+>]","",ptree$frame$var)),
                            threshold=as.numeric(ptree$frame$splits),
                            label=paste(ptree$frame$yval))
  
  neighbors <- get_neighbors(train, tree_matrix)
  candidate_index <- get_candidates_v5(train, ptree, tree_matrix, neighbors, beta)
  return(candidate_index)
}

rf_nodemap <- function(treemap, n_nodes) {
  map <- rep(0, n_nodes)
  map[1] <- 1
  for (i in 1:n_nodes) {
    if (treemap[i, 1] == 0) {next}
    map[treemap[i, 1]] <- map[i] * 2
    map[treemap[i, 2]] <- map[i] * 2 + 1
  }
  return(map)
}

assign_node <- function(train, tree_matrix) {
  where <- rep(0, nrow(train))
  for (i in 1:nrow(train)) {
    instance <- train[i, -1]
    node <- 1
    var <- tree_matrix$var[node]
    while(var != "leaf") {
      if (instance[var] <= tree_matrix$threshold[node]) {
        node <- tree_matrix$left[node]
      } else {
        node <- tree_matrix$right[node]
      }
      var <- tree_matrix$var[node]
    }
    where[i] <- node
  }
  return(where)
}


ensembleTRSD <- function(train, mincut=5, mindev=1e-5, beta) {
  ntree <- 10
  rf.m <- randomForest(y~., data=train, ntree=ntree, sampsize=ceiling(nrow(train)/10))
  relative_distance <- rep(0, nrow(train))
  for (i in 1:ntree) {
    n_nodes <- rf.m$forest$ndbigtree[i]
    nodeMap <- rf_nodemap(rf.m$forest$treemap[, , i], n_nodes)
    var <- rf.m$forest$bestvar[1:n_nodes, i]
    var[var == 0] <- "leaf"
    splits <- rf.m$forest$xbestsplit[1:n_nodes, i]
    y_val <- rf.m$forest$nodepred[1:n_nodes, i]
    tree_matrix <- data.frame(node=nodeMap,
                              var=var,
                              threshold=splits,
                              label=y_val,
                              left=rf.m$forest$treemap[, , i][1:n_nodes, 1],
                              right=rf.m$forest$treemap[, , i][1:n_nodes, 2])
    where <- assign_node(train, tree_matrix)
    neighbors <- get_neighbors(train, tree_matrix)
    distance <- get_relative_distance(train, where, tree_matrix, neighbors, beta)
    relative_distance <- relative_distance + distance
  }
  count <- round(nrow(train) * beta)
  candidate_index <- order(relative_distance, decreasing=TRUE)[1:count]
  return(candidate_index)
}

DTFLD <- function(train, mincut=5, mindev=1e-5, beta=0.1) {
  rtree <- rpart(y~., data=train, minsplit=mincut, method="class")
  ptree <- prune(rtree, cp=get_cp(rtree))
  ptree$frame$splits <- get_splits(ptree)
  tree_matrix <- data.frame(node=as.numeric(rownames(ptree$frame)),
                            var=paste(gsub("[X<+>]","",ptree$frame$var)),
                            threshold=as.numeric(ptree$frame$splits),
                            label=paste(ptree$frame$yval))
  
  neighbor <- get_pairwise_neighbors(train, tree_matrix)
  candidate_index <- FLD(train, ptree, tree_matrix, neighbor, beta)
  return(candidate_index)
}

pDTFLD <- function(train, mincut=5, mindev=1e-5, beta=0.1) {
  rtree <- rpart(y~., data=train, minsplit=mincut, method="class")
  ptree <- prune(rtree, cp=get_cp(rtree))
  ptree$frame$splits <- get_splits(ptree)
  tree_matrix <- data.frame(node=as.numeric(rownames(ptree$frame)),
                            var=paste(gsub("[X<+>]","",ptree$frame$var)),
                            threshold=as.numeric(ptree$frame$splits),
                            label=paste(ptree$frame$yval))
  
  neighbor <- get_neighbors(train, tree_matrix)
  candidate_index <- RFLD(train, ptree, tree_matrix, neighbor, beta)
  return(candidate_index)
}

DTIFLD_v1 <- function(train, mincut=5, mindev=1e-5, beta=0.1) {
  dtree <- tree(y~., data=train, mincut=mincut, mindev=mindev)
  
  tree_matrix <- data.frame(node=as.numeric(rownames(dtree$frame)),
                            var=paste(gsub("[X<+>]","",dtree$frame$var)),
                            threshold=as.numeric(gsub("[<+>]","",dtree$frame$splits[, 1])),
                            label=paste(dtree$frame$yval))
  
  neighbor <- get_neighbors(train, tree_matrix)
  candidate_index <- IFLD_v1(train, dtree, tree_matrix, neighbor, beta)
  return(candidate_index)
}

DTIFLD_v2 <- function(train, mincut=5, mindev=1e-5, beta=0.1) {
  dtree <- tree(y~., data=train, mincut=mincut, mindev=mindev)
  
  tree_matrix <- data.frame(node=as.numeric(rownames(dtree$frame)),
                            var=paste(gsub("[X<+>]","",dtree$frame$var)),
                            threshold=as.numeric(gsub("[<+>]","",dtree$frame$splits[, 1])),
                            label=paste(dtree$frame$yval))
  
  neighbor <- get_neighbors(train, tree_matrix)
  candidate_index <- IFLD_v2(train, dtree, tree_matrix, neighbor, beta)
  return(candidate_index)
}

CVS <- function(train, mincut=5, mindev=1e-5, beta=0.1) {
  candidate_index <- get_cvs(train, beta)
  return(candidate_index)
}

DTCVS <- function(train, mincut=5, mindev=1e-5, beta=0.1) {
  dtree <- tree(y~., data=train, mincut=mincut, mindev=mindev)
  
  tree_matrix <- data.frame(node=as.numeric(rownames(dtree$frame)),
                            var=paste(gsub("[X<+>]","",dtree$frame$var)),
                            threshold=as.numeric(gsub("[<+>]","",dtree$frame$splits[, 1])),
                            label=paste(dtree$frame$yval))
  
  neighbors <- get_neighbors(train, tree_matrix)
  candidate_index <- get_dtcvs(train, dtree, tree_matrix, neighbors, beta)
  return(candidate_index)
}

pDTCVS <- function(train, mincut=5, mindev=1e-5, beta=0.1) {
  rtree <- rpart(y~., data=train, minsplit=mincut, method="class")
  ptree <- prune(rtree, cp=get_cp(rtree))
  ptree$frame$splits <- get_splits(ptree)
  tree_matrix <- data.frame(node=as.numeric(rownames(ptree$frame)),
                            var=paste(gsub("[X<+>]","",ptree$frame$var)),
                            threshold=as.numeric(ptree$frame$splits),
                            label=paste(ptree$frame$yval))
  
  neighbors <- get_neighbors(train, tree_matrix)
  candidate_index <- get_dtcvs(train, ptree, tree_matrix, neighbors, beta)
  return(candidate_index)
}

remove_constant_features <- function(data) {
  redundant_index <- c()
  for (i in 2:ncol(data)) {
    if (length(unique((data[,i]))) == 1) {
      redundant_index <- c(redundant_index, i)
    }
  }
  if(length(redundant_index) != 0) {
    data <- data[, -redundant_index]
  }
  return(data)
}

SVDISC <- function(train, mincut=5, mindev=1e-5, beta=0.1) {
  sub_train <- get_initial_selection(train)
  sub_train <- remove_constant_features(sub_train)
  isvm <- svm(y~., data=sub_train)
  sv_index <- rep(0, nrow(sub_train))
  sv_index[isvm$index] <- 1
  sub_train$sv <- sv_index
  dtree <- tree(as.factor(sv)~.,data=sub_train[, -1], mincut=mincut)
  candidate_index <- which(predict(dtree, train[, -1], type="class") == 1)
  return(candidate_index)
}

SVIS <- function(train, mincut=5, mindev=1e-5, beta=0.1) {
  rf.m <- randomForest(y~., data=train, ntree=100)#, sampsize=ceiling(nrow(train)/10))
  margins <- as.numeric(abs(rf.m$votes[,1] - rf.m$votes[,2]))
  count <- round(length(margins) * beta)
  candidate_index <- order(margins)[1:count]
  return(candidate_index)
}

RVS <- function(train, mincut=5, mindev=1e-5, beta=0.1) {
  sample_index <- sample(nrow(train), round(nrow(train) * beta))
  return(sample_index)
}

FULL <- function(train, mincut=5, mindev=1e-5, beta=0.1) {
  all_index <- 1:nrow(train)
  return(all_index)
}
