library(MASS)
library(mvtnorm)
library(RSKC)
library(magrittr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(trimcluster)
library(GPArotation)


source("data-generator.R")
source("gap-alternative.R")

#distance
Eudist<-function(x,y){
  distance<-sqrt(sum((x-y)^2))
  return (distance)
}


#thresholding function
soft_threshold <- function(x, c = c) {
  sign(x) * pmax(abs(x) - c, 0)
}


scadThreshold <- function(entry, lambda, a) {
  e1 <- abs(entry) <= 2 * lambda
  e2 <- abs(entry) > 2 * lambda & abs(entry) <= a * lambda
  
  entry[e1] <- ifelse(
    abs(entry[e1]) - lambda > 0, sign(entry[e1]) * (abs(entry[e1]) - lambda), 0
  )
  
  entry[e2] <- ((a - 1) * entry[e2] - sign(entry[e2]) * a * lambda) / (a - 2)
  
  return(entry)
}


##function for OKM
theta_ipod_kmean_involve_lasso <- function(dataset, lambda, k, w){
  total_mean <- apply(dataset, 2 ,mean)
  all_dis <- c()
  for(i in 1:length(dataset[,1])){
    all_dis[i] <- Eudist(dataset[i,], total_mean)
  }
  
  a <- length(dataset[,1])*0.2
  most_far_dis <- sort(all_dis ,decreasing = T)[1:a] 
  most_far_pos <- which(all_dis%in%most_far_dis)
  
  E_ini <- matrix(0, ncol = length(dataset[1,]), nrow = length(dataset[,1]))
  E_ini <- data.frame(E_ini)
  E_n <- data.frame(matrix(0,ncol = length(dataset[1,]), nrow = length(dataset[,1])))
  
  E_ini[most_far_pos,] <- dataset[most_far_pos,]
  
  E <- E_ini
  
  ##soft thresholding of E
  new_E <- function(x, m, lambda){
    return((x - m)*max(0,1 - ((sqrt(sum((x - m)^2)))^-1)*lambda))
  }
  
  
  #scad
  scadd <- function(x_values, lambda1 ,a){
    return(x_values - scadThreshold(x_values, lambda1 ,a))
  }
  
  new_E <- function(x,m, lambda, gamma = 3.7) {
    z <- x - m
    z_norm <- sqrt(sum(z^2))
    return(z*(1-(scadd(z_norm, lambda1 ,a = gamma)/z_norm)))
  }
  
  
  
  ##initial parameter 
  kmean_dataset <- dataset - E
  kmean_dataset <- sweep(kmean_dataset, 2, w, "*")
  kmean_res <- kmeans(kmean_dataset, k)
  mu_new <- kmean_res$center
  mu_old <- mu_new
  iter = 0
  ##iteration part
  repeat{
    cluster_number <- sort(unique(kmean_res$cluster))
    G_res <- kmean_res$cluster
    G_conbin_number <- list()
    
    
    for(i in cluster_number){
      G_conbin_number[[i]] <- which(G_res == i)
    }
    
    for(i in 1:length(G_conbin_number)){
      for(j in G_conbin_number[[i]]){
        E_n[j,] <- new_E(sweep(dataset, 2, w, "*")[j,], mu_new[i,], lambda = lambda)
      }
    }
    
    
    kmean_dataset <- sweep(dataset, 2, w, "*") - E_n
    
    kmean_res <- kmeans(kmean_dataset, k)
    cluster_res_final <- kmean_res$cluster
    mu_new <- kmean_res$center
    if(Eudist(mu_new, mu_old) < 0.01 || iter == 200)break
    mu_old <- mu_new
    E <- E_n
    iter = iter + 1
    
    
  }
  
  
  
  ##judge which point is outliers
  E_2 <- diag(as.matrix(E)%*%t(as.matrix(E)))
  judge_outliner_idx <- which(E_2 != 0)
  
  
  # return part
  return(list(okm_iter = iter ,E = E, k = k, mu = mu_old, outlier_idx = judge_outliner_idx, cluster_res_final = cluster_res_final))
}






#variable selection
find_a <- function(w, E_res, k , c, dataset, cl_rest_a){
  tata <- c()
  for(i in 1:length(w)){
    tata[i] <- ifelse(w[i] == 0, 1, w[i])
  }
  
  E_return <- sweep(E_res, 2, tata,"/")
  
  
  new_dataset <- dataset - E_return
  data_with_predict <- as.data.frame(cbind(new_dataset, cl_rest_a))
  predict_cluster_category <- sort(unique(data_with_predict[,ncol(data_with_predict)]))
  
  ##all difference 
  all_mean <- apply(data_with_predict[,-ncol(data_with_predict)], 2, mean)
  
  all_distence <- (sweep(data_with_predict[,-ncol(data_with_predict)], 2, all_mean, "-"))^2
  
  ##each group difference
  each_dim_matrix <- matrix(data = NA, nrow = k, ncol = ncol(dataset))
  
  for(z in 1:length(predict_cluster_category)){
    group_individual_idx <-which(data_with_predict[,ncol(data_with_predict)] == z)
    each_group <- data_with_predict[,-ncol(data_with_predict)][group_individual_idx,]
    
    dim_dis_e <- c()
    for(p in 1:ncol(each_group)){
      distance_squ_all <- matrix(NA, ncol = length(each_group[,p]), nrow = length(each_group[,p]))
      for(i in 1:length(each_group[,p])){
        for(j in 1:length(each_group[,p])){
          distance_squ_all[i,j] <- (each_group[,p][i] - each_group[,p][j])^2
        }
      }
      dim_dis_e[p] <- sum(apply(distance_squ_all,2,mean))/2
    }
    each_dim_matrix[z,] <- dim_dis_e
  }
  
  
  
  
  
  
  
  a_j <- apply(all_distence, 2, sum) - apply(each_dim_matrix, 2, sum)
  
  new_a_j <- soft_threshold(a_j, c = c)
  new_a_j <- scadThreshold(a_j, c ,3.7)
  w <- new_a_j/norm(new_a_j, type = "2")
  
  return(w)
}



##algorithm
RSKOD <- function(k , c, lambda, dataset, w){
  w_path <- NULL
  iter_1 <- 0
  j_w_ser <- c()
  repeat{
    cat("==s")
    test_res <- theta_ipod_kmean_involve_lasso(dataset = dataset, lambda = lambda, k = k, w = w)
    w_a <- find_a(w = w, E_res =  test_res$E, k = k, c = c, dataset = dataset, cl_rest_a = test_res$cluster_res_final)
    judge_w <- sum(abs(w_a - w))/sum(abs(w))
    cat(judge_w)
    if(judge_w < 0.01 || iter_1 == 50)break
    w <- w_a
    j_w_ser <- c(judge_w, j_w_ser)
    w_path <- rbind(w, w_path)
    iter_1 = iter_1 + 1
  }
  
  w_f <- w_path[1,]
  
  return(list(w_path = w_path, okm_it = test_res$okm_iter,t_iter = iter_1, w_f = w_f, spare_okm_cluster = test_res$cluster_res_final ,
              spare_okm_E = test_res$E, spare_okm_outlier_idx = test_res$outlier_idx))
}




