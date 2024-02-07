library(MASS)
library(mvtnorm)
library(RSKC)
library(magrittr)
library(dplyr)
library(ggplot2)
library(tidyr)



# ###data generator
# make_culster_data <- function(n, dimen, noise_v, number, plt = T){
#   k_number <- length(n)
#   center <- NULL
#   zero_dim <- sample(dimen, noise_v)
#   for(i in 1:dimen){
#     center <- cbind(center, ifelse(runif(k_number) < 0.5, runif(k_number, min=-6, max=-3), runif(k_number, min=3, max=6)))
#   }
#   center[,zero_dim] <- 0
#   Sigmar <- diag(1,nrow = dimen)
#   data_without_oulier <- NULL
#   for(i in 1:k_number){
#     data_without_oulier <- rbind(data_without_oulier, cbind(mvrnorm(n = n[i], mu = center[i,], Sigma = Sigmar),i))
#   }
#   outliner_idx <- sample(length(data_without_oulier[,1]), number*length(data_without_oulier[,1]))
#   
#   
#   
#   r1 <- runif(100, min = 7,max = 13)
#   r2 <- runif(100, min = -13,max = -7)
#   r <-c(r1,r2)
#   gap_matrix <- NULL
#   for(i in 1:dimen){
#     gap_matrix <- cbind(gap_matrix, sample(r, length(outliner_idx), replace = T))
#   }
#   data_without_oulier <- data.frame(data_without_oulier)
#   
#   data_without_oulier[,c(1:dimen)][outliner_idx,] <- data_without_oulier[,c(1:dimen)][outliner_idx,]+ gap_matrix
#   realoutliner <- data_without_oulier[,c(1:dimen)][outliner_idx,]
#   if(plt == T){
#     plot(data_without_oulier[,c(1,2)], pch = 19, col = "navy",xlab = "X", ylab = "Y")
#     points(realoutliner, col = 5,pch =4)
#   }
#   return(list(gap = gap_matrix ,dataset = data_without_oulier, center = center, outliner = realoutliner, outliner_idx = outliner_idx))
#   
# }





#distance
Eudist<-function(x,y){
  distance<-sqrt(sum((x-y)^2))
  return (distance)
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
  
  ##hard thresholding of E
  new_E <- function(x, m, lambda){
    return((x - m)*max(0,1 - ((sqrt(sum((x - m)^2)))^-1)*lambda))
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




###lasso part
soft_threshold <- function(x, c = c) {
  sign(x) * pmax(abs(x) - c, 0)
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
  
  w <- new_a_j/norm(new_a_j, type = "2")
  
  return(w)
}



##algorithm_1
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




#algorithm_2_Gap_altern
GetWCSS <- function(x, Cs, ws=NULL){
  wcss.perfeature <- numeric(ncol(x))
  for(k in unique(Cs)){
    whichers <- (Cs==k)
    if(sum(whichers)>1) wcss.perfeature <- wcss.perfeature + apply(scale(x[whichers,],center=TRUE, scale=FALSE)^2, 2, sum)
  }
  bcss.perfeature <- apply(scale(x, center=TRUE, scale=FALSE)^2, 2, sum)-wcss.perfeature
  if(!is.null(ws)) return(list(wcss.perfeature=wcss.perfeature, wcss=sum(wcss.perfeature), wcss.ws=sum(wcss.perfeature*ws),
                               bcss.perfeature=bcss.perfeature))
  if(is.null(ws)) return(list(wcss.perfeature=wcss.perfeature, wcss=sum(wcss.perfeature), bcss.perfeature=bcss.perfeature))
}



cal_Gap <- function(dataset, k, c, lambda){
  w_init <- rep(1/sqrt(ncol(dataset)), ncol(dataset))

  sample_run <- RSKOD(k = k , c = c, lambda = lambda, dataset = dataset, w = w_init)


  x_star <- dataset - sample_run$spare_okm_E


  o_all <- sum(sample_run$w_f*GetWCSS(x_star, sample_run$spare_okm_cluster)$bcss.perfeature)


  nperms <- 45
  permx <- list()
  nnonzerows <- NULL
  for(i in 1:nperms){
    permx[[i]] <- matrix(NA, nrow=nrow(dataset), ncol=ncol(dataset))
    for(j in 1:ncol(dataset)){
      permx[[i]][,j] <- sample(dataset[,j])
    }
  }


  permtots <- c()
  for(K in 1:nperms){
    perm.out <- RSKOD(k = k , c = c, lambda = lambda, dataset = permx[[K]], w = w_init)
    perm.bcss <- GetWCSS(dataset - sample_run$spare_okm_E ,sample_run$spare_okm_cluster)$bcss.perfeature
    permtots[K] <- sum(perm.out$w_f*perm.bcss)
  }


  Gap <- ifelse((sample_run$t_iter >= 15)|(sample_run$okm_it >= 50), NA, log(o_all)- mean(log(permtots)))

  return(list(Gap = Gap, k = k, c = c, lambda = lambda, B_model = sample_run)) ## choose the s corresponding largest gap
}




gap_statistic_lam <- function(dataset, k ,lambda, c) {
  w_init <- rep(1/sqrt(ncol(dataset)), ncol(dataset))
  
  sample_run <- RSKOD(k = k , c = c, lambda = lambda, dataset = dataset, w = w_init)
  
  
  x_star <- dataset - sample_run$spare_okm_E
  x_star <- sweep(x_star, 2, sample_run$w_f,"*")
  
  
  B = 45
  wk <- sum(kmeans(x_star, centers = k)$withinss)
  
  rand_wks <- numeric(B)
  for(b in 1:B) {
    rand_data <- matrix(rnorm(nrow(x_star) * ncol(x_star)), nrow = nrow(x_star))
    rand_wks[b] <- sum(kmeans(rand_data, centers = k)$withinss)
  }
  
  gap <- log(mean(rand_wks)) - log(wk)
  Gap <- ifelse((sample_run$t_iter >= 15)|(sample_run$okm_it >= 50), NA, gap)
  
  return(c(gap = Gap))
}



c_var <- sort(2^seq(log2(10), log2(200), length = 10))
lam_var <- c(1.1:10.1)


##Alternating Optimization
select_hp <- function(c_var, lam_var, dataset,lambda_in){
  gap_serise_c <-c()
  gap_serise_l <-c()
  for(j in 1:length(c_var)){
    gap_serise_c[j] <- cal_Gap(dataset = dataset, k = 3, c = c_var[j], lambda = lambda_in)$Gap
  }
  
  b_c <- c_var[which.max(gap_serise_c)]
  
  for(j in 1:length(lam_var)){
    gap_serise_l[j] <- gap_statistic_lam(dataset = dataset, k = 3, c = b_c, lambda = lam_var[j])
  }
  
  b_lamda <- lam_var[which.min(gap_serise_l)]
  
  return(list(b_lamda = b_lamda, b_c = b_c,gap_serise_l,gap_serise_c))
  
}

#total function of algorithm_2
select_hp(c_var = c_var, lam_var = lam_var, dataset = dataset, lambda_in = 3.1)








#select_real_data <- NULL
#for(j in 1:100){ 
#  set.seed(12 + j*34)
#  data <- make_culster_data(n = c(50,50,50), dimen = 50, noise_v = 45, number = 0.1)
#  dataset <- data$dataset[,-ncol(data$dataset)]
#  w_init <- rep(1/sqrt(ncol(dataset)), ncol(dataset))


#  select_res <- select_hp(c_var = c_var, lam_var = lam_var, dataset = dataset, lambda_in = 3.1)


#  sample_run_1 <- RSKOD(k = 3, c = select_res$b_c, lambda = select_res$b_lamda, dataset = dataset, w = w_init)

#  real_idx <- data$dataset["i"]
# predicted_labels_1 <- as.factor(sample_run_1$spare_okm_cluster[-data$outliner_idx])
# actual_labels <- as.factor(real_idx[-data$outliner_idx,])
# CER_1 <- CER(predicted_labels_1,actual_labels)

# outlier_model_1 <- sample_run_1$spare_okm_outlier_idx
#number_o_1 <- length(outlier_model_1)
#o_TRP_1 <- length(which(data$outliner_idx %in% outlier_model_1))/length(data$outliner_idx)


#weigth_each_feature_1 <- as.vector(sample_run_1$w_f)
#o_feature_1 <- length(which(weigth_each_feature_1 != 0))


#real_z_v <- which(apply(data$center, 2, mean) != 0)
#pre_z_v_1 <- which(weigth_each_feature_1 != 0)
#v_TPR_1 <- length(which(real_z_v%in%pre_z_v_1))
#select_real_data <- rbind(select_real_data ,c(CER_1, number_o_1, o_TRP_1, o_feature_1, v_TPR_1))
#}



#data.var <- as.data.frame(expand.grid(c_var, lam_var))
#names(data.var) <- c("c", "lam")



######grip
#for(j in 1:nrow(data.var)){
#  cat("========== case", j)
#  aa <- data.var[j, "c"]
#  bb <- data.var[j, "lam"]

#  gap_serise[j] <- cal_Gap(dataset = dataset, k = 3, c = aa, lambda = bb)$Gap
#}





