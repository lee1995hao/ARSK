library(MASS)
library(mvtnorm)
library(RSKC)
library(magrittr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(trimcluster)
library(GPArotation)

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
  
  
  nperms <- 25
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
  
  
  B = 25
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



