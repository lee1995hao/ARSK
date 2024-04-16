library(MASS)
library(mvtnorm)
library(RSKC)
library(magrittr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(trimcluster)
library(GPArotation)
 
 ###data generator variable indenpent
 make_culster_data <- function(n, dimen, noise_v, number, plt = T){
   k_number <- length(n)
   center <- NULL
   zero_dim <- sample(dimen, noise_v)
   for(i in 1:dimen){
     center <- cbind(center, ifelse(runif(k_number) < 0.5, runif(k_number, min=-6, max=-3), runif(k_number, min=3, max=6)))
   }
   center[,zero_dim] <- 0
   Sigmar <- diag(1,nrow = dimen)
   data_without_oulier <- NULL
   for(i in 1:k_number){
     data_without_oulier <- rbind(data_without_oulier, cbind(mvrnorm(n = n[i], mu = center[i,], Sigma = Sigmar),i))
   }
   outliner_idx <- sample(length(data_without_oulier[,1]), number*length(data_without_oulier[,1]))



   r1 <- runif(100, min = 7,max = 13)
   r2 <- runif(100, min = -13,max = -7)
   r <-c(r1,r2)
   gap_matrix <- NULL
   for(i in 1:dimen){
     gap_matrix <- cbind(gap_matrix, sample(r, length(outliner_idx), replace = T))
   }
   data_without_oulier <- data.frame(data_without_oulier)

   data_without_oulier[,c(1:dimen)][outliner_idx,] <- data_without_oulier[,c(1:dimen)][outliner_idx,]+ gap_matrix
   realoutliner <- data_without_oulier[,c(1:dimen)][outliner_idx,]
   if(plt == T){
     plot(data_without_oulier[,c(1,2)], pch = 19, col = "navy",xlab = "X", ylab = "Y")
     points(realoutliner, col = 5,pch =4)
   }
   return(list(gap = gap_matrix ,dataset = data_without_oulier, center = center, outliner = realoutliner, outliner_idx = outliner_idx))

}

# data generator variable correlated
 make_culster_data <- function(n, dimen, noise_v, number, plt = T){
   k_number <- length(n)
   center <- NULL
   zero_dim <- sample(dimen, noise_v)
   for(i in 1:dimen){
     center <- cbind(center, ifelse(runif(k_number) < 0.5, runif(k_number, min=-6, max=-3), runif(k_number, min=3, max=6)))
   }
   center[,zero_dim] <- 0

   # groups covariance  structure
   size_grp <- n
   rho_grp_range <-list(min=0.1,max=0.5)
   n_grp <-  length(size_grp)
   ss_grp <- rep(1,n_grp)
   rrho_grp <- runif(n_grp,min=rho_grp_range$min,max=rho_grp_range$max)
   p_inf <- dimen


   s_grp <- list()
   for( j in 1: n_grp){
     S_g <- matrix(rrho_grp[j], ncol=p_inf, nrow=p_inf)
     diag(S_g) <- ss_grp[j]
     R <- Random.Start(p_inf)

     s_grp[[j]] <- R %*% S_g %*% t(R)
   }



   Sigmar <- s_grp
   data_without_oulier <- NULL
   for(i in 1:k_number){
     data_without_oulier <- rbind(data_without_oulier, cbind(mvrnorm(n = n[i], mu = center[i,], Sigma = Sigmar[[i]]),i))
   }
   outliner_idx <- sample(length(data_without_oulier[,1]), number*length(data_without_oulier[,1]))



   r1 <- runif(100, min = 7,max = 13)
   r2 <- runif(100, min = -13,max = -7)
   r <-c(r1,r2)
   gap_matrix <- NULL
   for(i in 1:dimen){
     gap_matrix <- cbind(gap_matrix, sample(r, length(outliner_idx), replace = T))
   }
   data_without_oulier <- data.frame(data_without_oulier)

   data_without_oulier[,c(1:dimen)][outliner_idx,] <- data_without_oulier[,c(1:dimen)][outliner_idx,]+ gap_matrix
   realoutliner <- data_without_oulier[,c(1:dimen)][outliner_idx,]
   if(plt == T){
     plot(data_without_oulier[,c(1,2)], pch = 19, col = "navy",xlab = "X", ylab = "Y")
     points(realoutliner, col = 5,pch =4)
   }
   return(list(gap = gap_matrix ,dataset = data_without_oulier, center = center, outliner = realoutliner, outliner_idx = outliner_idx))

}
