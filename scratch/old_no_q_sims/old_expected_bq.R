# This script simulates and analyzes data from the Study 3 task, focusing on trials in which the outcome of the 
# chosen cue was received. The script simulates choices and valence ratings under simplifying assumptions, 
# estimates the effects of reward and affect associations on these choices using methods similar to those used
# in the manuscript, and computes expected effect estimates based on these results. 

library(tidyverse)
library(mclust)
library(sigmoid)
library(future)
library(future.apply)
library(lme4)


# Define functions

# Adds columns to a data frame with values of cue associations (e.g., reward associations, affect associations)
# on each trial
add_assoc <- function(data,out_col,assoc_name,lrn_rate,choice_col="choice",
                      assoc_num=2,init_val=0, for_rate=1){
  # Initialize df with association values
  assoc <- data.frame(matrix(init_val, ncol = assoc_num, nrow = nrow(data)))
  colnames(assoc) <- paste0(assoc_name, 1:assoc_num)
  
  for(r in 1:(nrow(data)-1)){
    assoc[r+1,] <- assoc[r,]*for_rate # Default to same value as above
    # For the chosen option, update the association
    assoc[r+1,data[r,choice_col]] <- assoc[r,data[r,choice_col]] + 
                                     lrn_rate*(data[r,out_col] - assoc[r,data[r,choice_col]])
  }
  cbind(data,assoc)
}

# Adds columns to a data frame with values of cue associations (e.g., reward associations, affect associations)
# on each trial. In contrast to add_assoc, this allows for updating of unchosen options.
add_assoc_counter <- function(data,out_cols,assoc_name,lrn_rate,init_val=0){
  # Initialize df with association values
  outs <- data[out_cols]
  assoc <- data.frame(matrix(init_val, ncol = ncol(outs), nrow = nrow(data)))
  colnames(assoc) <- paste0(assoc_name, 1:ncol(outs))
  
  for(r in 1:(nrow(outs)-1)){
    assoc[r+1,] <- assoc[r,] + lrn_rate*(outs[r,] - assoc[r,])
  }
  cbind(data,assoc)
}

# The main function for estimating effects of Beta_Q and Beta_A based on simulated data
sim_bq_ba <- function(x,n_dsets, n_subs=323, n_blocks=3, n_tr_bl=24, alpha_m = -0.46, alpha_sd = 4.98, phi_m = 1.08, 
                      phi_sd = 2.86, frate_m = 0.264, frate_sd = .74, beta_0_m = 0.00, beta_0_sd = 0.06, 
                      beta_rew_m = 11.45, beta_rew_sd = 5.55, beta_bv_m = -0.07, beta_bv_sd = 0.08, beta_pr_m = 0.12, beta_pr_sd = 0.12, 
                      resid_sigma = 0.47, beta_a_m = 1.73, beta_a_sd = 0.94, beta_c_m = 0.047,beta_c_sd = 1.03, 
                      beta_q_m = 0, beta_q_sd = 0){
  # Estimated effects of reward and affect associations on choice
  est_beta_q <- c()
  est_beta_a <- c()
  sd_qdiff <- c()
  sd_adiff <- c()
  for(d in 1:n_dsets){
    dset <- data.frame(sub=c(),block = c(), trial = c(), choice = c(),val_rat = c(),
                       prev_rat = c(),rew = c(), bv = c())
    for(s in 1:n_subs){
      sub_data <- data.frame(sub=c(),block = c(), trial = c(), choice = c(),val_rat = c(),
                             prev_rat = c(),rew = c(), bv = c())
      # Get parameter values for the subject
      beta_rew <- rnorm(1,beta_rew_m,beta_rew_sd)
      beta_bv <- rnorm(1,beta_bv_m,beta_bv_sd)
      beta_pr <-  rnorm(1,beta_pr_m,beta_pr_sd)
      beta_0 <-  rnorm(1,beta_0_m,beta_0_sd)
      beta_q <- rnorm(1,beta_q_m,beta_q_sd)
      beta_a <- rnorm(1,beta_a_m,beta_a_sd)
      beta_c <- rnorm(1,beta_c_m,beta_c_sd)
      untrans_alpha <- rnorm(1,alpha_m,alpha_sd)
      alpha <- sigmoid(untrans_alpha)
      untrans_phi <- rnorm(1,phi_m,phi_sd)
      phi <- sigmoid(untrans_phi)
      untrans_frate <- rnorm(1,frate_m,frate_sd)
      frate <- sigmoid(untrans_frate)
      for(b in 1:n_blocks){
        block_data <- data.frame(matrix(ncol = 8, nrow = 24))
        colnames(block_data) <- c("sub","block","trial","choice","val_rat","prev_rat","rew","bv")
        q <- c(0,0)
        a <- c(0,0)
        c <- c(0,0)
        prev_rat <- 0
        for(t in 1:n_tr_bl){
          bv <- sample(c(-2,2),size=1) # Get box value for the trial
          choice_probs <- softmax(c(beta_q*q[1]+beta_a*a[1]+beta_c*c[1],beta_q*q[2]+beta_a*a[2]+beta_c*c[2]))
          choice <- sample(c(1,2),size=1,prob=choice_probs) # Get choice
          rew <- sample(c(-0.05,0,0.05),size=1) # Get reward
          true_val <- beta_0 + beta_rew*rew + beta_bv*bv # Get valence
          val_rat <- rnorm(1,true_val + beta_pr*prev_rat,sd=resid_sigma) # Get rated valence
          # Update Q, A, and C values
          q[choice] <- q[choice] + alpha*(rew - q[choice])
          a[choice] <- a[choice] + alpha*(true_val - a[choice])
          c[choice] <- c[choice] + phi*(1 - c[choice])
          if(choice == 1){
            c[2] <- c[2] + phi*(0 - c[2])
          } else{
            c[1] <- c[1] + phi*(0 - c[1])
          }
          # Fill out block_data
          block_data$choice[t] <- choice
          block_data$val_rat[t] <- val_rat
          block_data$prev_rat[t] <- prev_rat
          block_data$rew[t] <- rew
          block_data$bv[t] <- bv
          
          prev_rat <- val_rat # Set previous valence rating for the next trial
        }
        block_data$block <- b
        block_data$trial <- c(1:24)
        block_data$sub <- s
        sub_data <- rbind(sub_data,block_data)
      }
      valrat_fit <- lm(val_rat ~ rew + bv + prev_rat,sub_data) # Estimate effects of rew and bv on valence
      sub_data <- mutate(sub_data, mod_val = valrat_fit$coefficients[1] + valrat_fit$coefficients[2]*rew 
                         + valrat_fit$coefficients[3]*bv) # Model valence based on fit
      sub_data <- sub_data %>%
        mutate(ch1 = ifelse(choice == 1,1,0)) %>%
        mutate(ch2 = ifelse(choice == 2,1,0))
      # Compute Q, A, and C values, given mod_val, and add to the dataaset
      subd_list <- split(sub_data, sub_data$block)
      subd_list <- lapply(subd_list,add_assoc,out_col="rew",assoc_name="q",lrn_rate=alpha,for_rate=frate)
      subd_list <- lapply(subd_list,add_assoc,out_col="mod_val",assoc_name="a",lrn_rate=alpha,for_rate=frate)
      subd_list <- lapply(subd_list,add_assoc_counter,out_cols=c("ch1","ch2"),assoc_name="c",lrn_rate=phi)
      
      # Add to dataset
      sub_data <- bind_rows(subd_list)
      dset <- rbind(dset,sub_data)
    }
    # Calculate and log the effects of Q and A in this dataset
    dset <- dset %>% 
      mutate(q_diff = q2 - q1) %>%
      mutate(a_diff = a2 - a1) %>%
      mutate(c_diff = c2 - c1) %>%
      mutate(choice_binary = choice - 1)
    choice_fit <- glm(choice_binary ~ q_diff + a_diff + c_diff, dset,family=binomial)
    est_beta_q <- c(est_beta_q,choice_fit$coefficients[2])
    est_beta_a <- c(est_beta_a,choice_fit$coefficients[3])
    
    
    sd_qdiff <- c(sd_qdiff,sd(dset$q_diff))
    sd_adiff <- c(sd_adiff,sd(dset$a_diff))
  }
  return(data.frame(est_beta_q = est_beta_q, est_beta_a = est_beta_a, 
                    sd_qdiff = sd_qdiff, sd_adiff = sd_adiff))
}


plan(multisession)

est_eff_list <- future_lapply(1:10, sim_bq_ba, n_dsets = 100, future.seed=TRUE)
est_eff_df <- do.call(rbind,est_eff_list)

est_eff_df <- est_eff_df %>% 
                mutate(ba_sc = est_beta_a*sd_adiff) %>%
                mutate(bq_sc = est_beta_q*sd_qdiff)

#summary(lm(est_eff_df$est_beta_a ~ 1))
summary(lm(est_eff_df$est_beta_q ~ 1))
#summary(lm(est_eff_df$ba_sc ~ 1))
summary(lm(est_eff_df$bq_sc ~ 1))

hist(est_eff_df$est_bq_sc)
