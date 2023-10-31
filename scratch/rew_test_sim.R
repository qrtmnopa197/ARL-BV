library(MASS)
library(LDATS)
sim_rew_full <- function(id,rc_mu,rc_sigma,ra_mu,ra_sigma,r_ca,alpha_mu,alpha_sigma,resid){
  tr <- data.frame(id=rep(id,75),sub_index=rep(id,75),fA_img=c(rep(1,25),rep(2,25),rep(3,25)),
             fB_img=c(rep(4,25),rep(5,25),rep(6,25)),rew=sample(c(.05,0,-.05),75,replace=T),
             choice=rep(NA,75),valrat_z=rep(NA,75))
  
  #get parameter values
  cov_mat <-  matrix(c(ra_sigma^2, r_ca * ra_sigma * rc_sigma,
                       r_ca * ra_sigma * rc_sigma, rc_sigma^2), nrow = 2)
  r_w <- mvrnorm(1,c(ra_mu,rc_mu),cov_mat)
  alpha <- sigmoid(rnorm(1,alpha_mu,alpha_sigma))
  
  Q <- matrix(0, nrow = 75, ncol = 6)
  for(t in 1:75){
    #get choice
    probs <- softmax(c(r_w[2]*Q[t,tr$fA_img[t]], r_w[2]*Q[t,tr$fB_img[t]]))
    tr$choice[t] <- sample(c(1,2),1,prob=probs)
    
    #get affect rating
    tr$valrat_z[t] <- rnorm(1,r_w[1]*tr$rew[t],resid)
    
    #update Qs
    if(t < 75){
      Q[t+1,] = Q[t,]
      if(tr$choice[t] == 1){
        Q[t+1,tr$fA_img[t]] = Q[t,tr$fA_img[t]] + alpha*(tr$rew[t] - Q[t,tr$fA_img[t]])
      } else if(tr$choice[t] ==2){
        Q[t+1,tr$fB_img[t]] = Q[t,tr$fB_img[t]] + alpha*(tr$rew[t] - Q[t,tr$fB_img[t]])
      }
    }
  }
  return(tr)
}

#1:
samp1_list <- lapply(c(1:100),sim_rew_full,rc_mu=100,rc_sigma=1,ra_mu=15,ra_sigma=1,r_ca=-.8,
                     alpha_mu=-.05,alpha_sigma=1.7,resid=.3)
samp1_df <- do.call(rbind,samp1_list)
samp1_df <- add_probe_number(samp1_df,newcol="rat_number",val_col="valrat_z") 

#2
samp2_list <- lapply(c(1:100),sim_rew_full,rc_mu=30,rc_sigma=1,ra_mu=0,ra_sigma=1,r_ca=0,
                     alpha_mu=-.05,alpha_sigma=1.7,resid=.3)
samp2_df <- do.call(rbind,samp2_list)
samp2_df <- add_probe_number(samp2_df,newcol="rat_number",val_col="valrat_z") 

#3
samp3_list <- lapply(c(1:100),sim_rew_full,rc_mu=-100,rc_sigma=1,ra_mu=-15,ra_sigma=1,r_ca=.8,
                     alpha_mu=-.05,alpha_sigma=1.7,resid=.3)
samp3_df <- do.call(rbind,samp3_list)
samp3_df <- add_probe_number(samp3_df,newcol="rat_number",val_col="valrat_z") 

#4
samp4_list <- lapply(c(1:50),sim_rew_full,rc_mu=0,rc_sigma=100,ra_mu=0,ra_sigma=15,r_ca=-.9,
                     alpha_mu=-.05,alpha_sigma=1.7,resid=.2)
samp4_df <- do.call(rbind,samp4_list)
samp4_df <- add_probe_number(samp4_df,newcol="rat_number",val_col="valrat_z") 
samp4_df <- samp4_df %>% rename(choice_numeric = choice,out=rew)
samp4_df <- samp4_df %>% mutate(box_val=rep(1,3750),show_fres=rep(1,3750))




