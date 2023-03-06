simulate_s222 <- function(n_samp,ns,nb,n_pb,n_tp,pair_type="fixed",mb_size=4, alpha_r_mean,alpha_r_sd=0,beta_r_mean,beta_r_sd,
                          forget_alpha_r="match",forget_alpha_a="match", forget_r_mean=NULL,forget_r_sd=NULL,alpha_a_mean,
                          alpha_a_sd=0,beta_a_mean,beta_a_sd,forget_a_mean=NULL,forget_a_sd=NULL,phi_mean,phi_sd,tau_mean,tau_sd,
                          w_r_mean,w_r_sd=0,w_fo_mean,w_fo_sd=0,w_q_mean,w_q_sd=0,w_opv_mean,w_opv_sd=0,v_resid_var,
                          prob_r,prob_fo,r_seq=NULL,fo_seq=NULL,return_plot_dfs=FALSE){
  
  library(tidyverse)
  library(gridExtra)
  library(cowplot)
  library(matrixStats)
  
  return_list <- list() #initialize return list
  lcg_list <- list() #will be used later, to store learning curve grids temporarily before they're transferred
                     #to the return list
  if(return_plot_dfs){
    plot_df_list <- list() #will be used later to store plot dataframes before they're transferred to the return list
  }
  
  d_cor_df <- as.data.frame(matrix(ncol=2,nrow=n_samp))
  d_cor_names <- c("dQ_dA","dQ_dF")
  names(d_cor_df) <- d_cor_names #declare df for the corrrelations between dQ values - one row per sample, with
  #cols for each correlation
  pred_cor_df <- as.data.frame(matrix(ncol=6,nrow=n_samp))
  pred_cor_names <- c("r_fo","r_q","r_opv","fo_q","fo_opv","q_opv")
  names(pred_cor_df) <- pred_cor_names #declare df for the correlations between affect predictors
  
  samp_ds <- array(dim=c(n_tp,n_samp,2,2,n_pb))  #Create array for learning curves from each sample.
  #Think of this as a collection of 2D tables. Each table is specific
  #to a certain SD correction (1 is uncorrected), dQ/dA value (1 is dQ),
  #and pair, for a total of 8 tables. Each table has one row per trial (for each
  #trial per pair) and one column per sample. This will be filled out with
  #learning curves for every sample.
  
  for(samp_num in 1:n_samp){
    samp_df <- sim_samp_s222(ns,nb,n_pb,n_tp,pair_type,mb_size,
                             alpha_r_mean,alpha_r_sd,beta_r_mean,beta_r_sd,
                             forget_alpha_r,forget_alpha_a,
                             forget_r_mean,forget_r_sd,alpha_a_mean,
                             alpha_a_sd,beta_a_mean,beta_a_sd,forget_a_mean,
                             forget_a_sd,phi_mean,phi_sd,tau_mean,tau_sd,
                             w_r_mean,w_r_sd,w_fo_mean,w_fo_sd,w_q_mean,
                             w_q_sd,w_opv_mean,w_opv_sd,v_resid_var,
                             prob_r,prob_fo,r_seq,fo_seq)
    samp_df <- samp_df %>% mutate("samp_num" = samp_num)
    # if(samp_num == 1){
    #   sim_df <- samp_df #create trial-level dataset for the simulation
    # } else{
    #   sim_df <- rbind(sim_df,samp_df)
    # }
    #Calculating mean subject-level correlations for all relevant variable pairs
    d_cor_df$dQ_dA[samp_num] <- mean(by(samp_df,samp_df$subject,function(x) cor(x$dQ,x$dA)))
    d_cor_df$dQ_dF[samp_num] <- mean(by(samp_df,samp_df$subject,function(x) cor(x$dQ,x$dF)))
    pred_cor_df$r_fo[samp_num] <- mean(by(samp_df,samp_df$subject,function(x) cor(x$reward,x$fake_out)))
    pred_cor_df$r_q[samp_num] <- mean(by(samp_df,samp_df$subject,function(x) cor(x$reward,x$Q_chosen)))
    pred_cor_df$r_opv[samp_num] <- mean(by(samp_df,samp_df$subject,function(x) cor(x$reward,x$onpol_val)))
    pred_cor_df$fo_q[samp_num] <- mean(by(samp_df,samp_df$subject,function(x) cor(x$fake_out,x$Q_chosen)))
    pred_cor_df$fo_opv[samp_num] <- mean(by(samp_df,samp_df$subject,function(x) cor(x$fake_out,x$onpol_val)))
    pred_cor_df$q_opv[samp_num] <- mean(by(samp_df,samp_df$subject,function(x) cor(x$Q_chosen,x$onpol_val)))
    
    #loop over pairs...
    for(p in 1:n_pb){
      cur_fA <- p*2-1 #get the fA value of the current pair
      pair_df <- filter(samp_df,fA == cur_fA) #filter dataset by pair
      #loop over d values...
      for(d in c("dQ","dA")){
        if(d == "dQ"){
          dnum <- 1
        } else if(d=="dA"){
          dnum <- 2
        }
        d_df <- as.data.frame(matrix(nrow=n_tp,ncol=ns*nb)) #create df for subject-level learning curves: trials
        #per pair as rows, subjects*nb as columns. The reason to multiply the number of subjects is you'll
        #have one learning curve per pair per block per subject.
        d_sd_df <- d_df #duplicate df - this one for SD-corrected values
        #for each subject...
        for(s in 1:ns){
          sub_pair_df <- filter(pair_df,subject==s) #grab subject's trials (for the current pair)
          for(b in 1:nb){
            block_pair_df <- filter(sub_pair_df,block==1) #get the data for just one block
            d_df[(s-1)*nb+b] <- block_pair_df[d] #drop the d column in the subject-level learning curve df
            
            sub_df <- filter(samp_df,subject==s) #grab all of the subject's trials
            #figure out the SD of the relevant outcome for this subject
            if(d == "dQ"){
              sd <- sd(sub_df$reward)
            } else if(d == "dA"){
              sd <- sd(sub_df$valence)
            }
            d_sd_df[(s-1)*nb+b] <- block_pair_df[d]/sd #correct the d column by the subject's overall SD, and drop it in the corrected df
          }
        }
        samp_ds[,samp_num,1,dnum,p] <- rowMeans(d_df) #calculate the group average learning curve, and put in the
        #ds array column for this sample, correction (non-SD),
        #d value, and pair number
        samp_ds[,samp_num,2,dnum,p] <- rowMeans(d_sd_df) #diddo for the corrected df
      }
    }
  }
  
  d_hists <- lapply(d_cor_names, plot_hist, df=d_cor_df) #create histograms for the d-value correlations
  d_cor_grid <- plot_grid(plotlist=d_hists) #put in a grid
  #do the same for the affect predictor correlations
  pred_hists <- lapply(pred_cor_names, plot_hist, df=pred_cor_df)
  pred_cor_grid <- plot_grid(plotlist=pred_hists)
  
  #get plots of learning curves at the same level
  learning_curve_grids <- list()
  #for each pair...
  for(p in 1:n_pb){
    #declare a df that will contain the info for learning curves for this pair
    plot_df_names <- c("dQ_mean","dQ_sd","dA_mean","dA_sd","dQ_mean_z","dQ_sd_z","dA_mean_z","dA_sd_z")
    plot_df <- as.data.frame(matrix(ncol=length(plot_df_names),nrow=n_tp))
    colnames(plot_df) <- plot_df_names
    
    #create an array that will initially be filled out with the learning curve info - the info will later be
    #transferred to the above df. Create 4 tables with two columns and one row per trial. The four tables are
    #split along dimension of correction and d value.
    init_plot_df <- array(dim=c(n_tp,2,2,2)) 
    #loop through each of the four tables for each pair...
    for(correct in 1:2){
      for(d in 1:2){
        temp_mat <- samp_ds[,,correct,d,p] #assign the table to a temporary matrix
        init_plot_df[,1,d,correct] <- rowMeans(temp_mat) #assign the means to the first column in the relevant table
        init_plot_df[,2,d,correct] <- rowSds(temp_mat) #assign the SDs to the same table
      }
    }
    #collapse the four tables in to the plot_df
    plot_df[1:2] <- init_plot_df[,,1,1]
    plot_df[3:4] <- init_plot_df[,,2,1]
    plot_df[5:6] <- init_plot_df[,,1,2]
    plot_df[7:8] <- init_plot_df[,,2,2]
    
    lcg <- lc_grid(plot_df) #pass the df to a function that creates a grid of the desired learning curves
    lcg_list <- c(lcg_list,list(lcg)) #add this grid to the return list
    if(return_plot_dfs){
      plot_df_list <- c(plot_df_list,list(plot_df)) #add the plot_df too, if requested
    }
  }
  
  #Add all 
  return_list <- c(return_list,d_cor_grid=list(d_cor_grid))
  return_list <- c(return_list,pred_cor_grid=list(pred_cor_grid))
  return_list <- c(return_list,lcg_list=list(lcg_list))
  inputs <- as.list(match.call())[-1] #get function inputs
  if(return_plot_dfs){
    return_list <- c(return_list,plot_df_list=list(plot_df_list))
  }
  return_list <- c(return_list,func_inputs=list(inputs))
  
  return(return_list)
}









