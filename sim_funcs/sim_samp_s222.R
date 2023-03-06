# Simulates a sample from teh s22-2 experiment, returning a df with the trial-level data from the simulation
# ns: number subjects
# nb: number blocks
# n_pb: number of pairs per block
# n_tb: number of trials per pair
# pair_type: fixed or mixed
# mb_size:size of mini-blocks within which to shuffle pair presentation order (for fixed pair type only)
# alpha_r_mean: mean of distribution to draw the alpha for rewards from
# other mean/sd arguments: interepreted similarly to the above
# w_r_mean: mean of distribution to draw the valence weight for rewards from
# other w_ arguments: interpreted similarly to the above
# prob_r: a matrix with 1 column for each possible outcome (-$5 to +$5) and 1 row for each fractal. Values are 
#   the probability of receiving that reward when choosing that fractal
# prob_fo: similar to the above, but for fake outcomes

sim_samp_s222 <- function(ns,nb,n_pb,n_tp,pair_type="fixed",mb_size=4,
                          alpha_r_mean,alpha_r_sd=0,beta_r_mean,beta_r_sd,
                          forget_alpha_r="match",forget_alpha_a="match",
                          forget_r_mean=NULL,forget_r_sd=NULL,alpha_a_mean,
                          alpha_a_sd=0,beta_a_mean,beta_a_sd,forget_a_mean=NULL,
                          forget_a_sd=NULL,phi_mean,phi_sd,tau_mean,tau_sd,
                          w_r_mean,w_r_sd=0,w_fo_mean,w_fo_sd=0,w_q_mean,
                          w_q_sd=0,w_opv_mean,w_opv_sd=0,v_resid_var,
                          prob_r,prob_fo,r_seq=NULL,fo_seq=NULL){
  
  library(tidyverse)
  library(sigmoid)
  library(truncnorm)
  
  #add some variables based on the inputs
  n_tb <- n_pb*n_tp #number of trials per block
  n_fb <- n_pb*2 #number of fractals per block

  #Create df with the trial level data for the full sample (samp_df)
  frac_qs <- paste0("Q",1:n_fb)
  frac_as <- paste0("A",1:n_fb)
  frac_cs <- paste0("C",1:n_fb)
  frac_fs <- paste0("F",1:n_fb)
  frac_ps <- paste0("p",1:n_fb)
  
  df_names <- c("subject","block","trial","fA","fB","QfA","QfB","dQ","Q_chosen","AfA","AfB","dA",
                "CfA","CfB","dC","FfA","FfB","dF",frac_qs,frac_as,frac_cs,frac_ps,frac_fs,"onpol_val",
                "choice","reward","fake_out","valence",
                "alpha_r","beta_r","forget_r","alpha_a","beta_a","forget_a","phi","tau",
                "w_r","w_fo","w_q","w_opv")
  samp_df <- data.frame(matrix(ncol=length(df_names),nrow=0))
  names(samp_df) <- df_names
  template_df <- samp_df #this template will be used to create the block-level dfs
  
  for(s in 1:ns){
    #get subject-level empirical parameters
    cur_alpha_r <- rtruncnorm(n=1,a=0,b=1,alpha_r_mean,alpha_r_sd)
    cur_beta_r <- rtruncnorm(n=1,a=0,mean=beta_r_mean,sd=beta_r_sd)
    if(forget_alpha_r == "match"){
      cur_forget_r <- cur_alpha_r
    } else if(forget_alpha_f == "diff"){
      cur_forget_r <- rtruncnorm(n=1,a=0,b=1,forget_r_mean,forget_r_sd)
    } else{
      stop("acceptable values for forget_alpha_r are 'match' and 'diff'")
    }
    cur_alpha_a <- rtruncnorm(n=1,a=0,b=1,alpha_a_mean,alpha_a_sd)
    cur_beta_a <- rtruncnorm(n=1,a=0,mean=beta_a_mean,sd=beta_a_sd)
    if(forget_alpha_a == "match"){
      cur_forget_a <- cur_alpha_a
    } else if(forget_alpha_a == "diff"){
      cur_forget_a <- rtruncnorm(n=1,a=0,b=1,forget_a_mean,forget_a_sd)
    } else{
      stop("acceptable values for forget_alpha_a are 'match' and 'diff'")
    }
    cur_phi <- rtruncnorm(n=1,a=0,mean=phi_mean,sd=phi_sd)
    cur_tau <- rtruncnorm(n=1,a=0,b=1,tau_mean,tau_sd)
    cur_w_r <- rnorm1(w_r_mean,w_r_sd)
    cur_w_fo <- rnorm1(w_fo_mean,w_fo_sd)
    cur_w_q <- rnorm1(w_q_mean,w_q_sd)
    cur_w_opv <- rnorm1(w_opv_mean,w_opv_sd)
    
    for(b in 1:nb){
      #Create df with the trial level data for this block - same as samp_df, but with one row per trial in the block
      block_df <- template_df
      block_df[n_tb,] <- NA
      block_df$trial <- c(1:n_tb)
      block_df$block <- b
      
      #Add subject-specific variables to the df
      block_df$subject <- s
      block_df$alpha_r <- cur_alpha_r
      block_df$beta_r <- cur_beta_r
      block_df$forget_r <- cur_forget_r
      block_df$alpha_a <- cur_alpha_a
      block_df$beta_a <- cur_beta_a
      block_df$forget_a <- cur_forget_a
      block_df$phi <- cur_phi
      block_df$tau <- cur_tau
      block_df$w_r <- cur_w_r
      block_df$w_fo <- cur_w_fo
      block_df$w_q <- cur_w_q
      block_df$w_opv <- cur_w_opv
      
      #Initialize values to 0
      block_df$QfA <- 0
      block_df$QfB <- 0
      block_df$dQ <- 0
      block_df$AfA <- 0
      block_df$AfB <- 0
      block_df$dA <- 0
      block_df$CfA <- 0
      block_df$CfB <- 0
      block_df$dC <- 0
      block_df$FfA <- 0
      block_df$FfB <- 0
      block_df$dF <- 0
      block_df[c(frac_qs,frac_as,frac_cs,frac_fs)] <- 0
      block_df[grep("p\\d",names(block_df))] <- 0 #choice probabilities
      
      #Establish the fractals presented on each trial
      #If using a fixed pair task structure...
      if(pair_type == "fixed"){
        if(n_tb %% mb_size != 0){
          stop("n_tb must be a multiple of mb_size if using the fixed pair option")
        }
        pair_reps <- mb_size/n_pb #number of times each pair must appear in a mini-block to fill it
        unshuf_mb<- rep(1:n_pb,pair_reps) #create an unshuffled mini-block (e.g., c(1,2,1,2))
        pair_seq <- c() #vector with the pairs presented on each trial
        while(length(pair_seq) < n_tb){
          shuf_mb <- sample(unshuf_mb) #shuffle the mini-blcok
          pair_seq <- c(pair_seq,shuf_mb) #add it to the larger vector of pairs
        }
        block_df$fB <- pair_seq*2 #fractal B is the pair number*2, e.g. fractal B for the second pair should be fractal 4. 
        block_df$fA <- block_df$fB - 1 #fractal A is the preceding fractal, e.g., fractal A for the second pair should be 3
      }
      
      r_inds <- c(rep(0,n_fb)) #index for the number of the reward sequence we're on (for rand_r=FALSE only)
      f_inds <- c(rep(0,n_fb)) #ditto for fake outcomes
      
      for(t in 1:n_tb){
        #grab fA and fB for this trial
        fA <- block_df$fA[t]
        fB <- block_df$fB[t]
        #fill Q, A, C, and F values column
        block_df$QfA[t] <- block_df[t,paste0("Q",fA)]
        block_df$QfB[t] <- block_df[t,paste0("Q",fB)]
        block_df$AfA[t] <- block_df[t,paste0("A",fA)]
        block_df$AfB[t] <- block_df[t,paste0("A",fB)]
        block_df$CfA[t] <- block_df[t,paste0("C",fA)]
        block_df$CfB[t] <- block_df[t,paste0("C",fB)]
        block_df$FfA[t] <- block_df[t,paste0("F",fA)]
        block_df$FfB[t] <- block_df[t,paste0("F",fB)]
        #calculate dQ, dA, dC, and dF
        block_df$dQ[t] <- block_df$QfA[t] - block_df$QfB[t]
        block_df$dA[t] <- block_df$AfA[t] - block_df$AfB[t]
        block_df$dC[t] <- block_df$CfA[t] - block_df$CfB[t]
        block_df$dF[t] <- block_df$FfA[t] - block_df$FfB[t]
        #calculate choice probabilities, filling in 0 for fractals not presented
        prob_fA <- sigmoid(cur_beta_r*block_df$dQ[t] + cur_beta_a*block_df$dA[t]
                           + cur_phi*block_df$dC[t]) #calculate the prob of choosing fA
        prob_fB <- 1 - prob_fA #probability of choosing fB is the inverse
        #record probabilities
        block_df[t,paste0("p",fA)] <- prob_fA
        block_df[t,paste0("p",fB)] <- prob_fB
        #the other fractals' choice probabilities are left at 0
        
        #calculate on-policy value
        onpol_val <- as.numeric(block_df[t,grep("p\\d",names(block_df))]) %*% as.numeric(block_df[t,grep("Q\\d",names(block_df))])
        block_df$onpol_val[t] <- onpol_val #record it
        
        chosen_frac <- sample(c(fA,fB),prob=c(prob_fA,prob_fB),size=1) #sample choice
        block_df$choice[t] <- chosen_frac #record chosen fractal
        #get unchosen fractal (for C update later on)
        if(chosen_frac == fA){
          unchosen_frac <- fB
        } else{
          unchosen_frac <- fA
        }
        
        if(any(!is.na(prob_r[chosen_frac]))){
          reward <- sample(c(-5:5),prob=prob_r[chosen_frac,],size=1) #select reward, if there is a dist. for this fractal
        } else{
          #if there's no probability dist. for this fractal, that means there's a fixed reward sequence...
          r_inds[chosen_frac] <- r_inds[chosen_frac] + 1 # add 1 to the sequence index value of the chosen fractal,
                                                         # updating how many times you've selected the fractal
          r_ind <- r_inds[chosen_frac]          
          reward <- r_seq[chosen_frac,r_ind] #get the reward value for htis trial of the sequence
        }
        #the selected fractal
        block_df$reward[t] <- reward #log reward
        
        #do the same for fake outcomes...
        if(any(!is.na(prob_fo[chosen_frac]))){
          fake_out <- sample(c(-5:5),prob=prob_fo[chosen_frac,],size=1) 
        } else{
          f_inds[chosen_frac] <- f_inds[chosen_frac] + 1
          f_ind <- f_inds[chosen_frac]
          fake_out <- fo_seq[chosen_frac,f_ind]
        }
        block_df$fake_out[t] <- fake_out
        
        chosen_q <- block_df[t, paste0("Q",chosen_frac)] #get the Q value of the chosen fractal
        block_df$Q_chosen[t] <- chosen_q
        chosen_a <- block_df[t, paste0("A",chosen_frac)] #diddo A value
        chosen_c <- block_df[t, paste0("C",chosen_frac)] #diddo C value
        chosen_f <- block_df[t, paste0("F",chosen_frac)] #diddo F value
        unchosen_c <- block_df[t, paste0("C",unchosen_frac)] #useful for the C update
        
        #calculate valence based on predictors + random error
        val_noise <- rnorm1(0,v_resid_var)
        valence <- cur_w_r*reward + cur_w_fo*fake_out + cur_w_q*chosen_q + cur_w_opv*onpol_val + val_noise
        block_df$valence[t] <- valence
        
        #if this isn't the last trial in the block, update Q, A, and C values for the next trial
        if(t < n_tb){
          #The below takes too long, so I'm not using forgetting rates for now
          # #update Q, A, and F with forgetting factor, but don't with C
          # block_df[t+1,frac_qs] <- delta_update(block_df[t,frac_qs],cur_forget_r,0)
          # block_df[t+1,frac_as] <- delta_update(block_df[t,frac_as],cur_forget_a,0)
          # block_df[t+1,frac_fs] <- delta_update(block_df[t,frac_fs],cur_forget_a,0) #use same forgetting rate as for A
          block_df[t+1,frac_qs] <- block_df[t,frac_qs]
          block_df[t+1,frac_as] <- block_df[t,frac_as]
          block_df[t+1,frac_fs] <- block_df[t,frac_fs]
          block_df[t+1,frac_cs] <- block_df[t,frac_cs]
          #update chosen fractal's Q, A, C, and F values
          block_df[t+1,paste0("Q",chosen_frac)] <- delta_update(chosen_q,cur_alpha_r,reward)
          block_df[t+1,paste0("A",chosen_frac)] <- delta_update(chosen_a,cur_alpha_a,valence)
          block_df[t+1,paste0("F",chosen_frac)] <- delta_update(chosen_f,cur_alpha_a,fake_out) #use same learning rate as for A
          block_df[t+1,paste0("C",chosen_frac)] <- delta_update(chosen_c,cur_tau,1)
          block_df[t+1,paste0("C",unchosen_frac)] <- delta_update(unchosen_c,cur_tau,0)
        }
      }
      samp_df <- rbind(samp_df,block_df) #append the block df to the sample df
    }
  }
  return(samp_df)
}





