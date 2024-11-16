# Estimates effects of Q and A on choice using an approximate data analysis method
est_bq_ba_s3 <- function(data){
  # Get Q_fr values
  data_list <- by(data,data$sub_index, add_assoc, cue_cols = "chosen_frac", out_cols = "out", update_col = "show_fres", 
                  num_assoc = 6, assoc_name = "q_fr", lrn_rate = "lrn_fr", for_rate = "dcy_fr")
  data <- do.call(rbind,data_list)
  for(r in 1:nrow(data)){
    data$q_fr_ch[r] <- data[r,paste0("q_fr_",data$chosen_frac[r])]
    data$q_fr_unch[r] <- data[r,paste0("q_fr_",data$unchosen_frac[r])]
  }

  data_list <- list()
  for(s in 1:max(data$sub_index)){
    sub_data <- filter(data,sub_index == s)
    sub_data_fres <- filter(sub_data,show_fres==1)
    sub_data_nofres <- filter(sub_data,show_fres==0)
    fres_fit <- lm(valrat_z ~ out + box_val + q_fr_ch + q_fr_unch + prat, sub_data_fres)
    nofres_fit <- lm(valrat_z ~ out + q_fr_ch + q_fr_unch + prat, sub_data_nofres)
    sub_data <- sub_data %>% 
      mutate(mod_val = ifelse(
        show_fres == 1,
        fres_fit$coefficients[1] + fres_fit$coefficients[2]*out +
          fres_fit$coefficients[3]*box_val + fres_fit$coefficients[4]*q_fr_ch +
          fres_fit$coefficients[5]*q_fr_unch,
        nofres_fit$coefficients[1] + nofres_fit$coefficients[2]*out +
          nofres_fit$coefficients[3]*q_fr_ch + nofres_fit$coefficients[4]*q_fr_unch
      )
      ) %>%
      mutate(resid = ifelse(is.na(valrat_z),0,valrat_z - mod_val))
    data_list[[s]] <- sub_data
  }
  data <- do.call(rbind,data_list)

          
  # Get remaining associations used to predict choice
  data <- data %>%
            mutate(no_fres=ifelse(show_fres==1,0,1)) %>%
            mutate(fA_chosen=ifelse(choice==1,1,0)) %>%
            mutate(fB_chosen=ifelse(choice==2,1,0))
  
  data <- data %>%
            split(data$sub_index) %>%
              lapply(add_assoc, cue_cols = "chosen_frac", out_cols = "out", update_col = "no_fres", 
                     num_assoc = 6, assoc_name = "q_nf", lrn_rate = "lrn_nf", for_rate = "dcy_nf") %>%
              lapply(add_assoc, cue_cols = "chosen_frac", out_cols = "mod_val", update_col = "show_fres", 
                     num_assoc = 6, assoc_name = "a_fr", lrn_rate = "lrn_fr", for_rate = "dcy_fr") %>%
              lapply(add_assoc, cue_cols = "chosen_frac", out_cols = "mod_val", update_col = "no_fres", 
                     num_assoc = 6, assoc_name = "a_nf", lrn_rate = "lrn_nf", for_rate = "dcy_nf") %>%
              lapply(add_assoc, cue_cols = "chosen_frac", out_cols = "resid", update_col = "show_fres", 
                     num_assoc = 6, assoc_name = "r_fr", lrn_rate = "lrn_fr", for_rate = "dcy_fr") %>%
              lapply(add_assoc, cue_cols = "chosen_frac", out_cols = "resid", update_col = "no_fres", 
                     num_assoc = 6, assoc_name = "r_nf", lrn_rate = "lrn_nf", for_rate = "dcy_nf") %>%
              lapply(add_assoc, cue_cols = c("fA_ix","fB_ix"), out_cols = c("fA_chosen","fB_chosen"),
                     num_assoc = 6, assoc_name = "c", lrn_rate = "lrn_c") %>%
            bind_rows()
  
  # Get choice predictor differences
  for(r in 1:nrow(data)){
    data[r,"c_A"] <- data[r,paste0("c_",data[r,"fA_ix"])]
    data[r,"c_B"] <- data[r,paste0("c_",data[r,"fB_ix"])]
    data[r,"q_fr_A"] <- data[r,paste0("q_fr_",data[r,"fA_ix"])]
    data[r,"q_fr_B"] <- data[r,paste0("q_fr_",data[r,"fB_ix"])]
    data[r,"a_fr_A"] <- data[r,paste0("a_fr_",data[r,"fA_ix"])]
    data[r,"a_fr_B"] <- data[r,paste0("a_fr_",data[r,"fB_ix"])]
    data[r,"r_fr_A"] <- data[r,paste0("r_fr_",data[r,"fA_ix"])]
    data[r,"r_fr_B"] <- data[r,paste0("r_fr_",data[r,"fB_ix"])]
    data[r,"q_nf_A"] <- data[r,paste0("q_nf_",data[r,"fA_ix"])]
    data[r,"q_nf_B"] <- data[r,paste0("q_nf_",data[r,"fB_ix"])]
    data[r,"a_nf_A"] <- data[r,paste0("a_nf_",data[r,"fA_ix"])]
    data[r,"a_nf_B"] <- data[r,paste0("a_nf_",data[r,"fB_ix"])]
    data[r,"r_nf_A"] <- data[r,paste0("r_nf_",data[r,"fA_ix"])]
    data[r,"r_nf_B"] <- data[r,paste0("r_nf_",data[r,"fB_ix"])]
  }
  data <- data %>%
            mutate(c_diff = c_A - c_B) %>%
            mutate(q_fr_diff = q_fr_A - q_fr_B) %>%
            mutate(a_fr_diff = a_fr_A - a_fr_B) %>%
            mutate(r_fr_diff = r_fr_A - r_fr_B) %>%
            mutate(q_nf_diff = q_nf_A - q_nf_B) %>%
            mutate(a_nf_diff = a_nf_A - a_nf_B) %>%
            mutate(r_nf_diff = r_nf_A - r_nf_B)
            
  
  # Estimate standardized effects on choice
  choice_fit <- glm(fA_chosen ~ scale(q_fr_diff) + scale(a_fr_diff) + scale(r_fr_diff) + 
                                scale(q_nf_diff) + scale(a_nf_diff) + scale(r_nf_diff) +  scale(c_diff),
                                data, family = "binomial")
  
  list("bQ"=choice_fit$coefficients[2],"bA"=choice_fit$coefficients[3])
}

# Adds a column to a trial-level dataset with the estimated affective impact of the cue outcome on that trial
add_aff_imp <- function(){

}

# Adds columns to a data frame with values of cue associations (e.g., reward associations, affect associations)
# on each trial
# cue_cols: columns with cues for which an outcome was observed
# out_cols: columns with the outcomes of the cues (in the same order)
# update_col: columns indicating whehter to do any updating on the current trial. Defaults to assuming that 
#   updating should occur on every trial.
# num_assoc: number of associations to estimate
# assoc_name: how to name the associations
# lrn_rate: learning rate; a string for a column of the df containing the learning rate, a # otherwise
# for_rate: forgetting rate; a string for a column of the df containing the learning rate, a # otherwise
# init_val: initial values for the associations
add_assoc <- function(data, cue_cols, out_cols, update_col = NULL, num_assoc, assoc_name, lrn_rate, for_rate=0, 
                      init_val=0){
  # Initialize df with association values
  cues <- data[cue_cols]
  outs <- data[out_cols]
  if(is.null(update_col)){
    update <- rep(1,nrow(data))
  } else{
    update <- data[[update_col]]
  }

  assoc <- data.frame(matrix(init_val, ncol = num_assoc, nrow = nrow(data)))
  colnames(assoc) <- paste0(assoc_name, "_", 1:num_assoc)
  
  if(is.character(lrn_rate)){
    lrn_rate <- data[1,lrn_rate]
  } 
  if(is.character(for_rate)){
    for_rate <- data[1,for_rate]
  } 

  for(r in 1:(nrow(data)-1)){
    assoc[r+1,] <- assoc[r,]*(1-for_rate) # Default to same value as above
    if(update[r] == 1){
      # For cues with outcomes, update the association
      cues_r <- as.vector(unlist(cues[r,]))
      assoc[r+1,cues_r] <- assoc[r,cues_r] + 
        lrn_rate*(outs[r,] - assoc[r,cues_r])
    }
  }
  cbind(data,assoc)
}

# Adds a column to a trial-level dataset with the mean estimates of subject-level parameters
# sum: model summary object
# param: name of parameter
# trials: trial-level dataset
# sig_trans: if true, will run the parameter mean through an inverse logistic function
add_param_means <- function(sum,param,trials,sig_trans=T){
  means_list <- ncp_mean_hist(sum,param,print_hist=F)
  if(sig_trans){
    means_vec <- means_list$means %>% sigmoid()
  } else{
    means_vec <- means_list$means
  }
  means_df <- data.frame(means_vec,c(1:length(means_vec)))
  colnames(means_df) <- c(param,"sub_index")
  left_join(trials,means_df,by="sub_index")
}


# This function works like vec_optim, but it is tailored to approximate the effects of rewards and box values 
# on choices using vectors that represent reward associations and affect associations.
rew_bv_vec_optim <- function(eff_vec,init_pars = c(0,0)){
  target <- eff_vec[1:2] # Effects of reward and bv on choice
  rew_dvec <- c(1,0)
  aff_dvec <- eff_vec[3:4]/(abs(eff_vec[3]) + abs(eff_vec[4])) # Effects of reward and bv on valence, normalized
  dvec_mat <- cbind(rew_dvec,aff_dvec)
  
  optim_fit <- optim(par = init_pars, fn = vec_ss, dvec=dvec_mat, target=target, method = "L-BFGS-B", lower=0)
  return(optim_fit$par)
}

#accepts a subject-level trials df, and returns the same df with columns for the differences between the past 
#valence/reward/etc. values of the chosen and unchosen fractal
sub_past_diffs <- function(s_trials,columns=c("valrat_z","out"),past_trials=3){
  #initialize the columns being added
  new_cols <- c()
  for (col in columns) {
    for (trial in 1:past_trials) {
      new_cols <- c(new_cols, paste0(col,"_diff_",trial))
    }
  }
  s_trials[new_cols] <- NA
  
  for(t in 1:nrow(s_trials)){
    #get values for chosen and unchosen fractals
    past_chosen <- past_diffs(df=s_trials,frac_ix = s_trials[t,"chosen_frac"],columns,t,past_trials)
    past_unchosen <- past_diffs(df=s_trials,frac_ix = s_trials[t,"unchosen_frac"],columns,t,past_trials)
      
    s_trials[t,new_cols] <- past_chosen-past_unchosen #assign to the current trial
  }
  
  return(s_trials)
}

#Gets a vector of past values for given columns
past_diffs <- function(df,frac_ix, columns, t, past_trials){
  fres_df <- df %>% filter(show_fres == 1 & chosen_frac == frac_ix & trial <= t) 
  past_vals <- fres_df[columns]
  
  #if fewer than past_trials columns in this df, add rows with 0 in them
  if(nrow(past_vals) < past_trials){
    num_zero_rows <- past_trials - nrow(past_vals)
    zero_df <- data.frame(matrix(0,nrow=num_zero_rows,length(columns)))
    names(zero_df) <- columns
    past_vals <- rbind(zero_df,past_vals)
  }
  tail_past_vals <- tail(past_vals,past_trials)
  tpv_ordered <- tail_past_vals[past_trials:1, ]
  return(unlist(tpv_ordered)) #just get the values from the last three trials
}

#Creates a list of Prolific submissions to be approved/rejected based on attention checks
#sub: path to a subject-level CSV with columns for id and att_checks_passed
#date_min/max: the earliest/latest date from which to review participants. By default, every participant in the sub CSV is included.
#min_passed: the minimum number of attention checks that must be passed to be approved
#ids: Optional: a list of prolific IDs to review. If a value is passed in for this, date min/max is ignored
atcheck_review <- function(sub_path,date_min="1000-05-13_19h21.26.622",date_max ="4000-05-13_19h21.26.622",min_passed = NULL, ids=NULL){
  sub <- read.csv(sub_path)
  if(is.null(ids)){
    sub <- sub %>% filter(date >= date_min & date <= date_max) #grab subjects within the date range
  } else{
    sub <- sub %>% filter(id %in% ids) #if IDs were passed in, instead grab subjects with an acceptable ID
  }
  
  #get subjects who should be approved/rejected
  approve_df <- sub %>% filter(att_checks_passed >= min_passed)
  reject_df <- sub %>% filter(att_checks_passed < min_passed)
  # Print the result
  cat("Reject:\n")
  cat(paste(reject_df$id, collapse = "\n"))
  cat("Approve:\n")
  cat(paste(approve_df$id, collapse = "\n"))
  
  return(list(approve_ids = approve_df$id,reject_ids = reject_df$id))
}


#Creates a list with Prolific IDs and bonus payments to be copied and pasted onto Prolific
#sub: path to a subject-level df with columns for earnings,date, andid
#date_min/max: the earliest/latest date from which to pay participants. By default, there is no minimum or maximum limit.
#ids: Optional: a list of prolific IDs for which to get payments. If a value is passed in for this, date min/max is ignored
prolific_bp <- function(sub_path,date_min="1000-05-13_19h21.26.622",date_max ="4000-05-13_19h21.26.622",ids=NULL){
  sub <- read.csv(sub_path)
  if(is.null(ids)){
    sub <- sub %>% filter(date >= date_min & date <= date_max) #grab subjects within the date range
  } else{
    sub <- sub %>% filter(id %in% ids) #if IDs were passed in, instead grab subjects with an acceptable ID
  }

  # Convert data frame to list
  bp_list <- apply(sub, 1, function(row) paste(row["id"], row["earnings"], sep = ","))
  
  # Print the result
  cat(paste(bp_list, collapse = "\n"))
}

#Returns a list of data for input to stan, given trial-level data
#trials: trial-level data
#n_t: number of trials; must be set manually
#n_f: number of fractls. If NULL, it will be set to however many distinct fractal numbers are found in the data.
stan_data_arlbv <- function(trials,n_t){
  library(gdata)
  
  n_s <- length(unique(trials$id)) #get number of subjects
  
  frac_vec <- as.vector(as.matrix(trials[c("fA_img","fB_img")])) #get a big vector of fractals...
  n_f <- length(unique(frac_vec)) #and then get the number of fractals
  
  #Get matrices with subjects on the rows and trials in the columns, containing indices for fractal A/B
  fA <- sub_by_trial_matrix(trials,"fA_ix")   
  fB <- sub_by_trial_matrix(trials,"fB_ix")
  
  choice <- sub_by_trial_matrix(trials,"choice_numeric") #ditto whether fractal A (1) or B(2) was chosen  
  
  chosen_frac <- sub_by_trial_matrix(trials,"chosen_frac")
  unchosen_frac <- sub_by_trial_matrix(trials,"unchosen_frac")

  rew <- sub_by_trial_vec_list(trials,"out") #get a list of vectors, one per subject, each containing the reward for each trial
  bv <- sub_by_trial_vec_list(trials,"box_val")
  #resid <- sub_by_trial_vec_list(trials,"resids")
  
  fres <- sub_by_trial_matrix(trials,"show_fres")
  
  rat_num <- sub_by_trial_matrix(trials,"rat_number")
  n_rat <- max(trials$rat_number)
  
  #get ratings
  rat_trials <- filter(trials,rat_number != 0)
  rat <- rat_trials$valrat_z
  
  #get filled-in affect ratings for each trial
  #fi_rat <- sub_by_trial_vec_list(trials,"fi_rat") 
  
  prev_rat <- sub_by_trial_vec_list(trials,"prev_rate") #previous valence rating
  
  data <- list(
    n_t = n_t,
    n_s = n_s,
    n_f = n_f,
    fA = fA,
    fB = fB,
    rew = rew,
    choice = choice,
    chosen_frac = chosen_frac,
    unchosen_frac = unchosen_frac,
    bv = bv,
    fres = fres,
    rat = rat,
    rat_num = rat_num,
    n_rat = n_rat,
    prev_rat = prev_rat
  )
  return(data)
}