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
  cat("\nApprove:\n")
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

# Build population-level simulation anchors from a fitted Study 3 choice model summary.
# Returns:
#   anchor: data frame with one row per hierarchical parameter and columns for mu/sigma medians
#   resid_sigma: posterior median residual SD for simulated valence ratings
#   focal_scales: absolute fitted magnitudes used to set the true-mean simulation range for B_Q/B_A/B_R
s3_build_choice_recovery_anchor <- function(sum_df) {
  if (!exists("get_s1_median", mode = "function")) {
    stop("A summary-median helper is required (source s22_utilities.R).")
  }
  median_from_summary <- get_s1_median # Set alias

  # Hierarchical parameters with subject-level draws: param ~ N(mu, sigma).
  hier_params <- c(
    "rew_fr_sens", "aff_fr_sens", "resid_fr_sens",
    "rew_nf_sens", "aff_nf_sens", "resid_nf_sens",
    "ls_bias", "csens",
    "dcy_fr", "dcy_nf", "lrn_fr", "lrn_nf", "lrn_c",
    "B_0", "B_rew_fr", "B_rew_nf", "B_bv_fr",
    "B_q_fr", "B_q_nf", "B_pwqd_fr", "B_pwqd_nf", "B_auto"
  )

  # Learning/decay rates are constrained in Stan and represented on the logit scale.
  bounded_params <- c("dcy_fr", "dcy_nf", "lrn_fr", "lrn_nf", "lrn_c")

  anchor <- data.frame(
    param = hier_params,
    mu = NA_real_,
    sigma = NA_real_,
    bounded = hier_params %in% bounded_params,
    stringsAsFactors = FALSE
  )

  # Read posterior medians for each hierarchical mu/sigma pair.
  for (i in seq_len(nrow(anchor))) {
    p <- anchor$param[i]
    anchor$mu[i] <- median_from_summary(sum_df, paste0(p, "_mu"))
    anchor$sigma[i] <- median_from_summary(sum_df, paste0(p, "_sigma"))
  }

  # Residual SD is a single population parameter.
  resid_sigma <- median_from_summary(sum_df, "resid_sigma")

  # These fitted magnitudes define the simulation range for the three focal population means.
  focal_scales <- c(
    rew_fr_sens = abs(median_from_summary(sum_df, "rew_fr_sens_mu")),
    aff_fr_sens = abs(median_from_summary(sum_df, "aff_fr_sens_mu")),
    resid_fr_sens = abs(median_from_summary(sum_df, "resid_fr_sens_mu"))
  )

  list(anchor = anchor, resid_sigma = resid_sigma, focal_scales = focal_scales)
}

# Draw one set of true focal population means for one simulated dataset.
# lower_mult and upper_mult set the uniform sampling range around each fitted magnitude.
s3_draw_true_fr_means <- function(focal_scales, lower_mult = -1.5, upper_mult = 1.5) {
  c(
    rew_fr_sens = runif(1, min = lower_mult * focal_scales["rew_fr_sens"], max = upper_mult * focal_scales["rew_fr_sens"]),
    aff_fr_sens = runif(1, min = lower_mult * focal_scales["aff_fr_sens"], max = upper_mult * focal_scales["aff_fr_sens"]),
    resid_fr_sens = runif(1, min = lower_mult * focal_scales["resid_fr_sens"], max = upper_mult * focal_scales["resid_fr_sens"])
  )
}

# Draw subject-level parameters for one simulated dataset.
# Focal FR means are replaced with simulation-specific truth; all other means remain anchored.
s3_draw_subject_params <- function(anchor_list, sub_ids, true_fr_means) {
  anchor <- anchor_list$anchor
  n_s <- length(sub_ids)

  subj_params <- data.frame(sub_index = sub_ids, stringsAsFactors = FALSE)

  # Draw each subject-level parameter from its anchored population distribution.
  for (i in seq_len(nrow(anchor))) {
    p <- anchor$param[i]

    mu_val <- anchor$mu[i]
    # Replace only the focal FR means; all other means remain anchored to fitted medians.
    if (p %in% names(true_fr_means)) {
      mu_val <- true_fr_means[p]
    }

    draws <- rnorm(n_s, mean = mu_val, sd = anchor$sigma[i])
    if (isTRUE(anchor$bounded[i])) {
      draws <- plogis(draws)
    }
    subj_params[[p]] <- draws
  }

  subj_params
}

# Simulate one subject's full trial sequence from the fitted choice/learning equations.
# This function generates:
#   1) trial-wise choices from a softmax policy with side-bias and autocorrelation
#   2) valence ratings from the fitted valence equations
#   3) RL state trajectories used as predictors in recovery regressions
s3_simulate_subject_choice_model <- function(sub_design, par_row, resid_sigma, n_f) {
  sub_design <- sub_design[order(sub_design$overall_trial_nl), , drop = FALSE]

  # Allocate trial-level outputs and predictors used in recovery regressions.
  sub_design$choice_sim <- NA_integer_
  sub_design$chose_a_sim <- NA_integer_
  sub_design$chosen_frac_sim <- NA_integer_
  sub_design$unchosen_frac_sim <- NA_integer_
  sub_design$p_choose_a_sim <- NA_real_
  sub_design$p_chosen_sim <- NA_real_
  sub_design$valrat_sim <- NA_real_
  sub_design$q_fr_diff <- NA_real_
  sub_design$a_fr_diff <- NA_real_
  sub_design$r_fr_diff <- NA_real_
  sub_design$q_nf_diff <- NA_real_
  sub_design$a_nf_diff <- NA_real_
  sub_design$r_nf_diff <- NA_real_
  sub_design$c_diff <- NA_real_

  # Initialize latent state vectors.
  Q_fr <- rep(0, n_f)
  Q_nf <- rep(0, n_f)
  A_fr <- rep(0, n_f)
  A_nf <- rep(0, n_f)
  R_fr <- rep(0, n_f)
  R_nf <- rep(0, n_f)
  C <- rep(0, n_f)

  prev_sim <- 0

  # Step through trials in order, updating states after each simulated choice/outcome.
  for (t in seq_len(nrow(sub_design))) {
    cue_a <- sub_design$fA_ix[t]
    cue_b <- sub_design$fB_ix[t]

    # Store predictor differences before choice on this trial (the regressors used in recovery).
    sub_design$q_fr_diff[t] <- Q_fr[cue_a] - Q_fr[cue_b]
    sub_design$a_fr_diff[t] <- A_fr[cue_a] - A_fr[cue_b]
    sub_design$r_fr_diff[t] <- R_fr[cue_a] - R_fr[cue_b]
    sub_design$q_nf_diff[t] <- Q_nf[cue_a] - Q_nf[cue_b]
    sub_design$a_nf_diff[t] <- A_nf[cue_a] - A_nf[cue_b]
    sub_design$r_nf_diff[t] <- R_nf[cue_a] - R_nf[cue_b]
    sub_design$c_diff[t] <- C[cue_a] - C[cue_b]

    # Simulate choice under a two-option softmax with side bias on option A.
    value_a <- par_row$rew_fr_sens * Q_fr[cue_a] + par_row$rew_nf_sens * Q_nf[cue_a] +
      par_row$aff_fr_sens * A_fr[cue_a] + par_row$aff_nf_sens * A_nf[cue_a] +
      par_row$resid_fr_sens * R_fr[cue_a] + par_row$resid_nf_sens * R_nf[cue_a] +
      par_row$csens * C[cue_a] + par_row$ls_bias

    value_b <- par_row$rew_fr_sens * Q_fr[cue_b] + par_row$rew_nf_sens * Q_nf[cue_b] +
      par_row$aff_fr_sens * A_fr[cue_b] + par_row$aff_nf_sens * A_nf[cue_b] +
      par_row$resid_fr_sens * R_fr[cue_b] + par_row$resid_nf_sens * R_nf[cue_b] +
      par_row$csens * C[cue_b]

    p_choose_a <- 1 / (1 + exp(value_b - value_a))
    choice_sim <- ifelse(runif(1) < p_choose_a, 1L, 2L)
    chosen_frac <- ifelse(choice_sim == 1L, cue_a, cue_b)
    unchosen_frac <- ifelse(choice_sim == 1L, cue_b, cue_a)
    p_chosen <- ifelse(choice_sim == 1L, p_choose_a, 1 - p_choose_a)

    sub_design$choice_sim[t] <- choice_sim
    sub_design$chose_a_sim[t] <- ifelse(choice_sim == 1L, 1L, 0L)
    sub_design$chosen_frac_sim[t] <- chosen_frac
    sub_design$unchosen_frac_sim[t] <- unchosen_frac
    sub_design$p_choose_a_sim[t] <- p_choose_a
    sub_design$p_chosen_sim[t] <- p_chosen

    # Compute model-predicted valence for this trial before adding residual noise.
    if (sub_design$show_fres[t] == 1) {
      curr_pred <- par_row$B_0 + par_row$B_rew_fr * sub_design$out[t] +
        par_row$B_bv_fr * sub_design$box_val[t] +
        par_row$B_q_fr * Q_fr[chosen_frac] +
        par_row$B_pwqd_fr * (Q_fr[chosen_frac] * p_chosen + Q_fr[unchosen_frac] * (1 - p_chosen))
    } else {
      # The fitted choice model uses Q_fr terms in both FR and NF valence equations.
      curr_pred <- par_row$B_0 + par_row$B_rew_nf * sub_design$out[t] +
        par_row$B_q_nf * Q_fr[chosen_frac] +
        par_row$B_pwqd_nf * (Q_fr[chosen_frac] * p_chosen + Q_fr[unchosen_frac] * (1 - p_chosen))
    }

    nuis <- par_row$B_auto * prev_sim

    # Generate one valence rating on every trial.
    y <- curr_pred + nuis + rnorm(1, mean = 0, sd = resid_sigma)
    resid <- y - (curr_pred + nuis)
    sub_design$valrat_sim[t] <- y
    prev_sim <- y

    # Apply decay to all cue associations before outcome-based updates.
    Q_fr <- (1 - par_row$dcy_fr) * Q_fr
    Q_nf <- (1 - par_row$dcy_nf) * Q_nf
    A_fr <- (1 - par_row$dcy_fr) * A_fr
    A_nf <- (1 - par_row$dcy_nf) * A_nf
    R_fr <- (1 - par_row$dcy_fr) * R_fr
    R_nf <- (1 - par_row$dcy_nf) * R_nf

    # Update FR vs NF streams based on whether the fractal result was shown on this trial.
    if (sub_design$show_fres[t] == 1) {
      Q_fr[chosen_frac] <- Q_fr[chosen_frac] + par_row$lrn_fr * (sub_design$out[t] - Q_fr[chosen_frac])
      A_fr[chosen_frac] <- A_fr[chosen_frac] + par_row$lrn_fr * (curr_pred - A_fr[chosen_frac])
      R_fr[chosen_frac] <- R_fr[chosen_frac] + par_row$lrn_fr * (resid - R_fr[chosen_frac])
    } else {
      Q_nf[chosen_frac] <- Q_nf[chosen_frac] + par_row$lrn_nf * (sub_design$out[t] - Q_nf[chosen_frac])
      A_nf[chosen_frac] <- A_nf[chosen_frac] + par_row$lrn_nf * (curr_pred - A_nf[chosen_frac])
      R_nf[chosen_frac] <- R_nf[chosen_frac] + par_row$lrn_nf * (resid - R_nf[chosen_frac])
    }

    # Update trial-to-trial choice autocorrelation values for the displayed cues.
    choice_a <- ifelse(choice_sim == 1L, 1, 0)
    choice_b <- ifelse(choice_sim == 2L, 1, 0)
    C[cue_a] <- C[cue_a] + par_row$lrn_c * (choice_a - C[cue_a])
    C[cue_b] <- C[cue_b] + par_row$lrn_c * (choice_b - C[cue_b])
  }

  sub_design
}

# Simulate one complete Study 3 dataset on the fixed task scaffold.
# Runs subject-level simulation on each participant's trial structure, then stacks results.
s3_simulate_one_choice_dataset <- function(trial_template, subj_params, resid_sigma) {
  sub_ids <- sort(unique(trial_template$sub_index))
  n_f <- max(c(trial_template$fA_ix, trial_template$fB_ix), na.rm = TRUE)

  sim_list <- vector("list", length(sub_ids))
  for (i in seq_along(sub_ids)) {
    s <- sub_ids[i]
    sub_design <- trial_template[trial_template$sub_index == s, , drop = FALSE]
    par_row <- subj_params[subj_params$sub_index == s, , drop = FALSE]
    sim_list[[i]] <- s3_simulate_subject_choice_model(
      sub_design = sub_design,
      par_row = par_row,
      resid_sigma = resid_sigma,
      n_f = n_f
    )
  }

  dplyr::bind_rows(sim_list)
}

# Recover B_Q, B_A, and B_R from one simulated dataset using FE logistic regression.
# Failure handling marks simulations with model errors/non-convergence/separation warnings.
s3_fit_choice_recovery_fe <- function(sim_df) {
  warning_messages <- character(0)

  # Fit with all nuisance predictors plus subject fixed effects.
  fit <- withCallingHandlers(
    tryCatch(
      glm(
        chose_a_sim ~ q_fr_diff + a_fr_diff + r_fr_diff +
          q_nf_diff + a_nf_diff + r_nf_diff + c_diff + factor(sub_index),
        data = sim_df,
        family = "binomial"
      ),
      error = function(e) e
    ),
    warning = function(w) {
      warning_messages <<- c(warning_messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )

  if (inherits(fit, "error")) {
    return(list(success = FALSE, reason = "glm_error", warnings = warning_messages, coefs = NA))
  }

  warn_string <- paste(warning_messages, collapse = " | ")
  has_sep_warning <- grepl("fitted probabilities numerically 0 or 1 occurred", warn_string)
  if (!isTRUE(fit$converged) || has_sep_warning) {
    reason <- ifelse(has_sep_warning, "separation_warning", "glm_not_converged")
    return(list(success = FALSE, reason = reason, warnings = warning_messages, coefs = NA))
  }

  fit_coefs <- stats::coef(fit)
  needed <- c("q_fr_diff", "a_fr_diff", "r_fr_diff")
  if (!all(needed %in% names(fit_coefs))) {
    return(list(success = FALSE, reason = "missing_coef", warnings = warning_messages, coefs = NA))
  }

  if (!all(is.finite(fit_coefs[needed]))) {
    return(list(success = FALSE, reason = "nonfinite_coef", warnings = warning_messages, coefs = NA))
  }

  list(
    success = TRUE,
    reason = "ok",
    warnings = warning_messages,
    coefs = c(
      B_Q = unname(fit_coefs["q_fr_diff"]),
      B_A = unname(fit_coefs["a_fr_diff"]),
      B_R = unname(fit_coefs["r_fr_diff"])
    )
  )
}

# Run Study 3 parameter recovery across many simulated datasets.
# Workflow per simulation:
#   1) draw true focal means
#   2) draw subject-level parameters
#   3) simulate a full dataset
#   4) fit FE recovery model and save recovered coefficients
run_s3_choice_param_recovery <- function(
  trials,
  sum_df,
  n_sims = 250,
  seed = 19,
  lower_mult = -1.5,
  upper_mult = 1.5
) {
  set.seed(seed)

  # Keep only fields needed to reconstruct the trial sequence and outcomes.
  trial_template <- trials |>
    dplyr::select(sub_index, block, overall_trial_nl, fA_ix, fB_ix, show_fres, out, box_val, rat_number) |>
    dplyr::arrange(sub_index, overall_trial_nl)

  anchor_list <- s3_build_choice_recovery_anchor(sum_df = sum_df)
  sub_ids <- sort(unique(trial_template$sub_index))

  param_rows <- c("B_Q", "B_A", "B_R")
  out_rows <- vector("list", length = n_sims * length(param_rows))
  row_idx <- 1

  # Repeat the full simulate-and-recover pipeline n_sims times.
  for (sim in seq_len(n_sims)) {
    true_means <- s3_draw_true_fr_means(
      focal_scales = anchor_list$focal_scales,
      lower_mult = lower_mult,
      upper_mult = upper_mult
    )
    subj_params <- s3_draw_subject_params(
      anchor_list = anchor_list,
      sub_ids = sub_ids,
      true_fr_means = true_means
    )

    sim_df <- s3_simulate_one_choice_dataset(
      trial_template = trial_template,
      subj_params = subj_params,
      resid_sigma = anchor_list$resid_sigma
    )

    fit_out <- s3_fit_choice_recovery_fe(sim_df)

    # Save one row per focal parameter so downstream plotting can facet by parameter.
    for (p in param_rows) {
      if (p == "B_Q") {
        true_val <- true_means["rew_fr_sens"]
        rec_val <- if (isTRUE(fit_out$success)) fit_out$coefs["B_Q"] else NA_real_
      } else if (p == "B_A") {
        true_val <- true_means["aff_fr_sens"]
        rec_val <- if (isTRUE(fit_out$success)) fit_out$coefs["B_A"] else NA_real_
      } else {
        true_val <- true_means["resid_fr_sens"]
        rec_val <- if (isTRUE(fit_out$success)) fit_out$coefs["B_R"] else NA_real_
      }

      out_rows[[row_idx]] <- data.frame(
        sim = sim,
        parameter = p,
        true_value = as.numeric(true_val),
        recovered_value = as.numeric(rec_val),
        success = isTRUE(fit_out$success),
        fail_reason = fit_out$reason,
        stringsAsFactors = FALSE
      )
      row_idx <- row_idx + 1
    }
  }

  estimates <- dplyr::bind_rows(out_rows)

  fail_summary <- estimates |>
    dplyr::group_by(parameter, fail_reason) |>
    dplyr::summarise(n = dplyr::n_distinct(sim), .groups = "drop")

  list(
    estimates = estimates,
    fail_summary = fail_summary
  )
}

# Plot true-vs-recovered Study 3 choice parameters with a 45-degree reference line.
# Each panel corresponds to one focal parameter (B_Q, B_A, or B_R).
plot_s3_choice_param_recovery <- function(
  estimates_df,
  point_alpha = 0.35,
  point_size = 1.35
) {
  plot_df <- estimates_df[is.finite(estimates_df$recovered_value), , drop = FALSE]

  ggplot2::ggplot(plot_df, ggplot2::aes(x = true_value, y = recovered_value)) +
    ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dotted", linewidth = 0.6, color = "#4D4D4D") +
    ggplot2::geom_point(alpha = point_alpha, size = point_size, color = "#3E7CB1") +
    ggplot2::facet_wrap(~parameter, scales = "free") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      x = "True Population Mean",
      y = "Recovered Estimate"
    ) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "none"
    )
}

# Rebuild Study 3 recovery plots directly from a saved estimates CSV.
build_s3_choice_param_recovery_plot_from_file <- function(
  estimates_file,
  point_alpha = 0.35,
  point_size = 1.35
) {
  estimates_df <- read.csv(estimates_file, stringsAsFactors = FALSE)
  panel_plot <- plot_s3_choice_param_recovery(
    estimates_df = estimates_df,
    point_alpha = point_alpha,
    point_size = point_size
  )

  list(estimates = estimates_df, panel_plot = panel_plot)
}
