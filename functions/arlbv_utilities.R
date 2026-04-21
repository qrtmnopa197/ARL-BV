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
  point_size = 1.35,
  equal_axes = TRUE,
  n_major_breaks = 5
) {
  plot_df <- estimates_df[is.finite(estimates_df$recovered_value), , drop = FALSE]
  
  if (isTRUE(equal_axes)) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      stop("equal_axes=TRUE requires the patchwork package.")
    }

    # Build one square-coordinate panel per parameter so each panel has matching x/y scales,
    # while still allowing different limits across parameters.
    parameter_levels <- unique(as.character(plot_df$parameter))
    panel_list <- lapply(parameter_levels, function(param_name) {
      d <- plot_df[plot_df$parameter == param_name, , drop = FALSE]

      # Set panel-specific equal-axis limits directly from the data with a small cushion.
      # Round only displayed tick labels by parameter: nearest 1 for B_Q/B_A and nearest 0.1 for B_R.
      max_abs <- max(abs(c(d$true_value, d$recovered_value)), na.rm = TRUE)
      if (!is.finite(max_abs)) {
        max_abs <- 0
      }
      target_max <- max_abs * 1.08
      round_to <- ifelse(param_name == "B_R", 0.1, 1)

      if (target_max == 0) {
        # Handle degenerate no-variation panels with a small symmetric default span.
        axis_limits <- c(-round_to, round_to)
      } else {
        axis_limits <- c(-target_max, target_max)
      }

      axis_breaks <- seq(axis_limits[1], axis_limits[2], length.out = n_major_breaks)
      label_accuracy <- round_to

      ggplot2::ggplot(d, ggplot2::aes(x = true_value, y = recovered_value)) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dotted", linewidth = 0.6, color = "#4D4D4D") +
        ggplot2::geom_point(alpha = point_alpha, size = point_size, color = "#3E7CB1") +
        ggplot2::coord_equal(xlim = axis_limits, ylim = axis_limits, expand = FALSE) +
        ggplot2::scale_x_continuous(breaks = axis_breaks, labels = scales::label_number(accuracy = label_accuracy)) +
        ggplot2::scale_y_continuous(breaks = axis_breaks, labels = scales::label_number(accuracy = label_accuracy)) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = param_name,
          x = "True Population Mean",
          y = "Recovered Estimate"
        ) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          legend.position = "none",
          plot.title = ggplot2::element_text(hjust = 0.5)
        )
    })

    return(patchwork::wrap_plots(panel_list, nrow = 1))
  }

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
  point_size = 1.35,
  equal_axes = TRUE,
  n_major_breaks = 5
) {
  estimates_df <- read.csv(estimates_file, stringsAsFactors = FALSE)
  panel_plot <- plot_s3_choice_param_recovery(
    estimates_df = estimates_df,
    point_alpha = point_alpha,
    point_size = point_size,
    equal_axes = equal_axes,
    n_major_breaks = n_major_breaks
  )

  list(estimates = estimates_df, panel_plot = panel_plot)
}

# List final_ls posterior variables needed for the Study 3 closed-loop amputation.
# These are the hierarchical choice, learning, and valence parameters used to rebuild
# trial-wise latent states under the fitted model.
s3_amputation_required_vars <- function() {
  hier_params <- c(
    "rew_fr_sens", "aff_fr_sens", "resid_fr_sens",
    "rew_nf_sens", "aff_nf_sens", "resid_nf_sens",
    "ls_bias", "csens",
    "dcy_fr", "dcy_nf", "lrn_fr", "lrn_nf", "lrn_c",
    "B_0", "B_rew_fr", "B_rew_nf", "B_bv_fr",
    "B_q_fr", "B_q_nf", "B_pwqd_fr", "B_pwqd_nf", "B_auto"
  )

  c(
    as.vector(rbind(paste0(hier_params, "_mu"), paste0(hier_params, "_sigma"), paste0(hier_params, "_z"))),
    "resid_sigma"
  )
}

# Convert selected final_ls posterior draws into subject-level parameter rows.
# The returned data frame has one row per posterior draw per participant, which is
# the format used by the closed-loop simulation helpers below.
s3_posterior_draws_to_subject_params <- function(draws_array, sub_ids, n_draws = 250, seed = 19) {
  set.seed(seed)

  # Identify and sample posterior draws after flattening iterations and chains.
  draw_dims <- dim(draws_array)
  draw_names <- dimnames(draws_array)[[3]]
  flat_draws <- matrix(
    draws_array,
    nrow = draw_dims[1] * draw_dims[2],
    ncol = draw_dims[3],
    dimnames = list(NULL, draw_names)
  )
  draw_rows <- sample(seq_len(nrow(flat_draws)), size = min(n_draws, nrow(flat_draws)), replace = FALSE)

  # Define which hierarchical parameters need logit transformation after reconstruction.
  hier_params <- c(
    "rew_fr_sens", "aff_fr_sens", "resid_fr_sens",
    "rew_nf_sens", "aff_nf_sens", "resid_nf_sens",
    "ls_bias", "csens",
    "dcy_fr", "dcy_nf", "lrn_fr", "lrn_nf", "lrn_c",
    "B_0", "B_rew_fr", "B_rew_nf", "B_bv_fr",
    "B_q_fr", "B_q_nf", "B_pwqd_fr", "B_pwqd_nf", "B_auto"
  )
  bounded_params <- c("dcy_fr", "dcy_nf", "lrn_fr", "lrn_nf", "lrn_c")

  # Reconstruct subject-level parameters for each sampled posterior draw.
  out_list <- vector("list", length(draw_rows))
  for (d in seq_along(draw_rows)) {
    draw_vec <- flat_draws[draw_rows[d], ]
    draw_df <- data.frame(
      draw_id = d,
      posterior_row = draw_rows[d],
      sub_index = sub_ids,
      resid_sigma = as.numeric(draw_vec["resid_sigma"])
    )

    # Combine population means, population SDs, and subject z-scores.
    for (p in hier_params) {
      mu <- as.numeric(draw_vec[paste0(p, "_mu")])
      sigma <- as.numeric(draw_vec[paste0(p, "_sigma")])
      z_names <- paste0(p, "_z[", seq_along(sub_ids), "]")
      z <- as.numeric(draw_vec[z_names])
      vals <- mu + sigma * z

      # Learning and decay parameters live on the logit scale in Stan.
      if (p %in% bounded_params) {
        vals <- plogis(vals)
      }
      draw_df[[p]] <- vals
    }

    out_list[[d]] <- draw_df
  }

  dplyr::bind_rows(out_list)
}

# Read final_ls posterior draws needed for the Study 3 closed-loop amputation.
s3_read_amputation_draws <- function(
  model = "final_ls",
  model_out_dir = "~/projects/ARL_bv/output/results/stan_model_fits/final_samp_fits/"
) {
  fit_csvs <- list.files(paste0(model_out_dir, model), pattern = "\\.csv$", full.names = TRUE)
  if (length(fit_csvs) > 4) {
    stop("There is more than one set of model fit CSV files in the specified folder.")
  } else if (length(fit_csvs) < 2) {
    stop("There is not a full set of model fit CSV files in the specified folder.")
  }

  cmdstanr::read_cmdstan_csv(
    files = fit_csvs,
    variables = s3_amputation_required_vars()
  )
}

# Sample subject-level parameters from already-read final_ls posterior draws.
# Passing fit_draws keeps simulation reruns fast because the CmdStan CSVs do not
# need to be read again.
s3_sample_amputation_subject_params <- function(
  fit_draws = NULL,
  model = "final_ls",
  model_out_dir = "~/projects/ARL_bv/output/results/stan_model_fits/final_samp_fits/",
  sub_ids,
  n_draws = 250,
  seed = 19
) {
  # Fall back to reading draws here for older code paths, but prefer passing fit_draws.
  if (is.null(fit_draws)) {
    fit_draws <- s3_read_amputation_draws(
      model = model,
      model_out_dir = model_out_dir
    )
  }

  # Accept either the cmdstanr read object or the raw post_warmup_draws array.
  draws_array <- if (is.list(fit_draws) && !is.null(fit_draws$post_warmup_draws)) {
    fit_draws$post_warmup_draws
  } else {
    fit_draws
  }

  s3_posterior_draws_to_subject_params(
    draws_array = draws_array,
    sub_ids = sub_ids,
    n_draws = n_draws,
    seed = seed
  )
}

# Prepare the fixed trial scaffold used by the closed-loop amputation.
# The observed schedule, outcomes, and ratings are retained as task inputs while
# simulated choices determine the model's own future cue histories.
s3_prepare_amputation_template <- function(trials) {
  trials |>
    dplyr::filter(!is.na(choice_numeric)) |>
    dplyr::select(
      sub_index, overall_trial_nl, block, trial_nl, fA_ix, fB_ix, choice_numeric,
      chosen_frac, unchosen_frac, show_fres, out, box_val,
      rat_number, valrat_z, prev_rate
    ) |>
    dplyr::arrange(sub_index, overall_trial_nl)
}

# Simulate one participant in a closed-loop Study 3 amputation.
# Unlike replay, simulated choices determine subsequent Q/A/R/C histories. The task
# schedule and outcome stream are kept fixed, so this isolates feedback from the
# model's own choices without needing to regenerate task outcomes.
s3_closed_loop_subject_amputation <- function(sub_design, par_row, n_f, aff_fr_multiplier = 1, model_label = "Full closed-loop") {
  sub_design <- sub_design[order(sub_design$overall_trial_nl), , drop = FALSE]

  # Allocate simulated outputs and keep the observed trial fields for condition summaries.
  sub_design$model <- model_label
  sub_design$draw_id <- par_row$draw_id[1]
  sub_design$p_choose_a <- NA_real_
  sub_design$choice_sim <- NA_integer_
  sub_design$chosen_frac_sim <- NA_integer_
  sub_design$unchosen_frac_sim <- NA_integer_
  sub_design$valrat_sim <- NA_real_

  # Initialize the latent state vectors at the start of the task.
  Q_fr <- rep(0, n_f)
  Q_nf <- rep(0, n_f)
  A_fr <- rep(0, n_f)
  A_nf <- rep(0, n_f)
  R_fr <- rep(0, n_f)
  R_nf <- rep(0, n_f)
  C <- rep(0, n_f)
  prev_sim <- 0

  # Generate choices and update histories from the model's own simulated path.
  for (t in seq_len(nrow(sub_design))) {
    cue_a <- sub_design$fA_ix[t]
    cue_b <- sub_design$fB_ix[t]

    # Evaluate the full or A_fr-amputated policy at the current simulated state.
    value_a <- par_row$rew_fr_sens * Q_fr[cue_a] + par_row$rew_nf_sens * Q_nf[cue_a] +
      aff_fr_multiplier * par_row$aff_fr_sens * A_fr[cue_a] + par_row$aff_nf_sens * A_nf[cue_a] +
      par_row$resid_fr_sens * R_fr[cue_a] + par_row$resid_nf_sens * R_nf[cue_a] +
      par_row$csens * C[cue_a] + par_row$ls_bias

    value_b <- par_row$rew_fr_sens * Q_fr[cue_b] + par_row$rew_nf_sens * Q_nf[cue_b] +
      aff_fr_multiplier * par_row$aff_fr_sens * A_fr[cue_b] + par_row$aff_nf_sens * A_nf[cue_b] +
      par_row$resid_fr_sens * R_fr[cue_b] + par_row$resid_nf_sens * R_nf[cue_b] +
      par_row$csens * C[cue_b]

    p_choose_a <- 1 / (1 + exp(value_b - value_a))
    choice_sim <- ifelse(runif(1) < p_choose_a, 1L, 2L)
    chosen_frac <- ifelse(choice_sim == 1L, cue_a, cue_b)
    unchosen_frac <- ifelse(choice_sim == 1L, cue_b, cue_a)
    p_chosen <- ifelse(choice_sim == 1L, p_choose_a, 1 - p_choose_a)

    sub_design$p_choose_a[t] <- p_choose_a
    sub_design$choice_sim[t] <- choice_sim
    sub_design$chosen_frac_sim[t] <- chosen_frac
    sub_design$unchosen_frac_sim[t] <- unchosen_frac

    # Generate the valence prediction from the simulated choice path.
    if (sub_design$show_fres[t] == 1) {
      curr_pred <- par_row$B_0 + par_row$B_rew_fr * sub_design$out[t] +
        par_row$B_bv_fr * sub_design$box_val[t] +
        par_row$B_q_fr * Q_fr[chosen_frac] +
        par_row$B_pwqd_fr * (Q_fr[chosen_frac] * p_chosen + Q_fr[unchosen_frac] * (1 - p_chosen))
    } else {
      curr_pred <- par_row$B_0 + par_row$B_rew_nf * sub_design$out[t] +
        par_row$B_q_nf * Q_fr[chosen_frac] +
        par_row$B_pwqd_nf * (Q_fr[chosen_frac] * p_chosen + Q_fr[unchosen_frac] * (1 - p_chosen))
    }

    # Simulate a valence rating so residual associations can evolve in closed loop.
    nuis <- par_row$B_auto * prev_sim
    valrat_sim <- curr_pred + nuis + rnorm(1, mean = 0, sd = par_row$resid_sigma)
    resid <- valrat_sim - (curr_pred + nuis)
    sub_design$valrat_sim[t] <- valrat_sim
    prev_sim <- valrat_sim

    # Decay all streams before applying the outcome update for this trial.
    Q_fr <- (1 - par_row$dcy_fr) * Q_fr
    Q_nf <- (1 - par_row$dcy_nf) * Q_nf
    A_fr <- (1 - par_row$dcy_fr) * A_fr
    A_nf <- (1 - par_row$dcy_nf) * A_nf
    R_fr <- (1 - par_row$dcy_fr) * R_fr
    R_nf <- (1 - par_row$dcy_nf) * R_nf

    # Update the simulated chosen cue in the FR or NF stream.
    if (sub_design$show_fres[t] == 1) {
      Q_fr[chosen_frac] <- Q_fr[chosen_frac] + par_row$lrn_fr * (sub_design$out[t] - Q_fr[chosen_frac])
      A_fr[chosen_frac] <- A_fr[chosen_frac] + par_row$lrn_fr * (curr_pred - A_fr[chosen_frac])
      R_fr[chosen_frac] <- R_fr[chosen_frac] + par_row$lrn_fr * (resid - R_fr[chosen_frac])
    } else {
      Q_nf[chosen_frac] <- Q_nf[chosen_frac] + par_row$lrn_nf * (sub_design$out[t] - Q_nf[chosen_frac])
      A_nf[chosen_frac] <- A_nf[chosen_frac] + par_row$lrn_nf * (curr_pred - A_nf[chosen_frac])
      R_nf[chosen_frac] <- R_nf[chosen_frac] + par_row$lrn_nf * (resid - R_nf[chosen_frac])
    }

    # Let simulated choices, rather than observed choices, drive autocorrelation.
    choice_a <- ifelse(choice_sim == 1L, 1, 0)
    choice_b <- ifelse(choice_sim == 2L, 1, 0)
    C[cue_a] <- C[cue_a] + par_row$lrn_c * (choice_a - C[cue_a])
    C[cue_b] <- C[cue_b] + par_row$lrn_c * (choice_b - C[cue_b])
  }

  sub_design
}

# Simulate all full and A_fr-amputated closed-loop paths for the requested draws.
s3_run_closed_loop_amputation <- function(trial_template, subject_params) {
  sub_ids <- sort(unique(trial_template$sub_index))
  draw_ids <- sort(unique(subject_params$draw_id))
  n_f <- max(c(trial_template$fA_ix, trial_template$fB_ix), na.rm = TRUE)
  draw_sims <- vector("list", length(draw_ids))

  # Process one posterior draw at a time, then bind the draw-level trial paths.
  for (i in seq_along(draw_ids)) {
    d <- draw_ids[i]
    draw_params <- subject_params[subject_params$draw_id == d, , drop = FALSE]
    draw_sim_list <- vector("list", length(sub_ids) * 2)
    sim_idx <- 1

    # Simulate each participant under both the full and A_fr-amputated closed-loop policies.
    for (s in sub_ids) {
      sub_design <- trial_template[trial_template$sub_index == s, , drop = FALSE]
      par_row <- draw_params[draw_params$sub_index == s, , drop = FALSE]

      draw_sim_list[[sim_idx]] <- s3_closed_loop_subject_amputation(
        sub_design = sub_design,
        par_row = par_row,
        n_f = n_f,
        aff_fr_multiplier = 1,
        model_label = "Full closed-loop"
      )
      sim_idx <- sim_idx + 1

      draw_sim_list[[sim_idx]] <- s3_closed_loop_subject_amputation(
        sub_design = sub_design,
        par_row = par_row,
        n_f = n_f,
        aff_fr_multiplier = 0,
        model_label = "A_fr amputated closed-loop"
      )
      sim_idx <- sim_idx + 1
    }

    draw_sims[[i]] <- dplyr::bind_rows(draw_sim_list)
  }

  dplyr::bind_rows(draw_sims)
}

# Add simulated stay/switch outcomes for the next presentation of the same cue pair.
# The final presentation of each cue pair has no future comparison and is therefore
# assigned NA.
s3_add_simulated_stay <- function(closed_loop_sim, choice_col = "choice_sim", stay_col = "stay_sim") {
  required_cols <- c("draw_id", "model", "sub_index", "overall_trial_nl", "fA_ix", "fB_ix", choice_col)
  missing_cols <- setdiff(required_cols, names(closed_loop_sim))
  if (length(missing_cols) > 0) {
    stop("Missing columns for simulated stay construction: ", paste(missing_cols, collapse = ", "))
  }

  # Include block when available so identical cue labels cannot accidentally bridge task blocks.
  pair_group_cols <- intersect(c("draw_id", "model", "sub_index", "block", "fA_ix", "fB_ix"), names(closed_loop_sim))

  closed_loop_sim |>
    dplyr::group_by(dplyr::across(dplyr::all_of(pair_group_cols))) |>
    dplyr::arrange(overall_trial_nl, .by_group = TRUE) |>
    dplyr::mutate(
      next_choice_for_pair = dplyr::lead(.data[[choice_col]]),
      "{stay_col}" := as.numeric(.data[[choice_col]] == next_choice_for_pair)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-next_choice_for_pair)
}

# Fit the model-free stay regression and return the two focal log-odds coefficients.
# Subject fixed intercepts are included when more than one subject is present, giving
# the observed and simulated regressions the same fast, stable estimator.
s3_fit_stay_effect_glm <- function(data, stay_col = "stay", source = "Observed", draw_id = NA_integer_) {
  fit_df <- data |>
    dplyr::filter(show_fres == 1, !is.na(.data[[stay_col]]), !is.na(out), !is.na(box_val)) |>
    dplyr::mutate(
      stay_response = as.numeric(.data[[stay_col]]),
      out_z = as.numeric(scale(out)),
      box_val_z = as.numeric(scale(box_val))
    )

  use_subject_fe <- "sub_index" %in% names(fit_df) && dplyr::n_distinct(fit_df$sub_index) > 1
  fit_formula <- if (use_subject_fe) {
    stay_response ~ out_z + box_val_z + factor(sub_index)
  } else {
    stay_response ~ out_z + box_val_z
  }

  # Return explicit missing rows if a simulated draw produces an unstable fit.
  fit <- tryCatch(
    stats::glm(fit_formula, data = fit_df, family = stats::binomial()),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(data.frame(
      draw_id = draw_id,
      source = source,
      predictor = c("Cue outcome", "Box value"),
      estimate = NA_real_,
      se = NA_real_,
      ci_low = NA_real_,
      ci_high = NA_real_,
      n_trial = nrow(fit_df)
    ))
  }

  coef_df <- as.data.frame(summary(fit)$coefficients)
  coef_df$term <- rownames(coef_df)
  term_map <- data.frame(
    term = c("out_z", "box_val_z"),
    predictor = c("Cue outcome", "Box value"),
    stringsAsFactors = FALSE
  )

  coef_df |>
    dplyr::right_join(term_map, by = "term") |>
    dplyr::transmute(
      draw_id = draw_id,
      source = source,
      predictor = predictor,
      estimate = Estimate,
      se = `Std. Error`,
      ci_low = Estimate - stats::qnorm(.975) * `Std. Error`,
      ci_high = Estimate + stats::qnorm(.975) * `Std. Error`,
      n_trial = nrow(fit_df)
    )
}

# Fit the stay-effect GLM separately for each posterior draw and amputation condition.
# Each row of the output is a draw-level log-odds coefficient for cue outcome or box value.
s3_stay_effect_sim_draw_summary <- function(closed_loop_sim, stay_col = "stay_sim") {
  # Accept either the raw closed-loop simulation or a version with stay already added.
  sim_with_stay <- if (stay_col %in% names(closed_loop_sim)) {
    closed_loop_sim
  } else {
    s3_add_simulated_stay(closed_loop_sim, stay_col = stay_col)
  }

  # Fit one regression per posterior draw and amputation condition.
  split_cols <- c("draw_id", "model")
  split_index <- do.call(
    interaction,
    c(sim_with_stay[split_cols], list(drop = TRUE, lex.order = TRUE))
  )

  # Run the GLM on each split dataset and bind the two focal coefficients.
  split(sim_with_stay, split_index) |>
    lapply(function(draw_df) {
      s3_fit_stay_effect_glm(
        data = draw_df,
        stay_col = stay_col,
        source = draw_df$model[1],
        draw_id = draw_df$draw_id[1]
      )
    }) |>
    dplyr::bind_rows()
}

# Collapse draw-level simulated stay effects into posterior predictive intervals.
# The median is used as the point estimate so occasional unstable simulated fits do
# not dominate the plotted summary.
s3_stay_effect_sim_interval_summary <- function(sim_draw_summary) {
  sim_draw_summary |>
    dplyr::filter(!is.na(estimate)) |>
    dplyr::group_by(source, predictor) |>
    dplyr::summarise(
      # Use a temporary median name so the quantiles are computed from the
      # original draw-level coefficient vector, not the just-created median.
      estimate_median = stats::median(estimate),
      ci_low = stats::quantile(estimate, .025),
      ci_high = stats::quantile(estimate, .975),
      n_draw = dplyr::n(),
      .groups = "drop"
    ) |>
    dplyr::rename(estimate = estimate_median)
}

# Plot observed stay-effect coefficients against full and A_fr-amputated simulations.
# Observed effects are shown as bars with Wald CIs; simulated effects are shown as
# posterior predictive median dots.
plot_s3_stay_effect_amputation <- function(observed_summary, sim_summary) {
  # Keep the observed bars and simulated dots in a stable left-to-right order.
  predictor_levels <- c("Cue outcome", "Box value")
  sim_levels <- c("Full closed-loop", "A_fr amputated closed-loop")
  sim_colors <- c(
    "Full closed-loop" = "#2B6CB0",
    "A_fr amputated closed-loop" = "#B83227"
  )

  # Apply factor levels before plotting so both data sources align on the same x-axis.
  observed_plot <- observed_summary |>
    dplyr::mutate(predictor = factor(predictor, levels = predictor_levels))
  sim_plot <- sim_summary |>
    dplyr::mutate(
      predictor = factor(predictor, levels = predictor_levels),
      source = factor(source, levels = sim_levels)
    )

  # Offset the simulated model dots just enough that they are distinct from the bar tops.
  sim_dodge <- ggplot2::position_dodge(width = .15)

  # Draw observed coefficients as bars, then overlay the simulated posterior medians.
  ggplot2::ggplot() +
    ggplot2::geom_col(
      data = observed_plot,
      ggplot2::aes(x = predictor, y = estimate),
      fill = "grey82",
      color = "grey35",
      width = .56
    ) +
    ggplot2::geom_errorbar(
      data = observed_plot,
      ggplot2::aes(x = predictor, ymin = ci_low, ymax = ci_high),
      width = .16,
      color = "grey20"
    ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = .35, color = "black") +
    ggplot2::geom_point(
      data = sim_plot,
      ggplot2::aes(x = predictor, y = estimate, color = source, group = source),
      position = sim_dodge,
      size = 3
    ) +
    ggplot2::scale_color_manual(values = sim_colors) +
    ggplot2::labs(
      x = "",
      y = NULL,
      color = ""
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      legend.position = "bottom",
      axis.line.x = ggplot2::element_blank(),
      axis.ticks.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "#E6E6E6", linewidth = .35)
    )
}
