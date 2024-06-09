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