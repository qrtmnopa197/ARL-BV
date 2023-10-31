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
  
  #Create indices from 1:n_f for each fractal image. To do this, first create a mini-df with one column having all the fA_img values and the other two
  #columns having indices for fA and fB. This assumes that every fA_img is paired with a unique fB_img.
  f_index_df <- data.frame(fA_img = unique(trials$fA_img),fA_ix = 1:length(unique(trials$fA_img)),fB_ix = (1:length(unique(trials$fA_img))+length(unique(trials$fA_img))))
  trials <- left_join(trials,f_index_df,by="fA_img")
  
  #Get matrics with subjects on the rows and trials in the columns, containing indices for fractal A/B
  fA <- sub_by_trial_matrix(trials,"fA_ix")   
  fB <- sub_by_trial_matrix(trials,"fB_ix")
  
  choice <- sub_by_trial_matrix(trials,"choice_numeric") #ditto whether fractal A (1) or B(2) was chosen  

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
    bv = bv,
    fres = fres,
    rat = rat,
    rat_num = rat_num,
    n_rat = n_rat,
    prev_rat = prev_rat
  )
  return(data)
}