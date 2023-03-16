#takes a trial-level df with no lates, along with an fA image number that specifies a particular pair,
#and returns a two-column df with the pair presentation number to put in the relevant row indices
create_pair_pres_nums <- function(fa_im,df){
  df <- filter(df,fA_img == fa_im) #get trials with that pair
  df$pair_pres_num <- c(1:nrow(df)) #number the rows/trials
  df_join <- df %>% select(row_index, pair_pres_num)
  return(df_join)
}

#Takes a trial-level df and a block and pair type specifier.
#Returns the difference in average affect between the high-valued and low-valued fractal.
diff_av_val <- function(df,blocknum,type){
  df <- df %>% filter(block == blocknum) %>% filter(pair_type == type) #get the data for a particular block and pair type
  
  #get the mean valence rating on trials where the higher-valued fractal was chosen
  df_hi <- df %>% filter(hi_chosen == 1) 
  hi_av <- mean(df_hi$val_rat,na.rm=T)
  #if hi was never chosen, then the average affect is 0
  if(nrow(df_hi) == 0){
    hi_av <- 0
  }
  
  #get the mean valence rating on trials where the lower-valued fractal was chosen
  df_lo <- df %>% filter(hi_chosen == 0) 
  lo_av <- mean(df_lo$val_rat,na.rm=T)
  #if lo was never chosen, then the average lo affect is 0
  if(nrow(df_lo) == 0){
    lo_av = 0
  }
 
  to_ret <- hi_av - lo_av #get the difference
  
  return(to_ret)
}

#Returns a list of data for input to stan, given trial-level data
#trials: trial-level data
#n_t: number of trials; must be set manually
#n_f: number of fractls. If NULL, it will be set to however many distinct fractal numbers are found in the data.
stan_data_arl <- function(trials,n_t){
  library(gdata)
  
  n_s <- length(unique(trials$id)) #get number of subjects
  
  frac_vec <- as.vector(as.matrix(trials[c("fA_img","fB_img")])) #get a big vector of fractals...
  n_f <- length(unique(frac_vec)) #and then get the number of fractals
  
  #Create indices from 1:n_f for each fractal image. To do this, first create a mini-df with one column having all the fA_img values and the other two
  #columns having indices for fA and fB. This assumes that every fA_img is paired with a unique fB_img.
  f_index_df <- data.frame(fA_img = unique(trials$fA_img),fA_ix = 1:length(unique(trials$fA_img)),fB_ix = (1:length(unique(trials$fA_img))+length(unique(trials$fA_img))))
  trials <- left_join(trials,f_index_df,by="fA_img")
  #Get matrics with subjects on the rows and trials in the columns, containing indices for fractal A/B
  #NB: Because these fractal indices are based on the fractal image, they do not refer to the same fractal across different subjects
  
  fA <- sub_by_trial_matrix(trials,"fA_ix")   
  fB <- sub_by_trial_matrix(trials,"fB_ix")
  
  pres_num <- sub_by_trial_matrix(trials,"pair_pres_num") #diddo the number of times a particular pair has been presented
  
  choice <- sub_by_trial_matrix(trials,"choice_numeric") #diddo whether fractal A (1) or B(2) was chosen  
  

  reward <- sub_by_trial_vec_list(trials,"reward") #get a list of vectors, one per subject, each containing the reward for each trial
  fake_out <- sub_by_trial_vec_list(trials,"fake_out") #diddo fake_out
  affect <- sub_by_trial_vec_list(trials,"valrat_z") #diddo valence ratings
  aff_resid <- sub_by_trial_vec_list(trials,"resid")
  
  
  data <- list(
    n_t = n_t,
    n_s = n_s,
    n_f = n_f,
    fA = fA,
    fB = fB,
    reward = reward,
    affect = affect,
    choice = choice,
    fake_out = fake_out,
    aff_resid = aff_resid
  )
  return(data)
}