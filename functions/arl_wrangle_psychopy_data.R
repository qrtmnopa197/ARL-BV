#This function takes a single subject's raw psychopy data and reformats it into a long dataset usable for analysis. 
#It creates all variables of interest that can be created from the raw data alone.
wrangle_psychopy_data <- function(csv_path){
  library(ddpcr)
  library(tidyverse)
  library(sigmoid)
  if(!is.function(get_csvs_to_analyze)){
    stop("earn_rel_cum_earn() is not loaded into the environment of wrangle_psychopy_data(), suggesting that the functions in s22_utilities.R are not present in the function's environment.")
  }
  print(csv_path) #print the csv path for debugging purposes - so you know what subject you're on
  df_full <- read_csv(csv_path) #read in subject's data
  #if the subject never pressed the right arrow key on affect probes, the RT columns won't appear. Thus, you need to create fully empty RT columns for them.
  if(!("key_resp_2.rt" %in% names(df_full))){
    df_full$key_resp_2.rt <- NA
  }
  if(!("key_resp_5.rt" %in% names(df_full))){
    df_full$key_resp_5.rt <- NA
  }
  #make one df with all the subject-level info, and a second with all the trial-level info
  sub_info <- select(df_full, any_of(c(id = "participant", "age_gender_formresponse", "date", "earnings_before_ec", acheck = "acheck_pressq.keys", 
                                       fcheck = "feed_check_pressq.keys", total_experiment_time = "total_experiment_time", 
                                       choice_noprobe_pt = "mouse.clicked_name", choice_probe_pt = "mouse_4.clicked_name", 
                                       noprobe_qs = "mouse_3.clicked_name", probe_qs = "mouse_7.clicked_name")), 
                     contains("key_resp"))
  trials <- select(df_full, any_of(c(pair = "pair_desc", fractalA = "fa_desc", 
                                     fractalB = "fb_desc", fractalA_img = "fractalA", fractalB_img = "fractalB",
                                     x_fa = "TMPx_fa", dec_affect_probe = "dec_probe", feed_affect_probe = "feed_probe", pair_onset_time = "pair_onset_time",
                                     dec_rt = "mouse_2.rt", dec_time = "dec_choicetime", choice = "mouse_2.clicked_name",
                                     dec_probe_onset_time = "dec_probe_onset_time", valence_dec = "valence_slider_dec.response",
                                     arousal_dec = "arousal_slider_dec.response", dec_probe_rt = "key_resp_2.rt", too_slow = "too_slow",
                                     feed_onset_time = "feedback_onset_time", feed_probe_onset_time = "feed_probe_onset_time", outcome_a = "outcome_a",
                                     outcome_b = "outcome_b", chosen_outcome = "chosen_outcome", valence_feed = "valence_slider_feed.response",
                                     arousal_feed = "arousal_slider_feed.response", feed_probe_rt = "key_resp_5.rt", 
                                     trial_raw = "trials.thisN", block_raw = "blocks.thisN", makeup_repetition = "trials_and_makeups.thisN")),
                   ends_with(".ran"))
  
  #FIRST, REFORMAT THE TRIALS DF...
  
  #Find the rows containing the first and last trials, getting rid of everything above and below those...
  practice_q_runs <- which(trials$probe_practice_qs.ran == 1) #identify all rows on which probe_practice_qs - the last loop before the main task - ran
  main_task_start <- max(practice_q_runs) + 1 #the row after the last run of this loop should be the first trial
  block_runs <- which(trials$blocks.ran == 1) #identify all rows on which the blocks loop ran. 
  main_task_end <- max(block_runs) #the last run of this loop signifies the end of the main task
  trials <- trials[main_task_start:main_task_end,] #keep only rows within the main task
  
  #Create columns for block number and what makeup/repetition a trial represents
  #Block number and makeup/repetition number are output at the end of each block/repetition loop,
  #so you need to fill out all the rows above each loop end
  #with that value. You can use the custom "fill_vec" function to do this (see s22_functions.R in this folder)
  trials$block_raw <- fill_vec(trials$block_raw, bottom_up = TRUE)
  trials$makeup_repetition <- fill_vec(trials$makeup_repetition, bottom_up = TRUE)
  
  # You want one row per trial, meaning you want to ratchet down one row at the end of every trial. 
  # Fortunately, the end of each trial is marked by a loop end, so psychopy does indeed ratchet down one row
  # at the end of every trial. However, there are a few loops that sometimes end in between trials - 
  #  i.e., before the trial starts, with the first trial always starting on the right row - 
  # leading to an unnecessary ratchet-down before the start of the next trial, before any trial data has been collected.
  # To resolve this issue, you can simply delete all rows on which these loops have run/the pre-trial ratchet-down 
  # has occurred (since no trial data is collected before the ratchet-down).
  # A second issue is that there’s a loop which sometimes ends in the middle of a trial (remaining_trial), 
  # which results in a single trial’s data being spread across two lines. 
  # To address this, you should identify trials in which there was a ratchet-down mid-trial. 
  # In these cases, you know that the current row and the row below it represent a single trial’s data, 
  # and that you need to combine them into one row. The simplest way to do this is to copy the trial data
  # on the second row to the first row, thus ensuring that the first row contains the full trial’s data. 
  # Then, delete the second row, which contains only redundant information.
  trials$delete_row <- 0
  for(row in 1:nrow(trials)){
    #On rows where a pre-trial ratchet-down has occurred...
    if (!is.na(trials$trials_and_makeups.ran[row]) | !(is.na(trials$feed_check.ran[row])) | !(is.na(trials$affect_check.ran[row])) | !(is.na(trials$blocks.ran[row]))){
      trials$delete_row[row] <- 1 #mark row for deletion
    } else if (!(is.na(trials$remaining_trial.ran[row]))){
      #Mid-trial ratchet-down. First, copy all the valuable data from the row below to the row on which the trial started.
      row_below <- row + 1
      trials$trial_raw[row] <- trials$trial_raw[row_below]
      trials$feed_affect_probe[row] <- trials$feed_affect_probe[row_below]
      trials$dec_affect_probe[row] <- trials$dec_affect_probe[row_below]
      trials$fractalB_img[row] <- trials$fractalB_img[row_below]
      trials$fractalA_img[row] <- trials$fractalA_img[row_below]
      trials$fractalB[row]<- trials$fractalB[row_below]
      trials$fractalA[row] <- trials$fractalA[row_below]
      trials$pair[row] <- trials$pair[row_below]
      if(any(names(trials) == "too_slow")){
        trials$too_slow[row] <- trials$too_slow[row_below]
      }
      #Now, mark the row below (containing only redundant data) for deletion
      trials$delete_row[row_below] <- 1 
    }
  }
  trials <- filter(trials, delete_row == 0) #actually delete the marked rows
  trials <- select(trials, -delete_row,-ends_with(".ran")) #get rid of .ran and delete_row columns
  
  
  #Do some variable recoding...
  
  #X_fa is output every time it's reset, and is NA otherwise. So first, we need to fill out the NA values to match
  #the most recent reset - which can be done with the custom "fill_vec" function (see s22_functions.R in this folder)
  filled_xfa <- fill_vec(trials$x_fa)
  #recode x_fa values to left/right, assigning to a new column
  trials$fractalA_side <- ""
  #for each row in trials...
  for(row in 1:nrow(trials)){
    #if xfa is negative, fractal A is on the left; if positive, it's on the right
    if(filled_xfa[row] < 0){
      trials$fractalA_side[row] <- "left"
    } else if (filled_xfa[row] > 0){
      trials$fractalA_side[row] <- "right"
    } else{
      stop("x_fa has value of 0")
    }
  }
  trials <- select(trials,-x_fa) #don't need x_fa anymore
  
  #if the data has a "too slow" column in it...
  if(any(names(trials) == "too_slow")){
    #if the subject was "too slow" on a trial, mark the trial as late
    for(row in 1:nrow(trials)){
      if(trials$too_slow[row] == "TRUE"){
        trials$choice[row] <- "late"
      }
    }
    #this may leave some N/A choices, which are trials on which the subject was not late, but the choice wasn't recorded.
    #for these, get the choice based on the outcomes.
    for(row in 1:nrow(trials)){
      if(is.na(trials$choice[row])){
        if(trials$chosen_outcome[row] == trials$outcome_a[row]){
          trials$choice[row] <- "fractal_a_dec"
        } else if(trials$chosen_outcome[row] == trials$outcome_b[row]){
          trials$choice[row] <- "fractal_b_dec"
        }
      }
    }
    trials <- select(trials,-too_slow) #don't need this anymore
  } else{
    #Recode N/A choices as late. However, if a chosen outcome was registered, then this N/A was a glitch, and you should fill in the proper choice
    for(row in 1:nrow(trials)){
      if(is.na(trials$choice[row]) & is.na(trials$chosen_outcome[row])){
        trials$choice[row] <- "late"
      } else if(is.na(trials$choice[row]) & !is.na(trials$chosen_outcome[row])){
        if(trials$chosen_outcome[row] == trials$outcome_a[row]){
          trials$choice[row] <- "fractal_a_dec"
        } else if(trials$chosen_outcome[row] == trials$outcome_b[row]){
          trials$choice[row] <- "fractal_b_dec"
        }
      }
    }
    #for rows where a choice was recorded, but no outcomes were shown, mark the trial as late. On these occasions, the participant right-clicked
    #on a fractal, which typically causes the trial to immediately drop out to the "too-slow" screen, although a choice is registered.
    for(row in 1:nrow(trials)){
      if((trials$choice[row] == "fractal_a_dec" ||  trials$choice[row] == "fractal_b_dec") && (is.na(trials$outcome_a[row]) || is.na(trials$outcome_b[row]))){
        trials$choice[row] <- "late"
      }
    }
  }
  
  #python indexing starts at 0, so this rectifies that
  trials$trial <- trials$trial_raw + 1
  trials$block <- trials$block_raw + 1
  trials <- select(trials,-c(trial_raw,block_raw)) #don't need raw rows anymore
  
  #Trials actually give you the trial within the trials_and_makeups loop, not the trial within the block.
  #You can rectify this by going through the df top to bottom, making each repetition trial equal to the preceding trial plus one
  #NB:This only works assuming the trials are in ascending order and there are no gaps between trials. It's very difficult to do this
  #in a more principled way, so you'll just have to go with this for now.
  for(row in 1:nrow(trials)){
    #if you're on a repetition
    if(trials$makeup_repetition[row] != 0){
      trials$trial[row] <- trials$trial[row-1] + 1 #set the trial value to one more than the previous trial
    }
  }
  
  #Add columns for the number associated with each fractal image file. These will be used to label fractals in RL model-fitting
  trials <- mutate(trials, fractal_a_num = str_extract(fractalA_img,"\\d+"), fractal_b_num = str_extract(fractalB_img,"\\d+"))
  
  #create a column for the unchosen outcome
  trials$unchosen_outcome <- ""
  #go row by row...
  for(row in 1:nrow(trials)){
    #identify the choice - a, b, or late - and fill logically using the outcome_a/b value
    if(trials$choice[row] =="fractal_a_dec"){
      trials$unchosen_outcome[row] <- trials$outcome_b[row]
    } else if(trials$choice[row] =="fractal_b_dec"){
      trials$unchosen_outcome[row] <- trials$outcome_a[row]
    } else if(trials$choice[row] == "late"){
      trials$unchosen_outcome[row] <- NA
    }
  }
  trials$unchosen_outcome <- as.numeric(trials$unchosen_outcome) #make numeric
  
  trials$row_index <- c(1:nrow(trials)) #get the row numbers for this subject, which will be useful for anumber of things to follow
  
  trials <- earn_rel_cum_earn(trials,"sub") #get cumulative earnings and trial earnings relative to past earnings using a custom funciton
  blockwise_dfs <- by(trials,trials$block,earn_rel_cum_earn,data_type="block",simplify=FALSE) #do the same thing, but broken up by block
  trials <- do.call(rbind,blockwise_dfs) #combine the multiple dfs created by the above back into a single df
  
  #create a "stay" column, indicating whether the participant made the same choice the next time the pair was presented (if so, 1) or not (if so, 0)
  stay_data <- filter(trials,choice != "late") #just look at non-late choices (coding this is tricky if you look at everything)
  stay_data$stay <- NA #initialize stay column
  #go trial-by-trial...
  for(row in 1:nrow(stay_data)){
    pair_rows <- which(stay_data$pair == stay_data$pair[row]) #identify the rows on which the participant played with the same pair
    #if there are future rows on which they played with the same pair...
    if(any(pair_rows > row)){
      future_pres <- pair_rows[pair_rows > row] #identify all future trials/rows in which this pair was presented
      next_pres <- min(future_pres) #identify the next choice with this pair by taking the min of this vector
      if(stay_data$choice[next_pres] == stay_data$choice[row]){
        stay_data$stay[row] <- 1 #if the choice at the next presentation is the same as the current choice, the participant "stayed"
      } else if(stay_data$choice[next_pres] != stay_data$choice[row]){
        stay_data$stay[row] <- 0 #otherwise they switched
      }
    } 
  }
  stay_join <- select(stay_data,row_index,stay)
  trials <- left_join(trials,stay_join,by="row_index")
  
  
  #get a trial-by-trial measure of absolute regret - just the difference between the unchosen and chosen outcome
  trials$abs_regret <- ""
  for(row in 1:nrow(trials)){
    if(trials$choice[row] == "fractal_a_dec"){
      trials$abs_regret[row] <- trials$outcome_b[row] - trials$outcome_a[row]
    } else if(trials$choice[row] == "fractal_b_dec"){
      trials$abs_regret[row] <- trials$outcome_a[row] - trials$outcome_b[row]
    }
  }
  
  #Get the side that was chosen on each trial
  trials$side_chosen <- NA #initalize side chosen variable
  for(row in 1:nrow(trials)){
    #only update side chosen if a choice was made
    if(trials$choice[row] != "late"){
      #if fractal a was chosen, then assign side chosen with the side A is on
      if(trials$choice[row] == "fractal_a_dec"){
        if(trials$fractalA_side[row] == "left"){
          trials$side_chosen[row] <- "left"
        } else if(trials$fractalA_side[row] == "right"){
          trials$side_chosen[row] <- "right"
        }
        #if B was chosen, assign side chosen to the opposite of the side A is on
      } else if(trials$choice[row] == "fractal_b_dec"){
        if(trials$fractalA_side[row] == "left"){
          trials$side_chosen[row] <- "right"
        } else if(trials$fractalA_side[row] == "right"){
          trials$side_chosen[row] <- "left"
        }
      }
    }
  }
  trials$side_chosen_numeric <- ifelse(trials$side_chosen=="left",1,0) #convert the side column to a numeric variable
  
  #determine whether the correct choice was made on each trial
  dist_evs <- read.csv(paste0(path_to_project_directory,"code/dist_evs.csv")) #read in key identifying the EV of each fractal
  #create a column with EVs for fractal A, on each trial
  dist_evs_a <- rename(dist_evs,fractalA=f_desc,fractal_a_ev=ev) 
  trials <- left_join(trials,dist_evs_a,by="fractalA")
  #create a column with EVs for fractal B, on each trial
  dist_evs_b <- rename(dist_evs,fractalB=f_desc,fractal_b_ev=ev)
  trials <- left_join(trials,dist_evs_b,by="fractalB")
  trials$correct <- apply(trials, 1, mark_correct_incorrect) #create a column with 1s for all correct choices and 0s for all incorrect ones, based on the difference in EVs
  
  
  
  
  #NOW, THE SUBJECT INFO DF...
  age <- sub_info$age_gender_formresponse[grepl("\\d+",sub_info$age_gender_formresponse)] #identify row with the age by searching for a digit value
  gender <- sub_info$age_gender_formresponse[grepl("[a-zA-Z]+",sub_info$age_gender_formresponse)] #identify row with gender by searching for letter value
  if(length(age)==0){
    age <- 999
  }
  if(length(gender)==0){
    gender <- 999
  }
  earnings_before_ec <- sub_info$earnings_before_ec[grepl("\\d+",sub_info$earnings_before_ec)] #identify row with total earnings by searching for a digit value
  total_experiment_time <- sub_info$total_experiment_time[grepl("\\d+",sub_info$total_experiment_time)] #same strategy for total experiment time
  #determine whether the subject passed or failed the affect and feedback attention checks by determining whether they pressed "q" for each attention check - as they are supposed to
  #first replace NA values with blanks so the "any" function works more smoothly - see this: https://stackoverflow.com/questions/42226525/why-does-any-return-na-when-no-true-values
  acheck_column <- sub_info$acheck 
  fcheck_column <- sub_info$fcheck
  acheck_column[is.na(acheck_column)] <- ""
  fcheck_column[is.na(fcheck_column)] <- ""
  #now do the pass/fail assignment
  if(any(acheck_column == "q")){
    affect_attention_check <- "pass"
    affect_attention_check_numeric <- 1
  } else{
    affect_attention_check <- "fail"
    affect_attention_check_numeric <- 0
  }
  if(any(fcheck_column == "q")){
    feed_attention_check <- "pass"
    feed_attention_check_numeric <- 1
  } else{
    feed_attention_check <- "fail"
    feed_attention_check_numeric <- 0
  }
  att_checks_passed <- affect_attention_check_numeric + feed_attention_check_numeric #get total checks passed
  
  id <- sub_info$id[1] #get id from one of the rows
  date <- sub_info$date[1] #get date from one of the rows
  #check the longest string of late choices the subject had. Long strings suggest that the participant left the task
  value_runs <- rle(trials$choice) #get the runs of each value
  #if late is one of the choice values...
  if(any(value_runs$values == "late")){
    late_indices <- which(value_runs$values == "late") #get the indices that correspond to runs of "late"
    late_consecutive_runs <- max(value_runs$lengths[late_indices]) #get the longest "late" run
  } else{
    late_consecutive_runs <- 0 #if none of the runs are "late"s, then the subject had no late choices
  }
  #get the standard deviations of valence and arousal ratings at decision and feedback
  sd_valence_dec <- sd(trials$valence_dec,na.rm=TRUE) 
  sd_arousal_dec <- sd(trials$arousal_dec,na.rm=TRUE)
  sd_valence_feed <- sd(trials$valence_feed,na.rm=TRUE)
  sd_arousal_feed <- sd(trials$arousal_feed,na.rm=TRUE)
  
  #Get attention-checking metrics from the instructions slides
  #Figure out total number of automatic processions during the instructions slides, and consecutive number
  inst_slides <- select(sub_info,contains("key_resp")) %>% select(contains(".keys")) %>% select(-contains(c("resp_5.","resp_2.","resp_36.","resp_35.","resp_39.","resp_30.","resp_31."))) #select columns with instructions slides
  inst_resp <- as.vector(apply(inst_slides,2,function(x) any(!is.na(x)))) #get a vector saying whether each column is blank or not. Blank columns correspond to no key being pressed on the corresponding screen
  process_counts <- rle(inst_resp) #get counts of consecutive sequence of trues and falses
  if(any(process_counts$values == FALSE)){
    auto_process_indices <- which(process_counts$values == FALSE) #get the indices that correspond to runs of automatic processions
    consecutive_auto_process <- max(process_counts$lengths[auto_process_indices]) #get the longest run of automatic processions
  } else{
    consecutive_auto_process <- 0 #if none of the runs are "late"s, then the subject had no late choices
  }
  total_auto_process <- length(which(inst_resp==FALSE)) #get the total number of automatic processions during the instructions
  
  #Get the number of times the subject made a choice on the no-probe practice trials, so as to assess whether they completed that part
  if("noprobe_qs" %in% names(sub_info)){
    #count number of right and wrong answers
    noprobe_answers_correct <- length(which(sub_info$noprobe_qs=="poly_true_dec"))
    noprobe_answers_incorrect <- length(which(sub_info$noprobe_qs=="poly_false_dec"))
  } else{
    #if the column doesn't exist, it means they didn't provide answers
    noprobe_answers_correct <- 0
    noprobe_answers_incorrect <- 0
  }
  if("probe_qs" %in% names(sub_info)){
    #count number of right and wrong answers
    probe_answers_correct <- length(which(sub_info$probe_qs=="poly_true_probe"))
    probe_answers_incorrect <- length(which(sub_info$probe_qs=="poly_false_probe"))
    #total_auto_process <- table(inst_resp)[names(table(inst_resp))==FALSE] #get the total number of automatic processions during the instructions
  } else{
    #if the column doesn't exist, it means they didn't provide answers
    probe_answers_correct <- 0
    probe_answers_incorrect <- 0
  }
  answers_correct <- noprobe_answers_correct + probe_answers_correct
  answers_incorrect <- noprobe_answers_incorrect + probe_answers_incorrect
  
  if("choice_noprobe_pt" %in% names(sub_info)){
    noprobe_pt_choices <- sum(!is.na(sub_info$choice_noprobe_pt))
  } else{
    noprobe_pt_choices <- 0
  }
  #diddo the with-probe practice trials
  
  if("choice_probe_pt" %in% names(sub_info)){
    probe_pt_choices <- sum(!is.na(sub_info$choice_probe_pt))
  } else{
    probe_pt_choices <- 0
  }
  
  #get means of several trial-level variables
  mean_dec_rt <- mean(trials$dec_rt, na.rm=TRUE) #grab mean choice RT 
  mean_dec_probe_rt <- mean(trials$dec_probe_rt, na.rm=TRUE) #get mean RT to decision affect probes
  mean_feed_probe_rt <- mean(trials$feed_probe_rt, na.rm=TRUE) #get mean RT to feedback affect probes
  percent_correct <- mean(trials$correct, na.rm=TRUE) #get the percentage correct 
  percent_left <- mean(trials$side_chosen_numeric, na.rm=TRUE) #get the percentage of trials on which the left fractal was chosen
  
  #RUN MODEL-FREE PREDICTIONS OF CHOICE BASED ON PAST OUTCOMES. ADDITIONALLY, CREATE COLUMNS THAT WILL ALLOW YOU TO RUN MODEL-FREE PREDICTIONS OF CHOICE BASED ON AFFECT.
  filt_data <- filter(trials, choice != "late") #filter out late trials
  #get the fractals being used
  fanums <- unique(filt_data$fractal_a_num)
  fbnums <- unique(filt_data$fractal_b_num)
  fnums <- append(fanums,fbnums) %>% unique(.)
  for(num in fnums){
    #Create columns for fractal cumulative totals, presentation number, valence rating total following selection of that fractal (at feedback), and the number of feedback affect probes following selection of the fractal
    total_col <- paste0("f_c_total_",num)
    pres_num <- paste0("f_pres_num_",num)
    filt_data[total_col] <- 0
    filt_data[pres_num] <- 0
  }
  #initialize columns..
  filt_data$av_a_val <- 0 #for the average earning of fractal a on that trial
  filt_data$av_b_val <- 0 #"" of fractal b
  filt_data$a_pres_num <- 0 #the number of presentation of fractal a (and thus the pair) on that trial
  filt_data <- mutate(filt_data, numerical_choice = ifelse(choice == "fractal_a_dec", 1,0)) #make choice into a 1/0 variable, with 1 for fractal a chosen
  for(row in 1:(nrow(filt_data)-1)){
    #fill out the row below
    row_below <- row + 1
    filt_data[row_below, grep("f_",names(filt_data))] <- filt_data[row, grep("f_",names(filt_data))] #make the values for all fractal totals and presentation numbers the same as the current row
    filt_data[row_below, paste0("f_c_total_",filt_data$fractal_a_num[row])] <- filt_data[row, paste0("f_c_total_",filt_data$fractal_a_num[row])] + filt_data$outcome_a[row] #identify fa, and add outcome_a to that value
    filt_data[row_below, paste0("f_pres_num_",filt_data$fractal_a_num[row])] <- filt_data[row, paste0("f_pres_num_",filt_data$fractal_a_num[row])] + 1 #then add 1 to the trial presentation number
    #same for b
    filt_data[row_below, paste0("f_c_total_",filt_data$fractal_b_num[row])] <- filt_data[row, paste0("f_c_total_",filt_data$fractal_b_num[row])] + filt_data$outcome_b[row] 
    filt_data[row_below, paste0("f_pres_num_",filt_data$fractal_b_num[row])] <- filt_data[row, paste0("f_pres_num_",filt_data$fractal_b_num[row])] + 1
    
    #now, fill out the average a and b values for the row below by dividing the cumulative total for fractal a/b by the number of presentations for fractal a/b
    filt_data$av_a_val[row_below] <-  filt_data[row_below, paste0("f_c_total_",filt_data$fractal_a_num[row_below])]/filt_data[row_below, paste0("f_pres_num_",filt_data$fractal_a_num[row_below])]
    filt_data$av_b_val[row_below] <-  filt_data[row_below, paste0("f_c_total_",filt_data$fractal_b_num[row_below])]/filt_data[row_below, paste0("f_pres_num_",filt_data$fractal_b_num[row_below])]
    filt_data$a_pres_num[row_below] <- filt_data[row_below, paste0("f_pres_num_",filt_data$fractal_a_num[row_below])] #get the presentation number for fractal a
  }
  #get rid of NaNs produced by dividing by 0
  filt_data$av_a_val[is.na(filt_data$av_a_val)] <- 0 
  filt_data$av_b_val[is.na(filt_data$av_b_val)] <- 0 
  #make these columns numeric ahead of numerical operations
  filt_data$av_a_val <- as.numeric(filt_data$av_a_val)
  filt_data$av_b_val <- as.numeric(filt_data$av_b_val)
  filt_data$a_pres_num <- as.numeric(filt_data$a_pres_num)
  #get the difference between the averages
  filt_data <- mutate(filt_data,av_val_diff = av_a_val-av_b_val)
  #get coefficients from the model
  sub_fit <- glm(numerical_choice ~ av_val_diff, data=filt_data, family = "binomial")
  av_outcome_diff_p <- coef(summary(sub_fit))[2,4]
  av_outcome_diff_beta <- coef(summary(sub_fit))[2,1] #assign the av_val_diff main effect p to the appropriate place in the subject summary stats df
  #add the presentation number to the trials df. You generally want to avoid manipulating the trials df in the subject-level section, but there was no way to avoid it here.
  join_filt_data <- transmute(filt_data,row_index=row_index, pair_pres_num = a_pres_num + 1) #grab columns to join, adding 1 to pres_num so that it starts at 1
  trials <- left_join(trials,join_filt_data,by="row_index")
  
 
  #get the percentage of trials on which the subject didn't make a choice in time
  late_numeric <- ifelse(trials$choice == "late",1,0)
  late_percent <- mean(late_numeric)
  
  #get the total number of trials completed
  trials_completed <- sum(trials$choice != "late")
  
  #get the percentage of affect probes on which the subject provided a response
  affect_cols <- c("valence_dec","arousal_dec","valence_feed","arousal_feed") #cols with affect ratings
  affect_response_nums <- lapply(affect_cols,get_affect_response_numbers,trials) #pass each col, and the df they came from (w/o late trials), to a function that returns a column vector with the number
                                                                            #of probes they did and didn't respond to
  affect_response_df<- do.call(cbind,affect_response_nums) #reformat these as a single data frame
  probe_skipped_percent <- sum(affect_response_df[1,])/sum(affect_response_df[2,]) #get the overall percent of affect ratings missed
  
  
  sub_info_final <- data.frame(id,date,age,gender,earnings_before_ec,affect_attention_check,feed_attention_check,affect_attention_check_numeric,feed_attention_check_numeric,att_checks_passed,
                               consecutive_late_choices=late_consecutive_runs,probe_skipped_percent,trials_completed,
                               val_dec_skipped_percent = affect_response_df$valence_dec[3],val_feed_skipped_percent = affect_response_df$valence_feed[3],
                               arousal_dec_skipped_percent = affect_response_df$arousal_dec[3],arousal_feed_skipped_percent = affect_response_df$arousal_feed[3],
                               sd_valence_dec,sd_arousal_dec,sd_valence_feed,sd_arousal_feed,total_experiment_time,consecutive_auto_process,total_auto_process,
                               mean_dec_rt,mean_dec_probe_rt,mean_feed_probe_rt, percent_correct,percent_left,late_percent,
                               noprobe_pt_choices,probe_pt_choices,answers_correct,answers_incorrect,av_outcome_diff_p,av_outcome_diff_beta)
  
  base_data <- cbind(sub_info_final,trials) #merge the two dfs
  
  #rearrange the column order
  col_order <- c("id","date","age","gender","earnings_before_ec","total_experiment_time","consecutive_auto_process","total_auto_process",
                 "answers_correct","answers_incorrect","noprobe_pt_choices", "probe_pt_choices",
                 "percent_correct","percent_left","probe_skipped_percent","val_dec_skipped_percent","val_feed_skipped_percent",
                 "arousal_dec_skipped_percent","arousal_feed_skipped_percent",
                 "mean_dec_rt","mean_dec_probe_rt","mean_feed_probe_rt","block","makeup_repetition",
                 "trial","pair","fractalA","fractal_a_ev","fractalA_img","fractalB","fractal_b_ev","fractalB_img","fractal_a_num","fractal_b_num","fractalA_side","pair_pres_num",
                 "dec_affect_probe","feed_affect_probe","pair_onset_time","dec_rt","dec_time","choice", "side_chosen","side_chosen_numeric",
                 "stay","correct","dec_probe_onset_time","dec_probe_rt","valence_dec","arousal_dec","feed_onset_time", 
                 "outcome_a","outcome_b","chosen_outcome","unchosen_outcome","cumulative_earnings_sub","feed_probe_onset_time","cumulative_earnings_block",
                 "feed_probe_rt","valence_feed","arousal_feed","affect_attention_check","feed_attention_check","affect_attention_check_numeric",
                 "feed_attention_check_numeric","att_checks_passed","consecutive_late_choices","late_percent","trials_completed",
                 "sd_valence_dec","sd_arousal_dec","sd_valence_feed","sd_arousal_feed", "abs_regret","earn_rel_past_sub",
                 "earn_rel_past_block","av_outcome_diff_beta","av_outcome_diff_p","row_index")
  #manually check that col_order contains the same strings as the names of full_data
  if(length(setdiff(col_order,names(base_data))) > 0 || length(setdiff(names(base_data),col_order))){
    print(paste0("In col_order but not base_data: ",setdiff(col_order,names(base_data))))
    print(paste0("In base_data but not col_order: ",setdiff(names(base_data),col_order)))
    stop("Mismatch between names in base_data and col_order")
  }

  base_data <- base_data[,col_order] #rearrange column order
  
  return(list(base_data,sub_info_final)) #return the trial-level data and subject-level-data
}