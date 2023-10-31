#This function takes a single subject's raw psychopy data and reformats it into a long dataset usable for analysis. 
#It creates all variables of interest that can be created from the raw data alone.
arlbv_wrangle_psychopy_data <- function(csv_path){
  library(ddpcr)
  library(tidyverse)
  library(sigmoid)
  if(!is.function(get_csvs_to_analyze)){
    stop("get_csvs_to_analyze() is not loaded into the environment of arlbv_wrangle_psychopy_data(), suggesting that the functions in s22_utilities.R are not present in the function's environment.")
  }
  print(csv_path) #print the csv path for debugging purposes - so you know what subject you're on
  df_full <- read_csv(csv_path, show_col_types = F) #read in subject's data
  #if the subject never pressed the right arrow key on affect probes, the RT columns won't appear. Thus, you need to create fully empty RT columns for them.
  if(!("val_rat.rt" %in% names(df_full))){
    df_full$val_rat.rt <- NA
  }
  #make one df with all the subject-level info, and a second with all the trial-level info
  sub_info <- select(df_full, any_of(c(id = "participant", "age_gender_formresponse", "date", "earnings_before_ec", "earnings",total_experiment_time = "total_experiment_time", 
                                       comp_qs = "mouse_3.clicked_name",instruct_keypress = "instruct_key.keys",
                                       bv_check = "pressq_bv.keys", rev_check = "pressq_out_2.keys", val_check = "pressq_val_atchk.keys"))) 
  trials <- select(df_full, any_of(c("fA_img","fB_img", dec_rt = "mouse.rt", too_slow = "too_slow",choice = "click",
                                     "out",box_val="block_val","show_fres",
                                     val_rat_raw = "val_slider_5.response", val_rat_rt = "val_rat.rt", 
                                     trial_raw = "trials.thisN", block_raw = "blocks.thisN", makeup_repetition = "trials_and_makeups.thisN")),
                   ends_with(".ran"))
  
  #FIRST, REFORMAT THE TRIALS DF...
  
  #Find the rows containing the first and last trials, getting rid of everything above and below those...
  aff_qs_run <- which(trials$aff_qs.ran == 1) #identify the rows on which aff_qs - the last loop before the main task - ran
  main_task_start <- max(aff_qs_run) + 1 #the row after the last run of this loop should be the first trial
  block_runs <- which(trials$blocks.ran == 1) #identify all rows on which the blocks loop ran. 
  main_task_end <- max(block_runs) #the last run of this loop signifies the end of the main task
  trials <- trials[main_task_start:main_task_end,] #keep only rows within the main task
  
  #Create columns for block number and makeup/repetition a trial represents, as well as what the fA_img and fB_img were
  #These things are all are output at the end of each block/repetition loop,
  #so you need to fill out all the rows above each loop end
  #with that value. You can use the custom "fill_vec" function to do this (see s22_functions.R in this folder)
  trials$fA_img <- fill_vec(trials$fA_img, bottom_up = TRUE)
  trials$fB_img <- fill_vec(trials$fB_img, bottom_up = TRUE)
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
    if (!is.na(trials$trials_and_makeups.ran[row]) | !(is.na(trials$blocks.ran[row]))){
      trials$delete_row[row] <- 1 #mark row for deletion
    } else if (!(is.na(trials$remaining_trial.ran[row]))){
      #Mid-trial ratchet-down. First, copy all the valuable data from the row below to the row on which the trial started.
      row_below <- row + 1
      trials$trial_raw[row] <- trials$trial_raw[row_below]
      trials$fB_img[row] <- trials$fB_img[row_below]
      trials$fA_img[row] <- trials$fA_img[row_below]
      trials$show_fres[row] <- trials$show_fres[row_below]
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
  
  #if the data has a "too slow" column in it...
  if(any(names(trials) == "too_slow")){
    #if the subject was "too slow" on a trial, mark the trial as late
    for(row in 1:nrow(trials)){
      if(trials$too_slow[row] == "TRUE"){
        trials$choice[row] <- "late"
      }
    }
    trials <- select(trials,-too_slow) #don't need this anymore
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
  
  #Reduce the fA/B_img columns to just the fractal number - not the whole path to the JPEG
  trials <- mutate(trials, fA_img = str_extract(fA_img,"\\d+"), fB_img = str_extract(fB_img,"\\d+"))
  
  trials$row_index <- c(1:nrow(trials)) #get the row numbers for this subject, which will be useful for a number of things to follow
  
  #create choice numeric column, assigning a value of 1 for choosing fractal A and 2 for choosing fractal B
  trials$choice_numeric <- NA
  for(row in 1:nrow(trials)){
    if(trials$choice[row] == "fA"){
      trials$choice_numeric[row] <- 1
    } else if(trials$choice[row] == "fB"){
      trials$choice_numeric[row] <- 2
    }
  }
  
  #create a "stay" column, indicating whether the participant made the same choice the next time the pair was presented (if so, 1) or not (if so, 0)
  stay_data <- filter(trials,choice != "late") #just look at non-late choices (coding this is tricky if you look at everything)
  stay_data$stay <- NA #initialize stay column
  #go trial-by-trial...
  for(row in 1:nrow(stay_data)){
    pair_rows <- which(stay_data$fA_img == stay_data$fA_img[row]) #identify the rows on which the participant played with the same pair
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
  #join this stay column back to the main trials df
  stay_join <- select(stay_data,row_index,stay) 
  trials <- left_join(trials,stay_join,by="row_index")
  
  trials <- trials %>% mutate(val_rat = val_rat_raw*-1 + 2) %>% select(-val_rat_raw) #recode valence to be on a 0-1 scale
  trials <- trials %>% mutate(valrat_z = (val_rat - mean(val_rat,na.rm=T))/sd(val_rat,na.rm=T))
  
  block_list <- by(trials,trials$block,add_prevrate,rat_col_name="valrat_z") #add prev_rate column to the df for each block
  trials <- do.call(rbind,block_list)
  trials$prev_rate <- as.numeric(trials$prev_rate)
  
  #NOW, THE SUBJECT INFO DF...
  age <- sub_info$age_gender_formresponse[grepl("\\d+",sub_info$age_gender_formresponse)] #identify row with the age by searching for a digit value
  gender <- sub_info$age_gender_formresponse[grepl("[a-zA-Z]+",sub_info$age_gender_formresponse)] #identify row with gender by searching for letter value
  if(length(age)==0){
    age <- 999
  }
  if(length(gender)==0){
    gender <- 999
  }
  earnings_before_ec <- sub_info$earnings_before_ec[grepl("\\d+",sub_info$earnings_before_ec)] #identify row with earnings before ec by searching for a digit value
  earnings <- sub_info$earnings[grepl("\\d+",sub_info$earnings)] #ditto for total earnings
  total_experiment_time <- sub_info$total_experiment_time[grepl("\\d+",sub_info$total_experiment_time)] #same strategy for total experiment time
  
  #count the number of each type of attention check passed
  bv_check_passed <- length(which(sub_info$bv_check == "q"))
  rev_check_passed <- length(which(sub_info$rev_check == "q"))
  val_check_passed <- length(which(sub_info$val_check == "q"))
  
  att_checks_passed <- bv_check_passed + rev_check_passed + val_check_passed #get the total
  
  id <- sub_info$id[1] #get id from one of the rows
  date <- sub_info$date[1] #get date from one of the rows
  
  #check the longest string of late choices the subject had. Long strings suggest that the participant left the task
  value_runs <- rle(trials$choice) #get the runs of each value
  #if late is one of the choice values...
  if(any(value_runs$values == "late")){
    late_indices <- which(value_runs$values == "late") #get the indices that correspond to runs of "late"
    consecutive_late_choices <- max(value_runs$lengths[late_indices]) #get the longest "late" run
  } else{
    consecutive_late_choices <- 0 #if none of the runs are "late"s, then the subject had no late choices
  }
  #get the standard deviations of valence ratings
  sd_valrat <- sd(trials$val_rat,na.rm=TRUE) 
  
  instruct_keypress <- length(which(sub_info$instruct_keypress == "right")) #number of manual processions on instructions slides; 
                                                                            #a small number is a red flag they weren't paying attention
  
  answers_correct <- length(which(sub_info$comp_qs=="poly_true"))
  answers_incorrect <- length(which(sub_info$comp_qs=="poly_false"))
  
  mean_dec_rt <- mean(trials$dec_rt, na.rm=TRUE) #grab mean choice RT 
  mean_valrat_rt <- mean(trials$val_rat_rt, na.rm=TRUE) #get mean RT to decision affect probes
  
  trials_completed <- sum(trials$choice != "late") #get the total number of trials completed
  percent_left <- length(which(trials$choice == "fA"))/trials_completed  #times fA was chosen divided by total choices

  b1 <- filter(trials,choice != "late" & block == 1)
  percent_left_b1 <- length(which(b1$choice == "fA"))/nrow(b1)
  
  b2 <- filter(trials,choice != "late" & block == 2)
  percent_left_b2 <- length(which(b2$choice == "fA"))/nrow(b2)
  
  b3 <- filter(trials,choice != "late" & block == 3)
  percent_left_b3 <- length(which(b3$choice == "fA"))/nrow(b3)
    
  
  #get the percentage of trials on which the subject didn't make a choice in time
  late_numeric <- ifelse(trials$choice == "late",1,0)
  late_percent <- mean(late_numeric)
  
  #fit reward and bvs to valence ratings on trials with both, geting beta coefficients and r squared
  trials_aff_fres <- trials %>% filter(valrat_z != 0 & show_fres == 1)
  valreg_fit <- lm(valrat_z ~ out + box_val, trials_aff_fres)
  fres_b <- coef(valreg_fit)["out"]
  bv_b <- coef(valreg_fit)["box_val"]
  valreg_rsq <- summary(valreg_fit)$r.squared
  fres_bv_ratio <- abs(fres_b)/abs(bv_b)
  
  #get the percentage of trials on which subjects made an affect rating
  trials_nl <- trials %>% filter(choice != "late") #provided affect rating on all non-late trials
  valrat_skipped_percent <- length(which(is.na(trials_nl$val_rat)))/nrow(trials_nl)
  
  valrat_autoprocess <- length(which(!is.na(trials$val_rat) & is.na(trials$val_rat_rt))) #number of trials on which a val rating was made but the button wasn't clicked
  
  
  ### Calculate some of the above variables for each block separately
  
  #late percent
  block1 <- filter(trials,block == 1)
  late_numeric_b1 <- ifelse(block1$choice == "late",1,0)
  late_percent_b1 <- mean(late_numeric_b1)
  
  block2 <- filter(trials,block == 2)
  late_numeric_b2 <- ifelse(block2$choice == "late",1,0)
  late_percent_b2 <- mean(late_numeric_b2)
  
  block3 <- filter(trials,block == 3)
  late_numeric_b3 <- ifelse(block3$choice == "late",1,0)
  late_percent_b3 <- mean(late_numeric_b3)
  
  #skipped percent
  block1_nl <- block1 %>% filter(choice != "late") #provided affect rating on all non-late trials
  valrat_skipped_percent_b1 <- length(which(is.na(block1_nl$val_rat)))/nrow(block1_nl)
  
  block2_nl <- block2 %>% filter(choice != "late") #provided affect rating on all non-late trials
  valrat_skipped_percent_b2 <- length(which(is.na(block2_nl$val_rat)))/nrow(block2_nl)
  
  block3_nl <- block3 %>% filter(choice != "late") #provided affect rating on all non-late trials
  valrat_skipped_percent_b3 <- length(which(is.na(block3_nl$val_rat)))/nrow(block3_nl)

  #valence prediction
  trials_aff_fres1 <- block1 %>% filter(valrat_z != 0 & show_fres == 1)
  valreg_fit1 <- lm(valrat_z ~ out + box_val, trials_aff_fres1)
  fres_b1 <- coef(valreg_fit1)["out"]
  bv_b1 <- coef(valreg_fit1)["box_val"]
  valreg_rsq1 <- summary(valreg_fit1)$r.squared
  
  trials_aff_fres2 <- block2 %>% filter(valrat_z != 0 & show_fres == 1)
  valreg_fit2 <- lm(valrat_z ~ out + box_val, trials_aff_fres2)
  fres_b2 <- coef(valreg_fit2)["out"]
  bv_b2 <- coef(valreg_fit2)["box_val"]
  valreg_rsq2 <- summary(valreg_fit2)$r.squared
  
  trials_aff_fres3 <- block3 %>% filter(valrat_z != 0 & show_fres == 1)
  valreg_fit3 <- lm(valrat_z ~ out + box_val, trials_aff_fres3)
  fres_b3 <- coef(valreg_fit3)["out"]
  bv_b3 <- coef(valreg_fit3)["box_val"]
  valreg_rsq3 <- summary(valreg_fit3)$r.squared
  
  
  #AFTER DOING THIS, CHECK THAT YOU HAVE EVERYTHING YOU SAID YOU'D GET
  sub_info_final <- data.frame(id,date,age,gender,earnings,earnings_before_ec,bv_check_passed,rev_check_passed,val_check_passed,att_checks_passed,
                               consecutive_late_choices,valrat_skipped_percent,valrat_autoprocess,trials_completed,
                               sd_valrat,total_experiment_time,instruct_keypress,
                               mean_dec_rt,mean_valrat_rt,percent_left,late_percent,
                               answers_correct,answers_incorrect,fres_b,bv_b,valreg_rsq,fres_bv_ratio,
                               percent_left_b1,percent_left_b2,percent_left_b3,
                               late_percent_b1,late_percent_b2,late_percent_b3, valrat_skipped_percent_b1,valrat_skipped_percent_b2,valrat_skipped_percent_b3,
                               valreg_rsq1, valreg_rsq2, valreg_rsq3,fres_b1,fres_b2,fres_b3,bv_b1,bv_b2,bv_b3,row.names=NULL)
  
  base_data <- cbind(sub_info_final,trials) #merge the two dfs
  
  return(list(base_data,sub_info_final)) #return the trial-level data and subject-level-data
}