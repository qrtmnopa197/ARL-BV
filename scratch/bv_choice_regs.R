trials_fres <- trials %>% filter(show_fres == 1)
trials_recode <- trials_fres %>% mutate(out_rc = ifelse(choice=="fA",out,-out),
                                        box_val_rc = ifelse(choice=="fA",box_val,-box_val))
trials_nf_lag <- add_lag_cols(trials_recode,c("out_rc","box_val_rc"),lags=c(1:4)) 
trials_bin <- trials_nf_lag %>% mutate(fA_chosen = ifelse(choice=="fA",1,0))
lag_reg <- glm(fA_chosen ~ 
      box_val_rc_lag1 + box_val_rc_lag2 + box_val_rc_lag3 + box_val_rc_lag4, data=trials_bin,family=binomial)
