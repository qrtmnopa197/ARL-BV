// RL model that assumes only affect associations influence choices

data {
  int<lower=1> n_t; // number of trials
  int<lower=1> n_s; // number of subjects
  int<lower=1> n_f; //number of fractals
  
  array[n_s,n_t] int fA; // identities of fractal A (the lefthand fractal) on each of that subject's trials (fractals are identified with numbers 1:n_f)
  array[n_s,n_t] int fB; // identities of fractal B (the righthand fractal) on each of that subject's trials
  
  array[n_s] vector[n_t] rew; //the reward recieved - i.e., the monetary value of the outcome received
  array[n_s] vector[n_t] bv; //the box value
  array[n_s,n_t] int fres; //1/0 variable indicating whether the fractal result was received/not received
  array[n_s,n_t] int choice; // The choice made on each trial: 1 for fractal A (the lefthand fractal), 2 for fractal B (the righthand fractal)
  
  int n_rat; //number of valence ratings made
  array[n_s,n_t] int rat_num; //rating number for each trial
  vector[n_rat] rat; //the valence rating on each trial
  array[n_s] vector[n_t] prev_rat; //the previous valence rating (used for autoregressive term)
  
  array[n_s,n_t] int chosen_frac; //indices of the chosen fractal for each trial
  array[n_s,n_t] int unchosen_frac; //ditto unchosen fractals     
}

parameters {
  //RL parameters
  
  //effect of model-predicted valence - on trials where fractal result was received - on choice
  real aff_fr_sens_mu; 
  real<lower=0> aff_fr_sens_sigma;
  vector[n_s] aff_fr_sens_z;
  
  //effect of valence residual on choice - on trials where fractal result was received
  real resid_fr_sens_mu; 
  real<lower=0> resid_fr_sens_sigma;
  vector[n_s] resid_fr_sens_z;
  
  //effect of model-predicted valence - on trials where fractal result was not received - on choice
  real aff_nf_sens_mu; 
  real<lower=0> aff_nf_sens_sigma;
  vector[n_s] aff_nf_sens_z;
  
  //effect of residual on choice - on trials where fractal result was not received
  real resid_nf_sens_mu; 
  real<lower=0> resid_nf_sens_sigma;
  vector[n_s] resid_nf_sens_z;
  
  //strength of left-side bias - higher values mean a greater bias toward choosing the lefthand fractal, regardless of its value
  real ls_bias_mu;
  real<lower=0> ls_bias_sigma;
  vector[n_s] ls_bias_z;
  
  //strength of choice autocorrelation
  real csens_mu;
  real<lower=0> csens_sigma;
  vector[n_s] csens_z;
  
  //Decay rates; fractal values assumed to decay toward 0 when not updated
  
  //for trials where fractal result was shown
  real dcy_fr_mu;
  real<lower=0> dcy_fr_sigma;
  vector[n_s] dcy_fr_z;
  
  //for trials where fractal result was NOT shown
  real dcy_nf_mu;
  real<lower=0> dcy_nf_sigma;
  vector[n_s] dcy_nf_z;
  
  //Learning rates for Q-values and autocorrelation value C
  
  //for trials where fractal result was shown
  real lrn_fr_mu;
  real<lower=0> lrn_fr_sigma;
  vector[n_s] lrn_fr_z;
  
  //for trials where fractal result was NOT shown
  real lrn_nf_mu;
  real<lower=0> lrn_nf_sigma;
  vector[n_s] lrn_nf_z;
  
  //for autocorrelation
  real lrn_c_mu;
  real<lower=0> lrn_c_sigma;
  vector[n_s] lrn_c_z;
  
  
  //Valence regression parameters
  //intercept
  real B_0_mu;
  real<lower=0> B_0_sigma;
  vector[n_s] B_0_z;
  
  //effect of reward - on trials where fractal result was received
  real B_rew_fr_mu; 
  real<lower=0> B_rew_fr_sigma;
  vector[n_s] B_rew_fr_z;
  
  //effect of reward - on trials where fractal result was not received
  real B_rew_nf_mu; 
  real<lower=0> B_rew_nf_sigma;
  vector[n_s] B_rew_nf_z;
  
  //effect of box values - on trials where fractal result was received
  real B_bv_fr_mu; 
  real<lower=0> B_bv_fr_sigma;
  vector[n_s] B_bv_fr_z;
  
  //effect of Q value
  real B_q_fr_mu;
  real<lower=0> B_q_fr_sigma;
  vector[n_s] B_q_fr_z;
  
  real B_q_nf_mu; 
  real<lower=0> B_q_nf_sigma;
  vector[n_s] B_q_nf_z;
  
  //effect of PWQD - on trials where fractal result was received
  real B_pwqd_fr_mu; 
  real<lower=0> B_pwqd_fr_sigma;
  vector[n_s] B_pwqd_fr_z;
  
  real B_pwqd_nf_mu; 
  real<lower=0> B_pwqd_nf_sigma;
  vector[n_s] B_pwqd_nf_z;
  
  //effect of previous rating (autoregressive term)
  real B_auto_mu; 
  real<lower=0> B_auto_sigma;
  vector[n_s] B_auto_z;
  
  //residual standard deviation
  real<lower=0> resid_sigma;
}

transformed parameters {
  vector[n_s*n_t] choice_lik; //log likelihood of the choice on each trial
  vector[n_rat] rat_pred; //the predicted valence rating on each trial
  
  //SDs of expected values
  real A_fr_sd;
  real A_nf_sd;
  real R_fr_sd;
  real R_nf_sd;
  real C_sd;
  
  {//anonymous_scope_start
  
    array[n_s,n_t] vector[n_f] Q_fr; //Q-value for trials where fractal result is received
    array[n_s,n_t] vector[n_f] A_fr; //ditto A-value...
    array[n_s,n_t] vector[n_f] A_nf; 
    array[n_s,n_t] vector[n_f] R_fr; //ditto R...
    array[n_s,n_t] vector[n_f] R_nf; 
    array[n_s,n_t] vector[n_f] C; // Choice autocorrelation value
    
    //differences between Q, A, and C values of fractal A and fractal B
    matrix[n_s,n_t] A_fr_diff;
    matrix[n_s,n_t] A_nf_diff;
    matrix[n_s,n_t] R_fr_diff;
    matrix[n_s,n_t] R_nf_diff;
    matrix[n_s,n_t] C_diff;
    
    real choice_a; // 1 if fractal A was chosen, 0 otherwise - used for C update
    real choice_b; // 1 if fractal B was chosen, 0 otherwise - used for C update
    
    //Calculate subject-level values for hierarchical parameters using non-centered parameterization
    vector[n_s] csens = csens_mu + csens_sigma*csens_z;
    vector[n_s] aff_fr_sens = aff_fr_sens_mu + aff_fr_sens_sigma*aff_fr_sens_z;
    vector[n_s] aff_nf_sens = aff_nf_sens_mu + aff_nf_sens_sigma*aff_nf_sens_z;
    vector[n_s] resid_fr_sens = resid_fr_sens_mu + resid_fr_sens_sigma*resid_fr_sens_z;
    vector[n_s] resid_nf_sens = resid_nf_sens_mu + resid_nf_sens_sigma*resid_nf_sens_z;
    
    vector[n_s] ls_bias = ls_bias_mu + ls_bias_sigma*ls_bias_z;
    
    vector[n_s] dcy_fr = inv_logit(dcy_fr_mu + dcy_fr_sigma*dcy_fr_z);
    vector[n_s] dcy_nf = inv_logit(dcy_nf_mu + dcy_nf_sigma*dcy_nf_z);

    vector[n_s] lrn_fr = inv_logit(lrn_fr_mu + lrn_fr_sigma*lrn_fr_z);
    vector[n_s] lrn_nf = inv_logit(lrn_nf_mu + lrn_nf_sigma*lrn_nf_z);
    vector[n_s] lrn_c = inv_logit(lrn_c_mu + lrn_c_sigma*lrn_c_z);
    
    vector[n_s] B_0 = B_0_mu + B_0_sigma*B_0_z;
    vector[n_s] B_rew_fr = B_rew_fr_mu + B_rew_fr_sigma*B_rew_fr_z;
    vector[n_s] B_rew_nf = B_rew_nf_mu + B_rew_nf_sigma*B_rew_nf_z;
    vector[n_s] B_bv_fr = B_bv_fr_mu + B_bv_fr_sigma*B_bv_fr_z;
    vector[n_s] B_q_fr = B_q_fr_mu + B_q_fr_sigma*B_q_fr_z;
    vector[n_s] B_q_nf = B_q_nf_mu + B_q_nf_sigma*B_q_nf_z;
    vector[n_s] B_pwqd_fr = B_pwqd_fr_mu + B_pwqd_fr_sigma*B_pwqd_fr_z;
    vector[n_s] B_pwqd_nf = B_pwqd_nf_mu + B_pwqd_nf_sigma*B_pwqd_nf_z;
    vector[n_s] B_auto = B_auto_mu + B_auto_sigma*B_auto_z;
    
    real curr_pred; //the predicted valence for the current trial
    real nuis; //the nuisance variation in valence for the current trial
    real resid; //affect residual on the current trial
    
    vector[2] softmax_arg;

    // loop through each subject and trial...
    for (s in 1:n_s) {
      for (t in 1:n_t){
        if(t == 1){
          // for the first trial of each subject, set all expected values to 0
          Q_fr[s,t] = rep_vector(0,n_f); 
          A_fr[s,t] = rep_vector(0,n_f); 
          A_nf[s,t] = rep_vector(0,n_f); 
          R_fr[s,t] = rep_vector(0,n_f); 
          R_nf[s,t] = rep_vector(0,n_f); 
          C[s,t] = rep_vector(0,n_f);
        }
        
        softmax_arg = aff_fr_sens[s]*A_fr[s,t,{fA[s,t],fB[s,t]}] + aff_nf_sens[s]*A_nf[s,t,{fA[s,t],fB[s,t]}] +
                      resid_fr_sens[s]*R_fr[s,t,{fA[s,t],fB[s,t]}] + resid_nf_sens[s]*R_nf[s,t,{fA[s,t],fB[s,t]}] +
                      csens[s]*C[s,t,{fA[s,t],fB[s,t]}] + [ls_bias[s],0]'; //set argument for softmax decision function
                      
        choice_lik[(s-1)*n_t+t] = categorical_logit_lpmf(choice[s,t] | softmax_arg); //get the likelihood of the choice on this trial
        
        //get differences on this trial
        A_fr_diff[s,t] = A_fr[s,t,fA[s,t]] - A_fr[s,t,fB[s,t]];
        A_nf_diff[s,t] = A_nf[s,t,fA[s,t]] - A_nf[s,t,fB[s,t]];
        R_fr_diff[s,t] = R_fr[s,t,fA[s,t]] - R_fr[s,t,fB[s,t]];
        R_nf_diff[s,t] = R_nf[s,t,fA[s,t]] - R_nf[s,t,fB[s,t]];
        C_diff[s,t] = C[s,t,fA[s,t]] - C[s,t,fB[s,t]];
        
        //Generate valence rating prediction
        if(fres[s,t] == 1){
          //if the fractal result was received...
          curr_pred = B_0[s] + B_rew_fr[s]*rew[s,t] + B_bv_fr[s]*bv[s,t] + B_q_fr[s]*Q_fr[s,t,chosen_frac[s,t]] +
                      B_pwqd_fr[s]*(Q_fr[s,t,chosen_frac[s,t]]*exp(choice_lik[(s-1)*n_t+t]) + Q_fr[s,t,unchosen_frac[s,t]]*(1-exp(choice_lik[(s-1)*n_t+t])));
        } else if (fres[s,t] == 0){
          //if not...
          curr_pred = B_0[s] + B_rew_nf[s]*rew[s,t] + B_q_nf[s]*Q_fr[s,t,chosen_frac[s,t]] + 
                      B_pwqd_nf[s]*(Q_fr[s,t,chosen_frac[s,t]]*exp(choice_lik[(s-1)*n_t+t]) + Q_fr[s,t,unchosen_frac[s,t]]*(1-exp(choice_lik[(s-1)*n_t+t])));
        }
        nuis = B_auto[s]*prev_rat[s,t];
        
        if(rat_num[s,t] != 0){
          //if the participant made a valence rating... 
          rat_pred[rat_num[s,t]] = curr_pred + nuis; //add to the vector of predicted ratings
          resid = rat[rat_num[s,t]] - (curr_pred + nuis);
        } else{
          resid = 0;
        }
        
        
        // unlesss this is the subject's last trial, set the Q and C values for the next trial
        if(t < n_t){
          //decay Q values toward 0
          Q_fr[s,t+1] = (1-dcy_fr[s])*Q_fr[s,t];
          A_fr[s,t+1] = (1-dcy_fr[s])*A_fr[s,t];
          A_nf[s,t+1] = (1-dcy_nf[s])*A_nf[s,t];
          R_fr[s,t+1] = (1-dcy_fr[s])*R_fr[s,t];
          R_nf[s,t+1] = (1-dcy_nf[s])*R_nf[s,t];
          if(choice[s,t] == 1){
            // if fA was chosen...
            //update C using 1 for choice a and 0 for choice b
            choice_a = 1;
            choice_b = 0;
            // if the outcome was the fractal result...
            if(fres[s,t] == 1){
              //update fr values...
              Q_fr[s,t+1,fA[s,t]] = Q_fr[s,t,fA[s,t]] + lrn_fr[s]*(rew[s,t] - Q_fr[s,t,fA[s,t]]);
              A_fr[s,t+1,fA[s,t]] = A_fr[s,t,fA[s,t]] + lrn_fr[s]*(curr_pred - A_fr[s,t,fA[s,t]]);
              R_fr[s,t+1,fA[s,t]] = R_fr[s,t,fA[s,t]] + lrn_fr[s]*(resid - R_fr[s,t,fA[s,t]]);
            } else if(fres[s,t] == 0){
              //otherwise update nf values...
              A_nf[s,t+1,fA[s,t]] = A_nf[s,t,fA[s,t]] + lrn_nf[s]*(curr_pred - A_nf[s,t,fA[s,t]]);
              R_nf[s,t+1,fA[s,t]] = R_nf[s,t,fA[s,t]] + lrn_nf[s]*(resid - R_nf[s,t,fA[s,t]]);
            }
          } else if(choice[s,t] == 2){
            // vice-versa if fB was chosen...
            choice_a = 0;
            choice_b = 1;
            // if the outcome was the fractal result...
            if(fres[s,t] == 1){
              //update Q_fr...
              Q_fr[s,t+1,fB[s,t]] = Q_fr[s,t,fB[s,t]] + lrn_fr[s]*(rew[s,t] - Q_fr[s,t,fB[s,t]]);
              A_fr[s,t+1,fB[s,t]] = A_fr[s,t,fB[s,t]] + lrn_fr[s]*(curr_pred - A_fr[s,t,fB[s,t]]);
              R_fr[s,t+1,fB[s,t]] = R_fr[s,t,fB[s,t]] + lrn_fr[s]*(resid - R_fr[s,t,fB[s,t]]);
            } else if(fres[s,t] == 0){
              //otherwise update Q_nf
              A_nf[s,t+1,fB[s,t]] = A_nf[s,t,fB[s,t]] + lrn_nf[s]*(curr_pred - A_nf[s,t,fB[s,t]]);
              R_nf[s,t+1,fB[s,t]] = R_nf[s,t,fB[s,t]] + lrn_nf[s]*(resid - R_nf[s,t,fB[s,t]]);
            }
          }
          //update C
          C[s,t+1] = C[s,t]; // Initialize the next trial's C values to be the same as the current trial's (no forgetting)
          C[s,t+1,fA[s,t]] = C[s,t,fA[s,t]] + lrn_c[s]*(choice_a - C[s,t,fA[s,t]]); // update C value of fA with 1 if fA was chosen and 0 if it wasn't.
          C[s,t+1,fB[s,t]] = C[s,t,fB[s,t]] + lrn_c[s]*(choice_b - C[s,t,fB[s,t]]); // ditto fB
        }
      }
    }
    A_fr_sd = sd(A_fr_diff);
    A_nf_sd = sd(A_nf_diff);
    R_fr_sd = sd(R_fr_diff);
    R_nf_sd = sd(R_nf_diff);
    C_sd = sd(C_diff);
  }//anonymous_scope_end
}
model{
  //effects on choice
  aff_fr_sens_mu ~ normal(0,4); 
  aff_fr_sens_sigma ~ normal(0,6); 
  aff_nf_sens_mu ~ normal(0,6); 
  aff_nf_sens_sigma ~ normal(0,8); 
  
  resid_fr_sens_mu ~ normal(0,4); 
  resid_fr_sens_sigma ~ normal(0,6); 
  resid_nf_sens_mu ~ normal(0,6); 
  resid_nf_sens_sigma ~ normal(0,8); 

  csens_mu ~ normal(0,5);
  csens_sigma ~ normal(0,10);
  
  ls_bias_mu ~ normal(0,4);
  ls_bias_sigma ~ normal(0,8);
  
  //learning and decay rates
  dcy_fr_mu ~ normal(-.05,1.7);
  dcy_fr_sigma ~ normal(0,4);
  dcy_nf_mu ~ normal(-.05,1.7);
  dcy_nf_sigma ~ normal(0,4);
  lrn_fr_mu ~ normal(-.05,1.7);
  lrn_fr_sigma ~ normal(0,4);
  lrn_nf_mu ~ normal(-.05,1.7);
  lrn_nf_sigma ~ normal(0,4);
  lrn_c_mu ~ normal(-.05,1.7);
  lrn_c_sigma ~ normal(0,4);
  
  //valence regression parameters
  B_0_mu ~ normal(0,2);
  B_0_sigma ~ normal(0,4); 
  B_rew_fr_mu ~ normal(0,40);
  B_rew_fr_sigma ~ normal(0,60); 
  B_rew_nf_mu ~ normal(0,0.7);
  B_rew_nf_sigma ~ normal(0,1); 
  B_q_fr_mu ~ normal(0,40);
  B_q_fr_sigma ~ normal(0,60); 
  B_q_nf_mu ~ normal(0,40);
  B_q_nf_sigma ~ normal(0,60); 
  B_pwqd_fr_mu ~ normal(0,40);
  B_pwqd_fr_sigma ~ normal(0,60); 
  B_pwqd_nf_mu ~ normal(0,40);
  B_pwqd_nf_sigma ~ normal(0,60); 
  B_bv_fr_mu ~ normal(0,0.7);
  B_bv_fr_sigma ~ normal(0,1); 
  B_auto_mu ~ normal(0,1.5);
  B_auto_sigma ~ normal(0,2);
  
  //participant z-scores
  aff_fr_sens_z ~ std_normal();
  aff_nf_sens_z ~ std_normal();
  resid_fr_sens_z ~ std_normal();
  resid_nf_sens_z ~ std_normal();
  csens_z ~ std_normal();
  dcy_fr_z ~ std_normal();
  dcy_nf_z ~ std_normal();
  lrn_fr_z ~ std_normal();
  lrn_nf_z ~ std_normal();
  lrn_c_z ~ std_normal();
  B_0_z ~ std_normal();
  B_rew_fr_z ~ std_normal();
  B_rew_nf_z ~ std_normal();
  B_q_fr_z ~ std_normal();
  B_q_nf_z ~ std_normal();
  B_pwqd_fr_z ~ std_normal();
  B_pwqd_nf_z ~ std_normal();
  B_bv_fr_z ~ std_normal();
  B_auto_z ~ std_normal();
  ls_bias_z ~ std_normal();

  resid_sigma ~ normal(0,2);
  
  target += choice_lik; //increment target likelihood with the likelihoods of choices under the model
  rat ~ normal(rat_pred, resid_sigma); //increment target likelihood with the likelihood of valence ratings under the model
}
generated quantities{
  vector[n_rat] affect_lik; //log likelihoods of all affect rating
  real val_mod_rsq = variance(rat_pred)/(variance(rat_pred) + square(resid_sigma)); //valence model R^2

  //scaled effects of expected values
  real scd_aff_fr_sens_mu = A_fr_sd*aff_fr_sens_mu;
  real scd_aff_nf_sens_mu = A_nf_sd*aff_nf_sens_mu;
  real scd_resid_fr_sens_mu = R_fr_sd*resid_fr_sens_mu;
  real scd_resid_nf_sens_mu = R_nf_sd*resid_nf_sens_mu;
  real scd_csens_mu = C_sd*csens_mu;
  
  
  // Simulate data
  array[n_s*n_t] int sim_choice;
  vector[n_rat] sim_rat;
  vector[n_s*n_t] sim_prat;
  {
    array[n_s,n_t] vector[n_f] Q_fr; //Q-value for trials where fractal result is received
    array[n_s,n_t] vector[n_f] A_fr; //ditto A-value...
    array[n_s,n_t] vector[n_f] A_nf; 
    array[n_s,n_t] vector[n_f] R_fr; //ditto R...
    array[n_s,n_t] vector[n_f] R_nf; 
    array[n_s,n_t] vector[n_f] C; // Choice autocorrelation value
    
    real choice_a; // 1 if fractal A was chosen, 0 otherwise - used for C update
    real choice_b; // 1 if fractal B was chosen, 0 otherwise - used for C update
    
    //Calculate subject-level values for hierarchical parameters using non-centered parameterization
    vector[n_s] csens = csens_mu + csens_sigma*csens_z;
    vector[n_s] aff_fr_sens = aff_fr_sens_mu + aff_fr_sens_sigma*aff_fr_sens_z;
    vector[n_s] aff_nf_sens = aff_nf_sens_mu + aff_nf_sens_sigma*aff_nf_sens_z;
    vector[n_s] resid_fr_sens = resid_fr_sens_mu + resid_fr_sens_sigma*resid_fr_sens_z;
    vector[n_s] resid_nf_sens = resid_nf_sens_mu + resid_nf_sens_sigma*resid_nf_sens_z;
    
    vector[n_s] ls_bias = ls_bias_mu + ls_bias_sigma*ls_bias_z;
    
    vector[n_s] dcy_fr = inv_logit(dcy_fr_mu + dcy_fr_sigma*dcy_fr_z);
    vector[n_s] dcy_nf = inv_logit(dcy_nf_mu + dcy_nf_sigma*dcy_nf_z);

    vector[n_s] lrn_fr = inv_logit(lrn_fr_mu + lrn_fr_sigma*lrn_fr_z);
    vector[n_s] lrn_nf = inv_logit(lrn_nf_mu + lrn_nf_sigma*lrn_nf_z);
    vector[n_s] lrn_c = inv_logit(lrn_c_mu + lrn_c_sigma*lrn_c_z);
    
    vector[n_s] B_0 = B_0_mu + B_0_sigma*B_0_z;
    vector[n_s] B_rew_fr = B_rew_fr_mu + B_rew_fr_sigma*B_rew_fr_z;
    vector[n_s] B_rew_nf = B_rew_nf_mu + B_rew_nf_sigma*B_rew_nf_z;
    vector[n_s] B_bv_fr = B_bv_fr_mu + B_bv_fr_sigma*B_bv_fr_z;
    vector[n_s] B_q_fr = B_q_fr_mu + B_q_fr_sigma*B_q_fr_z;
    vector[n_s] B_q_nf = B_q_nf_mu + B_q_nf_sigma*B_q_nf_z;
    vector[n_s] B_pwqd_fr = B_pwqd_fr_mu + B_pwqd_fr_sigma*B_pwqd_fr_z;
    vector[n_s] B_pwqd_nf = B_pwqd_nf_mu + B_pwqd_nf_sigma*B_pwqd_nf_z;
    vector[n_s] B_auto = B_auto_mu + B_auto_sigma*B_auto_z;
    
    real curr_pred; //the predicted valence for the current trial
    real nuis; //the nuisance variation in valence for the current trial
    real resid; //affect residual on the current trial
    
    // Chosen and unchosen fractals on the current trial
    int ch_frac;
    int unch_frac;
    
    real ch_lik;
    
    real curr_prat; // Current previous rating
    
    vector[2] softmax_arg;
    
    for (s in 1:n_s) {
      for (t in 1:n_t){
        if(t == 1){
          // for the first trial of each subject, set all expected values to 0
          Q_fr[s,t] = rep_vector(0,n_f); 
          A_fr[s,t] = rep_vector(0,n_f); 
          A_nf[s,t] = rep_vector(0,n_f); 
          R_fr[s,t] = rep_vector(0,n_f); 
          R_nf[s,t] = rep_vector(0,n_f); 
          C[s,t] = rep_vector(0,n_f);
        }
        // If it's the first trial of the block, set previous rating to 0
        if(t % 48 == 1){
          curr_prat = 0;
        }
        
        softmax_arg = aff_fr_sens[s]*A_fr[s,t,{fA[s,t],fB[s,t]}] + aff_nf_sens[s]*A_nf[s,t,{fA[s,t],fB[s,t]}] +
                      resid_fr_sens[s]*R_fr[s,t,{fA[s,t],fB[s,t]}] + resid_nf_sens[s]*R_nf[s,t,{fA[s,t],fB[s,t]}] +
                      csens[s]*C[s,t,{fA[s,t],fB[s,t]}] + [ls_bias[s],0]'; //set argument for softmax decision function
                      
        sim_choice[(s-1)*n_t+t] = categorical_logit_rng(softmax_arg); 
        ch_lik = categorical_logit_lpmf(sim_choice[(s-1)*n_t+t] | softmax_arg);
        
        // Set the chosen and unchosen fractal for this trial
        if(sim_choice[(s-1)*n_t+t] == 1){
          ch_frac = fA[s,t];
          unch_frac = fB[s,t];
        } else{
          ch_frac = fB[s,t];
          unch_frac = fA[s,t];
        }
        
        
        sim_prat[(s-1)*n_t+t] = curr_prat; // set the previous rating for the current trial
        
        //Generate valence rating prediction
        if(fres[s,t] == 1){
          //if the fractal result was received...
          curr_pred = B_0[s] + B_rew_fr[s]*rew[s,t] + B_bv_fr[s]*bv[s,t] + B_q_fr[s]*Q_fr[s,t,ch_frac] +
                      B_pwqd_fr[s]*(Q_fr[s,t,ch_frac]*exp(ch_lik) + Q_fr[s,t,unch_frac]*(1-exp(ch_lik)));
        } else if (fres[s,t] == 0){
          //if not...
          curr_pred = B_0[s] + B_rew_nf[s]*rew[s,t] + B_q_nf[s]*Q_fr[s,t,ch_frac] + 
                      B_pwqd_nf[s]*(Q_fr[s,t,ch_frac]*exp(ch_lik) + Q_fr[s,t,unch_frac]*(1-exp(ch_lik)));
        }
        nuis = B_auto[s]*curr_prat;
        
        
        
        if(rat_num[s,t] != 0){
          //if the participant made a valence rating... 
          sim_rat[rat_num[s,t]] = normal_rng(curr_pred + nuis,resid_sigma); 
          resid = sim_rat[rat_num[s,t]] - (curr_pred + nuis);
          curr_prat = sim_rat[rat_num[s,t]];
        } else{
          resid = 0;
        }
        
        
        // unlesss this is the subject's last trial, set the Q and C values for the next trial
        if(t < n_t){
          //decay Q values toward 0
          Q_fr[s,t+1] = (1-dcy_fr[s])*Q_fr[s,t];
          A_fr[s,t+1] = (1-dcy_fr[s])*A_fr[s,t];
          A_nf[s,t+1] = (1-dcy_nf[s])*A_nf[s,t];
          R_fr[s,t+1] = (1-dcy_fr[s])*R_fr[s,t];
          R_nf[s,t+1] = (1-dcy_nf[s])*R_nf[s,t];
          if(sim_choice[(s-1)*n_t+t] == 1){
            // if fA was chosen...
            //update C using 1 for choice a and 0 for choice b
            choice_a = 1;
            choice_b = 0;
            // if the outcome was the fractal result...
            if(fres[s,t] == 1){
              //update fr values...
              Q_fr[s,t+1,fA[s,t]] = Q_fr[s,t,fA[s,t]] + lrn_fr[s]*(rew[s,t] - Q_fr[s,t,fA[s,t]]);
              A_fr[s,t+1,fA[s,t]] = A_fr[s,t,fA[s,t]] + lrn_fr[s]*(curr_pred - A_fr[s,t,fA[s,t]]);
              R_fr[s,t+1,fA[s,t]] = R_fr[s,t,fA[s,t]] + lrn_fr[s]*(resid - R_fr[s,t,fA[s,t]]);
            } else if(fres[s,t] == 0){
              //otherwise update nf values...
              A_nf[s,t+1,fA[s,t]] = A_nf[s,t,fA[s,t]] + lrn_nf[s]*(curr_pred - A_nf[s,t,fA[s,t]]);
              R_nf[s,t+1,fA[s,t]] = R_nf[s,t,fA[s,t]] + lrn_nf[s]*(resid - R_nf[s,t,fA[s,t]]);
            }
          } else if(sim_choice[(s-1)*n_t+t] == 2){
            // vice-versa if fB was chosen...
            choice_a = 0;
            choice_b = 1;
            // if the outcome was the fractal result...
            if(fres[s,t] == 1){
              //update Q_fr...
              Q_fr[s,t+1,fB[s,t]] = Q_fr[s,t,fB[s,t]] + lrn_fr[s]*(rew[s,t] - Q_fr[s,t,fB[s,t]]);
              A_fr[s,t+1,fB[s,t]] = A_fr[s,t,fB[s,t]] + lrn_fr[s]*(curr_pred - A_fr[s,t,fB[s,t]]);
              R_fr[s,t+1,fB[s,t]] = R_fr[s,t,fB[s,t]] + lrn_fr[s]*(resid - R_fr[s,t,fB[s,t]]);
            } else if(fres[s,t] == 0){
              //otherwise update Q_nf
              A_nf[s,t+1,fB[s,t]] = A_nf[s,t,fB[s,t]] + lrn_nf[s]*(curr_pred - A_nf[s,t,fB[s,t]]);
              R_nf[s,t+1,fB[s,t]] = R_nf[s,t,fB[s,t]] + lrn_nf[s]*(resid - R_nf[s,t,fB[s,t]]);
            }
          }
          //update C
          C[s,t+1] = C[s,t]; // Initialize the next trial's C values to be the same as the current trial's (no forgetting)
          C[s,t+1,fA[s,t]] = C[s,t,fA[s,t]] + lrn_c[s]*(choice_a - C[s,t,fA[s,t]]); // update C value of fA with 1 if fA was chosen and 0 if it wasn't.
          C[s,t+1,fB[s,t]] = C[s,t,fB[s,t]] + lrn_c[s]*(choice_b - C[s,t,fB[s,t]]); // ditto fB
        }
      }
    }
  }
  
  // Get rating likelihoods
  for(i in 1:n_rat){
    affect_lik[i] = normal_lpdf(rat[i] | rat_pred[i], resid_sigma);
  }
}



