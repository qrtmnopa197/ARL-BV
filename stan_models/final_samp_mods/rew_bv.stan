// Predicts choice using reward associations and box value associations

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
}

parameters{
  //RL parameters
  
  //effect of reward - on trials where fractal result was received - on choice
  real rew_fr_sens_mu; //population mean
  real<lower=0> rew_fr_sens_sigma; //population sd
  vector[n_s] rew_fr_sens_z; //z-scores for subjects within population
  
  //effect of affect/valence rating - on trials where fractal result was received - on choice
  real aff_fr_sens_mu; 
  real<lower=0> aff_fr_sens_sigma;
  vector[n_s] aff_fr_sens_z;
  
  //effect of reward - on trials where fractal result was not received - on choice
  real rew_nf_sens_mu; //population mean
  real<lower=0> rew_nf_sens_sigma; //population sd
  vector[n_s] rew_nf_sens_z; //z-scores for subjects within population
  
  //strength of choice autocorrelation
  real csens_mu;
  real<lower=0> csens_sigma;
  vector[n_s] csens_z;
  
  //strength of left-side bias - higher values mean a greater bias toward choosing the lefthand fractal, regardless of its value
  real ls_bias_mu;
  real<lower=0> ls_bias_sigma;
  vector[n_s] ls_bias_z;
  
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
}

transformed parameters {
  vector[n_s*n_t] choice_lik; //log likelihood of the choice on each trial

  //SDs of expected values
  real Q_fr_sd;
  real A_fr_sd;
  real C_sd;
  
  {//anonymous_scope_start
    array[n_s,n_t] vector[n_f] Q_fr; //Q-value for trials where fractal result is received
    array[n_s,n_t] vector[n_f] Q_nf; //Q-value for trials where fractal result is not received
    array[n_s,n_t] vector[n_f] A_fr; 
    array[n_s,n_t] vector[n_f] C; // Choice autocorrelation value
    real choice_a; // 1 if fractal A was chosen, 0 otherwise - used for C update
    real choice_b; // 1 if fractal B was chosen, 0 otherwise - used for C update
    
    //differences between Q, A, and C values of fractal A and fractal B
    matrix[n_s,n_t] Q_diff;
    matrix[n_s,n_t] A_diff;
    matrix[n_s,n_t] C_diff;
    
    
    vector[2] softmax_arg; //goes inside choice likelihood calculation
    
    //Calculate subject-level values for hierarchical parameters using non-centered parameterization
    vector[n_s] csens = csens_mu + csens_sigma*csens_z;
    vector[n_s] rew_fr_sens = rew_fr_sens_mu + rew_fr_sens_sigma*rew_fr_sens_z;
    vector[n_s] rew_nf_sens = rew_nf_sens_mu + rew_nf_sens_sigma*rew_nf_sens_z;
    vector[n_s] aff_fr_sens = aff_fr_sens_mu + aff_fr_sens_sigma*aff_fr_sens_z;

    vector[n_s] ls_bias = ls_bias_mu + ls_bias_sigma*ls_bias_z;
    
    vector[n_s] dcy_fr = inv_logit(dcy_fr_mu + dcy_fr_sigma*dcy_fr_z);
    vector[n_s] dcy_nf = inv_logit(dcy_nf_mu + dcy_nf_sigma*dcy_nf_z);

    vector[n_s] lrn_fr = inv_logit(lrn_fr_mu + lrn_fr_sigma*lrn_fr_z);
    vector[n_s] lrn_nf = inv_logit(lrn_nf_mu + lrn_nf_sigma*lrn_nf_z);
    vector[n_s] lrn_c = inv_logit(lrn_c_mu + lrn_c_sigma*lrn_c_z);
    

    // loop through each subject and trial...
    for (s in 1:n_s) {
      for (t in 1:n_t){
        if(t == 1){
          // for the first trial of each subject, set all Q/C values to 0
          Q_fr[s,t] = rep_vector(0,n_f); 
          Q_nf[s,t] = rep_vector(0,n_f); 
          A_fr[s,t] = rep_vector(0,n_f); 
          C[s,t] = rep_vector(0,n_f);
        }
        
        softmax_arg = rew_fr_sens[s]*Q_fr[s,t,{fA[s,t],fB[s,t]}] + rew_nf_sens[s]*Q_nf[s,t,{fA[s,t],fB[s,t]}] +
                      aff_fr_sens[s]*A_fr[s,t,{fA[s,t],fB[s,t]}] +
                      csens[s]*C[s,t,{fA[s,t],fB[s,t]}] + [ls_bias[s],0]'; //set argument for softmax decision function
                      
        choice_lik[(s-1)*n_t+t] = categorical_logit_lpmf(choice[s,t] | softmax_arg); //get the likelihood of the choice on this trial
        
        // Get differences between Q, A, C
        Q_diff[s,t] = Q_fr[s,t,fA[s,t]] - Q_fr[s,t,fB[s,t]];
        A_diff[s,t] = A_fr[s,t,fA[s,t]] - A_fr[s,t,fB[s,t]];
        C_diff[s,t] = C[s,t,fA[s,t]] - C[s,t,fB[s,t]];
      
        // unlesss this is the subject's last trial, set the Q and C values for the next trial
        if(t < n_t){
          //decay Q values toward 0
          Q_fr[s,t+1] = (1-dcy_fr[s])*Q_fr[s,t];
          Q_nf[s,t+1] = (1-dcy_nf[s])*Q_nf[s,t];
          A_fr[s,t+1] = (1-dcy_fr[s])*A_fr[s,t];
          if(choice[s,t] == 1){
            // if fA was chosen...
            //update C using 1 for choice a and 0 for choice b
            choice_a = 1;
            choice_b = 0;
            // if the outcome was the fractal result...
            if(fres[s,t] == 1){
              //update Q_fr...
              Q_fr[s,t+1,fA[s,t]] = Q_fr[s,t,fA[s,t]] + lrn_fr[s]*(rew[s,t] - Q_fr[s,t,fA[s,t]]);
              A_fr[s,t+1,fA[s,t]] = A_fr[s,t,fA[s,t]] + lrn_fr[s]*(bv[s,t] - A_fr[s,t,fA[s,t]]);
            } else if(fres[s,t] == 0){
              //otherwise update Q_nf
              Q_nf[s,t+1,fA[s,t]] = Q_nf[s,t,fA[s,t]] + lrn_nf[s]*(rew[s,t] - Q_nf[s,t,fA[s,t]]);
            }
          } else if(choice[s,t] == 2){
            // vice-versa if fB was chosen...
            choice_a = 0;
            choice_b = 1;
            // if the outcome was the fractal result...
            if(fres[s,t] == 1){
              //update Q_fr...
              Q_fr[s,t+1,fB[s,t]] = Q_fr[s,t,fB[s,t]] + lrn_fr[s]*(rew[s,t] - Q_fr[s,t,fB[s,t]]);
              A_fr[s,t+1,fB[s,t]] = A_fr[s,t,fB[s,t]] + lrn_fr[s]*(bv[s,t] - A_fr[s,t,fB[s,t]]);
            } else if(fres[s,t] == 0){
              //otherwise update Q_nf
              Q_nf[s,t+1,fB[s,t]] = Q_nf[s,t,fB[s,t]] + lrn_nf[s]*(rew[s,t] - Q_nf[s,t,fB[s,t]]);
            }
          }
          //update C
          C[s,t+1] = C[s,t]; // Initialize the next trial's C values to be the same as the current trial's (no forgetting)
          C[s,t+1,fA[s,t]] = C[s,t,fA[s,t]] + lrn_c[s]*(choice_a - C[s,t,fA[s,t]]); // update C value of fA with 1 if fA was chosen and 0 if it wasn't.
          C[s,t+1,fB[s,t]] = C[s,t,fB[s,t]] + lrn_c[s]*(choice_b - C[s,t,fB[s,t]]); // ditto fB
        }
      }
    }
    Q_fr_sd = sd(Q_diff);
    A_fr_sd = sd(A_diff);
    C_sd = sd(C_diff);
  }//anonymous_scope_end
}
model{
  //effects on choice
  rew_fr_sens_mu ~ normal(0,125); 
  rew_fr_sens_sigma ~ normal(0,175); 
  rew_nf_sens_mu ~ normal(0,4); 
  rew_nf_sens_sigma ~ normal(0,6); 
  
  aff_fr_sens_mu ~ normal(0,4); 
  aff_fr_sens_sigma ~ normal(0,6); 

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
  
  //participant z-scores
  rew_fr_sens_z ~ std_normal();
  rew_nf_sens_z ~ std_normal();
  aff_fr_sens_z ~ std_normal();
  csens_z ~ std_normal();
  dcy_fr_z ~ std_normal();
  dcy_nf_z ~ std_normal();
  lrn_fr_z ~ std_normal();
  lrn_nf_z ~ std_normal();
  lrn_c_z ~ std_normal();
  ls_bias_z ~ std_normal();


  target += choice_lik; //increment target likelihood with the likelihoods of choices under the model
}
generated quantities{
  //scaled effects of expected values
  real scd_rew_fr_sens_mu = Q_fr_sd*rew_fr_sens_mu;
  real scd_aff_fr_sens_mu = A_fr_sd*aff_fr_sens_mu;
  real scd_csens_mu = C_sd*csens_mu;
}
