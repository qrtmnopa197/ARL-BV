// This is a rough first-pass at an RL model that includes effects of both reward and fake outcome,
// fit to get a sense for whether effects existed based on pilot data.
// ~DP
data {
  int<lower=1> n_t; // number of trials
  int<lower=1> n_s; // number of subjects
  int<lower=1> n_f; //number of fractals
  array[n_s,n_t] int fA; // 2d array - one array for each subject, containing a sub-arry with the identity of
                          // fractal A on each of that subject's trials (fractals are identified with numbers 1:n_f)
  array[n_s,n_t] int fB; // diddo for fractal B identities
  array[n_s] vector[n_t] reward; // 1d array of vectors - one for each subject - that contains the reward on each trial
  array[n_s] vector[n_t] fake_out; // 1d array of vectors - one for each subject - that contains the valence rating z score on each trial
  array[n_s] vector[n_t] affect; 
  array[n_s,n_t] int choice; // 2d array - one integer array for each subject - containing subject's choices - coded
                             // 1 for fractal A and 2 for fractal B. 
}

parameters {
  real alpha_mu; // mean of group-level learning rate distribution
  real<lower=0> alpha_sigma; // sd of group-level learning rate distribution

  real flearn_mu;
  real<lower=0> flearn_sigma;

  real beta_mu; // mean of group-level inverse temperature distribution
  real<lower=0> beta_sigma; // sd of group-level inverse temperature distribution

  real fsens_int; // mean of group-level inverse temperature distribution
  real fob_b;
  real<lower=0> fsens_err_sigma; // sd of group-level inverse temperature distribution
  vector[n_s] fsens_dev;
  

  real phi_mu; //Mean of group-level autocorrelation weight distribution; in same units as beta - although C is
               //generally held within a tighter range, so phi would need to be a bit bigger (bit less than 2x prob) than
               //beta for autocorrelation to have the same effect as reward.
  real<lower=0> phi_sigma; //sd of group-level distribution

  real tau_mu; //mean of group-level learning rate for C - the autocorrelation value
  real<lower=0> tau_sigma; //sd

  vector[n_s] alpha_z; // z-score of subject-level alpha value (relative to group distribution)
  vector[n_s] flearn_z;
  vector[n_s] beta_z; // z-score of subject-level beta value (relative to group distribution)
  vector[n_s] phi_z; //z-score of subject-level phi
  vector[n_s] tau_z; //z-score of subject-level tau
  
  real b0_mu;
  real<lower=0> b0_sigma;
  vector[n_s] b0_z;
  
  real reward_b_mu;
  real<lower=0> reward_b_sigma;
  vector[n_s] reward_b_z;
  
  real fo_b_mu;
  real<lower=0> fo_b_sigma;
  vector[n_s] fo_b_z;
  
  real<lower=0> rgr_sigma;
}

transformed parameters {
  vector[n_s*n_t] log_lik;
  real rsq_fob_fsens;
  {//anonymous_scope_start
    array[n_s,n_t] vector[n_f] Q; // The Q values. These are a 2d array containing vectors: the first dimension is
                                  // the subject, the second dimension is the trial, and for each trial there's a
                                  // vector with a Q value for each fractal
    array[n_s,n_t] vector[n_f] C; // Choice autocorrelation value, of the same structure as Q
    array[n_s,n_t] vector[n_f] F;

    real choice_a; // 1 if fractal A was chosen, 0 otherwise - used for C update
    real choice_b; // 1 if fractal B was chosen, 0 otherwise - used for C update

    vector[n_s] beta = beta_mu + beta_sigma*beta_z;  // get subject-level beta values using n.c.p.
    vector[n_s] phi = phi_mu + phi_sigma*phi_z; //diddo phi
    vector[n_s] alpha = inv_logit(alpha_mu + alpha_sigma*alpha_z); //diddo alpha. Also, use inv_logit to get alpha between 0 and 1
    vector[n_s] flearn = inv_logit(flearn_mu + flearn_sigma*flearn_z);
    vector[n_s] tau = inv_logit(tau_mu + tau_sigma*tau_z); //likewise for tau
    
    vector[n_s] fsens_pred = fsens_int + fob_b*fo_b_z;
    vector[n_s] fsens = fsens_pred + fsens_err_sigma*fsens_dev;
    
    vector[2] softmax_arg;
    
    rsq_fob_fsens = variance(fsens_pred)/(variance(fsens_pred) + square(fsens_err_sigma));
    
    // loop through each subject and trial...
    for (s in 1:n_s) {
      for (t in 1:n_t){
        if(t == 1){
          Q[s,t] = rep_vector(0,n_f); // for the first trial of each subject, set all Q-values to 0
          C[s,t] = rep_vector(0,n_f); // diddo C values
          F[s,t] = rep_vector(0,n_f);
        }
        
        softmax_arg = [beta[s]*Q[s,t,fA[s,t]] + phi[s]*C[s,t,fA[s,t]] + fsens[s]*F[s,t,fA[s,t]],
                       beta[s]*Q[s,t,fB[s,t]] + phi[s]*C[s,t,fB[s,t]] + fsens[s]*F[s,t,fB[s,t]]]'; //set softmax arg outside of the log_lik calculation for efficiency
                                                                                                  //(you should avoid transposing inside a likelihood funciton, apparently)
        log_lik[(s-1)*n_t+t] = categorical_logit_lpmf(choice[s,t] | softmax_arg); //get the likelihood of the choice on this trial. Doing the log lik

        // unlesss this is the subject's last trial, set the Q,C, and A values for the next trial
        if(t < n_t){
          Q[s,t+1] = (1-alpha[s])*Q[s,t]; // Decay Q values toward their initial values - 0. The only values
                                          // that will not be decayed are the ones chosen on this trial -
                                          // updated below.
          F[s,t+1] = (1-flearn[s])*F[s,t]; // Decay A values toward their initial values - 0

          if(choice[s,t] == 1){
            // if fA was chosen...
            // update the Q and A values for fA
            Q[s,t+1,fA[s,t]] = Q[s,t,fA[s,t]] + alpha[s]*(reward[s,t]-Q[s,t,fA[s,t]]);
            F[s,t+1,fA[s,t]] = F[s,t,fA[s,t]] + flearn[s]*(fake_out[s,t]-F[s,t,fA[s,t]]);
            //update C using 1 for choice a and 0 for choice b
            choice_a = 1;
            choice_b = 0;
          } else if(choice[s,t] == 2){
            // vice-versa if fB was chosen
            Q[s,t+1,fB[s,t]] = Q[s,t,fB[s,t]] + alpha[s]*(reward[s,t]-Q[s,t,fB[s,t]]);
            F[s,t+1,fB[s,t]] = F[s,t,fB[s,t]] + flearn[s]*(fake_out[s,t]-F[s,t,fB[s,t]]);
            // C value update
            choice_a = 0;
            choice_b = 1;
          }
          C[s,t+1] = C[s,t]; // Iniitalize the next trial's C values to be the same as the current trial's, as with Q values
          C[s,t+1,fA[s,t]] = C[s,t,fA[s,t]] + tau[s]*(choice_a - C[s,t,fA[s,t]]); // update C value of fA with 1 if fA was chosen and 0 if it wasn't.
          C[s,t+1,fB[s,t]] = C[s,t,fB[s,t]] + tau[s]*(choice_b - C[s,t,fB[s,t]]); // diddo fB
        }
      }
    }
  }//anonymous_scope_end
}

model {
  vector[n_s] fo_b = fo_b_mu + fo_b_sigma*fo_b_z;
  vector[n_s] reward_b = reward_b_mu + reward_b_sigma*reward_b_z;
  vector[n_s] b0 = b0_mu + b0_sigma*b0_z;
  
   // add the hyperpriors on the group-level parameters to the target density
  alpha_mu ~ normal(-.05,1.7); //This yields a nearly uniform 0-1 prior when logit transformed
  alpha_sigma ~ normal(0,4); //allow for large inter-subject heterogeneity without enforcing marginal subject-level priors that are excessively horeshoe-shaped
  flearn_mu ~ normal(-.05,1.7);
  flearn_sigma ~ normal(0,4);
  tau_mu ~ normal(-.05,1.7); //setting tau priors same as alpha
  tau_sigma ~ normal(0,4); // ""
  beta_mu ~ normal(1,5); // Weakly informative prior
  beta_sigma ~ normal(0,5); //Setting the prior on the SD with similar logic
  fsens_int ~ normal(0,5);
  fob_b ~ normal(0,5);
  fsens_err_sigma ~ normal(0,5);
  phi_mu ~ normal(2,10); // Setting a reguarlizing prior. Matching the phi prior to the beta prior, assuming that uatocorrelation and reward are apt to have
                         // equivalent effects on choice. But differences in phi between fA and fB will be about .5 at asymptote (based on data in which
                         // subjects tend to repeat choices 75% of the time), whereas it's more like 1.2 for beta, so the prior for phi will be set a bit wider
  phi_sigma ~ normal(0,10);
  // add the "priors" on the subject-level parameters derived form the group-level parameters to the target density
  alpha_z ~ normal(0,1);
  flearn_z ~ normal(0,1);
  beta_z ~ normal(0,1);
  fsens_dev ~ normal(0,1);
  phi_z ~ normal(0,1);
  tau_z ~ normal(0,1);

  // Add the joint likelihood to the target density, mapping beta-adjusted Q values to choice probabilities
  // using a softmax function. Unfortunately, categorical_logit does not support vectorization.
  target += log_lik;
  
  rgr_sigma ~ normal(0,1);
  reward_b_mu ~ normal(0,0.5);
  reward_b_sigma ~ normal(0,1);
  fo_b_mu ~ normal(0,0.5);
  fo_b_sigma ~ normal(0,1);
  b0_mu ~ normal(0,2);
  b0_sigma ~ normal(0,4);
  
  reward_b_z ~ normal(0,1);
  fo_b_z ~ normal(0,1);
  b0_z ~ normal(0,1);

  
  for(s in 1:n_s){
    affect[s] ~ normal(b0[s] + reward_b[s]*reward[s] + fo_b[s]*fake_out[s], rgr_sigma);
  }
}
