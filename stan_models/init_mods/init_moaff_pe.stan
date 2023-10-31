// This is a rough first-pass at an RL model that includes effects of both reward and affect,
// fit to get a sense for whether effects existed based on pilot data.
// ~DP
data {
  int<lower=1> n_t; // number of trials
  int<lower=1> n_s; // number of subjects
  int<lower=1> n_f; //number of fractals
  array[n_s,n_t] int fA; // 2d array - one array for each subject, containing a sub-arry with the identity of
                          // fractal A on each of that subject's trials (fractals are identified with numbers 1:n_f)
  array[n_s,n_t] int fB; // diddo for fractal B identities
  array[n_s] vector[n_t] out; // 1d array of vectors - one for each subject - that contains the outcome of each trial
  array[n_s] vector[n_t] affect;
  array[n_s,n_t] int fres; //2d array (subjects by trials) indicating whether the fractal result was shown on that trial (1-0)
  array[n_s,n_t] int choice; // 2d array - one integer array for each subject - containing subject's choices - coded
                             // 1 for fractal A and 2 for fractal B. 
}

parameters {
  real alpha_mu; // mean of group-level learning rate distribution
  real<lower=0> alpha_sigma; // sd of group-level learning rate distribution
  vector[n_s] alpha_z; // z-score of subject-level alpha value (relative to group distribution)
 
  real beta_mu; // mean of group-level inverse temperature distribution
  real<lower=0> beta_sigma; // sd of group-level inverse temperature distribution
  vector[n_s] beta_z; // z-score of subject-level beta value (relative to group distribution)
  
  //repeat the above but for the box value trials - trials on which there was no fractal result shown
  real alpha_mu_nfr; 
  real<lower=0> alpha_sigma_nfr;
  vector[n_s] alpha_z_nfr; // z-score of subject-level alpha value (relative to group distribution)
  
  real beta_mu_nfr; // mean of group-level inverse temperature distribution
  real<lower=0> beta_sigma_nfr; // sd of group-level inverse temperature distribution
  vector[n_s] beta_z_nfr; // z-score of subject-level beta value (relative to group distribution)
  
  
  real alearn_mu; // mean of group-level learning rate distribution
  real<lower=0> alearn_sigma; // sd of group-level learning rate distribution
  vector[n_s] alearn_z; // z-score of subject-level alpha value (relative to group distribution)
 
  real asens_mu; // mean of group-level inverse temperature distribution
  real<lower=0> asens_sigma; // sd of group-level inverse temperature distribution
  vector[n_s] asens_z; // z-score of subject-level beta value (relative to group distribution)
  
  //repeat the above but for the box value trials - trials on which there was no fractal result shown
  real alearn_mu_nfr; 
  real<lower=0> alearn_sigma_nfr;
  vector[n_s] alearn_z_nfr; // z-score of subject-level alpha value (relative to group distribution)
  
  real asens_mu_nfr; // mean of group-level inverse temperature distribution
  real<lower=0> asens_sigma_nfr; // sd of group-level inverse temperature distribution
  vector[n_s] asens_z_nfr; // z-score of subject-level beta value (relative to group distribution)


  real phi_mu; //Mean of group-level autocorrelation weight distribution; in same units as beta - although C is
               //generally held within a tighter range, so phi would need to be a bit bigger (bit less than 2x prob) than
               //beta for autocorrelation to have the same effect as reward.
  real<lower=0> phi_sigma; //sd of group-level distribution
  vector[n_s] phi_z; //z-score of subject-level phi

  real tau_mu; //mean of group-level learning rate for C - the autocorrelation value
  real<lower=0> tau_sigma; //sd
  vector[n_s] tau_z; //z-score of subject-level tau
}

transformed parameters {
  vector[n_s*n_t] log_lik;
  {//anonymous_scope_start
    array[n_s,n_t] vector[n_f] Q; //Q value for fractal result outcome
    array[n_s,n_t] vector[n_f] Q_nfr; //Q value for box value trial outcomes
    
    array[n_s,n_t] vector[n_f] A; //Q value for fractal result outcome
    array[n_s,n_t] vector[n_f] A_nfr; //Q value for box value trial outcomes
    
    array[n_s,n_t] vector[n_f] C; // Choice autocorrelation value, of the same structure as Q

    real choice_a; // 1 if fractal A was chosen, 0 otherwise - used for C update
    real choice_b; // 1 if fractal B was chosen, 0 otherwise - used for C update

    vector[n_s] beta = beta_mu + beta_sigma*beta_z;  // get subject-level beta values using n.c.p.
    vector[n_s] beta_nfr = beta_mu_nfr + beta_sigma_nfr*beta_z_nfr; 
    vector[n_s] asens = asens_mu + asens_sigma*asens_z;  // get subject-level beta values using n.c.p.
    vector[n_s] asens_nfr = asens_mu_nfr + asens_sigma_nfr*asens_z_nfr; 
    vector[n_s] phi = phi_mu + phi_sigma*phi_z; //diddo phi
    vector[n_s] alpha = inv_logit(alpha_mu + alpha_sigma*alpha_z); //diddo alpha. Also, use inv_logit to get alpha between 0 and 1
    vector[n_s] alpha_nfr = inv_logit(alpha_mu_nfr + alpha_sigma_nfr*alpha_z_nfr); //diddo alpha. Also, use inv_logit to get alpha between 0 and 1
    vector[n_s] alearn = inv_logit(alearn_mu + alearn_sigma*alearn_z); //diddo alpha. Also, use inv_logit to get alpha between 0 and 1
    vector[n_s] alearn_nfr = inv_logit(alearn_mu_nfr + alearn_sigma_nfr*alearn_z_nfr); 
    vector[n_s] tau = inv_logit(tau_mu + tau_sigma*tau_z); //likewise for tau
    
    vector[2] softmax_arg;
    
    
    // loop through each subject and trial...
    for (s in 1:n_s) {
      for (t in 1:n_t){
        if(t == 1){
          Q[s,t] = rep_vector(0,n_f); // for the first trial of each subject, set all Q-values to 0
          Q_nfr[s,t] = rep_vector(0,n_f); 
          A[s,t] = rep_vector(0,n_f); // for the first trial of each subject, set all Q-values to 0
          A_nfr[s,t] = rep_vector(0,n_f); 
          C[s,t] = rep_vector(0,n_f);
        }
        
        softmax_arg = [beta[s]*Q[s,t,fA[s,t]] + beta_nfr[s]*Q_nfr[s,t,fA[s,t]] + asens[s]*A[s,t,fA[s,t]] + asens_nfr[s]*A_nfr[s,t,fA[s,t]] + phi[s]*C[s,t,fA[s,t]],
                       beta[s]*Q[s,t,fB[s,t]] + beta_nfr[s]*Q_nfr[s,t,fB[s,t]] + asens[s]*A[s,t,fB[s,t]] + asens_nfr[s]*A_nfr[s,t,fB[s,t]] + phi[s]*C[s,t,fB[s,t]]]'; //set softmax arg outside of the log_lik calculation for efficiency
                                                                         //(you should avoid transposing inside a likelihood funciton, apparently)
        log_lik[(s-1)*n_t+t] = categorical_logit_lpmf(choice[s,t] | softmax_arg); //get the likelihood of the choice on this trial. Doing the log lik

        // unlesss this is the subject's last trial, set the Q,C, and A values for the next trial
        if(t < n_t){
          Q[s,t+1] = (1-alpha[s])*Q[s,t]; // Decay Q values toward their initial values - 0. The only values
                                          // that will not be decayed are the ones chosen on this trial -
                                          // updated below.
          Q_nfr[s,t+1] = (1-alpha_nfr[s])*Q_nfr[s,t];
          A[s,t+1] = (1-alearn[s])*A[s,t]; 
          A_nfr[s,t+1] = (1-alearn_nfr[s])*A_nfr[s,t];
          if(choice[s,t] == 1){
            // if fA was chosen...
            //update C using 1 for choice a and 0 for choice b
            choice_a = 1;
            choice_b = 0;
            // if the outcome was the fractal result...
            if(fres[s,t] == 1){
              //update Q...
              Q[s,t+1,fA[s,t]] = Q[s,t,fA[s,t]] + alpha[s]*(out[s,t]-Q[s,t,fA[s,t]]);
              A[s,t+1,fA[s,t]] = A[s,t,fA[s,t]] + alearn[s]*affect[s,t];
            } else if(fres[s,t] == 0){
              //otherwise update Q_nfr
              Q_nfr[s,t+1,fA[s,t]] = Q_nfr[s,t,fA[s,t]] + alpha_nfr[s]*(out[s,t]-Q_nfr[s,t,fA[s,t]]);
              A_nfr[s,t+1,fA[s,t]] = A_nfr[s,t,fA[s,t]] + alearn_nfr[s]*affect[s,t];;
            }
          } else if(choice[s,t] == 2){
            // vice-versa if fB was chosen...
            choice_a = 0;
            choice_b = 1;
            // if the outcome was the fractal result...
            if(fres[s,t] == 1){
              //update Q...
              Q[s,t+1,fB[s,t]] = Q[s,t,fB[s,t]] + alpha[s]*(out[s,t]-Q[s,t,fB[s,t]]);
              A[s,t+1,fB[s,t]] = A[s,t,fB[s,t]] + alearn[s]*affect[s,t];
            } else if(fres[s,t] == 0){
              //otherwise update Q_nfr
              Q_nfr[s,t+1,fB[s,t]] = Q_nfr[s,t,fB[s,t]] + alpha_nfr[s]*(out[s,t]-Q_nfr[s,t,fB[s,t]]);
              A_nfr[s,t+1,fB[s,t]] = A_nfr[s,t,fB[s,t]] + alearn_nfr[s]*affect[s,t];
            }
          }
          //update C
          C[s,t+1] = C[s,t]; // Iniitalize the next trial's C values to be the same as the current trial's (no forgetting)
          C[s,t+1,fA[s,t]] = C[s,t,fA[s,t]] + tau[s]*(choice_a - C[s,t,fA[s,t]]); // update C value of fA with 1 if fA was chosen and 0 if it wasn't.
          C[s,t+1,fB[s,t]] = C[s,t,fB[s,t]] + tau[s]*(choice_b - C[s,t,fB[s,t]]); // diddo fB
        }
      }
    }
  }//anonymous_scope_end
}

model {
   // add the hyperpriors on the group-level parameters to the target density
  alpha_mu ~ normal(-.05,1.7); //This yields a nearly uniform 0-1 prior when logit transformed
  alpha_sigma ~ normal(0,4); //allow for large inter-subject heterogeneity without enforcing marginal subject-level priors that are excessively horeshoe-shaped
  alpha_mu_nfr ~ normal(-.05,1.7); //This yields a nearly uniform 0-1 prior when logit transformed
  alpha_sigma_nfr ~ normal(0,4);
  alearn_mu ~ normal(-.05,1.7); //This yields a nearly uniform 0-1 prior when logit transformed
  alearn_sigma ~ normal(0,4); //allow for large inter-subject heterogeneity without enforcing marginal subject-level priors that are excessively horeshoe-shaped
  alearn_mu_nfr ~ normal(-.05,1.7); //This yields a nearly uniform 0-1 prior when logit transformed
  alearn_sigma_nfr ~ normal(0,4);
  tau_mu ~ normal(-.05,1.7); //setting tau priors same as alpha
  tau_sigma ~ normal(0,4); // ""
  beta_mu ~ normal(1,5); // Weakly informative prior
  beta_sigma ~ normal(0,5); //Setting the prior on the SD with similar logic
  beta_mu_nfr ~ normal(1,5); // Weakly informative prior
  beta_sigma_nfr ~ normal(0,5); //Setting the prior on the SD with similar logic
  asens_mu ~ normal(1,5); // Weakly informative prior
  asens_sigma ~ normal(0,5); //Setting the prior on the SD with similar logic
  asens_mu_nfr ~ normal(1,5); // Weakly informative prior
  asens_sigma_nfr ~ normal(0,5); //Setting the prior on the SD with similar logic
  phi_mu ~ normal(2,10); // Setting a reguarlizing prior. Matching the phi prior to the beta prior, assuming that uatocorrelation and reward are apt to have
                         // equivalent effects on choice. But differences in phi between fA and fB will be about .5 at asymptote (based on data in which
                         // subjects tend to repeat choices 75% of the time), whereas it's more like 1.2 for beta, so the prior for phi will be set a bit wider.
  phi_sigma ~ normal(0,10);


  // add the "priors" on the subject-level parameters derived form the group-level parameters to the target density
  alpha_z ~ normal(0,1);
  alpha_z_nfr ~ normal(0,1);
  alearn_z ~ normal(0,1);
  alearn_z_nfr ~ normal(0,1);
  beta_z ~ normal(0,1);
  beta_z_nfr ~ normal(0,1);
  asens_z ~ normal(0,1);
  asens_z_nfr ~ normal(0,1);
  phi_z ~ normal(0,1);
  tau_z ~ normal(0,1);

  // Add the joint likelihood to the target density, mapping beta-adjusted Q values to choice probabilities
  // using a softmax function. Unfortunately, categorical_logit does not support vectorization.
  target += log_lik;
}
