data {
  int<lower=1> n_t; // number of trials
  int<lower=1> n_s; // number of subjects
  int<lower=1> n_f; //number of fractals
  array[n_s,n_t] int fA; // 2d array - one array for each subject, containing a sub-arry with the identity of
                          // fractal A on each of that subject's trials (fractals are identified with numbers 1:n_f)
  array[n_s,n_t] int fB; // diddo for fractal B identities
  array[n_s] vector[n_t] rew; 
  array[n_s] vector[n_t] fi_rat; 
  array[n_s] vector[n_t] resid; 
  array[n_s,n_t] int fres; //2d array (subjects by trials) indicating whether the fractal result was shown on that trial (1-0)
  array[n_s,n_t] int choice; // 2d array - one integer array for each subject - containing subject's choices - coded
                             // 1 for fractal A and 2 for fractal B. 
}

parameters{
  real rewnf_sens_mu; //mean effects of reward (fres) on affect and choice
  real<lower=0> rewnf_sens_sigma;
  vector[n_s] rewnf_sens_z;
  
  real af_sens_mu; //mean effects of reward (fres) on affect and choice
  real<lower=0> af_sens_sigma;
  vector[n_s] af_sens_z;
  
  real anf_sens_mu; //mean effects of reward (fres) on affect and choice
  real<lower=0> anf_sens_sigma;
  vector[n_s] anf_sens_z;
  
  //autocorrelation sens.
  real csens_mu;
  real<lower=0> csens_sigma;
  vector[n_s] csens_z;
  
  //decay rates
  real dcy_f_mu;
  real<lower=0> dcy_f_sigma;
  vector[n_s] dcy_f_z;
  
  real dcy_nf_mu;
  real<lower=0> dcy_nf_sigma;
  vector[n_s] dcy_nf_z;
  
  //learning rates
  real lrn_f_mu;
  real<lower=0> lrn_f_sigma;
  vector[n_s] lrn_f_z;
  
  real lrn_nf_mu;
  real<lower=0> lrn_nf_sigma;
  vector[n_s] lrn_nf_z;
  
  real lrn_c_mu;
  real<lower=0> lrn_c_sigma;
  vector[n_s] lrn_c_z;
}

transformed parameters {
  vector[n_s*n_t] log_lik;
  {//anonymous_scope_start
    array[n_s,n_t] vector[n_f] Q_nf; 
    array[n_s,n_t] vector[n_f] A; //Q value for fractal result outcome
    array[n_s,n_t] vector[n_f] A_nf; 
    array[n_s,n_t] vector[n_f] C; // Choice autocorrelation value, of the same structure as Q
    real choice_a; // 1 if fractal A was chosen, 0 otherwise - used for C update
    real choice_b; // 1 if fractal B was chosen, 0 otherwise - used for C update
    vector[2] softmax_arg; //goes inside choice likelihood calculation
    
    //NCP for normal priors
    vector[n_s] csens = csens_mu + csens_sigma*csens_z;
    vector[n_s] rewnf_sens = rewnf_sens_mu + rewnf_sens_sigma*rewnf_sens_z;
    vector[n_s] af_sens = af_sens_mu + af_sens_sigma*af_sens_z;
    vector[n_s] anf_sens = anf_sens_mu + anf_sens_sigma*anf_sens_z;
    
    vector[n_s] dcy_f = inv_logit(dcy_f_mu + dcy_f_sigma*dcy_f_z);
    vector[n_s] dcy_nf = inv_logit(dcy_nf_mu + dcy_nf_sigma*dcy_nf_z);

    vector[n_s] lrn_f = inv_logit(lrn_f_mu + lrn_f_sigma*lrn_f_z);
    vector[n_s] lrn_nf = inv_logit(lrn_nf_mu + lrn_nf_sigma*lrn_nf_z);
    vector[n_s] lrn_c = inv_logit(lrn_c_mu + lrn_c_sigma*lrn_c_z);
    
    // loop through each subject and trial...
    for (s in 1:n_s) {
      for (t in 1:n_t){
        if(t == 1){
          Q_nf[s,t] = rep_vector(0,n_f); 
          A[s,t] = rep_vector(0,n_f); // for the first trial of each subject, set all Q-values to 0
          A_nf[s,t] = rep_vector(0,n_f); 
          C[s,t] = rep_vector(0,n_f);
        }
        
        //set softmax arg outside of the log_lik calculation for efficiency (you should avoid transposing inside a likelihood funciton, apparently)
        softmax_arg = af_sens[s]*A[s,t,{fA[s,t],fB[s,t]}] 
                      + anf_sens[s]*A_nf[s,t,{fA[s,t],fB[s,t]}] + rewnf_sens[s]*Q_nf[s,t,{fA[s,t],fB[s,t]}] 
                      + csens[s]*C[s,t,{fA[s,t],fB[s,t]}];

        log_lik[(s-1)*n_t+t] = categorical_logit_lpmf(choice[s,t] | softmax_arg); //get the likelihood of the choice on this trial. Doing the log lik
        
        // unlesss this is the subject's last trial, set the Q,C, and A values for the next trial
        if(t < n_t){
          //decay all values toward 0
          Q_nf[s,t+1] = (1-dcy_nf[s])*Q_nf[s,t];
          A[s,t+1] = (1-dcy_f[s])*A[s,t];
          A_nf[s,t+1] = (1-dcy_nf[s])*A_nf[s,t];
          if(choice[s,t] == 1){
            // if fA was chosen...
            //update C using 1 for choice a and 0 for choice b
            choice_a = 1;
            choice_b = 0;
            // if the outcome was the fractal result...
            if(fres[s,t] == 1){
              //update Q...
              A[s,t+1,fA[s,t]] = A[s,t,fA[s,t]] + lrn_f[s]*(fi_rat[s,t] - A[s,t,fA[s,t]]);
            } else if(fres[s,t] == 0){
              //otherwise update Q_nf
              Q_nf[s,t+1,fA[s,t]] = Q_nf[s,t,fA[s,t]] + lrn_nf[s]*(rew[s,t] - Q_nf[s,t,fA[s,t]]);
              A_nf[s,t+1,fA[s,t]] = A_nf[s,t,fA[s,t]] + lrn_nf[s]*(resid[s,t] - A_nf[s,t,fA[s,t]]);
            }
          } else if(choice[s,t] == 2){
            // vice-versa if fB was chosen...
            choice_a = 0;
            choice_b = 1;
            // if the outcome was the fractal result...
            if(fres[s,t] == 1){
              //update Q...
              A[s,t+1,fB[s,t]] = A[s,t,fB[s,t]] + lrn_f[s]*(fi_rat[s,t] - A[s,t,fB[s,t]]);
            } else if(fres[s,t] == 0){
              //otherwise update Q_nf
              Q_nf[s,t+1,fB[s,t]] = Q_nf[s,t,fB[s,t]] + lrn_nf[s]*(rew[s,t] - Q_nf[s,t,fB[s,t]]);
              A_nf[s,t+1,fB[s,t]] = A_nf[s,t,fB[s,t]] + lrn_nf[s]*(resid[s,t] - A_nf[s,t,fB[s,t]]);
            }
          }
          //update C
          C[s,t+1] = C[s,t]; // Iniitalize the next trial's C values to be the same as the current trial's (no forgetting)
          C[s,t+1,fA[s,t]] = C[s,t,fA[s,t]] + lrn_c[s]*(choice_a - C[s,t,fA[s,t]]); // update C value of fA with 1 if fA was chosen and 0 if it wasn't.
          C[s,t+1,fB[s,t]] = C[s,t,fB[s,t]] + lrn_c[s]*(choice_b - C[s,t,fB[s,t]]); // diddo fB
        }
      }
    }
  }//anonymous_scope_end
}
model{
  //effects on choice
  rewnf_sens_mu ~ normal(0,4); 
  rewnf_sens_sigma ~ normal(0,6); 
  
  af_sens_mu ~ normal(0,4); 
  af_sens_sigma ~ normal(0,6); 
  anf_sens_mu ~ normal(0,6); 
  anf_sens_sigma ~ normal(0,8); 
  
  //sensitivity to affect residuals and uatocorrelation
  csens_mu ~ normal(0,5);
  csens_sigma ~ normal(0,10);
  
  //learning and decay rates
  dcy_f_mu ~ normal(-.05,1.7);
  dcy_f_sigma ~ normal(0,4);
  dcy_nf_mu ~ normal(-.05,1.7);
  dcy_nf_sigma ~ normal(0,4);
  lrn_f_mu ~ normal(-.05,1.7);
  lrn_f_sigma ~ normal(0,4);
  lrn_nf_mu ~ normal(-.05,1.7);
  lrn_nf_sigma ~ normal(0,4);
  lrn_c_mu ~ normal(-.05,1.7);
  lrn_c_sigma ~ normal(0,4);
  
  //zs
  rewnf_sens_z ~ std_normal();
  af_sens_z ~ std_normal();
  anf_sens_z ~ std_normal();
  csens_z ~ std_normal();
  dcy_f_z ~ std_normal();
  dcy_nf_z ~ std_normal();
  lrn_f_z ~ std_normal();
  lrn_nf_z ~ std_normal();
  lrn_c_z ~ std_normal();
  
  target += log_lik;
}