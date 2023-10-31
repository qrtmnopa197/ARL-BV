data {
  int<lower=1> n_t; // number of trials
  int<lower=1> n_s; // number of subjects
  int<lower=1> n_f; //number of fractals
  array[n_s,n_t] int fA; // 2d array - one array for each subject, containing a sub-arry with the identity of
                          // fractal A on each of that subject's trials (fractals are identified with numbers 1:n_f)
  array[n_s,n_t] int fB; // diddo for fractal B identities
  array[n_s] vector[n_t] rew; 
  array[n_s] vector[n_t] fi_rat; 
  array[n_s,n_t] int fres; //2d array (subjects by trials) indicating whether the fractal result was shown on that trial (1-0)
  array[n_s,n_t] int choice; // 2d array - one integer array for each subject - containing subject's choices - coded
                             // 1 for fractal A and 2 for fractal B. 
}

parameters{
  real rew_sens_mu; //mean effects of reward (fres) on affect and choice
  real<lower=0> rew_sens_sigma;
  vector[n_s] rew_sens_z;
  
  real rew_omega_mu;
  real<lower=0> rew_omega_sigma;
  vector[n_s] rew_omega_z;

  real aff_sens_mu; //mean effects of reward (fres) on affect and choice
  real<lower=0> aff_sens_sigma;
  vector[n_s] aff_sens_z;
  
  real aff_omega_mu;
  real<lower=0> aff_omega_sigma;
  vector[n_s] aff_omega_z;
  
  //autocorrelation sens.
  real csens_mu;
  real<lower=0> csens_sigma;
  vector[n_s] csens_z;
  
  //decay rates
  real dcy_mu;
  real<lower=0> dcy_sigma;
  vector[n_s] dcy_z;

  //learning rates
  real lrn_q_mu;
  real<lower=0> lrn_q_sigma;
  vector[n_s] lrn_q_z;
  
  real lrn_a_mu;
  real<lower=0> lrn_a_sigma;
  vector[n_s] lrn_a_z;
  
  real lrn_c_mu;
  real<lower=0> lrn_c_sigma;
  vector[n_s] lrn_c_z;
}

transformed parameters {
  vector[n_s*n_t] log_lik;
  {//anonymous_scope_start
    array[n_s,n_t] vector[n_f] Q; //Q value for fractal result outcome
    array[n_s,n_t] vector[n_f] A; //Q value for fractal result outcome
    array[n_s,n_t] vector[n_f] C; // Choice autocorrelation value, of the same structure as Q
    real choice_a; // 1 if fractal A was chosen, 0 otherwise - used for C update
    real choice_b; // 1 if fractal B was chosen, 0 otherwise - used for C update
    vector[2] softmax_arg; //goes inside choice likelihood calculation
    
    //NCP for normal priors
    vector[n_s] csens = csens_mu + csens_sigma*csens_z;
    vector[n_s] rew_sens = rew_sens_mu + rew_sens_sigma*rew_sens_z;
    vector[n_s] aff_sens = aff_sens_mu + aff_sens_sigma*aff_sens_z;

    //NCP for normal priors with sigmoid transform
    vector[n_s] dcy = inv_logit(dcy_mu + dcy_sigma*dcy_z);
    
    vector[n_s] lrn_q = inv_logit(lrn_q_mu + lrn_q_sigma*lrn_q_z);
    vector[n_s] lrn_a = inv_logit(lrn_a_mu + lrn_a_sigma*lrn_a_z);
    vector[n_s] lrn_c = inv_logit(lrn_c_mu + lrn_c_sigma*lrn_c_z);
    
    vector[n_s] rew_omega = inv_logit(rew_omega_mu + rew_omega_sigma*rew_omega_z);
    vector[n_s] aff_omega = inv_logit(aff_omega_mu + aff_omega_sigma*aff_omega_z);
    
    
    // loop through each subject and trial...
    for (s in 1:n_s) {
      for (t in 1:n_t){
        if(t == 1){
          Q[s,t] = rep_vector(0,n_f); // for the first trial of each subject, set all Q-values to 0
          A[s,t] = rep_vector(0,n_f); // for the first trial of each subject, set all Q-values to 0
          C[s,t] = rep_vector(0,n_f);
        }
        
        //set softmax arg outside of the log_lik calculation for efficiency (you should avoid transposing inside a likelihood funciton, apparently)
        softmax_arg = rew_sens[s]*Q[s,t,{fA[s,t],fB[s,t]}] + aff_sens[s]*A[s,t,{fA[s,t],fB[s,t]}] + csens[s]*C[s,t,{fA[s,t],fB[s,t]}];

        log_lik[(s-1)*n_t+t] = categorical_logit_lpmf(choice[s,t] | softmax_arg); //get the likelihood of the choice on this trial. Doing the log lik
        
        // unlesss this is the subject's last trial, set the Q,C, and A values for the next trial
        if(t < n_t){
          //decay all values toward 0
          Q[s,t+1] = (1-dcy[s])*Q[s,t];
          A[s,t+1] = (1-dcy[s])*A[s,t];
          if(choice[s,t] == 1){
            // if fA was chosen...
            //update C using 1 for choice a and 0 for choice b
            choice_a = 1;
            choice_b = 0;
            // if the outcome was the fractal result...
            if(fres[s,t] == 1){
              //update Q...
              Q[s,t+1,fA[s,t]] = Q[s,t,fA[s,t]] + lrn_q[s]*(rew[s,t]*rew_omega[s] - Q[s,t,fA[s,t]]);
              A[s,t+1,fA[s,t]] = A[s,t,fA[s,t]] + lrn_a[s]*(fi_rat[s,t]*aff_omega[s] - A[s,t,fA[s,t]]);
            } else if(fres[s,t] == 0){
              //otherwise update Q_nf
              Q[s,t+1,fA[s,t]] = Q[s,t,fA[s,t]] + lrn_q[s]*(rew[s,t]*(1-rew_omega[s]) - Q[s,t,fA[s,t]]);
              A[s,t+1,fA[s,t]] = A[s,t,fA[s,t]] + lrn_a[s]*(fi_rat[s,t]*(1-aff_omega[s]) - A[s,t,fA[s,t]]);
            }
          } else if(choice[s,t] == 2){
            // vice-versa if fB was chosen...
            choice_a = 0;
            choice_b = 1;
            // if the outcome was the fractal result...
            if(fres[s,t] == 1){
              //update Q...
              Q[s,t+1,fB[s,t]] = Q[s,t,fB[s,t]] + lrn_q[s]*(rew[s,t]*rew_omega[s] - Q[s,t,fB[s,t]]);
              A[s,t+1,fB[s,t]] = A[s,t,fB[s,t]] + lrn_a[s]*(fi_rat[s,t]*aff_omega[s] - A[s,t,fB[s,t]]);
            } else if(fres[s,t] == 0){
              //otherwise update Q_nf
              Q[s,t+1,fB[s,t]] = Q[s,t,fB[s,t]] + lrn_q[s]*(rew[s,t]*(1-rew_omega[s]) - Q[s,t,fB[s,t]]);
              A[s,t+1,fB[s,t]] = A[s,t,fB[s,t]] + lrn_a[s]*(fi_rat[s,t]*(1-aff_omega[s]) - A[s,t,fB[s,t]]);
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
  rew_sens_mu ~ normal(0,125); 
  rew_sens_sigma ~ normal(0,175); 
  
  aff_sens_mu ~ normal(0,6); 
  aff_sens_sigma ~ normal(0,8); 

  //sensitivity to affect residuals and uatocorrelation
  csens_mu ~ normal(0,5);
  csens_sigma ~ normal(0,10);
  
  //learning and decay rates
  dcy_mu ~ normal(-.05,1.7);
  dcy_sigma ~ normal(0,4);
  
  lrn_q_mu ~ normal(-.05,1.7);
  lrn_q_sigma ~ normal(0,4);
  lrn_a_mu ~ normal(-.05,1.7);
  lrn_a_sigma ~ normal(0,4);
  lrn_c_mu ~ normal(-.05,1.7);
  lrn_c_sigma ~ normal(0,4);
  rew_omega_mu ~ normal(-.05,1.7);
  rew_omega_sigma ~  normal(0,4);
  aff_omega_mu ~ normal(-.05,1.7);
  aff_omega_sigma ~ normal(0,4);
  
  //zs
  rew_sens_z ~ std_normal();
  aff_sens_z ~ std_normal();
  csens_z ~ std_normal();
  dcy_z ~ std_normal();
  lrn_q_z ~ std_normal();
  lrn_a_z ~ std_normal();
  lrn_c_z ~ std_normal();
  rew_omega_z ~ std_normal();
  aff_omega_z ~ std_normal();
  
  target += log_lik;
}


