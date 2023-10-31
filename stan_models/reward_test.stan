data {
  int<lower=1> n_t; // number of trials
  int<lower=1> n_s; // number of subjects
  int<lower=1> n_f; //number of fractals
  array[n_s,n_t] int fA; // 2d array - one array for each subject, containing a sub-arry with the identity of
                          // fractal A on each of that subject's trials (fractals are identified with numbers 1:n_f)
  array[n_s,n_t] int fB; // diddo for fractal B identities
  array[n_s] vector[n_t] rew; 
  array[n_s,n_t] int choice; // 2d array - one integer array for each subject - containing subject's choices - coded
                             // 1 for fractal A and 2 for fractal B. 
  array[n_s,n_t] int rat_num; //rating number for each trial
  int n_rat;
  vector[n_rat] rat; //affect ratings
}

parameters{
  real incp_mu; //intercept for raitng prediction
  real<lower=0> incp_sigma;
  vector[n_s] incp_z;
  
  vector[2] rewf_w_mu; //mean effects of reward (fres) on affect and choice
  vector<lower=0>[2] rewf_w_tau;
  cholesky_factor_corr[2] rewf_w_lOmega;
  matrix[2,n_s] rewf_w_z;
  
  //learning rates
  real lrn_q_mu;
  real<lower=0> lrn_q_sigma;
  vector[n_s] lrn_q_z;
  
  real<lower=0> resid_sigma;
}

transformed parameters {
  vector[n_s*n_t] choice_lik;
  vector[n_rat] rat_pred;
  {//anonymous_scope_start
    array[n_s,n_t] vector[n_f] Q; //Q value for fractal result outcome
    real resid; //affect residual on the current trial
    vector[2] softmax_arg; //goes inside choice likelihood calculation
    
    //NCP for normal priors
    vector[n_s] incp = incp_mu + incp_sigma*incp_z;
    vector[n_s] lrn_q = inv_logit(lrn_q_mu + lrn_q_sigma*lrn_q_z);
  
    //NCP for multivariate normal priors
    matrix[2,n_s] rewf_w = rep_matrix(rewf_w_mu,n_s) + diag_pre_multiply(rewf_w_tau, rewf_w_lOmega) * rewf_w_z; //get subject-level weights for rewf
    
    // loop through each subject and trial...
    for (s in 1:n_s) {
      for (t in 1:n_t){
        //Choice predictions
        if(t == 1){
          Q[s,t] = rep_vector(0,n_f); // for the first trial of each subject, set all Q-values to 0
        }
        
        //set softmax arg outside of the log_lik calculation for efficiency (you should avoid transposing inside a likelihood funciton, apparently)
        softmax_arg = rewf_w[2,s]*Q[s,t,{fA[s,t],fB[s,t]}];

        choice_lik[(s-1)*n_t+t] = categorical_logit_lpmf(choice[s,t] | softmax_arg); //get the likelihood of the choice on this trial. Doing the log lik
        
        // Generating rating prediction and residual, if a rating was made
        if(rat_num[s,t] != 0){
          rat_pred[rat_num[s,t]] = incp[s] +  rewf_w[1,s]*rew[s,t];
          resid = rat[rat_num[s,t]] - rat_pred[rat_num[s,t]];
        } else{
          resid = 0; //if there was no rating on this trial, set the residual to 0
        }

        // unlesss this is the subject's last trial, set the Q,C, and A values for the next trial
        if(t < n_t){
          Q[s,t+1] = Q[s,t];
          if(choice[s,t] == 1){
            // if the outcome was the fractal result...
            Q[s,t+1,fA[s,t]] = Q[s,t,fA[s,t]] + lrn_q[s]*(rew[s,t] - Q[s,t,fA[s,t]]);
          } else if(choice[s,t] == 2){
            Q[s,t+1,fB[s,t]] = Q[s,t,fB[s,t]] + lrn_q[s]*(rew[s,t] - Q[s,t,fB[s,t]]);
          }
        }
      }
    }
  }//anonymous_scope_end
}
model{
  incp_mu ~ normal(0,3);
  incp_sigma ~ normal(0,5);
  
  //effects on affect
  rewf_w_mu[2] ~ normal(0,125); 
  rewf_w_tau[2] ~ normal(0,175); 
  //effects on choice
  rewf_w_mu[1] ~ normal(0,20); 
  rewf_w_tau[1] ~ normal(0,30); 
  //nearly uniform prior on correlations (slighly biased toward 0)
  rewf_w_lOmega ~ lkj_corr_cholesky(1);
  
 
  lrn_q_mu ~ normal(-.05,1.7);
  lrn_q_sigma ~ normal(0,4);
 
  resid_sigma ~ normal(0,2);
  
  incp_z ~ std_normal();
  to_vector(rewf_w_z) ~ std_normal();
  lrn_q_z ~ std_normal();
  
  target += choice_lik;
  rat ~ normal(rat_pred, resid_sigma);
  
  //decay rates
}
generated quantities{
  //get correlations between effects on affect and choice
  corr_matrix[2] rewf_w_Omega = rewf_w_lOmega * rewf_w_lOmega';

  //get effects of variables on choice when their effect on affect is 0
  real rewf_zero = rewf_w_mu[2] - (rewf_w_Omega[1,2]*rewf_w_tau[2]*rewf_w_mu[1])/rewf_w_tau[1];
}



