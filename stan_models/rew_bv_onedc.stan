data {
  int<lower=1> n_t; // number of trials
  int<lower=1> n_s; // number of subjects
  int<lower=1> n_f; //number of fractals
  array[n_s,n_t] int fA; // 2d array - one array for each subject, containing a sub-arry with the identity of
                          // fractal A on each of that subject's trials (fractals are identified with numbers 1:n_f)
  array[n_s,n_t] int fB; // diddo for fractal B identities
  array[n_s] vector[n_t] rew; 
  array[n_s] vector[n_t] bv; 
  array[n_s,n_t] int fres; //2d array (subjects by trials) indicating whether the fractal result was shown on that trial (1-0)
  array[n_s,n_t] int choice; // 2d array - one integer array for each subject - containing subject's choices - coded
                             // 1 for fractal A and 2 for fractal B. 
  array[n_s,n_t] int rat_num; //rating number for each trial
  int n_rat;
  vector[n_rat] rat; //affect ratings
}

parameters{
  real lsbias_mu;
  real<lower=0> lsbias_sigma;
  vector[n_s] lsbias_z;
  
  real incp_mu; //intercept for raitng prediction
  real<lower=0> incp_sigma;
  vector[n_s] incp_z;
  
  vector[2] rewf_w_mu; //mean effects of reward (fres) on affect and choice
  vector<lower=0>[2] rewf_w_tau;
  cholesky_factor_corr[2] rewf_w_lOmega;
  matrix[2,n_s] rewf_w_z;
  
  vector[2] rewnf_w_mu; //mean effects of reward (fres) on affect and choice
  vector<lower=0>[2] rewnf_w_tau;
  cholesky_factor_corr[2] rewnf_w_lOmega;
  matrix[2,n_s] rewnf_w_z;
  
  vector[2] bvf_w_mu; //mean effects of reward (fres) on affect and choice
  vector<lower=0>[2] bvf_w_tau;
  cholesky_factor_corr[2] bvf_w_lOmega;
  matrix[2,n_s] bvf_w_z;
  
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
  
  real lrn_b_mu;
  real<lower=0> lrn_b_sigma;
  vector[n_s] lrn_b_z;
  
  real lrn_qnf_mu;
  real<lower=0> lrn_qnf_sigma;
  vector[n_s] lrn_qnf_z;
  
  real lrn_c_mu;
  real<lower=0> lrn_c_sigma;
  vector[n_s] lrn_c_z;
  
  real<lower=0> resid_sigma;
}

transformed parameters {
  vector[n_s*n_t] choice_lik;
  vector[n_rat] rat_pred;
  {//anonymous_scope_start
    array[n_s,n_t] vector[n_f] Q; //Q value for fractal result outcome
    array[n_s,n_t] vector[n_f] B; 
    array[n_s,n_t] vector[n_f] Q_nf; 
    array[n_s,n_t] vector[n_f] C; // Choice autocorrelation value, of the same structure as Q
    real choice_a; // 1 if fractal A was chosen, 0 otherwise - used for C update
    real choice_b; // 1 if fractal B was chosen, 0 otherwise - used for C update
    vector[2] softmax_arg; //goes inside choice likelihood calculation
    
    //NCP for normal priors
    vector[n_s] lsbias = lsbias_mu + lsbias_sigma*lsbias_z;
    vector[n_s] incp = incp_mu + incp_sigma*incp_z;
    vector[n_s] csens = csens_mu + csens_sigma*csens_z;
    
    //NCP for normal priors with sigmoid transform
    vector[n_s] dcy = inv_logit(dcy_mu + dcy_sigma*dcy_z);
    
    vector[n_s] lrn_q = inv_logit(lrn_q_mu + lrn_q_sigma*lrn_q_z);
    vector[n_s] lrn_b = inv_logit(lrn_b_mu + lrn_b_sigma*lrn_b_z);
    vector[n_s] lrn_qnf = inv_logit(lrn_qnf_mu + lrn_qnf_sigma*lrn_qnf_z);
    vector[n_s] lrn_c = inv_logit(lrn_c_mu + lrn_c_sigma*lrn_c_z);
  
    //NCP for multivariate normal priors
    matrix[2,n_s] rewf_w = rep_matrix(rewf_w_mu,n_s) + diag_pre_multiply(rewf_w_tau, rewf_w_lOmega) * rewf_w_z; //get subject-level weights for rewf
    matrix[2,n_s] rewnf_w = rep_matrix(rewnf_w_mu,n_s) + diag_pre_multiply(rewnf_w_tau, rewnf_w_lOmega) * rewnf_w_z; //get subject-level weights for rewf
    matrix[2,n_s] bvf_w = rep_matrix(bvf_w_mu,n_s) + diag_pre_multiply(bvf_w_tau, bvf_w_lOmega) * bvf_w_z; //get subject-level weights for rewf
    
    // loop through each subject and trial...
    for (s in 1:n_s) {
      for (t in 1:n_t){
        if(t == 1){
          Q[s,t] = rep_vector(0,n_f); // for the first trial of each subject, set all Q-values to 0
          B[s,t] = rep_vector(0,n_f); 
          Q_nf[s,t] = rep_vector(0,n_f); 
          C[s,t] = rep_vector(0,n_f);
        }
        
        //set softmax arg outside of the log_lik calculation for efficiency (you should avoid transposing inside a likelihood funciton, apparently)
        softmax_arg = rewf_w[2,s]*Q[s,t,{fA[s,t],fB[s,t]}] + bvf_w[2,s]*B[s,t,{fA[s,t],fB[s,t]}] + 
                      rewnf_w[2,s]*Q_nf[s,t,{fA[s,t],fB[s,t]}] + csens[s]*C[s,t,{fA[s,t],fB[s,t]}] + [lsbias[s],0]';

        choice_lik[(s-1)*n_t+t] = categorical_logit_lpmf(choice[s,t] | softmax_arg); //get the likelihood of the choice on this trial. Doing the log lik
        
        // Generating rating prediction and residual, if a rating was made
        if(rat_num[s,t] != 0){
          if(fres[s,t] == 1){
            rat_pred[rat_num[s,t]] = incp[s] + rewf_w[1,s]*rew[s,t] + bvf_w[1,s]*bv[s,t];
          }else if (fres[s,t] == 0){
            rat_pred[rat_num[s,t]] = incp[s] + rewnf_w[1,s]*rew[s,t];
          }
        }

        // unlesss this is the subject's last trial, set the Q,C, and A values for the next trial
        if(t < n_t){
          //decay all values toward 0
          Q[s,t+1] = (1-dcy[s])*Q[s,t];
          B[s,t+1] = (1-dcy[s])*B[s,t];
          Q_nf[s,t+1] = (1-dcy[s])*Q_nf[s,t];
          if(choice[s,t] == 1){
            // if fA was chosen...
            //update C using 1 for choice a and 0 for choice b
            choice_a = 1;
            choice_b = 0;
            // if the outcome was the fractal result...
            if(fres[s,t] == 1){
              //update Q...
              Q[s,t+1,fA[s,t]] = Q[s,t,fA[s,t]] + lrn_q[s]*(rew[s,t] - Q[s,t,fA[s,t]]);
              B[s,t+1,fA[s,t]] = B[s,t,fA[s,t]] + lrn_b[s]*(bv[s,t] - B[s,t,fA[s,t]]);
            } else if(fres[s,t] == 0){
              //otherwise update Q_nf
              Q_nf[s,t+1,fA[s,t]] = Q_nf[s,t,fA[s,t]] + lrn_qnf[s]*(rew[s,t] - Q_nf[s,t,fA[s,t]]);
            }
          } else if(choice[s,t] == 2){
            // vice-versa if fB was chosen...
            choice_a = 0;
            choice_b = 1;
            // if the outcome was the fractal result...
            if(fres[s,t] == 1){
              //update Q...
              Q[s,t+1,fB[s,t]] = Q[s,t,fB[s,t]] + lrn_q[s]*(rew[s,t] - Q[s,t,fB[s,t]]);
              B[s,t+1,fB[s,t]] = B[s,t,fB[s,t]] + lrn_b[s]*(bv[s,t] - B[s,t,fB[s,t]]);
            } else if(fres[s,t] == 0){
              //otherwise update Q_nf
              Q_nf[s,t+1,fB[s,t]] = Q_nf[s,t,fB[s,t]] + lrn_qnf[s]*(rew[s,t] - Q_nf[s,t,fB[s,t]]);
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
  lsbias_mu ~ normal(0,4);
  lsbias_sigma ~ normal(0,8);
  
  incp_mu ~ normal(0,2);
  incp_sigma ~ normal(0,4);
  
  //effects on affect
  rewf_w_mu[1] ~ normal(0,40); 
  rewf_w_tau[1] ~ normal(0,60); 
  rewnf_w_mu[1] ~ normal(0,.7); 
  rewnf_w_tau[1] ~ normal(0,1); 
  bvf_w_mu[1] ~ normal(0,.7);
  bvf_w_tau[1] ~ normal(0,1);
  //effects on choice
  rewf_w_mu[2] ~ normal(0,125); 
  rewf_w_tau[2] ~ normal(0,175); 
  rewnf_w_mu[2] ~ normal(0,4); 
  rewnf_w_tau[2] ~ normal(0,6); 
  bvf_w_mu[2] ~ normal(0,4);
  bvf_w_tau[2] ~ normal(0,6);
  
  //nearly uniform prior on correlations (slighly biased toward 0)
  rewf_w_lOmega ~ lkj_corr_cholesky(1);
  rewnf_w_lOmega ~ lkj_corr_cholesky(1);
  bvf_w_lOmega ~ lkj_corr_cholesky(1);
  
  //sensitivity to affect residuals and uatocorrelation
  csens_mu ~ normal(0,5);
  csens_sigma ~ normal(0,10);
  
  //learning and decay rates
  dcy_mu ~ normal(-.05,1.7);
  dcy_sigma ~ normal(0,4);
  
  lrn_q_mu ~ normal(-.05,1.7);
  lrn_q_sigma ~ normal(0,4);
  lrn_b_mu ~ normal(-.05,1.7);
  lrn_b_sigma ~ normal(0,4);
  lrn_qnf_mu ~ normal(-.05,1.7);
  lrn_qnf_sigma ~ normal(0,4);
  lrn_c_mu ~ normal(-.05,1.7);
  lrn_c_sigma ~ normal(0,4);
  
  resid_sigma ~ normal(0,2);
  
  //zs
  lsbias_z ~ std_normal();
  incp_z ~ std_normal();
  to_vector(rewf_w_z) ~ std_normal();
  to_vector(rewnf_w_z) ~ std_normal();
  to_vector(bvf_w_z) ~ std_normal();
  csens_z ~ std_normal();
  dcy_z ~ std_normal();
  lrn_q_z ~ std_normal();
  lrn_b_z ~ std_normal();
  lrn_qnf_z ~ std_normal();
  lrn_c_z ~ std_normal();
  
  target += choice_lik;
  rat ~ normal(rat_pred, resid_sigma);
}
generated quantities{
  //get effects of variables on choice when their effect on affect is 0
  real rewf_zero = rewf_w_mu[2] - (rewf_w_lOmega[2,1]*rewf_w_tau[2]*rewf_w_mu[1])/rewf_w_tau[1];
  real rewnf_zero = rewnf_w_mu[2] - (rewnf_w_lOmega[2,1]*rewnf_w_tau[2]*rewnf_w_mu[1])/rewnf_w_tau[1];
  real bvf_zero = bvf_w_mu[2] - (bvf_w_lOmega[2,1]*bvf_w_tau[2]*bvf_w_mu[1])/bvf_w_tau[1];
  
  vector[n_rat] affect_lik; //log likelihoods of all affect rating
  vector[n_s*n_t+n_rat] log_lik; //log likelihoods of all choices followed by log likelihoods of all affect ratings
  for(i in 1:n_rat){
    affect_lik[i] = normal_lpdf(rat[i] | rat_pred[i], resid_sigma);
  }
  log_lik = append_row(choice_lik,affect_lik);
}



