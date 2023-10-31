data {
  int<lower=1> n_t; // number of trials
  int<lower=1> n_s; // number of subjects
  int<lower=1> n_f; //number of fractals
  array[n_s,n_t] int fA; // 2d array - one array for each subject, containing a sub-arry with the identity of
                          // fractal A on each of that subject's trials (fractals are identified with numbers 1:n_f)
  array[n_s,n_t] int fB; // diddo for fractal B identities
  array[n_s] vector[n_t] bv; 
  array[n_s,n_t] int fres; //2d array (subjects by trials) indicating whether the fractal result was shown on that trial (1-0)
  array[n_s,n_t] int choice; // 2d array - one integer array for each subject - containing subject's choices - coded
                             // 1 for fractal A and 2 for fractal B. 
}

parameters{
  real lsbias_mu;
  real<lower=0> lsbias_sigma;
  vector[n_s] lsbias_z;
  
  real bvf_sens_mu; //mean effects of reward (fres) on affect and choice
  real<lower=0> bvf_sens_sigma;
  vector[n_s] bvf_sens_z;

  //decay rates
  real dcy_mu;
  real<lower=0> dcy_sigma;
  vector[n_s] dcy_z;
  
  real lrn_b_mu;
  real<lower=0> lrn_b_sigma;
  vector[n_s] lrn_b_z;
}

transformed parameters {
  vector[n_s*n_t] log_lik;
  {//anonymous_scope_start
    array[n_s,n_t] vector[n_f] B; 
    vector[2] softmax_arg; //goes inside choice likelihood calculation
    
    //NCP for normal priors
    vector[n_s] lsbias = lsbias_mu + lsbias_sigma*lsbias_z;
    vector[n_s] bvf_sens = bvf_sens_mu + bvf_sens_sigma*bvf_sens_z;
    
    //NCP for normal priors with sigmoid transform
    vector[n_s] dcy = inv_logit(dcy_mu + dcy_sigma*dcy_z);
    
    vector[n_s] lrn_b = inv_logit(lrn_b_mu + lrn_b_sigma*lrn_b_z);

    // loop through each subject and trial...
    for (s in 1:n_s) {
      for (t in 1:n_t){
        if(t == 1){
          B[s,t] = rep_vector(0,n_f); 
        }
        
        //set softmax arg outside of the log_lik calculation for efficiency (you should avoid transposing inside a likelihood funciton, apparently)
        softmax_arg = bvf_sens[s]*B[s,t,{fA[s,t],fB[s,t]}] + [lsbias[s],0]';

        log_lik[(s-1)*n_t+t] = categorical_logit_lpmf(choice[s,t] | softmax_arg); //get the likelihood of the choice on this trial. Doing the log lik
        
        // unlesss this is the subject's last trial, set the Q,C, and A values for the next trial
        if(t < n_t){
          //decay all values toward 0
          B[s,t+1] = (1-dcy[s])*B[s,t];
          if(choice[s,t] == 1){
            // if the outcome was the fractal result...
            if(fres[s,t] == 1){
              //update Q...
              B[s,t+1,fA[s,t]] = B[s,t,fA[s,t]] + lrn_b[s]*(bv[s,t] - B[s,t,fA[s,t]]);
            }
          } else if(choice[s,t] == 2){
            if(fres[s,t] == 1){
              //update Q...
              B[s,t+1,fB[s,t]] = B[s,t,fB[s,t]] + lrn_b[s]*(bv[s,t] - B[s,t,fB[s,t]]);
            }
          }
        }
      }
    }
  }//anonymous_scope_end
}
model{
  lsbias_mu ~ normal(0,4);
  lsbias_sigma ~ normal(0,8);
  //effects on choice
  bvf_sens_mu ~ normal(0,4);
  bvf_sens_sigma ~ normal(0,6);
  
  //learning and decay rates
  dcy_mu ~ normal(-.05,1.7);
  dcy_sigma ~ normal(0,4);
  
  lrn_b_mu ~ normal(-.05,1.7);
  lrn_b_sigma ~ normal(0,4);
  
  //zs
  lsbias_z ~ std_normal();
  bvf_sens_z ~ std_normal();
  dcy_z ~ std_normal();
  lrn_b_z ~ std_normal();
  
  target += log_lik;
}


