//This file contains stan code for hierarchical Bayesian prospect theory parameter estimation

data {
  int<lower=1> N;
  int<lower=1> T;
  int<lower=1, upper=T> Tsubj[N];
  real<lower=0> outcome1[N, T];
  real<lower=0> outcome2[N, T];
  real<lower=0> outcome3[N, T];
  real<lower=0> outcome4[N, T];
  real<lower=0, upper=1> probs[N, T]; //probability of best outcome
  int<lower=0, upper=1> choice[N, T];
  int<lower=0, upper=2> session[N,T];
  
  
}

parameters {
  // Hyper(group)-parameters
  vector[9] mu_pr;
  // subject-level raw parameters
  vector[N] rho_pr;
  vector[N] tau_pr;
  vector[N] prw_pr;
  vector[N] delta1_rho_pr;
  vector[N] delta1_tau_pr;
  vector[N] delta1_prw_pr;
  vector[N] delta2_rho_pr;
  vector[N] delta2_tau_pr;
  vector[N] delta2_prw_pr;
  vector<lower=0, upper=5>[9] sigma; 
  
}

transformed parameters {
   vector[N] rho0;
   vector[N] tau0;
   vector[N] prw0;
   vector[N] delta1_rho;
   vector[N] delta1_tau;
   vector[N] delta1_prw;
   vector[N] delta2_rho;
   vector[N] delta2_tau;
   vector[N] delta2_prw;

    for (i in 1:N){
      rho0[i] = 0+5*Phi_approx(mu_pr[1]+sigma[1]*rho_pr[i]);
      tau0[i]  = 0+10*Phi_approx(mu_pr[2]+sigma[2]*tau_pr[i]);
      prw0[i] = 0+5*Phi_approx(mu_pr[3]+sigma[3]*prw_pr[i]);
      
      delta1_rho[i] = -0.5+1*Phi_approx(mu_pr[4]+sigma[4]*delta1_rho_pr[i]);
      delta1_tau[i] = -5+10*Phi_approx(mu_pr[5]+sigma[5]*delta1_tau_pr[i]);
      delta1_prw[i] = -0.5+1*Phi_approx(mu_pr[6]+sigma[6]*delta1_prw_pr[i]);
      
      delta2_rho[i] = -0.5+1*Phi_approx(mu_pr[7]+sigma[7]*delta2_rho_pr[i]);
      delta2_tau[i] =-5+10*Phi_approx(mu_pr[8]+sigma[8]*delta2_tau_pr[i]);
      delta2_prw[i] = -0.5+1*Phi_approx(mu_pr[9]+sigma[9]*delta2_prw_pr[i]);
      
    }
 
}

model {
  // Hyper(group)-parameters
  mu_pr[1] ~normal(0, 1.0);
  mu_pr[2] ~normal(0, 1.0);
  mu_pr[3] ~normal(0, 1.0);
  mu_pr[4:9]  ~ normal(0, 1.0);
  sigma ~ uniform(0, 5);	  

  // Subject-level raw parameters (for Matt trick)
   rho_pr    ~ normal(0,1);
   tau_pr    ~ normal(0,1);
   prw_pr   ~ normal(0,1);
   delta1_rho_pr ~ normal(0,1);
   delta1_tau_pr ~ normal(0,1);
   delta1_prw_pr ~ normal(0,1);
   delta2_rho_pr ~ normal(0,1);
   delta2_tau_pr ~ normal(0,1);
   delta2_prw_pr ~ normal(0,1);

  for (i in 1:N) {
    for (t in 1:Tsubj[i]) {
   
        real EU_a_i;
        real EU_b_i;
        real pChooseA; //probability to choose lottery A over lottery B
        real prw;
        real rho;
        real tau;
        
        if (session[i,t]==0) {
        rho = rho0[i];
        tau = tau0[i];
        prw = prw0[i];
        
      } else if (session[i,t]==1){
        rho = rho0[i]+delta1_rho[i];
        tau = tau0[i]+delta1_tau[i];
        prw = prw0[i]+delta1_prw[i];
        
      } else {
        rho = rho0[i]+delta2_rho[i];
        tau = tau0[i]+delta2_tau[i];
        prw = prw0[i]+delta2_prw[i];
      }

      EU_a_i = ((probs[i, t]^prw)/(((probs[i, t])^prw+(1-probs[i, t])^prw)^(1/prw))) * ((outcome1[i, t])^rho)+(1 - ((probs[i, t]^prw)/(((probs[i, t])^prw+(1-probs[i,t])^prw)^(1/prw)))) *((outcome2[i, t])^rho);
        
      EU_b_i = ((probs[i, t]^prw)/(((probs[i, t])^prw+(1-probs[i, t])^prw)^(1/prw))) *((outcome3[i, t])^rho)+ (1 - ((probs[i, t]^prw)/(((probs[i, t])^prw+(1-probs[i, t])^prw)^(1/prw)))) * ((outcome4[i, t])^rho);
      
         pChooseA  = tau * (1/rho)*(log(EU_a_i)-log(EU_b_i));

         choice[i, t] ~ bernoulli_logit(pChooseA);

    }
  }
}
generated quantities {
  real<lower=0, upper=5> mu_rho;
  real<lower=0, upper=10> mu_tau;
  real<lower=0, upper=5> mu_prw;
  
   real<lower=-0.5, upper=0.5> mu_delta1_rho;
   real<lower=-5, upper=5> mu_delta1_tau;
   real<lower=-0.5, upper=0.5> mu_delta1_prw;
  
  real<lower=-0.5, upper=0.5> mu_delta2_rho;
  real<lower=-5, upper=5> mu_delta2_tau;
  real<lower=-0.5, upper=0.5> mu_delta2_prw;
  int j_t;

  real log_lik[N*T];

  // For posterior predictive check
  real y_pred[N, T];

  // Set all posterior predictions to 0 (avoids NULL values)
  for (i in 1:N) {
    for (t in 1:T) {
      y_pred[i, t] = -1;
    }
  }


  mu_rho    = Phi_approx(mu_pr[1]) * 5;
  mu_tau    = Phi_approx(mu_pr[2]) * 10;
  mu_prw    = Phi_approx(mu_pr[3]) * 5;
  
  mu_delta1_rho = Phi_approx(mu_pr[4]) * (1) -0.5 ;
  mu_delta1_tau = Phi_approx(mu_pr[5]) * (10) -5 ;
  mu_delta1_prw = Phi_approx(mu_pr[6]) * (1) -0.5 ;
  mu_delta2_rho = Phi_approx(mu_pr[7]) * (1) -0.5 ;
  mu_delta2_tau = Phi_approx(mu_pr[8]) * (10) -5 ;
  mu_delta2_prw = Phi_approx(mu_pr[9]) * (1) -0.5 ;
  
  j_t=1;
  { // local section, this saves time and space
    for (i in 1:N) {
      //log_lik[i] = 0;
      for (t in 1:Tsubj[i]) {
        real EU_a_i;
        real EU_b_i;
        real pChooseA; //probability to choose lottery A over lottery B
        real prw;
        real rho;
        real tau;
     
        if (session[i,t]==0) {
        rho = rho0[i];
        tau = tau0[i];
        prw = prw0[i];
        
      } else if (session[i,t]==1){
        rho = rho0[i]+delta1_rho[i];
        tau = tau0[i]+delta1_tau[i];
        prw = prw0[i]+delta1_prw[i];
        
      } else {
        rho = rho0[i]+delta2_rho[i];
        tau = tau0[i]+delta2_tau[i];
        prw = prw0[i]+delta2_prw[i];
      }
     
      
      EU_a_i = ((probs[i, t]^prw)/(((probs[i, t])^prw+(1-probs[i, t])^prw)^(1/prw)))* ((outcome1[i, t])^rho) + (1 - ((probs[i, t]^prw)/(((probs[i, t])^prw+(1-probs[i, t])^prw)^(1/prw)))) *((outcome2[i, t])^rho);
        
      EU_b_i = ((probs[i, t]^prw)/(((probs[i, t])^prw+(1-probs[i, t])^prw)^(1/prw))) *((outcome3[i, t])^rho)+ (1 - ((probs[i, t]^prw)/(((probs[i, t])^prw+(1-probs[i, t])^prw)^(1/prw)))) * ((outcome4[i, t])^rho);
	    pChooseA  = tau * (1/rho)*(log(EU_a_i)-log(EU_b_i));       
	

        log_lik[j_t] = bernoulli_logit_lpmf(choice[i, t] | pChooseA);
        j_t=j_t+1;
        
        // generate posterior prediction for current trial
        y_pred[i, t] = bernoulli_logit_rng(pChooseA);
      }
    }
  }
}


