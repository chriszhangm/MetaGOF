//
// This Stan program defines a simple model, with a
// vector of values 'y' modeled as normally distributed
// with mean 'mu' and standard deviation 'sigma'.
//
// Learn more about model development with Stan at:
//
//    http://mc-stan.org/users/interfaces/rstan.html
//    https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started
//

// The input data is a vector 'y' of length 'N'.
data {
    
      int K; /* Number of studies */
      int J; /*J=2*/
      int XC[K];
      int NXT[K];
      int NXC[K];
      int XT[K];
      real xcraw;
      real thetaraw;
  }
  
  transformed data {
    int NC[K];
    int NT[K];
    
    
    for(i in 1:K){
      NC[i] = XC[i] + NXC[i]; // num positive
      NT[i] = XT[i] + NXT[i]; // num negative
    }
  }
  
  parameters{
    vector[J] mu;
    vector[J] theta[K];
    cov_matrix[J] Q;
    vector<lower=0>[J] delta;
    //corr_matrix[J] Lcorr;
    
  }
  
  transformed parameters{
    vector[J] P[K] = inv_logit(theta);
    matrix[J,J] Vprior;
    cov_matrix[J] Sigma; //covariance matrix
    //define Sigma
    Sigma[1,1] = delta[1]*delta[1]*Q[1,1];
    Sigma[1,2] = delta[1]*delta[2]*Q[1,2];
    Sigma[2,2] = delta[2]*delta[2]*Q[2,2];
    Sigma[2,1] = Sigma[1,2];
    Vprior = (J+1)*diag_matrix(rep_vector(1,J));
  }
  
  model {
    
    for (i in 1:K){
      XC[i] ~ binomial(NC[i], P[i,1]);
      XT[i] ~ binomial(NT[i], P[i,2]);
      theta[i] ~ multi_normal(mu, Sigma);
    }
    
    
    //theta ~ multi_normal_cholesky(mu, L);
    
    
   // priors
    //xcraw:initial guess from moment of method.
    //thetaraw: initial guess from the Simple Average method.
    mu[1] ~ uniform(xcraw-5,xcraw+5);
    mu[2] ~ uniform(xcraw+thetaraw-5,xcraw+thetaraw+5);
    Q ~ inv_wishart((J+1),Vprior);
    delta ~ lognormal(0,1);
  }
  generated quantities {
    real<lower=-1,upper=1> rho=Sigma[2,1]/sqrt(Sigma[1,1]*Sigma[2,2]);
    real overeff = mu[2]-mu[1];
    real<lower=0> tauS = Sigma[1,1] + Sigma[2,2] - 2*Sigma[2,1];
  }