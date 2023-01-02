data {
  int<lower=0> U;       // Number of fixed p
  int<lower=0> N;       // Number of observations
  int<lower = 0, upper = 1> M;
  vector<lower=0>[N] y; // observations
  vector[N] x;          // regressor
  vector[2] p_beta;
  real lb_bound;
  real ub_bound;
}
transformed data {
  vector[N] lb;
  vector[N] ub;
  lb[1 : U] = rep_vector(lb_bound, U);
  lb[U + 1 : N] = rep_vector(0, N - U);
  ub[1 : N - U] = rep_vector(1, N - U);
  ub[N - U + 1 : N] = rep_vector(ub_bound, U);

  vector[N] x_trans = x - min(x);

  //real v_y = variance(y);
}
parameters {
  real<lower = 0, upper = 1> hyper_alpha_raw;
  real<lower = 0> alpha;               // Upper alpha parameter
  //real I_U;                            // Upper intercept
  //real<upper = I_U> I_L;               // Lower intercept
  real<lower = 0, upper = 1> S_U_raw;                 // Upper slope
  real<lower = 0, upper = 1> S_L_raw;                 // Upper slope
  vector<lower = lb, upper = ub>[N] p; // Fuzzy mixture paramater
  vector[2] eta;
}

transformed parameters {
  real<lower = 0> hyper_alpha = -log(hyper_alpha_raw);
  //real<lower = 0> alpha = -log(hyper_alpha_raw)*15;
  real S_U = -log(S_U_raw); // S_U ~ exponential(1)
  real S_L = -log(S_L_raw); // S_U ~ exponential(1)
  real I_L = -1 + eta[1]*.01; // I_L ~ normal(0, .1)
  real<lower = I_L> I_U = 1 + eta[2]*.05;
}

model {
  //hyper_alpha ~ exponential(100);
  alpha       ~ gamma(hyper_alpha, hyper_alpha);    // Prior on alpha
  //I_U         ~ normal(0, .1);             // Prior on lower intercept
  //I_L         ~ normal(-1, .05);
  eta         ~ std_normal();
  //S_L         ~ exponential(10);          // Prior on lower slope
  //S_U         ~ exponential(1);           // Prior on upper slope
  p           ~ beta(p_beta[1], p_beta[2]);     // Prior on mixture
  {
    vector[N] exp_beta;
    if (M) {
      exp_beta = exp(//
      p.*(I_L - S_L*x)//
      +//
      (1-p).*(I_U - S_U*x)//
      );         // exponential log-link
    }
    else {
      exp_beta = exp(//
        I_L - S_L*x//
        +//
        (1-p).*(I_U - S_U*x)//
      );
    }
    y ~ gamma(alpha, alpha ./ exp_beta); // Likelihood
  }
}

// generated quantities {
//   real nrmse;
//   {
//     vector[N] loop_var;
//     vector[N] exp_beta = exp(//
//       I_L - S_L*x//
//       +//
//       (1-p).*(I_U - S_U*x)//
//     );
//     array[N] real y_hat = gamma_rng(alpha, alpha ./ exp_beta);
//     for (i in  1:N) {
//       loop_var[i] = y_hat[i] - y[i];
//     }
//     loop_var = loop_var ^2;
//     nrmse = (sum(loop_var) / N) / v_y;
//   }
//   nrmse = sqrt(nrmse);
// }
