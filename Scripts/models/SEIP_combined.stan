functions {

  // SEIP_basic: constant beta and delta
  real[] SEIP_basic(real t, real[] y, real[] params, real[] x_r, int[] x_i) {
    real dydt[4];
    real S = y[1]; real E = y[2]; real I = y[3]; real P = y[4];
    real beta = params[1], alpha = params[2], m = params[3], gamma = params[4], delta = params[5];
    
    dydt[1] = beta * S * I + m * I;
    dydt[2] = beta * S * I - alpha * E;
    dydt[3] = alpha * E - gamma * I;
    dydt[4] = gamma * I - delta * P;
    
    return dydt;
  }
  
    // SEIP_temp_delta: temp-dependent delta
  real[] SEIP_temp_delta(real t, real[] y, real[] params, real[] x_r, int[] x_i) {
    real dydt[4];
    real S=y[1]; real E=y[2]; real I=y[3]; real P=y[4];
    real beta=params[1], alpha=params[2], m=params[3], gamma=params[4], a_delta=params[6], b_delta=params[7];
    real mean_temp=x_r[1], amplitude=x_r[2], phase_shift=x_r[3];
    real temp_t=mean_temp+amplitude*sin(2*pi()*t/365+phase_shift);
    real delta_t=exp(a_delta + b_delta*temp_t);
    real max_delta = 0.231;
    delta_t = fmin(delta_t, max_delta);
    
    dydt[1] = -beta * S * I + m * I;
    dydt[2] = beta * S * I - alpha * E;
    dydt[3] = alpha * E - gamma * I;
    dydt[4] = gamma * I - delta_t * P;
    
    return dydt;
  }

  // SEIP_temp_beta: temp-dependent progression probability * beta
  real[] SEIP_temp_beta(real t, real[] y, real[] params, real[] x_r, int[] x_i) {
    real dydt[4];
    real S=y[1]; real E=y[2]; real I=y[3]; real P=y[4];
    real beta=params[1], alpha=params[2], m=params[3], gamma=params[4], delta=params[5], b1=params[8], b0=params[9];
    real mean_temp=x_r[1], amplitude=x_r[2], phase_shift=x_r[3];
    real temp_t=mean_temp+amplitude*sin(2*pi()*t/365+phase_shift);
    real beta_prob=inv_logit(b0 + b1*temp_t);
    real eff_beta=beta * beta_prob;
    
    dydt[1] = -eff_beta * S * I + m * I;
    dydt[2] = eff_beta * S * I - alpha * E;
    dydt[3] = alpha * E - gamma * I;
    dydt[4] = gamma * I - delta * P;
    
    return dydt;
  }
  
    // SEIP_temp_beta: temp-dependent delta and beta
  real[] SEIP_temp_both(real t, real[] y, real[] params, real[] x_r, int[] x_i) {
    real dydt[4];
    real S=y[1]; real E=y[2]; real I=y[3]; real P=y[4];
    real beta=params[1], alpha=params[2], m=params[3], gamma=params[4], delta=params[5], b1=params[8], b0=params[9];
    real mean_temp=x_r[1], amplitude=x_r[2], phase_shift=x_r[3];
    real temp_t=mean_temp+amplitude*sin(2*pi()*t/365+phase_shift);
    real beta_prob=inv_logit(b0 + b1*temp_t);
    real eff_beta=beta * beta_prob;
    real delta_t=exp(a_delta + b_delta*temp_t);
    real max_delta = 0.231;
    delta_t = fmin(delta_t, max_delta)
    
    dydt[1] = -eff_beta * S * I + m * I;
    dydt[2] = eff_beta * S * I - alpha * E;
    dydt[3] = alpha * E - gamma * I;
    dydt[4] = gamma * I - delta_t * P;
    
    return dydt;
  }
}

}

data {
  int<lower = 1> n_obs;    // Number of days sampled
  int<lower = 1> n_difeq;  // Number of differential equations in the system
  int<lower = 1> n_fake;   // This is to generate "predicted"/"unsampled" data
  
  real<lower=0> min_ti;   // minimum time to infectious
  real<lower=0> alpha_m;  // death rate gamma shape
  real<lower=0> beta_m;   // death rate gamma rate
  real<lower=0> delta_est;
  real a_delta_mean;            // intercept
  real b_delta_mean;            // effect of temperature
  real<lower=0> shape_gamma;
  real<lower=0> rate_gamma;
  real b_temp_slope;
  real b_temp_inter;

  real mean_temp;
  real amplitude;
  real phase_shift;

  int<lower=0> y_obs[n_obs, 3];
  real y_edna[n_obs];

  real t0;        // Initial time point (zero)
  real ts[n_obs]; // Time points that were sampled
  
  real fake_ts[n_fake]; // Time points for "predicted"/"unsampled" data
}

transformed data {
  real x_r[3];
  int x_i[0];
  
  x_r[1] = mean_temp;
  x_r[2] = amplitude;
  x_r[3] = phase_shift;
}

parameters {
  real<lower=0> beta;  // rate of S->E

  // TODO: make this minimum value more flexible!
  real<lower=0, upper=1/min_ti> alpha; // T_infect > min_ti => T_infect = 1/alpha => alpha < 1/min_ti

  
  real<lower=0> m;     // rate of death from I
  
  real<lower=0> gamma; // rate of emitted particles from I (copy per ml / I / day)
  
  real a_delta;      // parameters for temp dependent decay rate (delta)
  real<upper=0> b_delta;

  real<lower=0> b1;      // parameters for temp dependent progression probability
  real b0;
  
  real<lower=0> delta;

  real<lower=0> sigma_edna;
  
  simplex[3] y0;  // Use a simplex for the initial state so that the compartments are proportions summing to 1.
  
  real<lower=0> P0;
}

transformed parameters{
  real params[8];
  params[1] = beta;
  params[2] = alpha;
  params[3] = m;
  params[4] = gamma;
  params[5] = delta;
  params[6] = a_delta;
  params[7] = b_delta;
  params[8] = b1;
  params[9] = b0;

  real y_init[4];
  y_init[1] = y0[1];
  y_init[2] = y0[2];
  y_init[3] = y0[3];
  y_init[4] = P0;
  
  real y_hat[n_obs, n_difeq]; // Output from the ODE solver for all four compartments
  { // creating temporary local environment
    real temp_y0[4];
    for (i in 1:4) {
      temp_y0[i] = y_init[i];
    }
    if (model_choice == 1) {
      y_hat = integrate_ode_rk45(SEIP_basic, temp_y0, t0, ts, params, x_r, x_i);
    } else if (model_choice == 2) {
      y_hat = integrate_ode_rk45(SEIP_temp_delta, temp_y0, t0, ts, params, x_r, x_i);
    } else if (model_choice == 3) {
      y_hat = integrate_ode_rk45(SEIP_temp_beta, temp_y0, t0, ts, params, x_r, x_i);
    } else if (model_choice == 4) {
      y_hat = integrate_ode_rk45(SEIP_temp_both, temp_y0, t0, ts, params, x_r, x_i);
    }
  }
}

// priors: P(params)
// sampling: P(obs|params) / likelihood function: L(params) given fixed data (obs))
// STAN uses gamma(shape, rate)
model {
  
  // alway relevant parameters
  beta ~ lognormal(log(0.02), 0.4); 
  // weakly informed, inferred from STAN
  // centered at a plausible value, but allows for a wide range, avoiding extreme tails
  // enforces positivity and infection rates are often right-skewed
  // lognormal over gamma because it allows heavier upper tail without excessive skew

  alpha ~ lognormal(log(0.01), 0.3);
  // Weakly informed, but there is an idea of what it should be
  // upper limit has been implemented (80.6)
  // I am thinking the time to infectious is probably around that limit (maybe around 90)
  // high spread to reflect high uncertainty
  
  m ~ inv_gamma(alpha_m, beta_m); // time to death rate, informed from lab data, 2.603128, 0.04443979
  // time to death ~ gamma(m_alpha, scale = 1/m_beta) => m = 1/time => m ~ inv_gamma(m_alpha, scale = m_beta)
  
  gamma ~ gamma(shape_gamma, rate_gamma); // I want this to be tighter than other variables are because the observed data is very weak
  
  sigma_edna ~ exponential(1);
  
  P0 ~ lognormal(log(0.005), 0.2);
  
  if (model_choice == 1 || model_choice == 3) {
    delta ~ lognormal(log(delta_est), 0.1)
  }
  
  if (model_choice == 2 || model_choice == 4) {
    a_delta ~ normal(a_delta_mean, 0.4);
    b_delta ~ normal(b_delta_mean, 0.05);
  }

  if (model_choice == 3 || model_choice == 4) {
    b0 ~ normal(b_temp_inter, 0.25);
    b1 ~ normal(b_temp_slope, 0.05) T[0, ];
  }
  
  // sampling distribution
  // likelihood: observe only the SEI compartments
  // y_t ~ Multinomial(n_t, proportions from ODE: [S,E,I]_t), where n_t = \sum (y[t]t)
  // comparing observed counts to model-predicted proportions
  for (t in 1:n_obs) {
    y_obs[t, ] ~ multinomial(to_vector(y_hat[t, 1:3]));
    
    // 0.1 is the LOD so values of 0.1 are left=-censored
    // using cdf at LOD, prob of being less than or equal to 0.1
    if (y_edna[t] == 0.1) {
      target += normal_lcdf(0.1 | y_hat[t, 4], sigma_edna); //+ 1e-6); // to remove instability
      // Add the log-probability that the observed eDNA value was less than or equal to 0.1, 
      // assuming it came from a normal distribution centered at y_hat[t,4]
    } else {
      y_edna[t] ~ normal(y_hat[t, 4], sigma_edna);
    }

  }
}

// Generate predicted data over a specified time frame iwth final parameter estimates:
generated quantities {
  real fake[n_fake, n_difeq];
    {
      real temp_y0[4];
      for (i in 1:4) {
        temp_y0[i] = y_init[i];
      }
      if (model_choice == 1) {
        fake = integrate_ode_rk45(SEIP_basic, temp_y0, t0, fake_ts, params, x_r, x_i);
      } else if (model_choice == 2) {
        fake = integrate_ode_rk45(SEIP_temp_delta, temp_y0, t0, fake_ts, params, x_r, x_i);
      } else if (model_choice == 3) {
        fake = integrate_ode_rk45(SEIP_temp_beta, temp_y0, t0, fake_ts, params, x_r, x_i);
      } else if (model_choice == 4) {
        fake = integrate_ode_rk45(SEIP_temp_both, temp_y0, t0, fake_ts, params, x_r, x_i);
      }
    }
}
