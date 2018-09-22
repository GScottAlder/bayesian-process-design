data {
	// Bypass runs to inform inlet peak distribution
	int<lower=1> n_bypasses;
	real z1b[n_bypasses]; // Isopropyl alcohol bypass peaks
	
	// Experimental trials
	int<lower=1> n_trials;
	real t[n_trials];  // Absolute temperature (K)
	real v[n_trials];  // Volumetric flow rate (m^3/s)
	real z1[n_trials]; // Isopropyl alcohol outlet peaks
	real z2[n_trials]; // Acetone outlet peaks
	real z3[n_trials]; // Propene outlet peaks
	
	// Hyperparameters
	real a1_mu;
	real a1_sd;
	real a2_mu;
	real a2_sd;
	real e1_mu;
	real e1_sd;
	real e2_mu;
	real e2_sd;
	
	// Infinitesimal for lower bounding gamma parameters and denominators
	real h;
}
parameters {
	// Dispersion Parameter (shared)
	real<lower=h> phi;
	
	// Inlet peak distribution
	real<lower=0> z10hat; // Expected value
	real<lower=0> z10;    // Samples
	
	// Arrhenius parameters
	real<lower=0> a1;
	real<lower=0> a2;
	real<lower=0> e1;
	real<lower=0> e2;
	real<lower=-5,upper=5> n1;
	real<lower=-5,upper=5> n2;
	
	real beta2[2];
	real beta3[2];
}
transformed parameters {
	// Thermal correction parameters
	real<lower=h> f2[n_trials]; // Acetone
	real<lower=h> f3[n_trials]; // Propene
	
	// Reaction rate constants
	real k1[n_trials]; // Reaction 1: Acetone
	real k2[n_trials]; // Reaction 2: Propene
	
	// Expected values of peaks
	real<lower=0> z1hat[n_trials]; // Isopropyl alcohol
	real<lower=0> z2hat[n_trials]; // Acetone
	real<lower=0> z3hat[n_trials]; // Propene
	
	// Thermodynamics
	real r_gas;
	r_gas = 8.3144598;
	
	for (i in 1:n_trials) {
		// Linear model of thermal response with temperature
		f2[i] = beta2[1] + beta2[2] * (t[i] - 298);
		f3[i] = beta3[1] + beta3[2] * (t[i] - 298);
		
		// Modified Arrhenius equation for reactions 1 and 2
		k1[i] = a1 * pow(t[i] / 298, n1) * exp(- e1 / r_gas / t[i]);
		k2[i] = a2 * pow(t[i] / 298, n2) * exp(- e2 / r_gas / t[i]);
			
		// Experimental reactor model with analytical solution
		z1hat[i] = z10 * exp( - r_gas * t[i] / v[i] * (k1[i] + k2[i]) * 0.0015);
		z2hat[i] = k1[i] / (k1[i] + k2[i]) * (z10 - z1hat[i]) / f2[i];
		z3hat[i] = k2[i] / (k1[i] + k2[i]) * (z10 - z1hat[i]) / f3[i];
	}
}
model {
	// Priors
	phi ~ gamma(1, 0.2);
	
	z10hat ~ gamma(25, 1.0/3.0);
	
	a1 ~ gamma(a1_mu ^ 2 / a1_sd ^ 2, a1_mu / a1_sd ^ 2);
	a2 ~ gamma(a2_mu ^ 2 / a2_sd ^ 2, a2_mu / a2_sd ^ 2);
	e1 ~ gamma(e1_mu ^ 2 / e1_sd ^ 2, e1_mu / e1_sd ^ 2);
	e2 ~ gamma(e2_mu ^ 2 / e2_sd ^ 2, e2_mu / e2_sd ^ 2);
	n1 ~ normal(1, 1);
	n2 ~ normal(1, 1);
	
	beta2[1] ~ gamma((85/86.0 / 0.05)^2, 85/86.0 / 0.05^2);
	beta3[1] ~ gamma((85/64.5 / 0.05)^2, 85/64.5 / 0.05^2);
	
	// Inlet/bypass peak distribution
	z1b ~ gamma(z10hat^2 / phi + h, z10hat / phi + h); // Observed
	z10 ~ gamma(z10hat^2 / phi + h, z10hat / phi + h); // Sampling 
	
	
    // Likelihoods
	for (i in 1:n_trials) {
		z1[i] ~ gamma(z1hat[i]^2 / phi + h, z1hat[i] / phi + h);
		z2[i] ~ gamma(z2hat[i]^2 / phi + h, z2hat[i] / phi + h);
		z3[i] ~ gamma(z3hat[i]^2 / phi + h, z3hat[i] / phi + h);
	}
}
generated quantities{
	real z1pred[n_trials]; // Predicted isopropyl alcohol peaks
	real z2pred[n_trials]; // Predicted acetone peaks
	real z3pred[n_trials]; // Predicted propene peaks
	
	for (i in 1:n_trials) {
		z1pred[i] = gamma_rng(z1hat[i]^2 / phi + h, z1hat[i] / phi + h);
		z2pred[i] = gamma_rng(z2hat[i]^2 / phi + h, z2hat[i] / phi + h);
		z3pred[i] = gamma_rng(z3hat[i]^2 / phi + h, z3hat[i] / phi + h);
	}
}
