data {
	int<lower=1> n_bypasses;
	real z0_a[n_bypasses];
}
parameters {
	real<lower=0.001> phi;
	real<lower=1> mu_z0_a;
	real<lower=0.001> f_a;
}
model {
	phi ~ gamma(1, 0.2);                         // Prior
	z0_a ~ gamma(mu_z0_a ^ 2 / phi, mu_z0_a / phi); // Likelihood
	f_a ~ gamma( (85 / 0.05) ^ 2, 85 / 0.05 ^ 2);
}
generated quantities {
	real z0pred_a;
	real t_inj;
	
	z0pred_a = gamma_rng(mu_z0_a ^ 2 / phi, mu_z0_a / phi); // Sampling Distribution
	t_inj = z0pred_a * 8.3144598 * 298 / (f_a * 5818.87 * 1.666667e-06);
}
