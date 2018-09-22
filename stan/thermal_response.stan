data {
	int<lower=1> n_bypasses;
	vector[n_bypasses] peaks;
	vector[3] hsd_intercept;
	vector[3] hsd_slope;
}

transformed data {
	vector[3] f = [85.0, 86.0, 64.5]'; // Literature Values
	vector[n_bypasses] scale = peaks ./ 60.0 ./ f[1];	
}

parameters {
	real<lower=0> sigma;
	
	real<lower=0> scale_factor;
	
	vector[3] alpha;
	vector[3] beta;
	
	
}

transformed parameters {
	matrix[16, 3] therm_corr_factor;
	for (i in 1:16) {
		for (c in 1:3) {
			// 298 to 673 K centered at 25 C
			therm_corr_factor[i, c] = alpha[c] * scale_factor + beta[c] * 25 * (i - 1);
		}
	}
 	
}

model {
	phi ~ gamma(2, 1.0 / 25.0);  // Prior
	scale_factor ~ normal(1, 0.5); // Prior
	
	scale ~ normal(scale_factor, sqrt(phi); // Likelihood
	
	alpha ~ normal(f, hsd_intercept); // Prior	
	 beta ~ normal(0, hsd_slope);     // Prior
} 
