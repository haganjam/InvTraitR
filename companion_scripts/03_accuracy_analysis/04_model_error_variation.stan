data{
     int<lower=0> N;
     int<lower=0> G;
     vector[N] abs_error;
     vector[N] hm;
    array[N] int bs;
     vector[N] td;
    array[N] int id;
}
parameters{
     matrix[4,G] Z;
     vector[4] abar;
     cholesky_factor_corr[4] L_Rho;
     vector<lower=0>[4] sigma_id;
     real<lower=0> sigma;
}
transformed parameters{
     vector[G] a_id;
     vector[G] b1_td;
     vector[G] b2_bs;
     vector[G] b3_hm;
     matrix[G,4] v;
    v = (diag_pre_multiply(sigma_id, L_Rho) * Z)';
    b3_hm = abar[4] + v[, 4];
    b2_bs = abar[3] + v[, 3];
    b1_td = abar[2] + v[, 2];
    a_id = abar[1] + v[, 1];
}
model{
     vector[N] mu;
    sigma ~ exponential( 1 );
    sigma_id ~ exponential( 1 );
    L_Rho ~ lkj_corr_cholesky( 4 );
    abar ~ normal( 0 , 1 );
    to_vector( Z ) ~ normal( 0 , 1 );
    for ( i in 1:N ) {
        mu[i] = a_id[id[i]] + b1_td[id[i]] * td[i] + b2_bs[id[i]] * bs[i] + b3_hm[id[i]] * hm[i];
    }
    abs_error ~ lognormal( mu , sigma );
}
generated quantities{
     matrix[4,4] Rho;
    Rho = multiply_lower_tri_self_transpose(L_Rho);
}