data{
     vector[1523] abs_error;
     vector[1523] hm;
    array[1523] int bs;
     vector[1523] td;
    array[1523] int id;
}
parameters{
     matrix[4,154] Z;
     vector[4] abar;
     cholesky_factor_corr[4] L_Rho;
     vector<lower=0>[4] sigma_id;
     real<lower=0> sigma;
}
transformed parameters{
     vector[154] a_id;
     vector[154] b1_td;
     vector[154] b2_bs;
     vector[154] b3_hm;
     matrix[154,4] v;
    v = (diag_pre_multiply(sigma_id, L_Rho) * Z)';
    b3_hm = abar[4] + v[, 4];
    b2_bs = abar[3] + v[, 3];
    b1_td = abar[2] + v[, 2];
    a_id = abar[1] + v[, 1];
}
model{
     vector[1523] mu;
    sigma ~ exponential( 1 );
    sigma_id ~ exponential( 1 );
    L_Rho ~ lkj_corr_cholesky( 4 );
    abar ~ normal( 0 , 1 );
    to_vector( Z ) ~ normal( 0 , 1 );
    for ( i in 1:1523 ) {
        mu[i] = a_id[id[i]] + b1_td[id[i]] * td[i] + b2_bs[id[i]] * bs[i] + b3_hm[id[i]] * hm[i];
    }
    abs_error ~ lognormal( mu , sigma );
}
generated quantities{
     matrix[4,4] Rho;
    Rho = multiply_lower_tri_self_transpose(L_Rho);
}