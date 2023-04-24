data{
     vector[2005] abs_error;
     vector[2005] hm;
     int bs[2005];
     vector[2005] td;
     int id[2005];
}
parameters{
     vector[160] b3_hm;
     vector[160] b2_bs;
     vector[160] b1_td;
     vector[160] a_id;
     real a;
     real b1;
     real b2;
     real b3;
     vector<lower=0>[4] sigma_id;
     real<lower=0> sigma;
     corr_matrix[4] Rho;
}
model{
     vector[2005] mu;
    Rho ~ lkj_corr( 4 );
    sigma ~ exponential( 1 );
    sigma_id ~ exponential( 1 );
    b3 ~ normal( 0 , 2 );
    b2 ~ normal( 0 , 2 );
    b1 ~ normal( 0 , 2 );
    a ~ normal( 0 , 2 );
    {
    real YY[160, 4];
    vector[4] MU;
    MU = [ a , b1 , b2 , b3 ]';
    for ( j in 1:160 ) YY[j] = [ a_id[j] , b1_td[j] , b2_bs[j] , b3_hm[j] ]';
    YY ~ multi_normal( MU , quad_form_diag(Rho , sigma_id) );
    }
    for ( i in 1:2005 ) {
        mu[i] = a_id[id[i]] + b1_td[id[i]] * td[i] + b2_bs[id[i]] * bs[i] + b3_hm[id[i]] * hm[i];
    }
    abs_error ~ normal( mu , sigma );
}
