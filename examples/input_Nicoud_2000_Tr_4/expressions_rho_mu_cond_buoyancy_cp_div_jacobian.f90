enth_arr(i,j,k)         = T_ext(i,j,k)
rho_arr(i,j,k)          = 1.d0/                           T_ext(i,j,k)
mu_arr(i,j,k)           = 1.d0/(180.0_dp)         / dsqrt(T_ext(i,j,k))
cond_arr(i,j,k)         = 1.d0/(180.0_dp*0.76_dp) / dsqrt(T_ext(i,j,k))
buoyancy_arr(i,j,k)     = 0.d0
cp_arr(i,j,k)           = 1.d0
div_jacobian_arr(i,j,k) = 1.d0

T_ext_ijk  =  T_ext(i,j,k)
enth_ijk   =  enth_arr(i,j,k)
rho_ijk    =  rho_arr(i,j,k)
mu_ijk     =  mu_arr(i,j,k)
cond_ijk   =  cond_arr(i,j,k)
cp_ijk     =  cp_arr(i,j,k)
