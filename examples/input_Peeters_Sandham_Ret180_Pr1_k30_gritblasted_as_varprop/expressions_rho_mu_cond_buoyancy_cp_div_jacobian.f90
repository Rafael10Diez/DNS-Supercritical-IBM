enth_arr(i,j,k)         = T_ext(i,j,k)
rho_arr(i,j,k)          = 1.
mu_arr(i,j,k)           = 0.005555555555555556d0 ! 1./180.
cond_arr(i,j,k)         = 0.005555555555555556d0 ! 1./180.
buoyancy_arr(i,j,k)     = 0.d0
cp_arr(i,j,k)           = 1.
div_jacobian_arr(i,j,k) = 0.

T_ext_ijk  =  T_ext(i,j,k)
enth_ijk   =  enth_arr(i,j,k)
rho_ijk    =  rho_arr(i,j,k)
mu_ijk     =  mu_arr(i,j,k)
cond_ijk   =  cond_arr(i,j,k)
cp_ijk     =  cp_arr(i,j,k)
