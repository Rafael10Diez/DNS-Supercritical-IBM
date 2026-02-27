module mod_cfd_ibm
  use  openacc
  use  mpi
  use, intrinsic :: iso_c_binding, only: C_INT, c_intptr_t, C_PTR, C_LOC
  use, intrinsic :: iso_fortran_env, only: i8 => int64, dp => real64
  use  mod_utils
  use  diezdecomp_api_generic
  implicit none
  contains
  subroutine ibm_part_1(obj_ibm_U, obj_ibm_V, obj_ibm_W, obj_ibm_T, obj_ibm_div, U, V, W, T, div_now, buffer_ibm, &
                        wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                        wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top, &
                        mpi_pos_y    , mpi_divs_y   , use_incompressible)
    use diezDecomp_api_ibm, only:  diezDecomp_ibm_type, diezdecomp_ibm_exec
    implicit none
    type(diezDecomp_ibm_type)  ::  obj_ibm_U, obj_ibm_V, obj_ibm_W, obj_ibm_T, obj_ibm_div
    real(dp)  ::  U(-1:,-1:,-1:), &
                  V(-1:,-1:,-1:), &
                  W(-1:,-1:,-1:), &
                  T(-1:,-1:,-1:), &
                  div_now(-1:,-1:,-1:), &
                  buffer_ibm(0:)   , &
                  wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                  wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top
    integer   ::  mpi_pos_y, mpi_divs_y
    logical   ::  use_incompressible
    call diezdecomp_ibm_exec(obj_ibm_U, U, buffer_ibm)
    call diezdecomp_ibm_exec(obj_ibm_V, V, buffer_ibm)
    call diezdecomp_ibm_exec(obj_ibm_W, W, buffer_ibm)
    call diezdecomp_ibm_exec(obj_ibm_T, T, buffer_ibm)
 
    if (mpi_pos_y ==  0            )  call set_val_slice(U, 0,  0, wall_BC_U_bot)
    if (mpi_pos_y ==  0            )  call set_val_slice(V, 0,  0, wall_BC_V_bot)
    if (mpi_pos_y ==  0            )  call set_val_slice(W, 0,  0, wall_BC_W_bot)
    if (mpi_pos_y ==  0            )  call set_val_slice(T, 0,  0, wall_BC_T_bot)
 
    if (mpi_pos_y == (mpi_divs_y-1))  call set_val_slice(U,-1, -1, wall_BC_U_top)
    if (mpi_pos_y == (mpi_divs_y-1))  call set_val_slice(V,-1, -1, wall_BC_V_top)
    if (mpi_pos_y == (mpi_divs_y-1))  call set_val_slice(W,-1, -1, wall_BC_W_top)
    if (mpi_pos_y == (mpi_divs_y-1))  call set_val_slice(T,-1, -1, wall_BC_T_top)
 
    if (.not.use_incompressible) then 
      call diezdecomp_ibm_exec(obj_ibm_div, div_now, buffer_ibm) ! div uses extrap coeffs
    end if
  end subroutine

  subroutine ibm_part_2(U_hat, V_hat, W_hat, T, &
                        rough_U_bot_inds, rough_U_top_inds, rough_W_bot_inds, rough_W_top_inds , &
                        rough_V_bot_inds, rough_V_top_inds, rough_T_bot_inds, rough_T_top_inds , &
                        wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                        wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top, &
                        nx_glob       , nz_loc       , ny_loc                      )
    real(dp)  ::  U_hat(-1:,-1:,-1:), &
                  V_hat(-1:,-1:,-1:), &
                  W_hat(-1:,-1:,-1:), &
                  T(-1:,-1:,-1:), &
                  wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                  wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top
    integer   ::  nx_glob, nz_loc, ny_loc, &
                  rough_U_bot_inds(0:,0:), rough_U_top_inds(0:,0:), rough_W_bot_inds(0:,0:), rough_W_top_inds(0:,0:), &
                  rough_V_bot_inds(0:,0:), rough_V_top_inds(0:,0:), rough_T_bot_inds(0:,0:), rough_T_top_inds(0:,0:)
    block 
     integer :: i,j,k
    
     ! ------------------------------ U_hat ------------------------------
     !$acc parallel loop collapse(2) default(present)
     do    j  = 0, nz_loc-1
     do    i  = 0, nx_glob-1
       do  k  = 0, rough_U_bot_inds(i,j), 1
         U_hat(i,j,k) = wall_BC_U_bot
       end do
     end do
     end do
    
     !$acc parallel loop collapse(2) default(present)
     do    j  = 0, nz_loc-1
     do    i  = 0, nx_glob-1
        do k  = rough_U_top_inds(i,j), ny_loc-1, 1
         U_hat(i,j,k) = wall_BC_U_top
        end do
     end do
     end do
    
     ! ------------------------------ W_hat ------------------------------
     !$acc parallel loop collapse(2) default(present)
     do    j  = 0, nz_loc-1
     do    i  = 0, nx_glob-1
        do k  = 0, rough_W_bot_inds(i,j), 1
         W_hat(i,j,k) = wall_BC_W_bot
        end do
     end do
     end do
    
     !$acc parallel loop collapse(2) default(present)
     do    j  = 0, nz_loc-1
     do    i  = 0, nx_glob-1
        do k  = rough_W_top_inds(i,j), ny_loc-1, 1
         W_hat(i,j,k) = wall_BC_W_top
        end do
     end do
     end do
    
     ! ------------------------------ V_hat ------------------------------
     !$acc parallel loop collapse(2) default(present)
     do    j  = 0, nz_loc-1
     do    i  = 0, nx_glob-1
        do k  = 0, rough_V_bot_inds(i,j), 1
         V_hat(i,j,k) = wall_BC_V_bot
        end do
     end do
     end do
    
     !$acc parallel loop collapse(2) default(present)
     do    j  = 0, nz_loc-1
     do    i  = 0, nx_glob-1
        do k  = rough_V_top_inds(i,j), ny_loc, 1
         V_hat(i,j,k) = wall_BC_V_top
        end do
     end do
     end do
    
     ! ------------------------------ T ------------------------------
     !$acc parallel loop collapse(2) default(present)
     do    j  = 0, nz_loc-1
     do    i  = 0, nx_glob-1
        do k  = 0, rough_T_bot_inds(i,j), 1
         T(i,j,k) = wall_BC_T_bot
        end do
     end do
     end do
    
     !$acc parallel loop collapse(2) default(present)
     do    j  = 0, nz_loc-1
     do    i  = 0, nx_glob-1
        do k  = rough_T_top_inds(i,j), ny_loc-1, 1
         T(i,j,k) = wall_BC_T_top
        end do
     end do
     end do
    end block
  end subroutine
end module