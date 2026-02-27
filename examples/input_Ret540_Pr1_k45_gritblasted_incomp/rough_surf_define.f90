module Mod_Height_Function
  use, intrinsic :: iso_fortran_env, only: i8 => int64, dp => real64
  implicit none

  ! --------------- public module API ---------------
  private
  public :: init_surf_obj, type_surf_obj, get_h

  ! --------------- surface object ---------------
  type type_surf_obj
    real(dp)               :: Lx, Lz, Ly
    integer              :: nreps_x, nreps_z
    ! additional private variables (for this module only)
    integer              :: private_L_mnab_coeffs
    real(dp), allocatable  :: private_mnab_coeffs(:,:)
  end type

  contains

  subroutine init_surf_obj(surf_obj)
    implicit none    
    type(type_surf_obj) :: surf_obj
    real(dp)              :: Lx, Lz, Ly, nu, cond, Ret, Pr, dtmax, Sf_x, Sf_z, Sq, & ! declared for compatibility with params.dat
                           Gb_x, Gb_z, Gb_y, bulk_target_Ub, bulk_target_Reb, &
                           wall_BC_T_bot, wall_BC_T_top, wall_BC_U_bot, wall_BC_U_top, wall_BC_V_bot, wall_BC_V_top, &
                           wall_BC_W_bot, wall_BC_W_top
    integer             :: nreps_x, nreps_z, nx, nz, ny, bulk_nprint, avg_iter_start, avg_nraw, avg_freq, avg_bin_nbins, &
                           snap_iter_start, snap_freq, filtering_steps, filtering_first, filtering_nstop, use_rough_surf, &
                           use_rough_surf_quady, &
                           bulk_use_target_Ub, bulk_use_target_Reb, bulk_use_rho_weighting, force_Sf_x_semilocal, &
                           force_Sf_x_pressure, force_Sf_x_accel, use_incompressible, use_utau_work, prof_1d_avg_use
     
    include "params.dat"
  
    surf_obj%Lx       =  2.815_dp
    surf_obj%Lz       =  1.4075_dp
    surf_obj%Ly       =  Ly
    surf_obj%nreps_x  =  nreps_x
    surf_obj%nreps_z  =  nreps_z
  
    open(19,file='./fortran_surf_files/mnab_coeffs.dat')
      read(19,*) surf_obj%private_L_mnab_coeffs
      allocate(surf_obj%private_mnab_coeffs(0:surf_obj%private_L_mnab_coeffs-1,0:3))
      block 
        integer :: i 
        do i=0,surf_obj%private_L_mnab_coeffs-1
          read(19,*) surf_obj%private_mnab_coeffs(i,:)
          surf_obj%private_mnab_coeffs(i,2) = 0.5_dp * surf_obj%private_mnab_coeffs(i,2)
          surf_obj%private_mnab_coeffs(i,3) = 0.5_dp * surf_obj%private_mnab_coeffs(i,3)
        end do
      end block
    close(19)
  end subroutine 

  function get_h(x, z, surf_obj, is_bottom) result(H)
    implicit none    
    type(type_surf_obj)  ::  surf_obj
    logical              ::  is_bottom
    real(dp)             ::  x,z,H, theta,m,n,a,b
    integer              ::  i
    real(dp), parameter  ::  const_2pi =  8.D0*DATAN(1.D0)

    H = 0.
    do i=0,surf_obj%private_L_mnab_coeffs-1
      m       =  surf_obj%private_mnab_coeffs(i,0)
      n       =  surf_obj%private_mnab_coeffs(i,1)
      a       =  surf_obj%private_mnab_coeffs(i,2)
      b       =  surf_obj%private_mnab_coeffs(i,3)
      theta   =  const_2pi * (m*(x/surf_obj%Lx) + n*(z/surf_obj%Lz))
      H       =  H + a*cos(theta) + b*sin(theta)
    end do

    if (.not.is_bottom) H = surf_obj%Ly - H

  end function
end module