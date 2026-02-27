module Mod_Height_Function
  implicit none

  ! --------------- public module API ---------------
  private
  public :: init_surf_obj, type_surf_obj, get_h

  ! --------------- surface object ---------------
  type type_surf_obj
    real*8               :: Lx, Lz, Ly
    integer              :: nreps_x, nreps_z
  end type

  contains

  subroutine init_surf_obj(surf_obj)
    implicit none    
    type(type_surf_obj) :: surf_obj
    real*8              :: Lx, Lz, Ly, nu, cond, Ret, Pr, dtmax, Sf_x, Sf_z, Sq, & ! declared for compatibility with params.dat
                           Gb_x, Gb_z, Gb_y, bulk_target_Ub, bulk_target_Reb, &
                           wall_BC_T_bot, wall_BC_T_top, wall_BC_U_bot, wall_BC_U_top, wall_BC_V_bot, wall_BC_V_top, &
                           wall_BC_W_bot, wall_BC_W_top
    integer             :: nreps_x, nreps_z, nx, nz, ny, bulk_nprint, avg_iter_start, avg_nraw, avg_freq, avg_bin_nbins, &
                           snap_iter_start, snap_freq, filtering_steps, filtering_first, filtering_nstop, use_rough_surf, &
                           bulk_use_target_Ub, bulk_use_target_Reb, bulk_use_rho_weighting, force_Sf_x_semilocal, &
                           force_Sf_x_pressure, force_Sf_x_accel, use_incompressible, use_utau_work, prof_1d_avg_use, &
                           use_rough_surf_quady
    
    include "params.dat"
    
    surf_obj%Lx                      =  Lx
    surf_obj%Lz                      =  Lz
    surf_obj%Ly                      =  Ly
    surf_obj%nreps_x                 =  nreps_x
    surf_obj%nreps_z                 =  nreps_z
  end subroutine

  function get_h(x, z, surf_obj, is_bottom) result(H)
    implicit none    
    type(type_surf_obj)  ::  surf_obj
    logical              ::  is_bottom
    real*8               ::  x,z,H
    
    if (is_bottom) then 
      H = 0.
    else 
      H = surf_obj%Ly
    end if
    
  end function
end module