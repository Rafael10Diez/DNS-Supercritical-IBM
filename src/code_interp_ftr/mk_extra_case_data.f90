
module mod_extended_params
  use, intrinsic :: iso_fortran_env, only: i8 => int64, dp => real64
  implicit none

  type type_extended_params_obj
    real(dp)  ::  Lx, Lz, Ly, Sf_x, Sf_z, Sq, &
                  wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                  wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top, &
                  wall_BC_rho_bot  , wall_BC_mu_bot   , wall_BC_cond_bot , &
                  wall_BC_rho_top  , wall_BC_mu_top   , wall_BC_cond_top
    integer   ::  nx, nz, ny, nreps_x, nreps_z
  end type

  contains

  subroutine write_extended_params_obj(params_obj)
    implicit none
    type(type_extended_params_obj)  ::  params_obj
    ! all parameters declared for compatibility with params.dat
    real(dp)              :: Lx, Lz, Ly, nu, cond, Ret, Pr, dtmax, Sf_x, Sf_z, Sq, &
                             Gb_x, Gb_z, Gb_y, bulk_target_Ub, bulk_target_Reb, &
                             wall_BC_T_bot, wall_BC_T_top, wall_BC_U_bot, wall_BC_U_top, wall_BC_V_bot, wall_BC_V_top, &
                             wall_BC_W_bot, wall_BC_W_top
    integer               :: nreps_x, nreps_z, nx, nz, ny, bulk_nprint, avg_iter_start, avg_nraw, avg_freq, avg_bin_nbins, &
                             snap_iter_start, snap_freq, filtering_steps, filtering_first, filtering_nstop, use_rough_surf, &
                             bulk_use_target_Ub, bulk_use_target_Reb, bulk_use_rho_weighting, force_Sf_x_semilocal, &
                             force_Sf_x_pressure, force_Sf_x_accel, use_incompressible, use_utau_work, prof_1d_avg_use, &
                             use_rough_surf_quady

    include "params.dat"

    params_obj%Lx             =  Lx
    params_obj%Lz             =  Lz
    params_obj%Ly             =  Ly
    params_obj%Sf_x           =  Sf_x
    params_obj%Sf_z           =  Sf_z
    params_obj%Sq             =  Sq
    params_obj%wall_BC_T_bot  =  wall_BC_T_bot
    params_obj%wall_BC_T_top  =  wall_BC_T_top
    params_obj%wall_BC_U_bot  =  wall_BC_U_bot
    params_obj%wall_BC_U_top  =  wall_BC_U_top
    params_obj%wall_BC_V_bot  =  wall_BC_V_bot
    params_obj%wall_BC_V_top  =  wall_BC_V_top
    params_obj%wall_BC_W_bot  =  wall_BC_W_bot
    params_obj%wall_BC_W_top  =  wall_BC_W_top
    params_obj%nx             =  nx
    params_obj%nz             =  nz
    params_obj%ny             =  ny
    params_obj%nreps_x        =  nreps_x
    params_obj%nreps_z        =  nreps_z
    if ((use_rough_surf_quady==1).and.(use_rough_surf.ne.1)) error stop '(use_rough_surf_quady.and.not.(use_rough_surf))'
    if (use_rough_surf.ne.1)                                 error stop '(.not.use_rough_surf)'

    block
      integer                ::  ii_side
      do ii_side=0,1
        block
          real(dp)  ::  T_ext_ijk, enth_ijk, rho_ijk, mu_ijk, cond_ijk, cp_ijk, drho_dT_ijk, &
                        T_ext(0:0,0:0,0:0)           , &
                        enth_arr(0:0,0:0,0:0)        , &
                        rho_arr(0:0,0:0,0:0)         , &
                        mu_arr(0:0,0:0,0:0)          , &
                        cond_arr(0:0,0:0,0:0)        , &
                        cp_arr(0:0,0:0,0:0)          , &
                        div_jacobian_arr(0:0,0:0,0:0), &
                        buoyancy_arr(0:0,0:0,0:0)
          integer   ::  i,j,k
          i = 0
          j = 0
          k = 0
          if (ii_side == 0) then
            T_ext     = params_obj%wall_BC_T_bot
            T_ext_ijk = params_obj%wall_BC_T_bot
          else
            T_ext     = params_obj%wall_BC_T_top
            T_ext_ijk = params_obj%wall_BC_T_top
          end if
          include 'expressions_rho_mu_cond_buoyancy_cp_div_jacobian.f90'
          if (abs(T_ext(i,j,k)-T_ext_ijk)>1e-10)   error stop '(abs(T_ext(i,j,k)-T_ext_ijk)>1e-10)'
          if (abs(rho_arr(i,j,k)-rho_ijk)>1e-10)   error stop '(abs(rho_arr(i,j,k)-rho_ijk)>1e-10)'
          if (abs(mu_arr(i,j,k)-mu_ijk)>1e-10)     error stop '(abs(mu_arr(i,j,k)-mu_ijk)>1e-10)'
          if (abs(cond_arr(i,j,k)-cond_ijk)>1e-10) error stop '(abs(cond_arr(i,j,k)-cond_ijk)>1e-10)'
          if (i.ne.0) error stop 'i.ne.0'
          if (j.ne.0) error stop 'j.ne.0'
          if (k.ne.0) error stop 'k.ne.0'
          if (ii_side == 0) then
            params_obj%wall_BC_rho_bot   = rho_arr(i,j,k)
            params_obj%wall_BC_mu_bot    = mu_arr(i,j,k)
            params_obj%wall_BC_cond_bot  = cond_arr(i,j,k)
          else
            params_obj%wall_BC_rho_top   = rho_arr(i,j,k)
            params_obj%wall_BC_mu_top    = mu_arr(i,j,k)
            params_obj%wall_BC_cond_top  = cond_arr(i,j,k)
          end if
        end block
      end do
    end block
  end subroutine

  subroutine read_extended_params_obj(params_obj)
    implicit none
    type(type_extended_params_obj)  ::  params_obj
    real(dp)  ::  Lx, Lz, Ly, Sf_x, Sf_z, Sq, &
                  wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                  wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top, &
                  wall_BC_rho_bot  , wall_BC_mu_bot   , wall_BC_cond_bot , &
                  wall_BC_rho_top  , wall_BC_mu_top   , wall_BC_cond_top
    integer   ::  nx, nz, ny, nreps_x, nreps_z

    include 'extended_params.dat'

    params_obj%Lx   = Lx
    params_obj%Lz   = Lz
    params_obj%Ly   = Ly
    params_obj%Sf_x = Sf_x
    params_obj%Sf_z = Sf_z
    params_obj%Sq   = Sq

    params_obj%wall_BC_U_bot = wall_BC_U_bot
    params_obj%wall_BC_V_bot = wall_BC_V_bot
    params_obj%wall_BC_W_bot = wall_BC_W_bot
    params_obj%wall_BC_T_bot = wall_BC_T_bot

    params_obj%wall_BC_U_top = wall_BC_U_top
    params_obj%wall_BC_V_top = wall_BC_V_top
    params_obj%wall_BC_W_top = wall_BC_W_top
    params_obj%wall_BC_T_top = wall_BC_T_top

    params_obj%wall_BC_rho_bot  = wall_BC_rho_bot
    params_obj%wall_BC_mu_bot   = wall_BC_mu_bot
    params_obj%wall_BC_cond_bot = wall_BC_cond_bot

    params_obj%wall_BC_rho_top  = wall_BC_rho_top
    params_obj%wall_BC_mu_top   = wall_BC_mu_top
    params_obj%wall_BC_cond_top = wall_BC_cond_top

    params_obj%nx = nx
    params_obj%nz = nz
    params_obj%ny = ny

    params_obj%nreps_x = nreps_x
    params_obj%nreps_z = nreps_z

    write(6,*) '---------------------- Read Parameters (extended_params.dat) ----------------------'
    write(6,'(A30,A3,E25.16)') 'Lx'              , ' = ' , params_obj%Lx
    write(6,'(A30,A3,E25.16)') 'Lz'              , ' = ' , params_obj%Lz
    write(6,'(A30,A3,E25.16)') 'Ly'              , ' = ' , params_obj%Ly
    write(6,'(A30,A3,E25.16)') 'Sf_x'            , ' = ' , params_obj%Sf_x
    write(6,'(A30,A3,E25.16)') 'Sf_z'            , ' = ' , params_obj%Sf_z
    write(6,'(A30,A3,E25.16)') 'Sq'              , ' = ' , params_obj%Sq

    write(6,'(A30,A3,E25.16)') 'wall_BC_U_bot'   , ' = ' , params_obj%wall_BC_U_bot
    write(6,'(A30,A3,E25.16)') 'wall_BC_V_bot'   , ' = ' , params_obj%wall_BC_V_bot
    write(6,'(A30,A3,E25.16)') 'wall_BC_W_bot'   , ' = ' , params_obj%wall_BC_W_bot
    write(6,'(A30,A3,E25.16)') 'wall_BC_T_bot'   , ' = ' , params_obj%wall_BC_T_bot

    write(6,'(A30,A3,E25.16)') 'wall_BC_U_top'   , ' = ' , params_obj%wall_BC_U_top
    write(6,'(A30,A3,E25.16)') 'wall_BC_V_top'   , ' = ' , params_obj%wall_BC_V_top
    write(6,'(A30,A3,E25.16)') 'wall_BC_W_top'   , ' = ' , params_obj%wall_BC_W_top
    write(6,'(A30,A3,E25.16)') 'wall_BC_T_top'   , ' = ' , params_obj%wall_BC_T_top

    write(6,'(A30,A3,E25.16)') 'wall_BC_rho_bot' , ' = ' , params_obj%wall_BC_rho_bot
    write(6,'(A30,A3,E25.16)') 'wall_BC_mu_bot'  , ' = ' , params_obj%wall_BC_mu_bot
    write(6,'(A30,A3,E25.16)') 'wall_BC_cond_bot', ' = ' , params_obj%wall_BC_cond_bot

    write(6,'(A30,A3,E25.16)') 'wall_BC_rho_top' , ' = ' , params_obj%wall_BC_rho_top
    write(6,'(A30,A3,E25.16)') 'wall_BC_mu_top'  , ' = ' , params_obj%wall_BC_mu_top
    write(6,'(A30,A3,E25.16)') 'wall_BC_cond_top', ' = ' , params_obj%wall_BC_cond_top

    write(6,'(A30,A3,I25)')    'nx'              , ' = ' , params_obj%nx
    write(6,'(A30,A3,I25)')    'nz'              , ' = ' , params_obj%nz
    write(6,'(A30,A3,I25)')    'ny'              , ' = ' , params_obj%ny
    write(6,'(A30,A3,I25)')    'nreps_x'         , ' = ' , params_obj%nreps_x
    write(6,'(A30,A3,I25)')    'nreps_z'         , ' = ' , params_obj%nreps_z
    flush(6)
  end subroutine
end module

program main
  use mod_extended_params
  use Mod_Height_Function
  implicit none
  type(type_extended_params_obj)  ::  params_obj
  type(type_surf_obj)             ::  surf_obj

  call  write_extended_params_obj(params_obj)

  open(19,file='extended_params.dat')
    write(19,'(A30,A3,E25.16)') 'Lx'              , ' = ' , params_obj%Lx
    write(19,'(A30,A3,E25.16)') 'Lz'              , ' = ' , params_obj%Lz
    write(19,'(A30,A3,E25.16)') 'Ly'              , ' = ' , params_obj%Ly
    write(19,'(A30,A3,E25.16)') 'Sf_x'            , ' = ' , params_obj%Sf_x
    write(19,'(A30,A3,E25.16)') 'Sf_z'            , ' = ' , params_obj%Sf_z
    write(19,'(A30,A3,E25.16)') 'Sq'              , ' = ' , params_obj%Sq

    write(19,'(A30,A3,E25.16)') 'wall_BC_U_bot'   , ' = ' , params_obj%wall_BC_U_bot
    write(19,'(A30,A3,E25.16)') 'wall_BC_V_bot'   , ' = ' , params_obj%wall_BC_V_bot
    write(19,'(A30,A3,E25.16)') 'wall_BC_W_bot'   , ' = ' , params_obj%wall_BC_W_bot
    write(19,'(A30,A3,E25.16)') 'wall_BC_T_bot'   , ' = ' , params_obj%wall_BC_T_bot

    write(19,'(A30,A3,E25.16)') 'wall_BC_U_top'   , ' = ' , params_obj%wall_BC_U_top
    write(19,'(A30,A3,E25.16)') 'wall_BC_V_top'   , ' = ' , params_obj%wall_BC_V_top
    write(19,'(A30,A3,E25.16)') 'wall_BC_W_top'   , ' = ' , params_obj%wall_BC_W_top
    write(19,'(A30,A3,E25.16)') 'wall_BC_T_top'   , ' = ' , params_obj%wall_BC_T_top

    write(19,'(A30,A3,E25.16)') 'wall_BC_rho_bot' , ' = ' , params_obj%wall_BC_rho_bot
    write(19,'(A30,A3,E25.16)') 'wall_BC_mu_bot'  , ' = ' , params_obj%wall_BC_mu_bot
    write(19,'(A30,A3,E25.16)') 'wall_BC_cond_bot', ' = ' , params_obj%wall_BC_cond_bot

    write(19,'(A30,A3,E25.16)') 'wall_BC_rho_top' , ' = ' , params_obj%wall_BC_rho_top
    write(19,'(A30,A3,E25.16)') 'wall_BC_mu_top'  , ' = ' , params_obj%wall_BC_mu_top
    write(19,'(A30,A3,E25.16)') 'wall_BC_cond_top', ' = ' , params_obj%wall_BC_cond_top

    write(19,'(A30,A3,I25)')    'nx'              , ' = ' , params_obj%nx
    write(19,'(A30,A3,I25)')    'nz'              , ' = ' , params_obj%nz
    write(19,'(A30,A3,I25)')    'ny'              , ' = ' , params_obj%ny

    write(19,'(A30,A3,I25)')    'nreps_x'         , ' = ' , params_obj%nreps_x
    write(19,'(A30,A3,I25)')    'nreps_z'         , ' = ' , params_obj%nreps_z
  close(19)

  block
    integer  ::  ii_side, ii_var
    call  init_surf_obj(surf_obj)
    do   ii_side = 0,1
      do ii_var  = 0,2
        block
          integer            ::  i, j
          real(dp)           ::  x_now, z_now, y_now, dx, dz, di, dj
          logical            ::  is_bottom
          character(len=2)   ::  tag
          character(len=3)   ::  side
          character(len=31)  ::  myfile

          if (ii_side == 0) then ; side='bot' ; is_bottom = .true.  ; end if
          if (ii_side == 1) then ; side='top' ; is_bottom = .false. ; end if

          if (ii_var  == 0) then ; tag='pp'; di = 0.5; dj = 0.5 ; end if
          if (ii_var  == 1) then ; tag='up'; di = 0  ; dj = 0.5 ; end if
          if (ii_var  == 2) then ; tag='pu'; di = 0.5; dj = 0   ; end if

          write(myfile, '(A21,A2,A1,A3,A4)') './surfs_hxz/surf_hxz_',tag,'_',side,'.dat'

          dx         =  params_obj%Lx / (params_obj%nx + 0.d0)
          dz         =  params_obj%Lz / (params_obj%nz + 0.d0)

          open(19,file=myfile, action='write', status='replace', form='unformatted', access='stream')
            write(19) params_obj%nx, params_obj%nz
            do    j=0,params_obj%nz-1
              do  i=0,params_obj%nx-1
                x_now  =  (i + di) * dx
                z_now  =  (j + dj) * dz
                y_now  =  get_h(x_now, z_now, surf_obj, is_bottom)
                write(19) y_now
              end do
            end do
          close(19)
        end block
      end do
    end do
  end block
end program