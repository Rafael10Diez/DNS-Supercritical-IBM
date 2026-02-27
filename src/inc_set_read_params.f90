
! read general parameters
block
  integer :: mpi_ierr, int_rough_surf, int_use_incompressible, int_use_utau_work, int_rough_surf_quady, &
             int_bulk_use_target_Ub, int_bulk_use_rho_weighting, int_bulk_use_target_Reb, &
             int_force_Sf_x_pressure, int_force_Sf_x_accel, int_force_Sf_x_semilocal
  if (irank_mpi == 0) then 
    block 
      character(len=16) :: aux_str_i
      call get_command_argument(1, aux_str_i) ; read(aux_str_i,*) mpi_divs_z
      call get_command_argument(2, aux_str_i) ; read(aux_str_i,*) mpi_divs_y
      if ((mpi_divs_y*mpi_divs_z).ne.nproc_mpi) error stop '(mpi_divs_y*mpi_divs_z).ne.nproc_mpi'
    end block 
    open(19,file="./geom_data/base_params.txt")
      read(19,*) nx_glob, nz_glob, ny_glob
      read(19,*) int_rough_surf, int_rough_surf_quady
      read(19,*) nstep, dtmax
      read(19,*) inv_dx, inv_dz ! first we read dx,dz and invert them later
      read(19,*) nu, cond
      read(19,*) Sf_x, Sf_z, Sq
      read(19,*) Gb_x, Gb_z, Gb_y
      read(19,*) avg_iter_start, avg_nraw, avg_freq, avg_bin_nbins
      read(19,*) snap_iter_start, snap_freq
      read(19,*) bulk_nprint
      read(19,*) filtering_steps, filtering_first, filtering_nstop
      if ((snap_iter_start<0).or.(snap_freq<0)) then 
          write(6,'(A68,2I10)') 'INFO: disabling intermediate snapshots [snap_iter_start,snap_freq]: ',&
          snap_iter_start, snap_freq;flush(6)
        snap_iter_start = nstep+100
        snap_freq       = nstep+100
      end if
      if ((filtering_steps<0).or.(filtering_first<0).or.(filtering_nstop<0)) then 
        write(6,'(A77,3I10)') 'INFO: disabling filtering [filtering_steps,filtering_first,filtering_nstop]: ',&
        filtering_steps,filtering_first,filtering_nstop;flush(6)
        filtering_steps = -1
        filtering_first = -1
        filtering_nstop = -1
      end if
      inv_dx = 1.d0/inv_dx ! invert dx
      inv_dz = 1.d0/inv_dz ! invert dz
    read(19,*) wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot
    read(19,*) wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top
    read(19,*) int_use_incompressible
    read(19,*) int_use_utau_work
    read(19,*) extrap_T_wall_0 , extrap_T_wall_1, extrap_T_wall_m1, extrap_T_wall_m2
    read(19,*) bulk_target_Ub, int_bulk_use_target_Ub, int_bulk_use_rho_weighting
    read(19,*) bulk_target_Reb, Ly, int_bulk_use_target_Reb
    read(19,*) int_force_Sf_x_pressure, int_force_Sf_x_accel, int_force_Sf_x_semilocal
    close(19)
  end if
  
  call MPI_BCAST(nx_glob              ,1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(nz_glob              ,1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(ny_glob              ,1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(mpi_divs_z           ,1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(mpi_divs_y           ,1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(int_rough_surf       ,1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(int_rough_surf_quady ,1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(nstep                ,1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(dtmax                ,1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(inv_dx               ,1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(inv_dz               ,1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(nu                   ,1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(cond                 ,1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(Sf_x                 ,1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(Sf_z                 ,1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(Sq                   ,1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(Gb_x                 ,1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(Gb_z                 ,1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(Gb_y                 ,1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
  
  call MPI_BCAST(avg_iter_start , 1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(avg_nraw       , 1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(avg_freq       , 1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(avg_bin_nbins  , 1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(snap_iter_start, 1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(snap_freq      , 1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(bulk_nprint    , 1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(filtering_steps, 1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(filtering_first, 1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(filtering_nstop, 1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)

  call MPI_BCAST(wall_BC_U_bot         , 1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(wall_BC_V_bot         , 1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(wall_BC_W_bot         , 1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(wall_BC_T_bot         , 1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(wall_BC_U_top         , 1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(wall_BC_V_top         , 1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(wall_BC_W_top         , 1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(wall_BC_T_top         , 1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(int_use_incompressible, 1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(int_use_utau_work     , 1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)

  call MPI_BCAST(extrap_T_wall_0       , 1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(extrap_T_wall_1       , 1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(extrap_T_wall_m1      , 1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(extrap_T_wall_m2      , 1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)

  call MPI_BCAST(bulk_target_Ub            ,1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(int_bulk_use_target_Ub    ,1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(int_bulk_use_rho_weighting,1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
  
  call MPI_BCAST(bulk_target_Reb           ,1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(Ly                        ,1,MPI_DOUBLE_PRECISION, 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(int_bulk_use_target_Reb   ,1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(int_force_Sf_x_pressure   ,1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(int_force_Sf_x_accel      ,1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
  call MPI_BCAST(int_force_Sf_x_semilocal       ,1,MPI_INTEGER         , 0, mpi_comm_world, mpi_ierr)
  
  if ((int_rough_surf<0).or.(int_rough_surf>1)) then 
    error stop '((int_rough_surf<0).or.(int_rough_surf>1))'
  end if 
  if ((int_rough_surf_quady<0).or.(int_rough_surf_quady>1)) then 
    error stop '((int_rough_surf_quady<0).or.(int_rough_surf_quady>1))'
  end if 
  if ((int_use_incompressible<0).or.(int_use_incompressible>1)) then 
    error stop '((int_use_incompressible<0).or.(int_use_incompressible>1))'
  end if 
  if ((int_use_utau_work<0).or.(int_use_utau_work>1)) then 
    error stop '((int_use_utau_work<0).or.(int_use_utau_work>1))'
  end if 
  if ((int_bulk_use_target_Ub<0).or.(int_bulk_use_target_Ub>1)) then 
    error stop '((int_bulk_use_target_Ub<0).or.(int_bulk_use_target_Ub>1))'
  end if 
  if ((int_bulk_use_rho_weighting<0).or.(int_bulk_use_rho_weighting>1)) then 
    error stop '((int_bulk_use_rho_weighting<0).or.(int_bulk_use_rho_weighting>1))'
  end if 
  if ((int_bulk_use_target_Reb<0).or.(int_bulk_use_target_Reb>1)) then 
    error stop '((int_bulk_use_target_Reb<0).or.(int_bulk_use_target_Reb>1))'
  end if 
  if ((int_force_Sf_x_pressure<0).or.(int_force_Sf_x_pressure>1)) then 
    error stop '((int_force_Sf_x_pressure<0).or.(int_force_Sf_x_pressure>1))'
  end if 
  if ((int_force_Sf_x_accel<0).or.(int_force_Sf_x_accel>1)) then 
    error stop '((int_force_Sf_x_accel<0).or.(int_force_Sf_x_accel>1))'
  end if 
  if ((int_force_Sf_x_semilocal<0).or.(int_force_Sf_x_semilocal>1)) then 
    error stop '((int_force_Sf_x_semilocal<0).or.(int_force_Sf_x_semilocal>1))'
  end if 

  use_rough_surf          =  (int_rough_surf             == 1)
  use_rough_surf_quady    =  (int_rough_surf_quady       == 1)
  use_incompressible      =  (int_use_incompressible     == 1)
  use_utau_work           =  (int_use_utau_work          == 1)
  bulk_use_target_Ub      =  (int_bulk_use_target_Ub     == 1)
  bulk_use_rho_weighting  =  (int_bulk_use_rho_weighting == 1)
  bulk_use_target_Reb     =  (int_bulk_use_target_Reb    == 1)
  force_Sf_x_pressure     =  (int_force_Sf_x_pressure    == 1)
  force_Sf_x_accel        =  (int_force_Sf_x_accel       == 1)
  force_Sf_x_semilocal         =  (int_force_Sf_x_semilocal        == 1)

  if (.not.bulk_use_target_Ub) then 
    if (abs(bulk_target_Ub)>1e-10) error stop 'abs(bulk_target_Ub)>1e-10'
  end if
  if (.not.bulk_use_target_Reb) then 
    if (abs(bulk_target_Reb)>1e-10) error stop 'abs(bulk_target_Reb)>1e-10'
  end if

  use_filtering      = ((filtering_steps>=0).and.(filtering_first>=0).and.(filtering_nstop>=0))

  if (use_incompressible.and.((abs(Gb_x)>0).or.(abs(Gb_y)>0).or.(abs(Gb_z)>0))) then 
      error stop 'ERROR: bouyancy not implemented for compressible flows'
  end if

end block
! expected conversions
  inv_dx__2  =  0.5*inv_dx
  inv_dz__2  =  0.5*inv_dz
  inv_dx2    =      inv_dx**2
  inv_dz2    =      inv_dz**2
  inv_dtmax  =  1.d0/dtmax
