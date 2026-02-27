module my_utils_prof_1d
  use, intrinsic :: iso_fortran_env, only: i8 => int64, dp => real64
  implicit none
  contains

  subroutine get_grid_params(nx, nz, ny, Ly, is_smooth)
    real(dp)   ::  Lx, Lz, Ly, nu, cond, Ret, Pr, dtmax, Sf_x, Sf_z, Sq, & ! declared for compatibility with params.dat
                   wall_BC_T_bot, wall_BC_T_top, &
                   wall_BC_U_bot, wall_BC_U_top, &
                   wall_BC_V_bot, wall_BC_V_top, &
                   wall_BC_W_bot, wall_BC_W_top,Gb_x,Gb_z,Gb_y, &
                   bulk_target_Reb, bulk_target_Ub
    integer  ::  nreps_x, nreps_z, nx, nz, ny, bulk_nprint, avg_iter_start, avg_nraw, avg_freq, avg_bin_nbins, &
                 snap_iter_start, snap_freq, filtering_steps, filtering_first, filtering_nstop, &
                 check_len_ref_angles, check_n_knots, i, j, k, ii, ni, nj, nk, use_rough_surf, &
                 use_incompressible, avg_apply_rms,use_utau_work, &
                 bulk_use_rho_weighting, bulk_use_target_Reb, bulk_use_target_Ub, &
                 force_Sf_x_accel, force_Sf_x_semilocal, force_Sf_x_pressure, prof_1d_avg_use, use_rough_surf_quady
    logical  ::  is_smooth
    include "../../input/params.dat"
    if ((use_rough_surf<0).or.(use_rough_surf>1)) error stop '((use_rough_surf<0).or.(use_rough_surf>1))'
    is_smooth = (use_rough_surf == 0)
  end subroutine 

  subroutine ibm_read_slabs(tag, side, A, ni_check, nj_check)
    implicit none
    integer       ::  A(0:,0:), iters, ni, nj, i,j, ni_check, nj_check
    character(len=1)  ::  tag,tag2
    character(len=3)  ::  side
    character(len=31)  ::  myfile
    tag2 = tag 
    if (tag2=='P') tag2='T'
    write(myfile, '(A22,A1,A1,A3,A4)') '../geom_data/ibm_inds_',tag2,'_',side,'.dat'
    open(19,file=myfile, action='read', status='old', form='unformatted', access='stream')
      read(19)  ni, nj
      if ((ni.ne.ni_check).or.(nj.ne.nj_check)) then 
        write(6,*) '(ni.ne.ni_check).or.(nj.ne.nj_check)', ni,nj,ni_check,nj_check ; flush(6)
        error stop '(ni.ne.ni_check).or.(nj.ne.nj_check)'
      end if
      do   j=0,nj-1
        do i=0,ni-1
          read(19)  A(i,j)
        end do
      end do 
    close(19)
  end subroutine

  subroutine read_y_dy_values(g_type, ny, y, dy, yu)
    implicit none
    character(len=1)       ::  g_type
    integer                ::  ny
    real(dp)               ::  y(0:), dy(0:), yu(0:)
    y  = 0
    dy = 0
    yu = 0
    open(19,file='../../input/yu.txt')
      block 
        integer :: i
        do i=0,ny
          read(19,*) yu(i)
        end do
      end block
    close(19)
    block 
      integer :: i
      do i=0,ny-1
        if      (g_type == 'p') then 
          y(i)    = (yu(i+1) + yu(i))*0.5_dp            !  yp(i)
          dy(i)   =  yu(i+1) - yu(i)                    !  dy_p(i)
        else if (g_type == 'u') then 
                        y(i)    =  yu(i)                              !  yu(i) (trimmed to [0:ny-1])
                        dy(i)   =  dy(i)   + (yu(i+1) - yu(i))*0.5_dp !  dy_u(i  ) += ...
          if ((i+1)<ny) dy(i+1) =  dy(i+1) + (yu(i+1) - yu(i))*0.5_dp !  dy_u(i+1) += ...
        else 
          error stop 'Unrecognized g_type'
        end if
      end do
    end block 
  end subroutine

  subroutine read_T_props(T_input, enth_ijk, rho_ijk, mu_ijk, cond_ijk, cp_ijk)
    real(dp), intent(in)   ::  T_input
    real(dp), intent(out)  ::  enth_ijk, rho_ijk, mu_ijk, cond_ijk, cp_ijk
    real(dp)               ::  T_ext(0:0,0:0,0:0)           , &
                               enth_arr(0:0,0:0,0:0)        , &
                               rho_arr(0:0,0:0,0:0)         , &
                               mu_arr(0:0,0:0,0:0)          , &
                               cond_arr(0:0,0:0,0:0)        , &
                               cp_arr(0:0,0:0,0:0)          , &
                               div_jacobian_arr(0:0,0:0,0:0), &
                               buoyancy_arr(0:0,0:0,0:0)    , drho_dT_ijk, T_ext_ijk
    integer                ::  i,j,k
    i          =  0
    j          =  0
    k          =  0
    T_ext_ijk  =  T_input
    T_ext      =  T_ext_ijk
    include '../thermophysical_properties/expressions_rho_mu_cond_buoyancy_cp_div_jacobian.f90'
    if (abs(T_ext(i,j,k)-T_ext_ijk)>1e-10)   error stop '(abs(T_ext(i,j,k)-T_ext_ijk)>1e-10)'
    if (abs(T_input-T_ext_ijk)>1e-10)        error stop '(abs(T_input-T_ext_ijk)>1e-10)'
    if (abs(enth_arr(i,j,k)-enth_ijk)>1e-10) error stop '(abs(enth_arr(i,j,k)-enth_ijk)>1e-10)'
    if (abs(rho_arr(i,j,k)-rho_ijk)>1e-10)   error stop '(abs(rho_arr(i,j,k)-rho_ijk)>1e-10)'
    if (abs(mu_arr(i,j,k)-mu_ijk)>1e-10)     error stop '(abs(mu_arr(i,j,k)-mu_ijk)>1e-10)'
    if (abs(cond_arr(i,j,k)-cond_ijk)>1e-10) error stop '(abs(cond_arr(i,j,k)-cond_ijk)>1e-10)'
    if (abs(cp_arr(i,j,k)-cp_ijk)>1e-10)     error stop '(abs(cp_arr(i,j,k)-cp_ijk)>1e-10)'
    if (i.ne.0)                              error stop 'i.ne.0'
    if (j.ne.0)                              error stop 'j.ne.0'
    if (k.ne.0)                              error stop 'k.ne.0'
  end subroutine
  
end module 

program main 
  use  my_utils_prof_1d
  implicit none 
  integer, allocatable    ::  slab_2d_bot(:,:), slab_2d_top(:,:)
  character(len=128)      ::  fname_avg
  character(len=128)      ::  fname_prof_1d
  real(dp), allocatable   ::  prof_1d(:,:), buffer_yu(:), y(:), dy(:)
  real(dp)                ::  Ly
  integer                 ::  nx, nz, ny, i_entry, n_entries
  logical                 ::  is_smooth
  character(len=1)        ::  tag_type, g_type
  character(len=8)        ::  tag_arr

  call get_grid_params(nx, nz, ny, Ly, is_smooth)
  allocate(slab_2d_bot(0:nx-1,0:nz-1) ,&
           slab_2d_top(0:nx-1,0:nz-1) ,&
           prof_1d(0:ny-1,0:3)        ,&
           buffer_yu(0:ny)            ,&
           y(0:ny-1), dy(0:ny-1)      )
           
  if (is_smooth) then
        slab_2d_bot = -1
        slab_2d_top =  ny
        write(6,*) 'INFO: Considering smooth surfaces. (balance_avg_1d_prof.f90)' ; flush(6)
 else
        write(6,*) 'INFO: Considering rough surfaces. (balance_avg_1d_prof.f90)' ; flush(6)
 end if
 
  open(39,file='all_interp_avg.dat')
    read(39,*) n_entries 
    do i_entry=0,n_entries-1
      read(39,*) tag_type, tag_arr
      read(39,'(A128)') fname_avg
      read(39,'(A128)') fname_prof_1d
      write(6,*) tag_type, ' ', tag_arr, ' ', fname_avg, ' ', fname_prof_1d;flush(6)

      if      ((tag_type == 'U').or.(tag_type == 'W').or.(tag_type == 'P').or.(tag_type == 'T')) then ; g_type = 'p'
      else if  (tag_type == 'V')                                                                 then ; g_type = 'u'
      else                                                                                     
        error stop 'Unrecognized tag_type'
      end if 
      call read_y_dy_values(g_type, ny, y, dy, buffer_yu)

      if (.not.(is_smooth)) then 
      call  ibm_read_slabs(tag_type, 'bot', slab_2d_bot, nx, nz)
      call  ibm_read_slabs(tag_type, 'top', slab_2d_top, nx, nz)
      end if
      block 
        real(dp)       ::  bulk, total_bulk
        integer        ::  nx_check,nz_check,ny_check
        open(19,file=trim(fname_avg), action='read', status='old', form='unformatted', access='stream')
          read(19) nx_check,nz_check,ny_check
          if ((nx.ne.nx_check).or.(nz.ne.nz_check).or.(ny.ne.ny_check)) then 
            write(6,*) 'Error: mismatch grid size', nx,nz,ny,nx_check,nz_check,ny_check;flush(6)
            error stop 'Error: mismatch grid size'
          end if 
          bulk       = 0
          if (is_smooth) then 
            total_bulk  =  Ly
            if (abs(buffer_yu(ny) - buffer_yu(0) - Ly)>1e-10) error stop '(abs(buffer_yu(ny) - buffer_yu(0) - Ly)>1e-10)'
          else
            total_bulk = 0
          end if
          block 
            real(dp)  ::  val, found, total
            integer   ::  i,j,k
            do     k=0,ny-1 
              total  =  0
              found  =  0
              do   j=0,nz-1
                do i=0,nx-1
                  read(19) val
                  if ((slab_2d_bot(i,j)<k).and.(k<slab_2d_top(i,j))) then 
                    total = total + val
                    found = found + 1
                  end if 
                end do 
              end do 
              prof_1d(k,0)  =  total/max(1e-10 + 0.d0, found)
              prof_1d(k,1)  =  found/((nx + 0.d0)*(nz + 0.d0))
              prof_1d(k,2)  =  y(k)
              prof_1d(k,3)  =  dy(k)
              bulk          =  bulk       + prof_1d(k,0)*prof_1d(k,1)*prof_1d(k,3)
              if (.not.is_smooth) then 
                total_bulk  =  total_bulk +              prof_1d(k,1)*prof_1d(k,3)
              else
                if (abs(prof_1d(k,1)-1)>1e-10) error stop 'abs(prof_1d(k,1)-1)>1e-10'
              end if 
            end do 
          end block
          if (abs(total_bulk)<1e-10) error stop '(abs(total_bulk)<1e-10)'
          bulk  =  bulk / total_bulk
        close(19)
        block 
          integer :: k
          open(19,file=trim(fname_prof_1d))
            write(19,'(I9,A46)') ny , ' ! next: val(i), frac_nonempty(i), y(i), dy(i)'
            do k=0,ny-1 
              write(19,'(4E25.16)') prof_1d(k,:)
            end do 
            write(19,'(E25.16,A3,A8,A46)') bulk,' ! ',tag_arr, '_bulk (weighted by frac_nonempty(:) and dy(:))'
            if (tag_arr == 'T_______') then
              block 
                real(dp) :: T_ext_bulk, enth_bulk, rho_bulk, mu_bulk, cond_bulk, cp_bulk
                T_ext_bulk = bulk
                call read_T_props(T_ext_bulk, enth_bulk, rho_bulk, mu_bulk, cond_bulk, cp_bulk)
                write(19,'(E25.16,A57)') enth_bulk ,' ! enth    _bulk (weighted by frac_nonempty(:) and dy(:))'
                write(19,'(E25.16,A57)') rho_bulk  ,' ! rho     _bulk (weighted by frac_nonempty(:) and dy(:))'
                write(19,'(E25.16,A57)') mu_bulk   ,' ! mu      _bulk (weighted by frac_nonempty(:) and dy(:))'
                write(19,'(E25.16,A57)') cond_bulk ,' ! cond    _bulk (weighted by frac_nonempty(:) and dy(:))'
                write(19,'(E25.16,A57)') cp_bulk   ,' ! cp      _bulk (weighted by frac_nonempty(:) and dy(:))'
              end block
            end if 
          close(19)
        end block 
      end block
    end do 
    close(39)
end program 
