
module mod_interp_routines
  use, intrinsic :: iso_fortran_env, only: i8 => int64, dp => real64

  implicit none

  type type_info
    real(dp)               ::  wall_BC, Lx, Lz, dx, dz, dy_ref
    integer                ::  nx, nz, ny, ny_fluid, ny_wall, iter_avg
    character(len=1)       ::  tag
    character(len=3)       ::  side
    real(dp), allocatable  ::  cen_slab_y  (:,:), var_slab_y(:,:), var_data_3d(:,:,:), var_y_3d(:,:,:), &
                               coords_yu(:), coords_yp(:), yy(:)
  end type

  contains

  subroutine init_info_obj(obj_info)
    implicit none
    type(type_info)                 ::  obj_info

    call get_command_argument(1, obj_info%tag)
    call get_command_argument(2, obj_info%side)

    block
      character(len=30) temp_chr
      call get_command_argument(3, temp_chr)
      read(temp_chr,*) obj_info%iter_avg
    end block

    call read_params_obj(obj_info)

    if       (obj_info%tag=='U') then ; obj_info%ny_fluid = 1 ; obj_info%ny_wall = 1 ;
    else if  (obj_info%tag=='V') then ; obj_info%ny_fluid = 1 ; obj_info%ny_wall = 1 ;
    else if  (obj_info%tag=='W') then ; obj_info%ny_fluid = 1 ; obj_info%ny_wall = 1 ;
    else if  (obj_info%tag=='T') then ; obj_info%ny_fluid = 1 ; obj_info%ny_wall = 1 ;
    else if  (obj_info%tag=='P') then ; obj_info%ny_fluid = 2 ; obj_info%ny_wall = 0 ;
    else ; error stop 'Unrecognized (tag)' ; end if

    write(6,*) ''
    write(6,*) '---------------------- Read Parameters (obj_info) ----------------------'
    write(6,'(A30,A3,A25)')    'tag'             , ' = ' , obj_info%tag
    write(6,'(A30,A3,A25)')    'side'            , ' = ' , obj_info%side
    write(6,'(A30,A3,I25)')    'nx'              , ' = ' , obj_info%nx
    write(6,'(A30,A3,I25)')    'nz'              , ' = ' , obj_info%nz
    write(6,'(A30,A3,I25)')    'ny'              , ' = ' , obj_info%ny
    write(6,'(A30,A3,I25)')    'ny_fluid'        , ' = ' , obj_info%ny_fluid
    write(6,'(A30,A3,I25)')    'ny_wall'         , ' = ' , obj_info%ny_wall
    write(6,'(A30,A3,E25.16)') 'Lx'              , ' = ' , obj_info%Lx
    write(6,'(A30,A3,E25.16)') 'Lz'              , ' = ' , obj_info%Lz
    write(6,'(A30,A3,E25.16)') 'dx'              , ' = ' , obj_info%dx
    write(6,'(A30,A3,E25.16)') 'dz'              , ' = ' , obj_info%dz
    write(6,'(A30,A3,E25.16)') 'wall_BC'         , ' = ' , obj_info%wall_BC
    flush(6)

    block
      integer :: k
      allocate(obj_info%coords_yu(0:obj_info%ny  ), &
               obj_info%coords_yp(0:obj_info%ny-1))
      open(19,file='yu.txt')
        do k=0,obj_info%ny
          read(19,*) obj_info%coords_yu(k)
        end do
      close(19)
    end block

    block
      integer :: k
      obj_info%dy_ref = obj_info%coords_yu(1) - obj_info%coords_yu(0)
      do k=1,obj_info%ny
        obj_info%coords_yp(k-1) = 0.5d0*( obj_info%coords_yu(k) + obj_info%coords_yu(k-1) )
        obj_info%dy_ref = min(obj_info%dy_ref, obj_info%coords_yu(k)- obj_info%coords_yu(k-1))
      end do
    end block

    if      (obj_info%tag=='U') then ; allocate(obj_info%yy, mold = obj_info%coords_yp) ; obj_info%yy = obj_info%coords_yp;
    else if (obj_info%tag=='W') then ; allocate(obj_info%yy, mold = obj_info%coords_yp) ; obj_info%yy = obj_info%coords_yp;
    else if (obj_info%tag=='V') then ; allocate(obj_info%yy, mold = obj_info%coords_yu) ; obj_info%yy = obj_info%coords_yu;
    else if (obj_info%tag=='T') then ; allocate(obj_info%yy, mold = obj_info%coords_yp) ; obj_info%yy = obj_info%coords_yp;
    else if (obj_info%tag=='P') then ; allocate(obj_info%yy, mold = obj_info%coords_yp) ; obj_info%yy = obj_info%coords_yp;
    else ; error stop 'Unrecognized tag' ; end if

  end subroutine

  subroutine read_2d_arr(tag, side, A)
    character(len=1)  ::  tag
    character(len=3)  ::  side
    character(len=2)  ::  mm
    character(len=31) ::  myfile
    real(dp)          ::  A(0:,0:)
    integer           ::  nx,nz,i,j

    if       (tag=='U') then ; mm = 'up' ;
    else if  (tag=='W') then ; mm = 'pu' ;
    else if  (tag=='V') then ; mm = 'pp' ;
    else if  (tag=='T') then ; mm = 'pp' ;
    else if  (tag=='P') then ; mm = 'pp' ;
    else ; error stop 'Unrecognized tag' ; end if
    write(myfile, '(A21,A2,A1,A3,A4)') './surfs_hxz/surf_hxz_',mm,'_',side,'.dat'

    open(19,file=myfile, action='read', status='old', form='unformatted', access='stream')
      read(19) nx, nz
      do   j=0,nz-1
        do i=0,nx-1
          read(19) A(i,j)
        end do
      end do
    close(19)
  end subroutine

  subroutine fill_3d_arrs(obj_info, tol_y)
    type(type_info)  ::  obj_info
    real(dp)         ::  tol_y

    allocate(obj_info%cen_slab_y(0:obj_info%nx-1, 0:obj_info%nz-1))
    allocate(obj_info%var_slab_y, mold=obj_info%cen_slab_y)

    allocate(obj_info%var_data_3d (0:obj_info%nx-1, 0:obj_info%nz-1, 0:(obj_info%ny_fluid+obj_info%ny_wall-1)))
    allocate(obj_info%var_y_3d     , mold=obj_info%var_data_3d)

    if ((obj_info%ny_fluid+obj_info%ny_wall).ne.2) error stop '(obj_info%ny_fluid+obj_info%ny_wall).ne.2'

    call read_2d_arr('T'         , obj_info%side, obj_info%cen_slab_y)
    call read_2d_arr(obj_info%tag, obj_info%side, obj_info%var_slab_y)

    if (obj_info%ny_wall>0) then
      if (obj_info%ny_fluid.ne.1) error stop 'obj_info%ny_fluid.ne.1'
      block
        integer :: k
        k = obj_info%ny_fluid+obj_info%ny_wall-1
        obj_info%var_y_3d   (:,:,k)  =  obj_info%var_slab_y(:,:)
        obj_info%var_data_3d(:,:,k)  =  obj_info%wall_BC
      end block
    end if

    block
      integer :: i,j,p,k0,K_dir
      logical :: is_bottom
      integer(kind=i8)   ::  offset_b, nx64, nz64, i64, j64, k64
      character(len=47)  ::  myfile
      nx64  =  obj_info%nx
      nz64  =  obj_info%nz

      if      (obj_info%side=='bot') then ; is_bottom=.true.  ; K_dir =  1;
      else if (obj_info%side=='top') then ; is_bottom=.false. ; K_dir = -1;
      else                                ; error stop 'Unrecognized side' ;
      end if

      p = obj_info%ny/2
      write(myfile,'(A23,A1,A10,I0.9,A4)') '../../avg/array_arr_3d_',obj_info%tag,'_avg_iter_',obj_info%iter_avg,'.dat'
      open(19,file=myfile, action='read', status='old', form='unformatted', access='stream')
        do      j=0,obj_info%nz-1
          do    i=0,obj_info%nx-1
            call move_pointer(p, obj_info%yy, obj_info%ny, obj_info%var_slab_y(i,j)+tol_y*K_dir, is_bottom)
            do k0=0,obj_info%ny_fluid-1
              i64 = i ; j64 = j ; k64 = p+k0*K_dir
              obj_info%var_y_3d(i,j,k0) = obj_info%yy(p+k0*K_dir)
              offset_b  =  1_i8 + 12_i8 + 8_i8*((k64*nz64 + j64)*nx64 + i64)
              read(19,pos=offset_b) obj_info%var_data_3d(i,j,k0)
            end do
          end do
        end do
      close(19)
    end block
  end subroutine

  subroutine move_pointer(p,A,n,y0,is_bottom)
    real(dp) :: A(0:), y0
    integer  :: p,n
    logical  :: is_bottom
    if (is_bottom) then
      do while (A(p)<=y0)
        p = p + 1
        if (p>=n) error stop 'p>=n'
      end do
      do while (A(p)>y0)
        p = p - 1
        if (p<0) error stop 'p<0'
      end do
      p = p + 1
    else
      do while (A(p)>=y0)
        p = p - 1
        if (p<0) error stop 'p<0'
      end do
      do while (A(p)<y0)
        p = p + 1
        if (p>=n) error stop 'p>=n'
      end do
      p = p - 1
    end if
  end subroutine


  subroutine read_params_obj(obj_info)
    implicit none
    type(type_info)  ::  obj_info
    real(dp)   ::  Lx, Lz, Ly, nu, cond, Ret, Pr, dtmax, Sf_x, Sf_z, Sq, &
                   wall_BC_T_bot, wall_BC_T_top, &
                   wall_BC_U_bot, wall_BC_U_top, &
                   wall_BC_V_bot, wall_BC_V_top, &
                   wall_BC_W_bot, wall_BC_W_top,Gb_x,Gb_z,Gb_y, &
                   bulk_target_Reb, bulk_target_Ub
    integer  ::  nreps_x, nreps_z, nx, nz, ny, bulk_nprint, avg_iter_start, avg_nraw, avg_freq, avg_bin_nbins, &
                 snap_iter_start, snap_freq, filtering_steps, filtering_first, filtering_nstop, &
                 check_len_ref_angles, check_n_knots, use_rough_surf, &
                 use_incompressible, avg_apply_rms,use_utau_work, &
                 bulk_use_rho_weighting, bulk_use_target_Reb, bulk_use_target_Ub, &
                 force_Sf_x_accel, force_Sf_x_semilocal, force_Sf_x_pressure, prof_1d_avg_use, use_rough_surf_quady

    include 'params.dat'

    obj_info%nx      = nx
    obj_info%nz      = nz
    obj_info%ny      = ny
    obj_info%Lx      = Lx
    obj_info%Lz      = Lz
    obj_info%dx      = Lx / (nx + 0.0d0)
    obj_info%dz      = Lz / (nz + 0.0d0)

    if       ((obj_info%tag=='U').and.(obj_info%side=='bot')) then ; obj_info%wall_BC = wall_BC_U_bot ;
    else if  ((obj_info%tag=='V').and.(obj_info%side=='bot')) then ; obj_info%wall_BC = wall_BC_V_bot ;
    else if  ((obj_info%tag=='W').and.(obj_info%side=='bot')) then ; obj_info%wall_BC = wall_BC_W_bot ;
    else if  ((obj_info%tag=='T').and.(obj_info%side=='bot')) then ; obj_info%wall_BC = wall_BC_T_bot ;

    else if  ((obj_info%tag=='U').and.(obj_info%side=='top')) then ; obj_info%wall_BC = wall_BC_U_top ;
    else if  ((obj_info%tag=='V').and.(obj_info%side=='top')) then ; obj_info%wall_BC = wall_BC_V_top ;
    else if  ((obj_info%tag=='W').and.(obj_info%side=='top')) then ; obj_info%wall_BC = wall_BC_W_top ;
    else if  ((obj_info%tag=='T').and.(obj_info%side=='top')) then ; obj_info%wall_BC = wall_BC_T_top ;
    else if   (obj_info%tag=='P')                             then ; obj_info%wall_BC = -1e20                    ;
    else ; error stop 'Unrecognized (tag,side)' ; end if

  end subroutine

  function reduce_stencil(x_cen, z_cen, i, j, y_cen, obj_info, x_adim, x_shifts, z_adim, z_shifts, &
                         is_fit_h, deriv, dx, dz) result(val)
    real(dp)         ::  x_cen, z_cen, y_cen, val, x_adim(0:), z_adim(0:)
    integer          ::  x_shifts(0:), z_shifts(0:), i,j,i0,j0,k0, deriv ! deriv: [0: ddx, 1: ddy, 2: ddz, 3: val]
    type(type_info)  ::  obj_info
    real(dp)         ::  A(0:2,0:2), dx, dz  ! stencil, not always filled
    logical          ::  is_fit_h

    block
      do     i0=0,size(x_shifts,1)-1
        do   j0=0,size(z_shifts,1)-1
          if (is_fit_h) then
            A(i0,j0) = obj_info%var_slab_y(modulo(i+x_shifts(i0),obj_info%nx), modulo(j+z_shifts(j0),obj_info%nz))
          else
            block
              real(dp) :: y_01(0:1), p_01(0:1), coeffs_f(0:2)
              do k0=0,1
                p_01(k0) = obj_info%var_data_3d(modulo(i+x_shifts(i0),obj_info%nx), modulo(j+z_shifts(j0),obj_info%nz), k0)
                y_01(k0) = obj_info%var_y_3d   (modulo(i+x_shifts(i0),obj_info%nx), modulo(j+z_shifts(j0),obj_info%nz), k0)
              end do
              if (deriv ==1) then
                if (abs(y_01(1) - y_01(0))<1e-10) error stop 'abs(y_01(1) - y_01(0))<1e-10'
                A(i0,j0) = (p_01(1) - p_01(0))/(y_01(1) - y_01(0))
              else
                call interp_coeffs_2_or_3(y_cen, y_01, coeffs_f, 2, .false.)
                A(i0,j0) = p_01(0)*coeffs_f(0) + p_01(1)*coeffs_f(1)
              end if
            end block
          end if
        end do
      end do
      call quick_reduction(A, x_cen, x_adim, 1, (deriv == 0))
      call quick_reduction(A, z_cen, z_adim, 2, (deriv == 2))
      val = A(0,0)
      if((deriv == 0)) val = val/dx
      if((deriv == 2)) val = val/dz
    end block
  end function

  subroutine quick_reduction(A, x_cen, x_adim, ax, is_normal)
    real(dp)  ::  x_cen, x_adim(0:), A(0:2,0:2), B(0:2,0:2), coeffs_f(0:2) ! stencil, not always filled
    integer   ::  ax,p,n
    logical   ::  is_normal
    B        = 0  ;  coeffs_f = 0 ; n = size(x_adim,1)
    call interp_coeffs_2_or_3(x_cen, x_adim, coeffs_f, n, is_normal)
    do p=0,n-1
      if      (ax == 1) then ; B(0,:) = B(0,:) + coeffs_f(p)*A(p,:) ;
      else if (ax == 2) then ; B(:,0) = B(:,0) + coeffs_f(p)*A(:,p) ;
      else ; error stop 'Unrecognized ax' ; end if
    end do
    A = B
  end subroutine

  subroutine interp_coeffs_2_or_3(xc, x_coords, coeffs_f, n, is_normal)
    real(dp), intent(in) :: xc, x_coords(0:*)
    real(dp)             :: coeffs_f(0:2), dx, x0, x1, x2
    integer              :: n
    logical              :: is_normal
    coeffs_f = 0
    if ((n<2).or.(n>3)) error stop '((n<2).or.(n>3))'
              dx = (x_coords(n-1) - x_coords(0))/(n - 1.0d0)
              x0 = (x_coords(0)   - xc)/dx
              x1 = (x_coords(1)   - xc)/dx
    if (n==3) x2 = (x_coords(2)   - xc)/dx
    if (n==2) then
      if (.not.is_normal) then
        coeffs_f(0) = -x1/(x0 - x1)
        coeffs_f(1) =  x0/(x0 - x1)
      else
        coeffs_f(0) =  1_dp/(x0 - x1)
        coeffs_f(1) = -1_dp/(x0 - x1)
      end if
    else if (n==3) then
      block
        real(dp) :: denom
        denom = (x0**2*x1 - x0**2*x2 - x0*x1**2 + x0*x2**2 + x1**2*x2 - x1*x2**2)
        if (.not.is_normal) then
          coeffs_f(0) = (x1**2*x2 - x1*x2**2)/denom
          coeffs_f(1) = (-x0**2*x2 + x0*x2**2)/denom
          coeffs_f(2) = (x0**2*x1 - x0*x1**2)/denom
        else
          coeffs_f(0) = (-x1**2 + x2**2)/denom
          coeffs_f(1) = (x0**2 - x2**2)/denom
          coeffs_f(2) = (-x0**2 + x1**2)/denom
        end if
      end block
    end if
    if (is_normal) coeffs_f = coeffs_f/dx
  end subroutine

 ! import sympy as sp
 ! for is_quad in [True, False]:
 !     dofs          =  3 if is_quad else 2
 !     all_x, all_y  =  [sp.symbols(f'{c}0,{c}1,{c}2',real=True,positive=True,nonzero=True)[:dofs] for c in 'xy']
 !     a,b,c         =   sp.symbols('a,b,c',real=True,positive=True,nonzero=True)
 !     if not is_quad: a = 0
 !     f    =  lambda x: a*(x**2) + b*x + c
 !     eqs  =  [(f(x_) - y_) for x_,y_ in zip(all_x,all_y)]
 !     if is_quad:
 !         sol  =  sp.solve(eqs,a,b,c)
 !     else:
 !         sol  =  sp.solve(eqs,b,c)
 !     coeffs  = [sol[c].diff(y_) for y_ in all_y]
 !     normals = [sol[b].diff(y_) for y_ in all_y]
 !     print(dofs,coeffs,normals)
 ! #
 ! # 3 [(x1**2*x2 - x1*x2**2)/(x0**2*x1 - x0**2*x2 - x0*x1**2 + x0*x2**2 + x1**2*x2 - x1*x2**2)  ,
 ! #    (-x0**2*x2 + x0*x2**2)/(x0**2*x1 - x0**2*x2 - x0*x1**2 + x0*x2**2 + x1**2*x2 - x1*x2**2) ,
 ! #    (x0**2*x1 - x0*x1**2)/(x0**2*x1 - x0**2*x2 - x0*x1**2 + x0*x2**2 + x1**2*x2 - x1*x2**2)]
 ! #
 ! #   [(-x1**2 + x2**2)/(x0**2*x1 - x0**2*x2 - x0*x1**2 + x0*x2**2 + x1**2*x2 - x1*x2**2),
 ! #    (x0**2 - x2**2)/(x0**2*x1 - x0**2*x2 - x0*x1**2 + x0*x2**2 + x1**2*x2 - x1*x2**2),
 ! #    (-x0**2 + x1**2)/(x0**2*x1 - x0**2*x2 - x0*x1**2 + x0*x2**2 + x1**2*x2 - x1*x2**2)]
 ! #
 ! #
 ! # 2 [-x1/(x0 - x1),
 ! #     x0/(x0 - x1)]
 ! #
 ! #   [ 1/(x0 - x1),
 ! #    -1/(x0 - x1)]

  subroutine quick_write(iter_avg, A, nx, nz, tag, side, deriv_val, deriv_h)
    real(dp)          :: A(0:,0:)
    integer           :: nx, nz, iter_avg, deriv_val, deriv_h
    character(len=3)  :: key_val, key_h
    character(len=1)  :: tag
    character(len=3)  :: side
    character(len=45) :: my_file

    if      (deriv_val== 0) then ; key_val = 'ddx';
    else if (deriv_val== 1) then ; key_val = 'ddy';
    else if (deriv_val== 2) then ; key_val = 'ddz';
    else if (deriv_val== 3) then ; key_val = 'val';
    else if (deriv_val==-1) then ; key_val = '___';
    else ; error stop 'deriv_val out of range'; end if 

    if      (deriv_h==0) then ; key_h = 'n_x';
    else if (deriv_h==1) then ; key_h = 'n_y';
    else if (deriv_h==2) then ; key_h = 'n_z';
    else if (deriv_h==3) then ; key_h = 'avg';
    else ; error stop 'deriv_h out of range'; end if 
    
    write(my_file,'(A18,A1,A1,A3,A1,A3,A1,A3,A1,I0.9,A4)') './interp_2d/array_',tag,'_',side,&
                                                     '_',key_val,'_',key_h,'_',iter_avg,'.dat'
    open(19,file=my_file,action='write',status='replace',form='unformatted',access='stream')
      write(19) nx, nz
      block
        integer :: i,j
        do   j=0,nz-1
          do i=0,nx-1
            write(19) A(i,j)
          end do
        end do
      end block
    close(19)
    write(6,'(A21,A41,A1)') 'Done! (writing file: ', my_file, ')'
    flush(6)
  end subroutine
end module

program main
  use               mod_interp_routines
  use, intrinsic :: iso_fortran_env, only: i8 => int64, dp => real64

  implicit none

  integer , allocatable ::  x_shifts(:), z_shifts(:), all_deriv_val(:), all_deriv_h(:)
  real(dp), allocatable ::  x_adim(:)  , z_adim(:), arr_val(:,:), arr_ddx(:,:), arr_ddz(:,:), arr_ddy(:,:), my_arr(:,:)
  real(dp)              ::  tol_y
  type(type_info)       ::  obj_info
  integer               :: numel_gl, n_fits
  real(dp), allocatable :: c_gl(:), w_gl(:)

  call init_info_obj(obj_info)

  tol_y    =  0.9*obj_info%dy_ref
  numel_gl =  5

  call fill_3d_arrs(obj_info, tol_y)

  if      (obj_info%tag=='U') then ; allocate(x_adim(0:1),x_shifts(0:1)) ; x_adim = [-0.5, 0.5  ] ; x_shifts = [ 0, 1   ] ;
  else if (obj_info%tag=='W') then ; allocate(x_adim(0:2),x_shifts(0:2)) ; x_adim = [-1  , 0  ,1] ; x_shifts = [-1, 0, 1] ;
  else if (obj_info%tag=='V') then ; allocate(x_adim(0:2),x_shifts(0:2)) ; x_adim = [-1  , 0  ,1] ; x_shifts = [-1, 0, 1] ;
  else if (obj_info%tag=='T') then ; allocate(x_adim(0:2),x_shifts(0:2)) ; x_adim = [-1  , 0  ,1] ; x_shifts = [-1, 0, 1] ;
  else if (obj_info%tag=='P') then ; allocate(x_adim(0:2),x_shifts(0:2)) ; x_adim = [-1  , 0  ,1] ; x_shifts = [-1, 0, 1] ;
  else ; error stop 'Unrecognized (tag)' ; end if

  if      (obj_info%tag=='U') then ; allocate(z_adim(0:2),z_shifts(0:2)) ; z_adim = [-1  , 0  ,1] ; z_shifts = [-1, 0, 1] ;
  else if (obj_info%tag=='W') then ; allocate(z_adim(0:1),z_shifts(0:1)) ; z_adim = [-0.5, 0.5  ] ; z_shifts = [ 0, 1   ] ;
  else if (obj_info%tag=='V') then ; allocate(z_adim(0:2),z_shifts(0:2)) ; z_adim = [-1  , 0  ,1] ; z_shifts = [-1, 0, 1] ;
  else if (obj_info%tag=='T') then ; allocate(z_adim(0:2),z_shifts(0:2)) ; z_adim = [-1  , 0  ,1] ; z_shifts = [-1, 0, 1] ;
  else if (obj_info%tag=='P') then ; allocate(z_adim(0:2),z_shifts(0:2)) ; z_adim = [-1  , 0  ,1] ; z_shifts = [-1, 0, 1] ;
  else ; error stop 'Unrecognized (tag)' ; end if

  allocate(my_arr(0:obj_info%nx-1,0:obj_info%nz-1))

  if      (obj_info%tag=='U') then ; n_fits =  9 ;
  else if (obj_info%tag=='V') then ; n_fits =  9 ;
  else if (obj_info%tag=='W') then ; n_fits =  9 ;
  else if (obj_info%tag=='T') then ; n_fits =  7 ;
  else if (obj_info%tag=='P') then ; n_fits = 10 ;
  else ; error stop 'Unrecognized (tag)' ; end if

  allocate(all_deriv_val(0:n_fits-1), all_deriv_h(0:n_fits-1))

  if      (obj_info%tag=='U') then; all_deriv_val(0:n_fits-1) = [            0, 1, 2, 3, 0, 1, 1, 2, 2];
  else if (obj_info%tag=='V') then; all_deriv_val(0:n_fits-1) = [            0, 1, 2, 3, 0, 0, 1, 2, 2];
  else if (obj_info%tag=='W') then; all_deriv_val(0:n_fits-1) = [            0, 1, 2, 3, 0, 0, 1, 1, 2];
  else if (obj_info%tag=='T') then; all_deriv_val(0:n_fits-1) = [            0, 1, 2, 3, 0, 1, 2]      ;
  else if (obj_info%tag=='P') then; all_deriv_val(0:n_fits-1) = [-1, -1, -1, 0, 1, 2, 3, 3, 3, 3]      ;
  else ; error stop 'Unrecognized (tag)' ; end if
  
  if      (obj_info%tag=='U') then; all_deriv_h(0:n_fits-1) = [         3, 3, 3, 3, 0, 0, 1, 0, 2];
  else if (obj_info%tag=='V') then; all_deriv_h(0:n_fits-1) = [         3, 3, 3, 3, 0, 1, 1, 1, 2];
  else if (obj_info%tag=='W') then; all_deriv_h(0:n_fits-1) = [         3, 3, 3, 3, 0, 2, 1, 2, 2];
  else if (obj_info%tag=='T') then; all_deriv_h(0:n_fits-1) = [         3, 3, 3, 3, 0, 1, 2];
  else if (obj_info%tag=='P') then; all_deriv_h(0:n_fits-1) = [0, 1, 2, 3, 3, 3, 3, 0, 1, 2];
  else ; error stop 'Unrecognized (tag)' ; end if 

  ! U [3, 0, 1, 1, 2, 2] [3, 0, 0, 1, 0, 2]
  ! V [3, 0, 0, 1, 2, 2] [3, 0, 1, 1, 1, 2]
  ! W [3, 0, 0, 1, 1, 2] [3, 0, 2, 1, 2, 2]
  ! T [3, 0, 1, 2] [3, 0, 1, 2]
  ! P [3, 3, 3, 3] [3, 0, 1, 2]

  ! import sympy as sp
  ! U = [[sp.symbols(f"U_{i}_{j}",real=True,positive=True,nonzero=True) for j in range(3)] for i in range(3)]
  ! N =  [sp.symbols(f"N_{j}"    ,real=True,positive=True,nonzero=True) for j in range(3)]
  ! T =  [sp.symbols(f"T_{j}"    ,real=True,positive=True,nonzero=True) for j in range(3)]
  ! P =  sp.symbols(f"P"    ,real=True,positive=True,nonzero=True)
  ! tensor = [[(U[i][j]+U[j][i]) for j in range(3)] for i in range(3)]
  ! F      = [(sum(N[j]*tensor[i][j] for j in range(3))+N[i]*P) for i in range(3)]
  ! Q      = [sum(N[j]*T[j] for j in range(3))]
  ! total = sum(F) + sum(Q)
  ! for tag in 'UVWTP':
  !     derivs = [3]
  !     hgrads = [3]
  !     if tag in 'UVW':
  !         var_ = U['UVW'.index(tag)] + [0]
  !     elif tag == 'T':
  !         var_ = T + [0]
  !     elif tag == 'P':
  !         var_ = [0,0,0] + [P]
  !     for i_grad in range(4):
  !         for i_normal in range(3):
  !             result = total + 0
  !             if var_[i_grad]:
  !                 result = result.diff(var_[i_grad] ).diff(N[i_normal])
  !                 if result:
  !                     derivs.append(i_grad)
  !                     hgrads.append(i_normal)
  !     print(tag,derivs,hgrads)

  allocate(c_gl(0:numel_gl-1), w_gl(0:numel_gl-1))
  if      ((numel_gl) == 1) then 
    c_gl(0:numel_gl-1) = [0.0]
    w_gl(0:numel_gl-1) = [2.0]
  else if ((numel_gl) == 2) then
    c_gl(0:numel_gl-1) = [-0.5773502691896257, 0.5773502691896257]
    w_gl(0:numel_gl-1) = [1.0, 1.0]   
  else if ((numel_gl) == 3) then
    c_gl(0:numel_gl-1) = [-0.7745966692414834, 0.0, 0.7745966692414834]
    w_gl(0:numel_gl-1) = [0.5555555555555557, 0.8888888888888888, 0.5555555555555557]
  else if ((numel_gl) == 4) then
    c_gl(0:numel_gl-1) = [-0.8611363115940526, -0.33998104358485626, 0.33998104358485626, 0.8611363115940526]
    w_gl(0:numel_gl-1) = [0.3478548451374537, 0.6521451548625462, 0.6521451548625462, 0.3478548451374537]
  else if ((numel_gl) == 5) then
    c_gl(0:numel_gl-1) = [-0.906179845938664, -0.5384693101056831, 0.0, 0.5384693101056831, 0.906179845938664]
    w_gl(0:numel_gl-1) = [0.23692688505618942, 0.4786286704993662, 0.568888888888889, 0.4786286704993662, 0.23692688505618942]
  else if ((numel_gl) == 6) then
    c_gl(0:numel_gl-1) = [-0.932469514203152, -0.6612093864662645, -0.23861918608319693, 0.23861918608319693, &
                           0.6612093864662645, 0.932469514203152]
    w_gl(0:numel_gl-1) = [0.17132449237916975, 0.36076157304813894, 0.46791393457269137, 0.46791393457269137, &
                          0.36076157304813894, 0.17132449237916975]
  else if ((numel_gl) == 7) then
    c_gl(0:numel_gl-1) = [-0.9491079123427585, -0.7415311855993945, -0.4058451513773972, 0.0, 0.4058451513773972, &
                           0.7415311855993945, 0.9491079123427585]
    w_gl(0:numel_gl-1) = [0.12948496616887065, 0.2797053914892766, 0.3818300505051183, 0.41795918367346896, &
                          0.3818300505051183, 0.2797053914892766, 0.12948496616887065]
  else if ((numel_gl) == 8) then
    c_gl(0:numel_gl-1) = [-0.9602898564975362, -0.7966664774136267, -0.525532409916329, -0.18343464249564978, &
                           0.18343464249564978, 0.525532409916329, 0.7966664774136267, 0.9602898564975362]
    w_gl(0:numel_gl-1) = [0.10122853629037669, 0.22238103445337434, 0.31370664587788705, 0.36268378337836177, &
                          0.36268378337836177, 0.31370664587788705, 0.22238103445337434, 0.10122853629037669]
  else if ((numel_gl) == 9) then
    c_gl(0:numel_gl-1) = [-0.9681602395076261, -0.8360311073266358, -0.6133714327005904, -0.3242534234038089, 0.0, &
                           0.3242534234038089, 0.6133714327005904, 0.8360311073266358, 0.9681602395076261]
    w_gl(0:numel_gl-1) = [0.08127438836157472, 0.18064816069485712, 0.26061069640293566, 0.3123470770400028, &
                          0.33023935500125967, 0.3123470770400028, 0.26061069640293566, 0.18064816069485712, 0.08127438836157472]
  else if ((numel_gl) == 10) then
    c_gl(0:numel_gl-1) = [-0.9739065285171717, -0.8650633666889845, -0.6794095682990244, -0.4333953941292472, &
                          -0.14887433898163122, 0.14887433898163122, 0.4333953941292472, 0.6794095682990244, &
                           0.8650633666889845, 0.9739065285171717]
    w_gl(0:numel_gl-1) = [0.06667134430868807, 0.14945134915058036, 0.219086362515982, 0.2692667193099965, &
                          0.295524224714753, 0.295524224714753, 0.2692667193099965, 0.219086362515982, &
                          0.14945134915058036, 0.06667134430868807]
  else 
    error stop 'numel_gl out of range'
  end if 

  block 
    integer :: ii_fit,deriv_val,deriv_h,K_dir
    
    if      (obj_info%side == 'bot') then ; K_dir =  1 ; 
    else if (obj_info%side == 'top') then ; K_dir = -1 ; 
    else ; error stop 'Unrecognized side' ; end if

    do ii_fit=0,n_fits-1
      deriv_val = all_deriv_val(ii_fit)
      deriv_h   = all_deriv_h  (ii_fit)
      block
        integer  :: i,j
        do   j=0,obj_info%nz-1
          do i=0,obj_info%nx-1
            block 
              integer   ::  i_gl, j_gl 
              real(dp)  ::  total_loc
              total_loc  = 0
              do   i_gl=0,numel_gl-1
                do j_gl=0,numel_gl-1
                  block 
                    real(dp) :: val_loc, n_comp_loc, h_loc
                    h_loc = reduce_stencil( c_gl(i_gl)/2. , c_gl(j_gl)/2. , i, j, 0. , obj_info, &
                                           x_adim, x_shifts, z_adim, z_shifts, .true., 3, obj_info%dx, obj_info%dz)
                    if (deriv_val==(-1)) then
                        val_loc = 1
                    else
                        val_loc = reduce_stencil( c_gl(i_gl)/2. , c_gl(j_gl)/2. , i, j, h_loc, obj_info, &
                                         x_adim, x_shifts, z_adim, z_shifts, .false., deriv_val, obj_info%dx, obj_info%dz)
                    end if 
                    if (deriv_h==1) then 
                      n_comp_loc   = K_dir
                    else if (deriv_h == 3) then 
                      n_comp_loc   = 1
                    else if ((deriv_h == 0).or.(deriv_h == 2)) then 
                      n_comp_loc   = -K_dir*reduce_stencil( c_gl(i_gl)/2. , c_gl(j_gl)/2. , i, j, 0. , obj_info, &
                                      x_adim, x_shifts, z_adim, z_shifts, .true., deriv_h, obj_info%dx, obj_info%dz)
                    else 
                      error stop 'deriv_h out of range'
                    end if
                    total_loc = total_loc + val_loc*n_comp_loc*w_gl(i_gl)*w_gl(j_gl)/4_dp
                  end block 
                end do 
              end do
              my_arr(i,j) =  total_loc
            end block 
          end do
        end do
      end block
      call quick_write(obj_info%iter_avg, my_arr, obj_info%nx, obj_info%nz, obj_info%tag, obj_info%side, deriv_val, deriv_h)
    end do 
  end block 
end program

