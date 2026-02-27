
module mod_interp_routines
  use, intrinsic :: iso_fortran_env, only: i8 => int64, dp => real64

  implicit none

  type type_extended_obj_info
    real(dp)  ::  Lx, Lz, Ly, dx, dz, &
                  wall_BC_mu_bot, wall_BC_cond_bot , wall_BC_T_bot, &
                  wall_BC_mu_top, wall_BC_cond_top , wall_BC_T_top
    integer   ::  nx, nz, ny, iter_avg
  end type

  type type_bulk
    real(dp)  ::  U_bulk, V_bulk, W_bulk, T_bulk, enth_bulk, rho_bulk, mu_bulk, cond_bulk, cp_bulk
  end type

  contains

  subroutine init_extended_obj_info(obj_info)
    implicit none
    type(type_extended_obj_info)  ::  obj_info
    real(dp)  ::  Lx, Lz, Ly, Sf_x, Sf_z, Sq, &
                  wall_BC_U_bot, wall_BC_V_bot, wall_BC_W_bot, wall_BC_T_bot, &
                  wall_BC_U_top, wall_BC_V_top, wall_BC_W_top, wall_BC_T_top, &
                  wall_BC_rho_bot  , wall_BC_mu_bot   , wall_BC_cond_bot , &
                  wall_BC_rho_top  , wall_BC_mu_top   , wall_BC_cond_top
    integer   ::  nx, nz, ny, nreps_x, nreps_z

    include 'extended_params.dat'

    obj_info%Lx   = Lx
    obj_info%Lz   = Lz
    obj_info%Ly   = Ly

    obj_info%wall_BC_mu_bot   = wall_BC_mu_bot
    obj_info%wall_BC_cond_bot = wall_BC_cond_bot
    obj_info%wall_BC_T_bot    = wall_BC_T_bot

    obj_info%wall_BC_mu_top   = wall_BC_mu_top
    obj_info%wall_BC_cond_top = wall_BC_cond_top
    obj_info%wall_BC_T_top    = wall_BC_T_top

    obj_info%nx = nx
    obj_info%nz = nz
    obj_info%ny = ny

    obj_info%dx = obj_info%Lx / (obj_info%nx + 0.d0)
    obj_info%dz = obj_info%Lz / (obj_info%nz + 0.d0)

    write(6,*) ''
    write(6,*) '---------------------- Read Parameters (extended_params.dat) ----------------------'
    write(6,'(A30,A3,E25.16)') 'Lx'              , ' = ' , obj_info%Lx
    write(6,'(A30,A3,E25.16)') 'Lz'              , ' = ' , obj_info%Lz
    write(6,'(A30,A3,E25.16)') 'Ly'              , ' = ' , obj_info%Ly

    write(6,'(A30,A3,E25.16)') 'dx'              , ' = ' , obj_info%dx
    write(6,'(A30,A3,E25.16)') 'dz'              , ' = ' , obj_info%dz

    write(6,'(A30,A3,E25.16)') 'wall_BC_mu_bot'  , ' = ' , obj_info%wall_BC_mu_bot
    write(6,'(A30,A3,E25.16)') 'wall_BC_cond_bot', ' = ' , obj_info%wall_BC_cond_bot

    write(6,'(A30,A3,E25.16)') 'wall_BC_mu_top'  , ' = ' , obj_info%wall_BC_mu_top
    write(6,'(A30,A3,E25.16)') 'wall_BC_cond_top', ' = ' , obj_info%wall_BC_cond_top

    write(6,'(A30,A3,I25)')    'nx'              , ' = ' , obj_info%nx
    write(6,'(A30,A3,I25)')    'nz'              , ' = ' , obj_info%nz
    write(6,'(A30,A3,I25)')    'ny'              , ' = ' , obj_info%ny
    flush(6)
  end subroutine

  subroutine quick_write(A, key, side, obj_info)
    real(dp)                      ::  A(0:,0:)
    type(type_extended_obj_info)  ::  obj_info
    character(len=15)             ::  key
    character(len=3)              ::  side
    character(len=64)             ::  my_file
    write(my_file,'(A17,A15,A6,A3,A10,I0.9,A4)') './interp_results/',key,'_side_',side,'_iter_avg_',obj_info%iter_avg,'.dat'
    open(19,file=my_file,action='write',status='replace',form='unformatted',access='stream')
      write(19) obj_info%nx, obj_info%nz
      block 
        integer :: i,j
        do   j=0,obj_info%nz-1 
          do i=0,obj_info%nx-1
            write(19) A(i,j)
          end do
        end do
      end block
    close(19)
    write(6,*) 'Done! (writing file: ', my_file, ')'
    flush(6)
  end subroutine
 
  subroutine load_arr(A,tag,i_val,i_h,side,obj_info, all_A, ncols, order)
    real(dp)                      ::  A(0:,0:), all_A(0:,0:,0:)
    character(len=128)            ::  order(0:)
    type(type_extended_obj_info)  ::  obj_info
    integer                       ::  i, j, nx_check, nz_check, i_val, i_h, ncols
    character(len=3)              ::  key_val, key_h, side
    character(len=1)              ::  tag
    character(len=45)             ::  my_file
    character(len=9)              ::  order_k
    
    A = 0 ! will be re-filled

    if      (i_val== 0) then ; key_val = 'ddx';
    else if (i_val== 1) then ; key_val = 'ddy';
    else if (i_val== 2) then ; key_val = 'ddz';
    else if (i_val== 3) then ; key_val = 'val';
    else if (i_val==-1) then ; key_val = '___';
    else ; error stop 'i_val out of range'; end if 

    if      (i_h==0) then ; key_h = 'n_x';
    else if (i_h==1) then ; key_h = 'n_y';
    else if (i_h==2) then ; key_h = 'n_z';
    else if (i_h==3) then ; key_h = 'avg';
    else ; error stop 'i_h out of range'; end if 

    write(my_file,'(A18,A1,A1,A3,A1,A3,A1,A3,A1,I0.9,A4)') './interp_2d/array_',tag,'_',side,'_',&
                                                            key_val,'_',key_h,'_',obj_info%iter_avg,'.dat'
    
    write(6,*) '    Reading array (load_arr): ', my_file ; flush(6)
    block 
      write(order_k,'(A1,A1,A3,A1,A3)') tag, '_', key_val, '_', key_h
      order(ncols) = trim(order_k)
    end block 
    open(19,file=my_file, action='read', status='old', form='unformatted', access='stream')
      read(19) nx_check, nz_check
      if (nx_check.ne.obj_info%nx) error stop 'nx_check.ne.obj_info%nx'
      if (nz_check.ne.obj_info%nz) error stop 'nz_check.ne.obj_info%nz'
      do   j=0,obj_info%nz-1
        do i=0,obj_info%nx-1
          read(19) A(i,j)
          all_A(i,j,ncols)  =  A(i,j)
        end do 
      end do 
    close(19)
    ncols = ncols + 1
    write(6,'(A37,A9,A1,E25.16)') '        -> Variable: (per unit area) ', order_k, ' ', &
                                  sum(A)/(obj_info%nx*obj_info%nz+0.d0) ; flush(6)
  end subroutine

  subroutine write_tec(fname, all_A, nx, nz, ncols, order)
    character(len=*)   ::  fname
    real(dp)           ::  all_A(0:,0:,0:)
    integer            ::  nx, nz, ncols, i, j, k
    character(len=128) ::  order(0:)
    character(len=128) ::  order_k
    character(len=7)   ::  fmt_order_7
    character(len=8)   ::  fmt_order_8
    character(len=9)   ::  fmt_order_9
    character(len=3)   ::  aux_end_3
    open(19,file=trim(fname))
      write(19,'(A13)', advance="no") 'VARIABLES = "'
      do k=0,ncols-1
        order_k = trim(order(k))
        aux_end_3                  = '"  '
        if (k<(ncols-1)) aux_end_3 = '","'
        if      (len_trim(order_k)<10) then 
          write(fmt_order_7,'(A2,I1,A4)') '(A', len_trim(order_k), ',A3)'
          write(19,fmt_order_7, advance="no") trim(order_k),aux_end_3
        else if (len_trim(order_k)<100) then 
          write(fmt_order_8,'(A2,I2,A4)') '(A', len_trim(order_k), ',A3)'
          write(19,fmt_order_8, advance="no") trim(order_k),aux_end_3
        else if (len_trim(order_k)<1000) then 
          write(fmt_order_9,'(A2,I3,A4)') '(A', len_trim(order_k), ',A3)'
          write(19,fmt_order_9, advance="no") trim(order_k),aux_end_3
        else  ; error stop 'len_trim(order_k)>=1000' ; end if
      end do
      write(19,'(A1)') ' '
      write(19,'(A7,I9,A3,I9,A8)') 'ZONE I=',nz,' K=',nx,' F=POINT'
      do     i=0,nx-1
        do   j=0,nz-1
          do k=0,ncols-1
            write(19,'(E25.16,A1)', advance="no") all_A(i,j,k),' '
          end do
          write(19,'(A1)') ' '
        end do 
      end do 
    close(19)
  end subroutine

  subroutine fill_bulk_obj(obj_bulk, iter_avg)
    type(type_bulk)  ::  obj_bulk
    integer          ::  iter_avg
    call read_bulk_value(obj_bulk%U_bulk   , 'U', iter_avg, 0, ' ! U________bulk (weighted by frac_nonempty(:) and dy(:))')
    call read_bulk_value(obj_bulk%V_bulk   , 'V', iter_avg, 0, ' ! V________bulk (weighted by frac_nonempty(:) and dy(:))')
    call read_bulk_value(obj_bulk%W_bulk   , 'W', iter_avg, 0, ' ! W________bulk (weighted by frac_nonempty(:) and dy(:))')
    call read_bulk_value(obj_bulk%T_bulk   , 'T', iter_avg, 0, ' ! T________bulk (weighted by frac_nonempty(:) and dy(:))')
    call read_bulk_value(obj_bulk%enth_bulk, 'T', iter_avg, 1, ' ! enth    _bulk (weighted by frac_nonempty(:) and dy(:))')
    call read_bulk_value(obj_bulk%rho_bulk , 'T', iter_avg, 2, ' ! rho     _bulk (weighted by frac_nonempty(:) and dy(:))')
    call read_bulk_value(obj_bulk%mu_bulk  , 'T', iter_avg, 3, ' ! mu      _bulk (weighted by frac_nonempty(:) and dy(:))')
    call read_bulk_value(obj_bulk%cond_bulk, 'T', iter_avg, 4, ' ! cond    _bulk (weighted by frac_nonempty(:) and dy(:))')
    call read_bulk_value(obj_bulk%cp_bulk  , 'T', iter_avg, 5, ' ! cp      _bulk (weighted by frac_nonempty(:) and dy(:))')
  end subroutine 

  subroutine read_bulk_value(bulk_val, tag, iter_avg, i_read, check_tag)
    real(dp)           ::  bulk_val
    character(len=1)   ::  tag
    integer            ::  iter_avg, i_read, i, n
    character(len=57)  ::  fname_avg, check_tag, check_tag2
    write(fname_avg,'(A37,A1,A6,I0.9,A4)') '../../prof_1d_from_avg3d/prof_1d_avg_',tag,'_iter_',iter_avg,'.dat'
    open(19,file=fname_avg)
      read(19,*) n
      do i=1,(n+i_read) ; read(19,*) ; end do 
      read(19,'(E25.16,A57)') bulk_val, check_tag2
      if (check_tag.ne.check_tag2) error stop 'check_tag.ne.check_tag2'
    close(19)
  end subroutine

end module

program main
  use               mod_interp_routines
  use, intrinsic :: iso_fortran_env, only: i8 => int64, dp => real64

  implicit none

  real(dp), allocatable        ::  total_arr(:,:)  , temp_arr(:,:), all_A(:,:,:)
  type(type_extended_obj_info) ::  obj_info
  type(type_bulk)              ::  obj_bulk
  character(len=128), allocatable ::  order(:)
  integer :: ncols_max 
  
  ! -> Initialize object
  call init_extended_obj_info(obj_info)

  ! -> Fill iter_avg
  block 
    character(len=30) temp_chr
    call get_command_argument(1, temp_chr)
    read(temp_chr,*) obj_info%iter_avg
  end block 

  ! -> Allocate arrays
  allocate(total_arr (0:obj_info%nx-1,0:obj_info%nz-1))
  allocate(temp_arr , mold = total_arr)

  ncols_max = 20
  allocate(all_A (0:obj_info%nx-1,0:obj_info%nz-1,0:ncols_max-1), &
           order(0:ncols_max-1))
  order(:) = ' '
  order(0) = 'X' ; order(1) = 'Y' ; order(2) = 'Z'


  ! -> Write coordinates xp/zp
  open(19,file='coords_xp.dat')
    block 
      integer :: i
      write(19,'(I9)') obj_info%nx
      do  i=0,obj_info%nx-1
        write(19,'(E25.16)') (i+0.5d0)*obj_info%dx
        all_A(i,:,0) = (i+0.5d0)*obj_info%dx
      end do
    end block
  close(19)

  open(19,file='coords_zp.dat')
    block 
      integer :: j
      write(19,'(I9)') obj_info%nz
      do  j=0,obj_info%nz-1
        write(19,'(E25.16)') (j+0.5d0)*obj_info%dz
        all_A(:,j,2) = (j+0.5d0)*obj_info%dz
      end do
    end block
  close(19)

  call fill_bulk_obj(obj_bulk, obj_info%iter_avg)

  ! -> Loop over sides
  block 
    integer          ::  ii_side
    character(len=3) ::  side
    real(dp)         ::  mu_side, cond_side, T_side

  do ii_side=0,1
    ! -> Get properties [mu,cond]
    if      (ii_side == 0) then; side = 'bot'; mu_side = obj_info%wall_BC_mu_bot ; cond_side = obj_info%wall_BC_cond_bot ;
    else if (ii_side == 1) then; side = 'top'; mu_side = obj_info%wall_BC_mu_top ; cond_side = obj_info%wall_BC_cond_top ;
    else ; error stop 'Unrecognized ii_side' ; end if 

    if      (ii_side == 0) then; T_side = obj_info%wall_BC_T_bot ;
    else if (ii_side == 1) then; T_side = obj_info%wall_BC_T_top ;
    else ; error stop 'Unrecognized ii_side' ; end if 
    
    block
      character(len=31)  ::  myfile
      integer            :: i,j
      write(myfile, '(A24,A3,A4)') './surfs_hxz/surf_hxz_pp_',side,'.dat'
      open(19,file=myfile, action='read', status='old', form='unformatted', access='stream')
        read(19) obj_info%nx, obj_info%nz
        do   j=0,obj_info%nz-1
          do i=0,obj_info%nx-1
            read(19) all_A(i,j,1)
          end do 
        end do 
      close(19)
    end block

    write(6,'(A24,A3,A13,E25.16,A11,E25.16,A10,E25.16)')  '---> Interpolating Side ', side, &
                                                          ' ; cond_side ', cond_side, ' ; mu_side ', mu_side, ' ; T_side ', T_side

    block 
      integer :: i,j,ncols
      character(len=1) :: Ui, Uj
      do i=0,2
        if (i==0) Ui = 'U'
        if (i==1) Ui = 'V'
        if (i==2) Ui = 'W'
        total_arr = 0
        ncols    = 4
        order(:) = ' '
        order(0) = 'X' ; order(1) = 'Y' ; order(2) = 'Z'
        ! ------------ static pressure correction ------------
        block
          real(dp) :: res_normal, P_avg
          call load_arr(temp_arr ,'P',-1,i, side, obj_info, all_A, ncols, order)
          call load_arr(total_arr,'P',3,3, side, obj_info, all_A, ncols, order)
          res_normal =  sum(temp_arr)/(obj_info%nx*obj_info%nz + 0.d0)
          P_avg      =  sum(total_arr)/(obj_info%nx*obj_info%nz + 0.d0)
          temp_arr   =  temp_arr - res_normal 
          total_arr  =  temp_arr*P_avg
          write(6,'(A35,I1,A18,E25.16,A11,E25.16,A1)')  '--> Static Pressure Correction: (F[',i,&
                                                        ']): (res_normal = ', res_normal, ') (P_avg = ', P_avg, ')'
          flush(6)
        end block
        do j=0,2
          if (j==0) Uj = 'U'
          if (j==1) Uj = 'V'
          if (j==2) Uj = 'W'
          if (i==j) then 
            call load_arr(temp_arr,'P',3,i, side, obj_info, all_A, ncols, order) ; total_arr = total_arr           - temp_arr
            call load_arr(temp_arr,Ui ,i,i, side, obj_info, all_A, ncols, order) ; total_arr = total_arr + 2*mu_side*temp_arr
          else 
            call load_arr(temp_arr,Ui ,j,j, side, obj_info, all_A, ncols, order) ; total_arr = total_arr +   mu_side*temp_arr
            call load_arr(temp_arr,Uj ,i,j, side, obj_info, all_A, ncols, order) ; total_arr = total_arr +   mu_side*temp_arr
          end if 
        end do
        write(6,'(A21,I1,A19,A3,A1,E25.16)')  '--> Summary: Total F[',i,']: (per unit area) ', side, ' ', &
                                       sum(total_arr)/(obj_info%nx*obj_info%nz+0.d0)
        block 
          character(len= 1) :: xi
          character(len=15) :: fname_total_force
          character(len=59) :: fname_tec
          if (i==0) xi = 'x'
          if (i==1) xi = 'y'
          if (i==2) xi = 'z'
          write(fname_total_force,'(A12,A1,A2)') 'total_force_',xi,'__'
          call quick_write(total_arr, fname_total_force, side, obj_info)
          write(fname_tec,'(A31,A1,A1,A3,A10,I0.9,A4)') './interp_results/summary_force_', xi, '_', side, &
                                                        '_iter_avg_', obj_info%iter_avg,'.tec'
          if (ncols>=ncols_max) error stop 'ncols>=ncols_max'
          all_A(:,:,3)   =  total_arr
          order(3)  =  trim(fname_total_force)

          total_arr = total_arr/(0.5_dp*obj_bulk%rho_bulk*(obj_bulk%U_bulk**2)) ! Cf
          write(fname_total_force,'(A8,A1,A1,A3,A2)') 'Cf_bulk_',xi,'_',side,'__'
          all_A(:,:,ncols)  =  total_arr
          order(ncols)      =  trim(fname_total_force)
          ncols             =  ncols + 1
          write(6,'(A22,I1,A19,A3,A1,E25.16)')  '--> Summary: Total Cf[',i,']: (per unit area) ', side, ' ', &
                                                       sum(total_arr)/(obj_info%nx*obj_info%nz+0.d0)
          call write_tec(fname_tec, all_A, obj_info%nx, obj_info%nz, ncols, order)
        end block
      end do 
    end block 

    block 
      integer :: i,ncols
      total_arr = 0
      ncols    = 4
      order(:) = ' '
      order(0) = 'X' ; order(1) = 'Y' ; order(2) = 'Z'
      do i=0,2
        call load_arr(temp_arr,'T',i,i, side, obj_info, all_A, ncols, order) ; total_arr = total_arr - cond_side*temp_arr
      end do 
      call quick_write(total_arr, 'total_heat_flux', side, obj_info)
      write(6,'(A38,A3,A1,E25.16)')  '--> Summary: Total Q  (per unit area) ',side,' ', &
                                     sum(total_arr)/(obj_info%nx*obj_info%nz+0.d0)
      block 
        character(len=61) :: fname_tec
        character(len=6)  :: fname_Nu
        real(dp)          :: sign_side
        write(fname_tec,'(A35,A3,A10,I0.9,A4)') './interp_results/summary_heat_flux_', side, '_iter_avg_', &
                                                obj_info%iter_avg,'.tec'
        if (ncols>=ncols_max) error stop 'ncols>=ncols_max'
        all_A(:,:,3)  =  total_arr
        order(3)  =  'heat_flux' 
        sign_side = 1
        if (sum(total_arr)<0) sign_side = -1
        total_arr = total_arr/(obj_bulk%cond_bulk*abs(obj_bulk%T_bulk-T_side)/obj_info%Ly)*sign_side ! Nu
        write(fname_Nu,'(A3,A3)') 'Nu_',side
        all_A(:,:,ncols)  =  total_arr
        order(ncols)      =  fname_Nu
        ncols             =  ncols + 1
        write(6,'(A38,A3,A1,E25.16)')  '--> Summary: Total Nu (per unit area) ', side, ' ', &
                                       sum(total_arr)/(obj_info%nx*obj_info%nz+0.d0)
            
        call write_tec(fname_tec, all_A, obj_info%nx, obj_info%nz, ncols, order)
      end block
    end block 
  end do 
  end block
end program

