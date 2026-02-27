program main 
  use Mod_Height_Function, only: init_surf_obj, type_surf_obj, get_h
  use ibm_coeffs_engine
  use, intrinsic :: iso_fortran_env, only: i8 => int64, dp => real64
  implicit none
  type(type_surf_obj)  ::  obj
  integer              ::  i, j, k, i2, ii, jj, Lmax, nx, nz, ny, inow, cc_x, cc_z
  real(dp) , allocatable ::  x(:), z(:), y(:), all_ca(:), all_cb(:), all_cw(:), hxz_bot(:,:), hxz_top(:,:), buffer_c(:)
  real(dp)               ::  time_start, time_finish, bc_w_bot, bc_w_top
  integer, allocatable ::  all_ijk_t(:,:), all_ijk_p1(:,:), all_ijk_p2(:,:), slab_bot(:,:), slab_top(:,:),buffer_ijk(:,:)
  logical, allocatable ::      is_ok(:,:,:)
  logical              ::  use_mode_ghost
  character(len=1)          :: tag
  call cpu_time(time_start)
  call init_surf_obj(obj)

  open(19,file='ibm_params.txt')
    block 
     integer :: int_mode_quady
      read(19,*) tag, nx, nz, ny, cc_x, cc_z, int_mode_quady
      read(19,*) bc_w_bot, bc_w_top
      if ((int_mode_quady<0).or.(int_mode_quady>1)) error stop '(int_mode_quady<0).or.(int_mode_quady>1)'
      use_mode_ghost = (int_mode_quady==0)
    end block 
    allocate(x(0:nx-1), z(0:nz-1), y(0:ny-1)   ,&
             hxz_bot(0:nx-1,0:nz-1)            ,&
             hxz_top(0:nx-1,0:nz-1)            ,&
             slab_bot(0:nx*cc_x-1,0:nz*cc_z-1) ,&
             slab_top(0:nx*cc_x-1,0:nz*cc_z-1) )
    if (use_mode_ghost) allocate(is_ok(0:nx-1,0:nz-1,0:ny-1))
    do i=0,nx-1 ; read(19,*) x(i) ; end do
    do j=0,nz-1 ; read(19,*) z(j) ; end do
    do k=0,ny-1 ; read(19,*) y(k) ; end do
  close(19)

  Lmax = nx*nz*ny
  do     i=0,nx-1 
    do   j=0,nz-1 
      hxz_bot(i,j) = get_h(x(i),z(j),obj,.true. )
      hxz_top(i,j) = get_h(x(i),z(j),obj,.false.)
      ! do k=0,ny-1
      !   if (((hxz_bot(i,j)-1e-10)<=y(k)).and.(y(k)<=(hxz_top(i,j)+1e-10))) Lmax = Lmax + 1 ! upper bound
      ! end do 
    end do 
  end do

  allocate(all_ijk_t (0:Lmax*cc_x*cc_z-1,0:3),&
           all_ijk_p1(0:Lmax*cc_x*cc_z-1,0:2),&
           all_ijk_p2(0:Lmax*cc_x*cc_z-1,0:2),&
               all_ca(0:Lmax*cc_x*cc_z-1    ),&
               all_cb(0:Lmax*cc_x*cc_z-1    ),&
               all_cw(0:Lmax*cc_x*cc_z-1    ) )
  inow = 0
  if (use_mode_ghost) then 
    call make_rough_surface_ghost(x, z, y, slab_bot, all_ijk_t, all_ijk_p1, all_ijk_p2, all_ca, all_cb, all_cw, hxz_bot, &
                                  is_ok, obj, .true. , inow, Lmax, tag, bc_w_bot)
    call make_rough_surface_ghost(x, z, y, slab_top, all_ijk_t, all_ijk_p1, all_ijk_p2, all_ca, all_cb, all_cw, hxz_top, &
                                  is_ok, obj, .false., inow, Lmax, tag, bc_w_top)
  else 
    call make_rough_surface_quady(x, z, y, slab_bot, all_ijk_t, all_ijk_p1, all_ijk_p2, all_ca, all_cb, all_cw, hxz_bot, &
                                  .true. , inow, Lmax, tag, bc_w_bot)
    call make_rough_surface_quady(x, z, y, slab_top, all_ijk_t, all_ijk_p1, all_ijk_p2, all_ca, all_cb, all_cw, hxz_top, &
                                  .false., inow, Lmax, tag, bc_w_top)
  end if
  Lmax = inow*cc_x*cc_z ! update Lmax
  
  !progagate slabs
  do     ii=0,cc_x-1
    do    i=0,nx-1
      do  j=0,nz-1
        slab_bot(i+ nx*ii,j) = slab_bot(i,j)
        slab_top(i+ nx*ii,j) = slab_top(i,j)
      end do 
    end do 
  end do 
  do     jj=0,cc_z-1
    do    j=0,nz-1
      do  i=0,nx*cc_x-1
        slab_bot(i,j+nz*jj) = slab_bot(i,j)
        slab_top(i,j+nz*jj) = slab_top(i,j)
      end do 
    end do 
  end do 

  ! propagate 
  do    ii=0,cc_x-1
    do  jj=0,cc_z-1
      do i=0,inow-1
        i2                =   i + inow*(ii*cc_z + jj)
        all_ijk_t (i2,0)  =  all_ijk_t (i,0) + nx*ii
        all_ijk_p1(i2,0)  =  all_ijk_p1(i,0) + nx*ii
        all_ijk_p2(i2,0)  =  all_ijk_p2(i,0) + nx*ii

        all_ijk_t (i2,1)  =  all_ijk_t (i,1) + nz*jj
        all_ijk_p1(i2,1)  =  all_ijk_p1(i,1) + nz*jj
        all_ijk_p2(i2,1)  =  all_ijk_p2(i,1) + nz*jj

        all_ijk_t (i2,2)  =  all_ijk_t (i,2)
        all_ijk_p1(i2,2)  =  all_ijk_p1(i,2)
        all_ijk_p2(i2,2)  =  all_ijk_p2(i,2)

        all_ca(i2)        =  all_ca(i)
        all_cb(i2)        =  all_cb(i)
        all_cw(i2)        =  all_cw(i)
      end do 
    end do 
  end do


  allocate(buffer_ijk(0:Lmax-1,0:3))
  do i=0,Lmax-1
    j              = all_ijk_t(i,0) 
    all_ijk_t(i,0) = all_ijk_t(i,2)
    all_ijk_t(i,2) = j
    all_ijk_t(i,3) = i
  end do 
  call mergesort_rows(all_ijk_t,0,Lmax-1,buffer_ijk)
  do i=0,Lmax-1
    j              = all_ijk_t(i,0) 
    all_ijk_t(i,0) = all_ijk_t(i,2)
    all_ijk_t(i,2) = j
  end do 

  buffer_ijk(0:Lmax-1,0:2) = all_ijk_p1(0:Lmax-1,0:2)
  do j=0,Lmax-1
    i                 =  all_ijk_t(j,3)
    all_ijk_p1(j,0:2) =  buffer_ijk(i,0:2)
  end do 

  buffer_ijk(0:Lmax-1,0:2) = all_ijk_p2(0:Lmax-1,0:2)
  do j=0,Lmax-1
    i                 =  all_ijk_t(j,3)
    all_ijk_p2(j,0:2) =  buffer_ijk(i,0:2)
  end do 
  deallocate(buffer_ijk)
  allocate(buffer_c  (0:Lmax-1))

  buffer_c(0:Lmax-1) = all_ca(0:Lmax-1)
  do j=0,Lmax-1
    i                 =  all_ijk_t(j,3)
    all_ca(j)         =  buffer_c(i)
  end do 

  buffer_c(0:Lmax-1) = all_cb(0:Lmax-1)
  do j=0,Lmax-1
    i                 =  all_ijk_t(j,3)
    all_cb(j)         =  buffer_c(i)
  end do 

  buffer_c(0:Lmax-1) = all_cw(0:Lmax-1)
  do j=0,Lmax-1
    i                 =  all_ijk_t(j,3)
    all_cw(j)         =  buffer_c(i)
  end do 
  deallocate(buffer_c)

  block 
    character(len=32) :: fname
    write(fname,'(A27,A1,A4)') '../../geom_data/ibm_coeffs_',tag,'.dat'
    open(19,file=fname, action='write', status='replace', form='unformatted', access='stream')
      write(19) Lmax
      do i=0,Lmax-1
        write(19) all_ijk_t(i,0:2), all_ca(i), all_ijk_p1(i,0:2), all_cb(i), all_ijk_p2(i,0:2), all_cw(i)
      end do 
    close(19)
  end block

  block 
    character(len=34) :: fname
    write(fname,'(A25,A1,A8)') '../../geom_data/ibm_inds_',tag,'_top.dat'
    open(19,file=fname, action='write', status='replace', form='unformatted', access='stream')
    write(19) nx*cc_x, nz*cc_z
    do    j=0,nz*cc_z-1
      do  i=0,nx*cc_x-1
        write(19) slab_top(i,j)
      end do 
    end do 
    close(19)
    write(fname,'(A25,A1,A8)') '../../geom_data/ibm_inds_',tag,'_bot.dat'
    open(19,file=fname, action='write', status='replace', form='unformatted', access='stream')
    write(19) nx*cc_x, nz*cc_z
    do    j=0,nz*cc_z-1
      do  i=0,nx*cc_x-1
        write(19) slab_bot(i,j)
      end do 
    end do 
    close(19)
  end block
  call cpu_time(time_finish)
  write(6,'(A6,A2,A2,F12.6)') 'Done: ',tag,' :',(time_finish-time_start); flush(6)
end program
