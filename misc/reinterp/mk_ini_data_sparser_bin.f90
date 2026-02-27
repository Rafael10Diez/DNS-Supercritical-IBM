module my_utils
  implicit none 

  type  gridInfo
    integer             :: nx, nz, ny
    character(len=1)         :: tag
    real*8, allocatable :: y(:), bot(:,:), top(:,:)
  end type

  contains

  subroutine read_data_g_info(this, fname_header)
    implicit none
    type(gridInfo)              ::  this
    character(len=*)            ::  fname_header
    integer                     ::  nx, nz, ny

    open(19,file=fname_header)
      read(19,*) nx, nz, ny
      ! allocate arrays
      this%nx  =  nx
      this%nz  =  nz
      this%ny  =  ny
      allocate(this%y  (0:ny-1)        , &
               this%bot(0:nx-1,0:nz-1) , &
               this%top(0:nx-1,0:nz-1) )

      this%y   = -1
      this%bot = -1
      this%top = -1

      block 
        integer :: i
        do   i=0,ny-1 
          read(19,*) this%y(i)
        end do
      end block

      block 
        integer :: i,j, iter
        real*8  :: val0, val1
        do  iter=1,nx*nz
          read(19,*) i,j, val0, val1
          this%bot(i,j) = val0
          this%top(i,j) = val1
        end do
      end block
    close(19)
  end subroutine

  subroutine read_coeffs_xz(ii, cc, fname)
    integer ::  ii(0:,0:), i, n
    real*8  ::  cc(0:,0:)
    character(len=*) :: fname
    open(19,file=fname)
      read(19,*) n
      if ((n.ne.size(ii,1)).or.(n.ne.size(cc,1))) write(6,*) 'ERROR: size mismatch', fname; flush(6)
      do i=0,n-1
        read(19,*) ii(i,0), cc(i,0), ii(i,1), cc(i,1), ii(i,2), cc(i,2)
      end do 
    close(19)
  end subroutine

  subroutine reinterp_arr_xz(source_arr, buffer_2d, gs, gt, ix,cx,jz,cz)
    implicit none
    real*8 :: source_arr(0:,0:,0:), buffer_2d(0:,0:), cx(0:,0:), cz(0:,0:)
    integer :: i,j,k, nx_t, nz_t, nz_s, ny_s,         ix(0:,0:), jz(0:,0:)
    type(gridInfo)       ::  gs, gt
    do      k=0, gs%ny-1
      do    i=0, gt%nx-1
        do  j=0, gs%nz-1
          buffer_2d(i,j) = source_arr(ix(i,0),j,k)*cx(i,0) + &
                           source_arr(ix(i,1),j,k)*cx(i,1) + &
                           source_arr(ix(i,2),j,k)*cx(i,2) 
        end do
      end do
      do    i=0, gt%nx-1
        do  j=0, gt%nz-1
          source_arr(i,j,k) = buffer_2d(i,jz(j,0))*cz(j,0) + &
                              buffer_2d(i,jz(j,1))*cz(j,1) + &
                              buffer_2d(i,jz(j,2))*cz(j,2)
        end do
      end do
    end do
  end subroutine

  function find_closest(A,i,j,lo_orig,hi_orig,val) result(k)
  integer  ::  i,j,k,lo,hi,mid,lo_orig,hi_orig
  real*8   ::  A(0:,0:,0:), val
  lo = lo_orig
  hi = hi_orig
    do while ((hi-lo)>1)
      mid = (lo+hi)/2
      if (val>=A(i,j,mid)) lo = mid
      if (val<=A(i,j,mid)) hi = mid
    end do
    k = lo
    if (abs(A(i,j,hi)-val)<abs(A(i,j,lo)-val)) k = hi
  end function

end module my_utils

program main
    use my_utils
    implicit none

    type(gridInfo)       ::  gs, gt
    real*8 , allocatable ::  ys_adim(:,:,:), buffer_2d(:,:), cx(:,:), cz(:,:), gs_A(:,:,:)
    integer, allocatable ::  ix(:,:), jz(:,:)
    integer              ::  nx_max, nz_max, ny_max
    character(len=1)          ::  tag

    call get_command_argument(1, tag)

    block 
      character(len=22) :: fname22
      write(fname22,'(A7,A1,A9,A1,A4)') './data_',tag,'/gs_info_',tag,'.dat'
      call read_data_g_info(gs, fname22)
      write(fname22,'(A7,A1,A9,A1,A4)') './data_',tag,'/gt_info_',tag,'.dat'
      call read_data_g_info(gt, fname22)
    end block
    
    ! allocate arrays
    nx_max = max(gs%nx, gt%nx)
    nz_max = max(gs%nz, gt%nz)
    ny_max = max(gs%ny, gt%ny)

    allocate(gs_A(0:nx_max-1, 0:nz_max-1, 0:ny_max-1))
    allocate(  ix(0:gt%nx -1, 0:2), &
               cx(0:gt%nx -1, 0:2))
    allocate(  jz(0:gt%nz -1, 0:2), &
               cz(0:gt%nz -1, 0:2))

    allocate(ys_adim  , mold=gs_A)
    allocate(buffer_2d(0:nx_max-1,&
                       0:nz_max-1))

    ! read gs_A
    block 
      integer  ::  i,j,k, iter, n0_check, n1_check, n2_check
      real*8   ::  t0, t1, t2, val
      character(len=27) :: fname27
      write(fname27,'(A7,A1,A14,A1,A4)') './data_',tag,'/source_array_',tag,'.bin'
      open(19,file=fname27,action='read',status='old',form='unformatted',access='stream')
        read(19) n0_check, n1_check, n2_check
        if ((n0_check.ne.gs%nx).or.&
            (n1_check.ne.gs%nz).or.&
            (n2_check.ne.gs%ny)) then 
            write(6,*) 'ERROR: grid mismatch', n0_check, n1_check, n2_check, gs%nx, gs%nz, gs%ny
            go to 4096
        end if 
        do     k=0, gs%ny-1
          do   j=0, gs%nz-1
            do i=0, gs%nx-1
              read(19) gs_A(i,j,k)
            end do 
          end do
        end do 
      close(19)
    end block

    ! read_coeffs_xz
    block 
      character(len=26) :: fname26
      write(fname26,'(A7,A1,A8,A1,A9)') './data_',tag,'/coeffs_',tag,'_ic_x.dat'
      call read_coeffs_xz(ix, cx, fname26)
      write(fname26,'(A7,A1,A8,A1,A9)') './data_',tag,'/coeffs_',tag,'_jc_z.dat'
      call read_coeffs_xz(jz, cz, fname26)
    end block 

    ! ys_adim
    block 
      integer :: i,j,k
      do     k=0,gs%ny-1
        do   j=0,gs%nz-1
          do i=0,gs%nx-1
            ys_adim(i,j,k) = (gs%y(k) - gs%bot(i,j))/(gs%top(i,j) - gs%bot(i,j))
          end do 
        end do 
      end do 
    end block

    ! reinterp_arr_xz
    call reinterp_arr_xz(ys_adim, buffer_2d, gs, gt, ix,cx,jz,cz)
    call reinterp_arr_xz(gs_A   , buffer_2d, gs, gt, ix,cx,jz,cz)
    
    ! process and write
    block 
      integer ::  i,j,k, ks
      real*8  ::  yy, y0, y1, y2, val, c0,c1,c2, deno
      character(len=27) :: fname27
      write(fname27,'(A7,A1,A14,A1,A4)') './data_',tag,'/target_array_',tag,'.bin'
      open(19,file=fname27,action='write',status='replace',form='unformatted',access='stream')
      write(19) gt%nx, gt%nz, gt%ny
        do      k=0,gt%ny-1
          do    j=0,gt%nz-1
            do  i=0,gt%nx-1
              yy    =  (gt%y(k) - gt%bot(i,j))/(gt%top(i,j) - gt%bot(i,j))
              yy    =  max(0.d0,min(1.d0,yy))
              ks    =  find_closest(ys_adim,i,j,1,gs%ny-2,yy)
              y0    =  ys_adim(i,j,ks-1) - yy
              y1    =  ys_adim(i,j,ks  ) - yy
              y2    =  ys_adim(i,j,ks+1) - yy
              deno  =  y0**2*y1 - y0**2*y2 - y0*y1**2 + y0*y2**2 + y1**2*y2 - y1*y2**2
              c0    =  (y1**2*y2 - y1*y2**2)/deno
              c1    =  (-y0**2*y2 + y0*y2**2)/deno
              c2    =  (y0**2*y1 - y0*y1**2)/deno
              val   =  gs_A(i,j,ks-1)*c0 + gs_A(i,j,ks  )*c1 + gs_A(i,j,ks+1)*c2
              write(19) val
            end do 
          end do 
        end do 
      close(19)
    end block
    
    deallocate(gs_A, ix, cx, jz, cz, ys_adim, buffer_2d, &
               gs%y, gs%bot, gs%top, &
               gt%y, gt%bot, gt%top, )
4096 continue
end program

! import sympy as sp 
! x,y0,y1,y2,a,b,c,f0,f1,f2  =  sp.symbols('x,y0,y1,y2,a,b,c,f0,f1,f2',real=True)
! f                          =  a + b*x + c*(x**2)
! eqs                        =  [f.subs(x,y0) - f0,
!                                f.subs(x,y1) - f1,
!                                f.subs(x,y2) - f2]
! sol                        =  sp.solve(eqs,a,b,c)
! f                          =  f.subs(a,sol[a]).subs(b,sol[b]).subs(c,sol[c]).subs(x,0)
! c0                         =  f.diff(f0)
! c1                         =  f.diff(f1)
! c2                         =  f.diff(f2)
! assert (f - (c0*f0 + c1*f1 + c2*f2)).expand() == 0 
! assert (c0 + c1 + c2).simplify() == 1
! print(f"{c0 = }")
! print(f"{c1 = }")
! print(f"{c2 = }")
! ## c0   = (y1**2*y2 - y1*y2**2)/(y0**2*y1 - y0**2*y2 - y0*y1**2 + y0*y2**2 + y1**2*y2 - y1*y2**2)
! ## c1   = (-y0**2*y2 + y0*y2**2)/(y0**2*y1 - y0**2*y2 - y0*y1**2 + y0*y2**2 + y1**2*y2 - y1*y2**2)
! ## c2   = (y0**2*y1 - y0*y1**2)/(y0**2*y1 - y0**2*y2 - y0*y1**2 + y0*y2**2 + y1**2*y2 - y1*y2**2)
! ## deno = y0**2*y1 - y0**2*y2 - y0*y1**2 + y0*y2**2 + y1**2*y2 - y1*y2**2
! ## c0   = (y1**2*y2 - y1*y2**2)/deno
! ## c1   = (-y0**2*y2 + y0*y2**2)/deno
! ## c2   = (y0**2*y1 - y0*y1**2)/deno

