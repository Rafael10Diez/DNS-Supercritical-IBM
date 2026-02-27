program main
    implicit none
    integer ::  nx, nz, ny, i, j, k
    real*8  :: myzero
    nx     =  560
    nz     =  280
    ny     =  280
    myzero = 0.d0

    open(19,file='./array_ini_P.dat', action='write', status='replace', form='unformatted',access='stream')
      write(19) nx, nz, ny
      do      k=0,ny-1
        do    j=0,nz-1
          do  i=0,nx-1
            write(19) myzero
          end do
        end do
      end do
    close(19)
end program

! cd /home/rafael/Downloads/orig_input/scanned_data_uttiya/dns_ref/reinterp_vels ; gfortran mk_fix_zero.f90 ; ./a.out