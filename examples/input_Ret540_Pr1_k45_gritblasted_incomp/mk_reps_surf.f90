program main
    implicit none
    integer              ::  nx, nz, ny, ii, jj, i, j, k, nreps_x, nreps_z
    real*8, allocatable  ::  U(:,:,:) , V(:,:,:) , W(:,:,:) , T(:,:,:)

    nreps_x = 2
    nreps_z = 2

    open(19,file='../orig_input_1x1tile/array_ini_U.dat',action='read',status='old',form='unformatted',access='stream')
      read(19) nx, nz, ny
      allocate(U(0:nx-1,0:nz-1,0:ny-1))
      do      k=0,ny-1
        do    j=0,nz-1
          do  i=0,nx-1
            read(19) U(i,j,k)
          end do
        end do
      end do
    close(19)

    open(19,file='../orig_input_1x1tile/array_ini_V.dat',action='read',status='old',form='unformatted',access='stream')
      read(19) nx, nz, ny
      allocate(V(0:nx-1,0:nz-1,0:ny-1))
      do      k=0,ny-1
        do    j=0,nz-1
          do  i=0,nx-1
            read(19) V(i,j,k)
          end do
        end do
      end do
    close(19)

    open(19,file='../orig_input_1x1tile/array_ini_W.dat',action='read',status='old',form='unformatted',access='stream')
      read(19) nx, nz, ny
      allocate(W(0:nx-1,0:nz-1,0:ny-1))
      do      k=0,ny-1
        do    j=0,nz-1
          do  i=0,nx-1
            read(19) W(i,j,k)
          end do
        end do
      end do
    close(19)

    open(19,file='../orig_input_1x1tile/array_ini_T.dat',action='read',status='old',form='unformatted',access='stream')
      read(19) nx, nz, ny
      allocate(T(0:nx-1,0:nz-1,0:ny-1))
      do      k=0,ny-1
        do    j=0,nz-1
          do  i=0,nx-1
            read(19) T(i,j,k)
          end do
        end do
      end do
    close(19)

    open(19,file='array_ini_U.dat', action='write', status='replace', form='unformatted',access='stream')
      write(19) nreps_x*nx, nreps_z*nz, ny
      do      k=0,ny-1
        do   jj=0,nreps_z-1
        do    j=0,nz-1
          do  ii=0,nreps_x-1
          do  i=0,nx-1
            write(19) U(i,j,k)
          end do
          end do
        end do
        end do
      end do
    close(19)

    open(19,file='array_ini_V.dat', action='write', status='replace', form='unformatted',access='stream')
      write(19) nreps_x*nx, nreps_z*nz, ny
      do      k=0,ny-1
        do   jj=0,nreps_z-1
        do    j=0,nz-1
          do  ii=0,nreps_x-1
          do  i=0,nx-1
            write(19) V(i,j,k)
          end do
          end do
        end do
        end do
      end do
    close(19)

    open(19,file='array_ini_W.dat', action='write', status='replace', form='unformatted',access='stream')
      write(19) nreps_x*nx, nreps_z*nz, ny
      do      k=0,ny-1
        do   jj=0,nreps_z-1
        do    j=0,nz-1
          do  ii=0,nreps_x-1
          do  i=0,nx-1
            write(19) W(i,j,k)
          end do
          end do
        end do
        end do
      end do
    close(19)

    open(19,file='array_ini_T.dat', action='write', status='replace', form='unformatted',access='stream')
      write(19) nreps_x*nx, nreps_z*nz, ny
      do      k=0,ny-1
        do   jj=0,nreps_z-1
        do    j=0,nz-1
          do  ii=0,nreps_x-1
          do  i=0,nx-1
            write(19) T(i,j,k)
          end do
          end do
        end do
        end do
      end do
    close(19)

    deallocate(U,V,W,T)
end program

! cd /media/rafael/DATA/Dropbox/scanned_data_uttiya/dns_ref ; gfortran mk_linear_temp.f90 ; ./a.out