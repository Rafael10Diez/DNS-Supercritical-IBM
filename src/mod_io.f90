module mod_io
  use mpi
  use openacc
  use, intrinsic :: iso_fortran_env, only: i8 => int64, dp => real64
  implicit none
  contains
  subroutine cfd_read_arr(my_file, A, ih, tr_io_bwd, slab_00k, temp_arr_n_xzy, work_tr, &
                          irank_mpi, nproc_mpi, ny_loc, nz_loc, ny_loc_00k, nx_glob, nz_glob, ny_glob)
    use diezdecomp_api_generic, only: diezdecomp_props_transp, diezdecomp_transp_execute_generic_buf
    implicit none
    integer, intent(in)            :: ih, irank_mpi, nproc_mpi, ny_loc, nz_loc, ny_loc_00k, nx_glob, nz_glob, ny_glob
    character(len=*)               ::  my_file
    type(diezdecomp_props_transp)  ::  tr_io_bwd
    real(dp)                       ::  A(ih:,ih:,ih:), slab_00k(0:,0:,0:), temp_arr_n_xzy(0:,0:,0:), work_tr(0:)
    block
      integer ::  iloc,i,j,k,mpi_ierr, &
                  nx_check, nz_check, ny_check, ny_here
      !$acc wait
      call MPI_BARRIER(mpi_comm_world,mpi_ierr)
      if (irank_mpi == (nproc_mpi-1)) then
          open(19,file=my_file, action='read', status='old', form='unformatted', access='stream')
            read(19) nx_check, nz_check, ny_check
      end if

      call MPI_BCAST(nx_check, 1, MPI_INTEGER, nproc_mpi-1, mpi_comm_world, mpi_ierr)
      call MPI_BCAST(nz_check, 1, MPI_INTEGER, nproc_mpi-1, mpi_comm_world, mpi_ierr)
      call MPI_BCAST(ny_check, 1, MPI_INTEGER, nproc_mpi-1, mpi_comm_world, mpi_ierr)

      if ((nx_check.ne. nx_glob ).or. &
          (nz_check.ne. nz_glob ).or. &
          (ny_check.ne.(ny_glob))) then
          write(6,*) 'ERROR: Mismatch in grid params (boilerplate_read_arr.f90)', &
          nx_glob, nz_glob, ny_glob, nproc_mpi, nx_check, nz_check, ny_check, my_file ; flush(6)
          error stop 'ERROR: Mismatch in grid params (boilerplate_read_arr.f90)'
      end if

    if (irank_mpi == (nproc_mpi-1)) then
      do iloc=0,nproc_mpi-1
        if (iloc==(nproc_mpi-1)) then
          ny_here = ny_loc_00k
        else
          call MPI_Recv(ny_here, 1, MPI_INT, iloc, 1, mpi_comm_world, mpi_status_ignore, mpi_ierr)
        end if
        do      k=0,ny_here-1
          do    j=0,nz_glob-1
            do  i=0,nx_glob-1
              read(19) slab_00k(i,j,k)
            end do
          end do
        end do
        if (iloc<(nproc_mpi-1)) then
          call MPI_Send(slab_00k, nx_glob*nz_glob*ny_here, MPI_DOUBLE_PRECISION, iloc, 1, mpi_comm_world, mpi_ierr)
        end if
      end do
      close(19)
    else
      call MPI_Send(ny_loc_00k, 1, MPI_INT, nproc_mpi-1, 1, mpi_comm_world, mpi_ierr)
      call MPI_Recv(slab_00k, nx_glob*nz_glob*ny_loc_00k, MPI_DOUBLE_PRECISION, nproc_mpi-1, 1, mpi_comm_world, &
                    mpi_status_ignore, mpi_ierr)
    end if
    !$acc wait
    call MPI_BARRIER(mpi_comm_world,mpi_ierr)
    !$acc update device(slab_00k)
    !$acc wait
    call MPI_BARRIER(mpi_comm_world,mpi_ierr)
    call diezdecomp_transp_execute_generic_buf(tr_io_bwd,slab_00k,temp_arr_n_xzy,work_tr)
    !$acc wait
    call MPI_BARRIER(mpi_comm_world,mpi_ierr)
    !$acc parallel loop collapse(3) private(i,j,k) default(present)
    do      k=0,ny_loc-1
      do    j=0,nz_loc-1
        do  i=0,nx_glob-1
          A(i,j,k) = temp_arr_n_xzy(i,j,k)
        end do
      end do
    end do
    !$acc wait
    !$acc update self(A)
    !$acc wait
    call MPI_BARRIER(mpi_comm_world,mpi_ierr)
    if (irank_mpi == 0) then
      write(6,*) 'INFO: read file: ', my_file, nx_glob, nz_glob, ny_glob, nproc_mpi, nx_check, nz_check, ny_check ; flush(6)
    end if
    end block
  end subroutine

  subroutine cfd_write_arr(my_file, A, ih, tr_io_fwd, slab_00k, temp_arr_n_xzy, work_tr, &
                          irank_mpi, nproc_mpi, ny_loc, nz_loc, ny_loc_00k, nx_glob, nz_glob, ny_glob)
    use diezdecomp_api_generic, only: diezdecomp_props_transp, diezdecomp_transp_execute_generic_buf
    implicit none
    integer, intent(in)            :: ih, irank_mpi, nproc_mpi, ny_loc, nz_loc, ny_loc_00k, nx_glob, nz_glob, ny_glob
    character(len=*)               ::  my_file
    type(diezdecomp_props_transp)  ::  tr_io_fwd
    real(dp)                       ::  A(ih:,ih:,ih:), slab_00k(0:,0:,0:), temp_arr_n_xzy(0:,0:,0:), work_tr(0:)

    !associate: my_file, A
    block
      integer :: i,j,k,iloc,mpi_ierr,ny_here
      !$acc wait
      call MPI_BARRIER(mpi_comm_world,mpi_ierr)
    !$acc parallel loop collapse(3) private(i,j,k) default(present)
    do      k=0,ny_loc-1
      do    j=0,nz_loc-1
        do  i=0,nx_glob-1
          temp_arr_n_xzy(i,j,k) = A(i,j,k)
        end do
      end do
    end do
    !$acc wait
    call MPI_BARRIER(mpi_comm_world,mpi_ierr)
    call diezdecomp_transp_execute_generic_buf(tr_io_fwd,temp_arr_n_xzy,slab_00k,work_tr)
    !$acc wait
    call MPI_BARRIER(mpi_comm_world,mpi_ierr)
    !$acc wait
    !$acc update self(slab_00k)
    !$acc wait
    call MPI_BARRIER(mpi_comm_world,mpi_ierr)

      if (irank_mpi == 0) then
        open(19,file=my_file, action='write', status='replace', form='unformatted', access='stream')
          write(19) nx_glob, nz_glob, ny_glob
          do iloc=0,nproc_mpi-1
          if (iloc==0) then
            ny_here = ny_loc_00k
          else
      call MPI_Recv(ny_here ,                       1, MPI_INT             , iloc, 1, mpi_comm_world, mpi_status_ignore, mpi_ierr)
      call MPI_Recv(slab_00k, nx_glob*nz_glob*ny_here, MPI_DOUBLE_PRECISION, iloc, 1, mpi_comm_world, mpi_status_ignore, mpi_ierr)
          end if
            do  k=0,ny_here-1
            do  j=0,nz_glob-1
            do  i=0,nx_glob-1
              write(19) slab_00k(i,j,k)
            end do
            end do
            end do
          end do
        close(19)
      else
        call MPI_Send(ny_loc_00k, 1                         , MPI_INT             , 0, 1, mpi_comm_world, mpi_ierr)
        call MPI_Send(slab_00k  , nx_glob*nz_glob*ny_loc_00k, MPI_DOUBLE_PRECISION, 0, 1, mpi_comm_world, mpi_ierr)
      end if
      !$acc wait
      call MPI_BARRIER(mpi_comm_world,mpi_ierr)
      if (irank_mpi == 0) then
        write(6,*) 'INFO: wrote file: ', my_file; flush(6)
      end if
    end block
  end subroutine
end module
