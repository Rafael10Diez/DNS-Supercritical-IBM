module mod_avg_apply
  use, intrinsic :: iso_fortran_env, only: i8 => int64, dp => real64
  implicit none
  contains
  subroutine avg_apply(nx_glob, nz_glob, ny_glob, nz_loc, ny_loc, irank_mpi, nproc_mpi, &
        !__[PYTHON_HOLDER_AVG_VARS]
        !__[PYTHON_HOLDER_PROF1D_COMPBULK_VARS]
    avg_bin_nbins, prof_1d_avg_N, avg_N_now, avg_nraw, iters_full, avg_iter_start, avg_freq, avg_inv_nraw, binary_cycle,&
    U,V,W,T,P_now,P_old,P_per,P,slab_00k,&
    mu_arr           , cond_arr         , rho_arr          , buoyancy_arr     , &
    enth_arr         , cp_arr           , div_jacobian_arr , &
    temp_arr_n_xzy, prof_1d_temp, tr_red1d_fwd, tr_red1d_bwd, tr_io_fwd, &
    mpi_divs_z, ny_loc_00k, tr_red1d_B, tr_red1d_B_1d, work_tr)
    use diezdecomp_api_generic, only: diezdecomp_props_transp
    use mod_utils, only: batched_binary_reduction
    use mod_io, only: cfd_write_arr
    implicit none
    !__[PYTHON_HOLDER_AVG_DECL_FUNC_VARS]
    !__[PYTHON_HOLDER_PROF1D_COMPBULK_FUNC_DECL_VARS]
    real(dp)                       ::  temp_arr_n_xzy(0:,0:,0:), prof_1d_temp(0:), work_tr(0:), &
                                       U(-1:,-1:,-1:), V(-1:,-1:,-1:), W(-1:,-1:,-1:), T(-1:,-1:,-1:), slab_00k(0:,0:,0:),&
                                       P_now(-1:,-1:,-1:), P_old(-1:,-1:,-1:), P_per(-1:,-1:,-1:), P(0:,0:,0:), avg_inv_nraw,&
                                       mu_arr          (-1:,-1:,-1:) , &
                                       cond_arr        (-1:,-1:,-1:) , &
                                       rho_arr         (-1:,-1:,-1:) , &
                                       buoyancy_arr    (-1:,-1:,-1:) , &
                                       enth_arr        (-1:,-1:,-1:) , &
                                       cp_arr          (-1:,-1:,-1:) , &
                                       div_jacobian_arr(-1:,-1:,-1:)
    integer                        ::  nx_glob, nz_glob, ny_glob, nz_loc, ny_loc, irank_mpi, nproc_mpi, &
                                       avg_bin_nbins, prof_1d_avg_N, avg_N_now, avg_nraw, iters_full, avg_iter_start, &
                                       avg_freq, mpi_divs_z, ny_loc_00k, binary_cycle(0:)
    real(dp)                       ::  tr_red1d_B(0:,0:,0:), tr_red1d_B_1d(0:)
    type(diezdecomp_props_transp)  ::  tr_red1d_fwd, tr_red1d_bwd, tr_io_fwd

    prof_1d_avg_N = prof_1d_avg_N + 1

    !__[PYTHON_HOLDER_PROF1D_APPLY_0]
    !block
    !  integer :: i,j,k
    !  !$acc parallel loop  collapse(3) default(present)
    !  do     k = 0, ny_loc-1
    !    do   j = 0, nz_loc-1
    !      do i = 0, nx_glob-1
    !        temp_arr_n_xzy(i,j,k) = {y}
    !      end do
    !    end do
    !  end do
    !  call batched_binary_reduction(temp_arr_n_xzy, prof_1d_temp, nx_loc*nz_loc, ny_loc, .true., .true.,&
    !                                tr_red1d_fwd,tr_red1d_bwd,mpi_divs_z,ny_loc_00k,tr_red1d_B,tr_red1d_B_1d,work_tr)
    !  !$acc parallel loop  collapse(1) default(present)
    !  do k = 0, ny_loc-1
    !    {x}(k) = {x}(k) + prof_1d_temp(k)
    !  end do
    !end block

    if ((iters_full >= avg_iter_start).and.(mod(iters_full-avg_iter_start,avg_freq)==0)) then
      block
        logical :: do_write

        ! ------------------------ report parameters ------------------------
        avg_N_now = avg_N_now + 1
        block
          integer :: i,j,k
          !$acc parallel loop  collapse(3) default(present)
          do  k  = 0,   ny_loc-1
          do  j  = 0,   nz_loc-1
          do  i  = 0,   nx_glob-1
           !__[PYTHON_HOLDER_AVG_APPLY_0]
           !{x}(i,j,k) = {x}(i,j,k) + ({y})*avg_inv_nraw
          end do
          end do
          end do
        end block
        if (avg_N_now==avg_nraw) then
          ! write file, before or after binary_tree average
          if (avg_bin_nbins>0) then
            block
              integer :: avg_bin_reach
              call avg_binary_tree(avg_bin_reach, avg_bin_nbins, &
                !__[PYTHON_HOLDER_AVG_VARS]
                !__[PYTHON_HOLDER_PROF1D_COMPBULK_VARS]
                binary_cycle)
              do_write = (avg_bin_reach == (2**avg_bin_nbins))
            end block
          else
            do_write = .true.
          end if

        if (do_write) then
          !__[PYTHON_HOLDER_AVG_APPLY_1]
          ! !$acc update host({x})
          ! !$acc wait
          ! block
          !   character(len={(len(x)+31)})  ::  my_file
          !     associate(A => {x})
          !       write(my_file,'(A{(len(x)+18)},I0.9,A4)') './avg/array_{x}_iter_',iters_full,'.dat'
          ! boilerplate_write_arr
          !     end associate
          ! end block
        end if

        ! set average variables to zero
        avg_N_now = 0
        block
          integer :: i,j,k
          !$acc parallel loop  collapse(3) default(present)
          do  k  = 0, ny_loc-1
          do  j  = 0, nz_loc-1
          do  i  = 0, nx_glob-1
            !__[PYTHON_HOLDER_AVG_APPLY_2]
            ! {x}(i,j,k) = 0.
          end do
          end do
          end do
        end block
        end if
      end block
    end if
  end subroutine

  subroutine avg_binary_tree(avg_bin_reach, avg_bin_nbins, &
    !__[PYTHON_HOLDER_AVG_VARS]
    !__[PYTHON_HOLDER_PROF1D_COMPBULK_VARS]
      binary_cycle)
    integer :: avg_bin_reach, avg_bin_nbins, binary_cycle(0:)
    !__[PYTHON_HOLDER_AVG_DECL_FUNC_VARS]
    !__[PYTHON_HOLDER_PROF1D_COMPBULK_FUNC_DECL_VARS]

    block
      integer :: ibin
      logical :: do_loop
      avg_bin_reach = 1
      ibin          = 0
      do_loop       = .true.
      do while ((ibin<avg_bin_nbins).and.(do_loop))
          if (binary_cycle(ibin)==0) then
    !__[PYTHON_HOLDER_BINAVG_0]
    ! binary_{x}(:,:,:, ibin)  =  {x}
          else
    !__[PYTHON_HOLDER_BINAVG_1]
    ! {x} =  0.5*({x} + binary_{x}(:,:,:, ibin))
              avg_bin_reach                   =  avg_bin_reach * 2
          end if
          do_loop                             =  (binary_cycle(ibin)==1)
          binary_cycle(ibin)                  =  mod(binary_cycle(ibin) + 1, 2)
          ibin                                =  ibin + 1
      end do
    end block

    ! import random
    ! import numpy as np
    ! nbins   =  6
    ! shape   =  3,4
    ! all_x   =  [ 20*(np.random.rand(*shape)-0.5) for _ in range(2**nbins)]
    ! random.shuffle(all_x)
    ! assert len(all_x) == 2**nbins
    ! cycle  =  np.zeros(nbins).astype(int)
    ! tree   =  np.zeros((nbins,) + all_x[0].shape)
    ! for x in all_x:
    !     reach = 1
    !     new   = x
    !     loop  = True
    !     i     = 0
    !     while (i<nbins) and loop:
    !         if cycle[i] == 0:
    !             tree[i] = new
    !         else:
    !             new    = 0.5*(new + tree[i])
    !             reach *= 2
    !         loop     = (cycle[i] == 1)
    !         cycle[i] = (cycle[i]+1)%2
    !         i       += 1
    !     print(reach)
    !     print(new)
    ! assert reach       == 2**nbins
    ! assert np.fabs(new-sum(all_x)/len(all_x)).max() <  1e-12
    ! print('Tests passed!')
  end subroutine
end module