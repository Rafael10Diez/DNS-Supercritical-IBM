module diezDecomp_api_ibm
  use mpi
  use, intrinsic                    ::  iso_fortran_env, only: i8 => int64, sp => real32, dp => real64
  use, intrinsic                    ::  iso_c_binding  , only: c_loc, c_size_t
  use diezdecomp_api_generic, only   :  diezdecomp_track_mpi_decomp, diezdecomp_parsed_mpi_ranks
#if defined(_OPENACC)
  use openacc, only: acc_handle_kind
#else
  use, intrinsic :: iso_fortran_env, only: acc_handle_kind => int64
#endif
  implicit none
! #if defined(_DIEZDECOMP_SINGLE) (single precision should not be enabled)
!   integer, parameter :: rp          = sp
!   integer, parameter :: MPI_REAL_RP = MPI_REAL
! #else
  integer, parameter :: rp          = dp
  integer, parameter :: MPI_REAL_RP = MPI_DOUBLE_PRECISION
! #endif
  real(rp), parameter            :: LARGE = huge(1._rp)
  integer(acc_handle_kind), save :: diezdecomp_stream_default = 1

  private

  public :: diezDecomp_ibm_type , &
            diezdecomp_ibm_init , &
            diezdecomp_ibm_exec

  ! -------------------- data types --------------------
  type diezDecomp_ibm_type_sendrecv
    integer                ::  n_recv, n_send, nreq, &
                               loc_recv_i0, loc_send_i0, loc_sendrecv_s, n_pairs
    integer(i8)            ::  wsize_buf_mpi
    logical                ::  initialized = .false.
    logical                ::  mode_forcing
    integer , allocatable  ::  recv_arr(:,:), send_arr(:,:), &
                               count_recv(:)  , count_send(:)  , &
                               rdispls(:)  , sdispls(:)  , &
                               rcounts(:)  , scounts(:)  , &
                               all_request(:), all_status (:,:), &
                               send_A_pos(:) , send_buff_pos(:), &
                               recv_A_pos(:) , recv_buff_pos(:,:), &
                               aux_ii_to_ibuf(:), &
                               i0(:), s(:), other(:)
    real(dp), allocatable  ::  recv_buff_cw12(:,:)
    logical , allocatable  ::  is_send(:)
  end type

  type diezDecomp_ibm_type
    integer                             ::  n_forcings, n_fluids, nx_glob, nz_glob, ny_glob, irank_mpi, nproc_mpi, &
                                            di_rank, dj_rank, dk_rank, nhalo, nskip
    integer(i8)                         ::  wsize_buf_mpi
    logical                             ::  initialized = .false.
    type (diezDecomp_ibm_type_sendrecv) ::  obj_sendrecv_forcing, obj_sendrecv_fluids
    integer , allocatable               ::  i2b(:), j2b(:), k2b(:)
    type(diezdecomp_parsed_mpi_ranks)   ::  obj_ranks
  end type

  contains

  ! -------------------- init wrappers --------------------
  subroutine diezdecomp_ibm_init(this, all_ijk, all_c, n_forcings, n_fluids, A_shape, &
                                 blocks_i0, blocks_j0, blocks_k0, blocks_nx, blocks_nz, blocks_ny, pos_mpi_ref, &
                                 nx_glob, nz_glob, ny_glob, nhalo, nskip, irank_mpi, nproc_mpi)
    implicit none
    type(diezDecomp_ibm_type)  ::  this
    real(rp)                   ::  all_c(0:,0:)
    integer                    ::  all_ijk(0:,0:,0:), n_forcings, n_fluids, A_shape(0:2), &
                                   blocks_i0(0:), blocks_j0(0:), blocks_k0(0:), &
                                   blocks_nx(0:), blocks_nz(0:), blocks_ny(0:), pos_mpi_ref(0:2), &
                                   nx_glob, nz_glob, ny_glob, &
                                   nhalo, nskip, irank_mpi, nproc_mpi

    if (this%initialized) then ; error stop 'this%initialized' ; end if
    call diezdecomp_track_mpi_decomp(pos_mpi_ref, this%obj_ranks, irank_mpi, nproc_mpi)

    allocate(this%i2b(0:nx_glob-1), this%j2b(0:nz_glob-1), this%k2b(0:ny_glob-1))
    call fill_i2b(this%i2b, blocks_i0, blocks_nx, this%obj_ranks%shape_mpi_ranks(0), nx_glob)
    call fill_i2b(this%j2b, blocks_j0, blocks_nz, this%obj_ranks%shape_mpi_ranks(1), nz_glob)
    call fill_i2b(this%k2b, blocks_k0, blocks_ny, this%obj_ranks%shape_mpi_ranks(2), ny_glob)

    this%n_forcings =  n_forcings
    this%n_fluids   =  n_fluids

    this%nx_glob    =  nx_glob
    this%nz_glob    =  nz_glob
    this%ny_glob    =  ny_glob
    this%irank_mpi  =  irank_mpi
    this%nproc_mpi  =  nproc_mpi
    this%nhalo      =  nhalo
    this%nskip      =  nskip

    this%di_rank    =  blocks_i0(this%obj_ranks%flat_mpi_ranks(irank_mpi,0))
    this%dj_rank    =  blocks_j0(this%obj_ranks%flat_mpi_ranks(irank_mpi,1))
    this%dk_rank    =  blocks_k0(this%obj_ranks%flat_mpi_ranks(irank_mpi,2))

    if (this%obj_ranks%mpi_ranks(this%i2b(this%di_rank), &
                                 this%j2b(this%dj_rank), &
                                 this%k2b(this%dk_rank)).ne.irank_mpi) error stop 'mpi_ranks(i2b(di),j2b(dj),k2b(dk)).ne.irank'

      call diezdecomp_ibm_sendrecv_init(this%obj_sendrecv_forcing , this, all_ijk, all_c, A_shape, nhalo, nskip, .true. )
      call diezdecomp_ibm_sendrecv_init(this%obj_sendrecv_fluids, this, all_ijk, all_c, A_shape, nhalo, nskip, .false.)

    this%wsize_buf_mpi  =  1
    this%wsize_buf_mpi  =  max(this%wsize_buf_mpi, &
                               max(this%obj_sendrecv_forcing%wsize_buf_mpi, this%obj_sendrecv_fluids%wsize_buf_mpi))
    this%initialized    =  .true.
  end subroutine

  subroutine diezdecomp_ibm_sendrecv_init(this, par, all_ijk, all_c, A_shape, nhalo, nskip, mode_forcing)
    type(diezDecomp_ibm_type_sendrecv)  ::  this
    type(diezDecomp_ibm_type)           ::  par
    logical                             ::  mode_forcing
    integer                             ::  A_shape(0:2), nhalo, nskip,s_inner(0:2), all_ijk(0:,0:,0:)
    real(rp)                            ::  all_c(0:,0:)
    s_inner = A_shape - 2*nskip - 2*nhalo

    if (this%initialized) then ; error stop 'this%initialized' ; end if
    this%mode_forcing = mode_forcing
    allocate(this%count_recv(0:par%nproc_mpi-1)    , this%count_send(0:par%nproc_mpi-1),&
             this%sdispls   (0:par%nproc_mpi-1)    , this%scounts   (0:par%nproc_mpi-1),&
             this%rdispls   (0:par%nproc_mpi-1)    , this%rcounts   (0:par%nproc_mpi-1))
    associate(n_recv       => this%n_recv    , n_send      => this%n_send     ,&
              count_recv   => this%count_recv, count_send  => this%count_send ,&
              sdispls      => this%sdispls   , scounts     => this%scounts    ,&
              rdispls      => this%rdispls   , rcounts     => this%rcounts    ,&
              di_rank      => par%di_rank    , dj_rank     => par%dj_rank     , dk_rank     => par%dk_rank   ,&
              nx_glob      => par%nx_glob    , nz_glob     => par%nz_glob     , ny_glob     => par%ny_glob   ,&
              irank_mpi    => par%irank_mpi  , nproc_mpi   => par%nproc_mpi   ,&
              i2b          => par%i2b        , &
              j2b          => par%j2b        , &
              k2b          => par%k2b        , &
              rank_3d_grid => par%obj_ranks%mpi_ranks)

      associate(pos => n_recv)
        block
          integer :: iter_outer,ii,jj,ig,jg,kg,p
          ! fill recv points
          do iter_outer=0,1
            pos = 0
            do   ii=0,par%n_forcings-1
              ! recv: point in 1st/2nd layer area (leave switch for it)
              ig = modulo(all_ijk(0,0,ii), nx_glob)
              jg = modulo(all_ijk(1,0,ii), nz_glob)
              kg = modulo(all_ijk(2,0,ii), ny_glob)
              if (mode_forcing) then
                ! identify missing halo area (to be fetched from another process)
                if (get_is_halo(ig, jg, kg, s_inner, di_rank,dj_rank,dk_rank, nx_glob, nz_glob, ny_glob, nhalo)) then
                  p                 =  rank_3d_grid(i2b(ig), j2b(jg), k2b(kg))
                  if (iter_outer == 1) this%recv_arr(:,pos)  =  (/ p, ii, 0 /)
                  pos                                        =  pos + 1
                end if
              else
                ! identify fluid points
                block
                  integer :: i,j,k
                  if (get_is_inner(ig, jg, kg, s_inner, di_rank,dj_rank,dk_rank, nx_glob, nz_glob, ny_glob)) then
                    ! forcing point in irank_mpi?
                    do jj=1,2
                      i  = modulo(all_ijk(0,jj,ii), nx_glob)
                      j  = modulo(all_ijk(1,jj,ii), nz_glob)
                      k  = modulo(all_ijk(2,jj,ii), ny_glob)
                      if (get_is_halo(i, j, k, s_inner, di_rank,dj_rank,dk_rank, nx_glob, nz_glob, ny_glob, nhalo)) then
                        ! if halo, copy locally
                        p  =  rank_3d_grid(i2b(ig), j2b(jg), k2b(kg))
                        if (p.ne.irank_mpi) error stop 'p.ne.irank_mpi'
                      else ! otherwise, fetch normally
                        p  =  rank_3d_grid(i2b(i), j2b(j), k2b(k))
                      end if
                      if (iter_outer == 1) this%recv_arr(:,pos)  =  (/ p, ii, jj /)
                      pos                                          =  pos + 1
                    end do
                  end if
                end block
              end if
            end do
            if (this%n_recv.ne.pos) error stop 'this%n_recv.ne.pos'
            if (iter_outer == 0) allocate(this%recv_arr(0:2,0:max(1,this%n_recv)-1))
          end do
        end block
      end associate
      block
        integer :: pos,p, mpi_ierr
        count_recv = 0
        count_send = 0
        do pos=0,n_recv-1
          p              =  this%recv_arr(0,pos)
          count_recv(p)  =  count_recv(p) + 1
        end do
        call mpi_alltoall(count_recv, 1, mpi_int, count_send, 1, mpi_int, mpi_comm_world, mpi_ierr)
        n_send   =  sum(count_send)

        if (this%n_send.ne.n_send) error stop 'this%n_send.ne.n_send'
        allocate(this%send_arr(0:2,0:max(1,max(this%n_recv,this%n_send))-1))

        call mergesort_cols(this%recv_arr,0,n_recv-1,this%send_arr) ! sorted by (proc,ii,jj) send_arr is buffer
        scounts  =  3*count_recv
        rcounts  =  3*count_send
        call fill_displs(sdispls, scounts, nproc_mpi)
        call fill_displs(rdispls, rcounts, nproc_mpi)
        call mpi_alltoallv(this%recv_arr, scounts, sdispls, mpi_int, &
                           this%send_arr, rcounts, rdispls, mpi_int, mpi_comm_world, mpi_ierr)
      end block
      block
        integer :: i
        this%nreq = 0
        do i=0,nproc_mpi-1
          if ((count_recv(i)>0).and.(irank_mpi.ne.i)) this%nreq = this%nreq + 1
          if ((count_send(i)>0).and.(irank_mpi.ne.i)) this%nreq = this%nreq + 1
        end do
      end block
      allocate(this%all_request(                 0:max(1,this%nreq)-1) , &
               this%all_status (MPI_STATUS_SIZE, 0:max(1,this%nreq)-1) , &
               this%i0         (                 0:max(1,this%nreq)-1) , &
               this%s          (                 0:max(1,this%nreq)-1) , &
               this%other      (                 0:max(1,this%nreq)-1) , &
               this%is_send    (                 0:max(1,this%nreq)-1) )
      block
        integer :: i,prev,ireq
        this%loc_send_i0     = -1
        this%loc_recv_i0     = -1
        this%loc_sendrecv_s  = -1
        prev = 0
        ireq = 0
        do i=0,nproc_mpi-1
          if (count_recv(i)>0) then
            if (irank_mpi.ne.i) then
              this%i0     (ireq) = prev
              this%s      (ireq) = count_recv(i)
              this%other  (ireq) = i
              this%is_send(ireq) = .false.
              ireq               = ireq + 1
            else
              this%loc_recv_i0     = prev
              this%loc_sendrecv_s  = count_recv(i)
            end if
            prev  =  prev + count_recv(i)
          end if
        end do
        do i=0,nproc_mpi-1
          if (count_send(i)>0) then
            if (irank_mpi.ne.i) then
              this%i0     (ireq) = prev
              this%s      (ireq) = count_send(i)
              this%other  (ireq) = i
              this%is_send(ireq) = .true.
              ireq               = ireq + 1
            else
              this%loc_send_i0 = prev
              if (this%loc_sendrecv_s.ne.count_send(i)) error stop 'this%loc_sendrecv_s.ne.count_send(i)'
            end if
            prev  =  prev + count_send(i)
          end if
        end do
        if (ireq.ne.this%nreq) error stop 'ireq.ne.this%nreq'
        if (prev.ne.(n_send+n_recv)) error stop 'prev.ne.(n_send+n_recv)'
      end block
      ! buffers:
      ! recv - forcing: A(i) =    buffer(j)
      ! recv - fluid: A(i) = c1*buffer(j1) + c2*buffer(j2) + cw
      ! send - forcing/fluid: buffer(j) = A(i)
      allocate(this%send_A_pos(0:max(1,n_send)-1), this%send_buff_pos(0:max(1,n_send)-1))
      block
        integer :: di,dj,dk
        di   = di_rank-nhalo-nskip
        dj   = dj_rank-nhalo-nskip
        dk   = dk_rank-nhalo-nskip
        if (mode_forcing) then
          ! recv - forcing: A(i) =    buffer(j)
          block
            integer :: pos,i,j,k,ii,jj,ibuf,si,sj,sk,iter_outer
            do iter_outer=0,1
              ibuf = 0
              do pos=0,n_recv-1
                ii                        =  this%recv_arr(1,pos)
                jj                        =  this%recv_arr(2,pos)
                if (jj.ne.0) error stop 'jj.ne.0'
                do si=0,2
                  do sj=0,2
                    do sk=0,2
                      i = all_ijk(0,jj,ii) + nx_glob*(si-1)
                      j = all_ijk(1,jj,ii) + nz_glob*(sj-1)
                      k = all_ijk(2,jj,ii) + ny_glob*(sk-1)
                      if (((i-di)>=0).and.((i-di)<A_shape(0)).and.&
                          ((j-dj)>=0).and.((j-dj)<A_shape(1)).and.&
                          ((k-dk)>=0).and.((k-dk)<A_shape(2))) then 
                        if (iter_outer==1) then 
                          this%recv_A_pos(ibuf)      =  (k-dk)*A_shape(1)*A_shape(0) + (j-dj)*A_shape(0) + (i-di)
                          this%recv_buff_pos(0,ibuf) =  pos ! (+ 0) because recv is first
                        end if
                        ibuf                       =  ibuf + 1
                      end if 
                    end do 
                  end do 
                end do 
              end do
              if (iter_outer==0) allocate(this%recv_buff_pos (0:0,0:max(1,ibuf)-1), this%recv_A_pos(0:max(1,ibuf)-1))
              this%n_pairs = ibuf
            end do 
          end block
        else
          block
            integer :: iter_outer,pos,ibuf,max_ibuf,ii,jj
            allocate(this%aux_ii_to_ibuf(0:max(1,par%n_forcings)-1))
            do iter_outer=0,1
              this%aux_ii_to_ibuf = -1
              if (iter_outer == 1) this%recv_buff_pos = -1
              ! recv - fluid: A(i) = c1*buffer(j1) + c2*buffer(j2) + cw
              ibuf = -1
              this%n_pairs = 0
              do pos=0,n_recv-1
                ii                             =  this%recv_arr(1,pos)
                jj                             =  this%recv_arr(2,pos)
                if (this%aux_ii_to_ibuf(ii)<0) then
                  ibuf                         =  this%n_pairs
                  this%n_pairs                 =  this%n_pairs + 1
                  this%aux_ii_to_ibuf(ii)      =  ibuf
                  if (iter_outer == 1) this%recv_buff_cw12(:,ibuf)  =  all_c(:,ii)
                  block
                    integer  ::  i,j,k
                    i                      =  pfix(all_ijk(0, 0,ii), di_rank-nhalo, di_rank + s_inner(0)-1+nhalo, nx_glob)
                    j                      =  pfix(all_ijk(1, 0,ii), dj_rank-nhalo, dj_rank + s_inner(1)-1+nhalo, nz_glob)
                    k                      =  pfix(all_ijk(2, 0,ii), dk_rank-nhalo, dk_rank + s_inner(2)-1+nhalo, ny_glob)
                    if ((i.ne.all_ijk(0, 0,ii)).or.&
                        (j.ne.all_ijk(1, 0,ii)).or.&
                        (k.ne.all_ijk(2, 0,ii))) error stop 'all_ijk(:,0,ii) modified by pfix'
                    if (((i-di)<0).or.((i-di)>=A_shape(0)).or.&
                        ((j-dj)<0).or.((j-dj)>=A_shape(1)).or.&
                        ((k-dk)<0).or.((k-dk)>=A_shape(2))) error stop 'error out-of-bounds A_shape'
                    if (iter_outer == 1) this%recv_A_pos(ibuf)  =  (k-dk)*A_shape(1)*A_shape(0) + (j-dj)*A_shape(0) + (i-di)
                  end block
                else
                  ibuf                     =  this%aux_ii_to_ibuf(ii)
                end if
                if (iter_outer == 1) this%recv_buff_pos(jj,ibuf)  =  pos
              end do
              if (iter_outer==0) allocate(this%recv_buff_pos (1:2,0:max(1,this%n_pairs)-1) , &
                                          this%recv_buff_cw12(0:2,0:max(1,this%n_pairs)-1) , &
                                          this%recv_A_pos(        0:max(1,this%n_pairs)-1) )
            end do
          end block
          if (minval(this%recv_buff_pos(:,0:this%n_pairs-1))<0) error stop 'minval(this%recv_buff_pos(:,0:this%n_pairs-1))<0'
        end if
        block
          integer :: pos,i,j,k,ii,jj
          ! send - forcing/fluid: buffer(j) = A(i)
          do pos=0,n_send-1
            ii                      =  this%send_arr(1,pos)
            jj                      =  this%send_arr(2,pos)
            if (      this%send_arr(0,pos).ne.irank_mpi    ) error stop '      this%send_arr(0,pos).ne.irank_mpi    '
            if (      mode_forcing .and.(jj.ne.0)         ) error stop '      mode_forcing .and.(jj.ne.0)         '
            if ((.not.mode_forcing).and.((jj<1).or.(jj>2))) error stop '(.not.mode_forcing).and.((jj<1).or.(jj>2))'
            i                       =  pfix(all_ijk(0,jj,ii), di_rank-nhalo, di_rank + s_inner(0)-1+nhalo, nx_glob)
            j                       =  pfix(all_ijk(1,jj,ii), dj_rank-nhalo, dj_rank + s_inner(1)-1+nhalo, nz_glob)
            k                       =  pfix(all_ijk(2,jj,ii), dk_rank-nhalo, dk_rank + s_inner(2)-1+nhalo, ny_glob)
            if (((i-di)<0).or.((i-di)>=A_shape(0)).or.&
                ((j-dj)<0).or.((j-dj)>=A_shape(1)).or.&
                ((k-dk)<0).or.((k-dk)>=A_shape(2))) error stop 'error out-of-bounds A_shape'
            this%send_A_pos(pos)    =  (k-dk)*A_shape(1)*A_shape(0) + (j-dj)*A_shape(0) + (i-di)
            this%send_buff_pos(pos) =  pos + n_recv ! send goes after recv
          end do
        end block
      end block
    end associate
    !$acc enter data copyin(this)
    !$acc enter data copyin(this%send_A_pos, this%send_buff_pos, this%recv_A_pos, this%recv_buff_pos)
    if (.not.mode_forcing) then
      !$acc enter data copyin(this%recv_buff_cw12)
    end if
    this%wsize_buf_mpi  =  max(1,this%n_recv) + max(1,this%n_send)
    this%initialized    =  .true.
  end subroutine

  ! -------------------- exec wrappers --------------------
  subroutine diezdecomp_ibm_exec(this, A, buffer)
    type(diezDecomp_ibm_type)  ::  this
    real(rp)                   ::  A(0:*), buffer(0:*)
    call diezdecomp_ibm_sendrecv_exec(this%obj_sendrecv_fluids, A, buffer)
    call diezdecomp_ibm_sendrecv_exec(this%obj_sendrecv_forcing , A, buffer)
  end subroutine

  subroutine diezdecomp_ibm_sendrecv_exec(this, A, buffer)
    type(diezDecomp_ibm_type_sendrecv) ::  this
    real(rp)                           ::  A(0:*), buffer(0:*)
    call diezdecomp_ibm_buffer_pack(A, buffer, this%send_A_pos, this%send_buff_pos, this%n_send, .true.)
    call copy_1d_buffer(buffer, buffer, this%loc_send_i0, this%loc_recv_i0, this%loc_sendrecv_s) ! if loc_sendrecv_s<=0, copy_1d_buffer will do nothing
    if (this%nreq>0) then
      !$acc wait
      block
        integer :: i, i0, s, other, mpi_ierr
        do i=0,this%nreq-1
          i0   = this%i0(i)
          s    = this%s(i)
          other= this%other(i)
          !$acc host_data use_device(buffer)
          if (this%is_send(i)) then
                 call MPI_ISend(buffer(i0), s, MPI_REAL_RP, other, 1, mpi_comm_world, this%all_request(i), mpi_ierr)
          else ; call MPI_IRecv(buffer(i0), s, MPI_REAL_RP, other, 1, mpi_comm_world, this%all_request(i), mpi_ierr)
          end if
          !$acc end host_data
        end do
        call mpi_waitall(this%nreq, this%all_request, this%all_status, mpi_ierr)
      end block
      !$acc wait
    end if
    if (this%mode_forcing) then
      call diezdecomp_ibm_buffer_pack(A, buffer, this%recv_A_pos, this%recv_buff_pos, this%n_pairs, .false.) ! 0:0 dimension will be flattened
    else
      call diezdecomp_ibm_apply_forcing(A, buffer, this%recv_A_pos, this%recv_buff_pos, this%recv_buff_cw12, this%n_pairs)
    end if
  end subroutine

  ! -------------------- openacc kernels --------------------
  subroutine diezdecomp_ibm_apply_forcing(A, buffer, pos_A, pos_b, cw12, n)
    real(rp) :: A(0:*), buffer(0:*), cw12(0:,0:)
    integer  :: pos_A(0:), pos_b(1:,0:), n, i
    if (n>0) then
      !$acc wait
      !$acc parallel loop collapse(1) default(present) private(i)
      do i=0,n-1
        A(pos_A(i)) = cw12(1,i)*buffer(pos_b(1,i)) + cw12(2,i)*buffer(pos_b(2,i)) +  cw12(0,i)
      end do
      !$acc wait
    end if
  end subroutine

  subroutine diezdecomp_ibm_buffer_pack(A, buffer, all_pos, all_b, n, mode_pack)
    real(rp) :: A(0:*), buffer(0:*)
    integer  :: all_pos(0:), all_b(0:*), n, i
    logical  :: mode_pack
    if (n>0) then
      !$acc wait
      if (mode_pack) then
        !$acc parallel loop collapse(1) default(present) private(i)
        do i=0,n-1
          buffer(all_b(i)) = A(all_pos(i))
        end do
      else
        !$acc parallel loop collapse(1) default(present) private(i)
        do i=0,n-1
          A(all_pos(i))   = buffer(all_b(i))
        end do
      end if
      !$acc wait
    end if
  end subroutine

  ! -------------------- utilities --------------------
  function fits_bbox(i, j, k, ia, ib, ja, jb, ka, kb, ni, nj, nk) result(is_fit)
    logical :: is_fit
    integer :: i, j, k, ia, ib, ja, jb, ka, kb, ni, nj, nk,&
               di,dj,dk,i2,j2,k2
    is_fit = .false.
    do     di=0,2
      do   dj=0,2
        do dk=0,2
          i2 = i + (di-1)*ni
          j2 = j + (dj-1)*nj
          k2 = k + (dk-1)*nk
          if ((ia<=i2).and.(i2<=ib).and.(ja<=j2).and.(j2<=jb).and.(ka<=k2).and.(k2<=kb)) is_fit=.true.
        end do
      end do
    end do
  end function

  function get_is_halo(i, j, k, s_inner, di,dj,dk, nx, nz, ny, nhalo) result(is_halo)
    logical :: is_halo
    integer :: i, j, k, i0a, i0b, j0a, j0b, k0a, k0b, i1a, i1b, j1a, j1b, k1a, k1b,&
               di,dj,dk,nhalo, nx, nz, ny, s_inner(0:2)
    i0a  =  -nhalo + di ; i0b  =  -1 + di ; i1a  =  s_inner(0) + di; i1b  =  s_inner(0)+nhalo-1 + di
    j0a  =  -nhalo + dj ; j0b  =  -1 + dj ; j1a  =  s_inner(1) + dj; j1b  =  s_inner(1)+nhalo-1 + dj
    k0a  =  -nhalo + dk ; k0b  =  -1 + dk ; k1a  =  s_inner(2) + dk; k1b  =  s_inner(2)+nhalo-1 + dk
    is_halo = (fits_bbox(i,j,k, i0a, i0b, j0a, j1b, k0a, k1b, nx, nz, ny).or.&
               fits_bbox(i,j,k, i1a, i1b, j0a, j1b, k0a, k1b, nx, nz, ny).or.&
               fits_bbox(i,j,k, i0a, i1b, j0a, j0b, k0a, k1b, nx, nz, ny).or.&
               fits_bbox(i,j,k, i0a, i1b, j1a, j1b, k0a, k1b, nx, nz, ny).or.&
               fits_bbox(i,j,k, i0a, i1b, j0a, j1b, k0a, k0b, nx, nz, ny).or.&
               fits_bbox(i,j,k, i0a, i1b, j0a, j1b, k1a, k1b, nx, nz, ny))
  end function

  function get_is_inner(i, j, k, s_inner, di,dj,dk, nx, nz, ny) result(is_inner)
    logical :: is_inner
    integer :: i, j, k, ia, ib, ja, jb, ka, kb,&
               di,dj,dk, nx, nz, ny, s_inner(0:2)
    ia = di ; ib = s_inner(0)-1+di
    ja = dj ; jb = s_inner(1)-1+dj
    ka = dk ; kb = s_inner(2)-1+dk
    is_inner = fits_bbox(i,j,k, ia, ib, ja, jb, ka, kb, nx, nz, ny)
  end function

  recursive subroutine mergesort_cols(A,i,j,buffer)
    implicit none
    integer :: i,j,L,mid,i0,j0,k,k2,pos,buffer(:,0:),A(:,0:)
    logical :: aux_bool
    L = j-i+1
    if (L>1) then
      if (L==2) then
        if (is_greater(A(:,i),A(:,j))) then
          buffer(:,0) = A(:,i)
          A(:,i)      = A(:,j)
          A(:,j)      = buffer(:,0)
        end if
      else
        mid = (i+j)/2
        call mergesort_cols(A,     i, mid, buffer)
        call mergesort_cols(A, mid+1,   j, buffer)
        i0 = i
        j0 = mid+1
        k  = 0
        do while (((i0<=mid).or.(j0<=j)).and.(k<=L))
          aux_bool = (j0>j)
          if (.not.aux_bool) then
            aux_bool = (i0<=mid)
            if (aux_bool) aux_bool=aux_bool.and.(.not.is_greater(A(:,i0),A(:,j0)))
          end if
          if (aux_bool) then !((j0>j).or.((i0<=mid).and.(.not.is_greater(A(:,i0),A(:,j0))))) then
            buffer(:,k) = A(:,i0)
            i0          = i0 + 1
          else
            buffer(:,k) = A(:,j0)
            j0          = j0 + 1
          end if
          k = k + 1
        end do
        pos = i
        do k2=0,k-1
          A(:,pos) = buffer(:,k2)
          pos      = pos + 1
        end do
      end if
    end if
  end subroutine

  function is_greater(a, b) result(res)
    implicit none
    logical :: res
    integer :: a(:), b(:), i
    res = .false.
    do i=1,max(size(a,1),size(b,1))
      if (a(i) > b(i)) then
        res = .true.
        return
      else if (a(i) < b(i)) then
        res = .false.
        return
      end if
    end do
    ! reached this line ---> all equal ---> is_greater=false --> ok
  end function

  subroutine copy_1d_buffer(A, B, ia, ib, sa)
    implicit none
    real(rp) :: A(0:*), B(0:*)
    integer  :: i, ia, ib, sa, di
    integer(acc_handle_kind) :: stream
    if (sa>0) then
      di = ia - ib
      !$acc parallel loop  collapse(1) default(present)
      do i = ib, ib+sa-1
        B(i) = A(i+di)
      end do
    end if
  end subroutine

  subroutine fill_i2b(i2b, b_i0, b_n0, sb, n)
    integer :: i2b(0:), b_i0(0:), b_n0(0:), sb, n, i,j
    i2b = -1
    if (sb>0) then
      do   i=0,sb-1
        do j=b_i0(i),(b_i0(i)+b_n0(i)-1)
          i2b(j) = i
        end do
      end do
    end if
    if (minval(i2b)<0) error stop 'minval(i2b)<0'
  end subroutine

  subroutine fill_displs(displs, counts, n)
    integer :: displs(0:), counts(0:), n, i, prev
    prev = 0
    if (n>0) then 
      do i=0,n-1
        displs(i) = prev
        prev      = prev + counts(i)
      end do
    end if
  end subroutine

  function pfix(i, i0,i1, n) result(new_i)
    implicit none
    integer :: i, i0,i1, n, i2, shift_i, new_i
    logical :: found
    new_i = i
    if (.not.(((i0<=i).and.(i<=i1)))) then
      found = .false.
      do shift_i=0,2
        i2 = i + n*(shift_i-1)
        if ((i0<=i2).and.(i2<=i1)) then
          new_i = i2
          found = .true.
        end if
      end do
      if (.not.found) error stop '.not.found'
    end if
  end function


end module
