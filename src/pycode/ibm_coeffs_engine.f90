module ibm_coeffs_engine
  use Mod_Height_Function, only : type_surf_obj, get_h
  use, intrinsic :: iso_fortran_env, only: i8 => int64, dp => real64
  implicit none
  contains
  subroutine search_mid_h(orig_xzy_a, orig_xzy_b, xzy_n, orig_h_a, orig_h_b, obj, is_bottom)
    implicit none
    type(type_surf_obj) ::  obj
    logical             ::  is_bottom
    real(dp)            ::       xzy_a(0:2),      xzy_b(0:2), xzy_n(0:2), &
                            orig_xzy_a(0:2), orig_xzy_b(0:2),&
                            x0,z0,y0,dx,dz,dy,lo,hi,mid,new_h, h_a, h_b, orig_h_a, orig_h_b
    integer             ::  iters
  
    xzy_a(0:2)  =  orig_xzy_a(0:2)
    xzy_b(0:2)  =  orig_xzy_b(0:2)
    h_a         =  orig_h_a
    h_b         =  orig_h_b
    
    if ((xzy_a(2)-h_a) > &
        (xzy_b(2)-h_b)) then 
      xzy_a(0:2)  =  orig_xzy_b(0:2) 
      xzy_b(0:2)  =  orig_xzy_a(0:2)
      h_a         =  orig_h_b
      h_b         =  orig_h_a
    end if
  
    dx = xzy_b(0) - xzy_a(0)
    dz = xzy_b(1) - xzy_a(1)
    dy = xzy_b(2) - xzy_a(2)
  
    ! eval_h  =  lambda mid: get_h(dx*mid + x0, dz*mid + z0) - (dy*mid + y0)
      
    if ((abs(dx) < 1e-10).and.(abs(dz) < 1e-10)) then ! edge case
      xzy_n(0)  =  0.5*(xzy_a(0)+xzy_b(0))
      xzy_n(1)  =  0.5*(xzy_a(1)+xzy_b(1))
      xzy_n(2)  =  0.5*(  h_a   +  h_b)
    else
      x0 = xzy_a(0)
      z0 = xzy_a(1)
      y0 = xzy_a(2)
      if ((h_b-xzy_b(2)-1e-10)> & ! eval_h(1)-1e-10
          (h_a-xzy_a(2))) then    ! eval_h(0)
        write(6,*) 'assert np.all(eval_h(1) <= (eval_h(0)+1e-10))' ; flush(6)
        error stop 'assert np.all(eval_h(1) <= (eval_h(0)+1e-10))'
      end if
      if ((xzy_b(2)+1e-10) < h_b) then 
        write(6,*) 'assert np.all(         (xzy_b[2] - get_h(*xzy_b[:2])) >= -1e-10) # "b" is above surface' ; flush(6)
        error stop 'assert np.all(         (xzy_b[2] - get_h(*xzy_b[:2])) >= -1e-10) # "b" is above surface'
      end if
      if ((xzy_a(2)-1e-10) > h_a) then 
        write(6,*) 'assert np.all(1e-10 >= (xzy_a[2] - get_h(*xzy_a[:2]))) # "a" is below surface' ; flush(6)
        error stop 'assert np.all(1e-10 >= (xzy_a[2] - get_h(*xzy_a[:2]))) # "a" is below surface'
      end if
  
      lo    = 0
      hi    = 1
      iters = 0
  
      do while ((hi-lo)>1e-10)
        iters = iters + 1
        if (iters>=10000) then
          write(6,*) 'assert iters < 10_000, (lo, hi, xzy_a, xzy_b)'; flush(6)
          error stop 'assert iters < 10_000, (lo, hi, xzy_a, xzy_b)'
        end if 
        mid   = 0.5*(lo+hi)
        new_h = get_h(dx*mid + x0, dz*mid + z0,obj,is_bottom)-(dy*mid + y0)
        if ((new_h-1e-10) <= 0) hi = mid
        if ((new_h+1e-10) >= 0) lo = mid
      end do 
      mid = 0.5*(lo+hi)
      xzy_n(0)  =  dx*mid + x0
      xzy_n(1)  =  dz*mid + z0
      xzy_n(2)  =  get_h(xzy_n(0),xzy_n(1),obj,is_bottom)
    end if 
  end subroutine
  
      
  subroutine get_contained_alpha(contained, alpha, T, A, B, ee, xzy_p1)
    implicit none
    real(dp)  :: T(0:2), A(0:2), B(0:2), ee(0:2), xzy_p1(0:2), alpha, norm_A, norm_B, norm_T
    logical :: contained
    norm_A     =  sum((A-xzy_p1)*ee)
    norm_B     =  sum((B-xzy_p1)*ee)
    norm_T     =  sum((T-xzy_p1)*ee)
    contained  =  (((norm_A-1e-10) <= norm_T).and.&
                                     (norm_T <= (norm_B+1e-10)))
    alpha      =  1-(norm_T-norm_A)/(norm_B-norm_A)
  end subroutine
  
  subroutine apply_Fadlun(ca, cb, cw, xzy_t, xzy_p1, xzy_p2, ht, h2, obj, is_bottom)      
    implicit none
    type(type_surf_obj) :: obj
    real(dp)              ::  ca, cb, cw, xzy_t(0:2), xzy_n(0:2), xzy_p1(0:2), xzy_p2(0:2), ee(0:2), xzy_g(0:2), &
                            alpha_n, alpha_p1, ht, h2
    logical             ::  contain_n1, contain_12, is_bottom

    call search_mid_h(xzy_t, xzy_p2, xzy_n, ht, h2, obj, is_bottom)

    xzy_g(0)  =  2*xzy_n(0) - xzy_t(0)
    xzy_g(1)  =  2*xzy_n(1) - xzy_t(1)
    xzy_g(2)  =  2*xzy_n(2) - xzy_t(2)
    
    ee        =  (xzy_p2-xzy_p1)/sum((xzy_p2-xzy_p1)**2)
    
    call get_contained_alpha(contain_n1, alpha_n , xzy_g, xzy_n , xzy_p1, ee, xzy_p1)
    call get_contained_alpha(contain_12, alpha_p1, xzy_g, xzy_p1, xzy_p2, ee, xzy_p1)
    
    if (contain_n1.and.contain_12) then
        if (.not.((abs(alpha_n)<1e-10).and.(abs(alpha_p1-1)<1e-10))) then 
          write(6,*) 'assert abs(alpha_n)<1e-10 and abs(alpha_p1-1)<1e-10' ; flush(6)
          error stop 'assert abs(alpha_n)<1e-10 and abs(alpha_p1-1)<1e-10'
        end if
    end if
    
    if (contain_n1) then 
      ca  =   -(1-alpha_n) ! negative signs
      cb  = 0
      cw  = 2 -(  alpha_n)
    else
      if (.not.contain_12) then 
        write(6,*) 'assert contain_12' ; flush(6)
        error stop 'assert contain_12'
      end if
      ca  =    -(    alpha_p1) ! negative signs
      cb  =    -(1 - alpha_p1)
      cw  =  2
    end if
  end subroutine

    subroutine make_rough_surface_quady(x, z, y, slab, all_ijk_t, all_ijk_p1, all_ijk_p2, all_ca, all_cb, all_cw, &
                                      hxz, is_bottom, inow, Lmax, tag, bc_w)
    implicit none
    integer             :: nx,nz,ny,i,j,k,inow,ncount_old,Lmax, &
                           all_ijk_t(0:,0:), all_ijk_p1(0:,0:), all_ijk_p2(0:,0:), slab(0:,0:), it,jt
    real(dp)            :: all_ca(0:), all_cb(0:), all_cw(0:), hxz(0:,0:), bc_w,&
                           ca, cb, cw, x(0:), z(0:), y(0:),&
                           time_start, time_finish
    character(len=1)    :: tag
    logical             :: is_bottom
    nx     =  size(x,1)
    nz     =  size(z,1)
    ny     =  size(y,1)

    if (is_bottom) then ; slab = -1
    else                ; slab = ny ; end if
    
    call cpu_time(time_start)
    ncount_old =  inow

    do     i=0,nx-1 
      do   j=0,nz-1 
        do k=0,ny-1
          if  (is_bottom.and.(          hxz(i,j) < y(k))) then 
            if               (y(k-1) <= hxz(i,j)        ) then 
              all_ijk_t (inow,0:2) =  (/   i ,   j ,   k      /)
              all_ijk_p1(inow,0:2) =  (/   i ,   j ,   k + 1  /)
              all_ijk_p2(inow,0:2) =  (/   i ,   j ,   k + 2  /)
              call apply_quady(ca, cb, cw, y(k), y(k+1), y(k+2), hxz(i,j))
              all_ca(inow) =  ca
              all_cb(inow) =  cb
              all_cw(inow) =  cw*bc_w
              it           =  all_ijk_t(inow,0)
              jt           =  all_ijk_t(inow,1)
              slab(it,jt)  =  max(k-1, slab(it,jt))
              inow = inow + 1
              if (inow>Lmax) error stop 'ERROR: inow>Lmax'
            end if 
          end if
          if  ((.not.is_bottom).and.(         y(k) < hxz(i,j)          )) then 
            if                      (                hxz(i,j)  <= y(k+1)) then 
              all_ijk_t (inow,0:2) =  (/   i ,   j ,   k      /)
              all_ijk_p1(inow,0:2) =  (/   i ,   j ,   k - 1  /)
              all_ijk_p2(inow,0:2) =  (/   i ,   j ,   k - 2  /)
              call apply_quady(ca, cb, cw, y(k), y(k-1), y(k-2), hxz(i,j))
              all_ca(inow) =  ca
              all_cb(inow) =  cb
              all_cw(inow) =  cw*bc_w
              it           =  all_ijk_t(inow,0)
              jt           =  all_ijk_t(inow,1)
              slab(it,jt)  =  min(k+1, slab(it,jt))
              inow = inow + 1
              if (inow>Lmax) error stop 'ERROR: inow>Lmax'
            end if 
          end if
        end do
      end do 
    end do 
    call cpu_time(time_finish)
    write(6,'(A4,A2,A2,L2,I12,F12.6)') '    ',tag,' :',is_bottom,(inow - ncount_old), &
    (time_finish-time_start); flush(6)
  end subroutine

  subroutine apply_quady(ca, cb, cw, y0, ya, yb, hw)
    real(dp) :: ca, cb, cw, y0, ya, yb, hw, xw, xa, xb
    ! w : (xa**2*xb - xa*xb**2)/(xa**2*xb - xa**2*xw - xa*xb**2 + xa*xw**2 + xb**2*xw - xb*xw**2)
    ! a : (xb**2*xw - xb*xw**2)/(xa**2*xb - xa**2*xw - xa*xb**2 + xa*xw**2 + xb**2*xw - xb*xw**2)
    ! b : (-xa**2*xw + xa*xw**2)/(xa**2*xb - xa**2*xw - xa*xb**2 + xa*xw**2 + xb**2*xw - xb*xw**2)
    xw  =  hw - y0
    xa  =  ya - y0
    xb  =  yb - y0
    cw  =  (xa**2*xb - xa*xb**2)/(xa**2*xb - xa**2*xw - xa*xb**2 + xa*xw**2 + xb**2*xw - xb*xw**2)
    ca  =  (xb**2*xw - xb*xw**2)/(xa**2*xb - xa**2*xw - xa*xb**2 + xa*xw**2 + xb**2*xw - xb*xw**2)
    cb  =  (-xa**2*xw + xa*xw**2)/(xa**2*xb - xa**2*xw - xa*xb**2 + xa*xw**2 + xb**2*xw - xb*xw**2)
    ! import sympy as sp 
    ! a,b,c,xw,xa,xb,yw,ya,yb  =  sp.symbols('a,b,c,xw,xa,xb,yw,ya,yb',real=True,positive=True) 
    ! f                        =  lambda x,a,b,c: a*x**2 + b*x + c
    ! eqs = [(f(x,a,b,c)-y) for x,y in zip([xw,xa,xb], [yw,ya,yb])]
    ! sol = sp.solve(eqs,a,b,c)
    ! yc = f(0,sol[a],sol[b],sol[c])
    ! weights = [yc.diff(y) for y in [yw,ya,yb]]
    ! assert (sum(w*y for w,y in zip(weights,[yw,ya,yb])) - yc).expand() == 0
    ! assert sum(weights).expand().simplify() == 1
    ! for w,c in zip(weights,'wab'): print(c,':', w)
  end subroutine

  subroutine make_rough_surface_ghost(x, z, y, slab, all_ijk_t, all_ijk_p1, all_ijk_p2, all_ca, all_cb, all_cw, &
                                      hxz, is_ok, obj, is_bottom, inow, Lmax, tag, bc_w)
    implicit none
    type(type_surf_obj) :: obj
    integer             :: nx,nz,ny,ii_d,dx,dz,dy,i,j,k,i2,j2,k2,inow,ncount_old,Lmax, &
                           all_ijk_t(0:,0:), all_ijk_p1(0:,0:), all_ijk_p2(0:,0:), slab(0:,0:), it,jt, sx,sz
    real(dp)            :: mag_x, mag_z, all_ca(0:), all_cb(0:), all_cw(0:), hxz(0:,0:), bc_w,&
                           ca, cb, cw, xzy_t(0:2), xzy_p1(0:2), xzy_p2(0:2),ht,h2, x(0:), z(0:), y(0:),&
                           time_start, time_finish
    character(len=1)    :: tag
    logical             :: is_fluid, is_bottom, is_ok(0:,0:,0:)
    integer, dimension(26, 3) :: all_dxzy = reshape((/ 0,  0, -1 ,&
                                                       0,  0,  1 ,&
                                                      -1,  0,  0 ,&
                                                       1,  0,  0 ,&
                                                       0, -1,  0 ,&
                                                       0,  1,  0 ,&
                                                      -1,  0, -1 ,&
                                                      -1,  0,  1 ,&
                                                       1,  0, -1 ,&
                                                       1,  0,  1 ,&
                                                       0, -1, -1 ,&
                                                       0, -1,  1 ,&
                                                       0,  1, -1 ,&
                                                       0,  1,  1 ,&
                                                      -1, -1,  0 ,&
                                                      -1,  1,  0 ,&
                                                       1, -1,  0 ,&
                                                       1,  1,  0 ,&
                                                      -1, -1, -1 ,&
                                                      -1, -1,  1 ,&
                                                      -1,  1, -1 ,&
                                                      -1,  1,  1 ,&
                                                       1, -1, -1 ,&
                                                       1, -1,  1 ,&
                                                       1,  1, -1 ,&
                                                       1,  1,  1 /), shape(all_dxzy), order=(/2,1/))
    nx     =  size(x,1)
    nz     =  size(z,1)
    ny     =  size(y,1)
    mag_x  =  (x(nx-1)-x(0))/(nx-1.d0)
    mag_z  =  (z(nz-1)-z(0))/(nz-1.d0)
    if ((abs(mag_x-(x(1)-x(0)))>1e-10).or.(abs(mag_z-(z(1)-z(0)))>1e-10)) then 
      write(6,*) 'ERROR: Mismatch in mag_x, mag_z', mag_x, mag_z, (x(1)-x(0)), (z(1)-z(0))
      error stop 'ERROR: Mismatch in mag_x, mag_z'
    end if 
    do     i=0,nx-1 
      do   j=0,nz-1 
        do k=0,ny-1
          if (is_bottom) then ; is_ok(i,j,k) = (hxz(i,j) < y(k))
          else                ; is_ok(i,j,k) = (hxz(i,j) > y(k)) ; end if
        end do
      end do 
    end do 

    if (is_bottom) then ; slab = -1
    else                ; slab = ny ; end if

    do ii_d=1,size(all_dxzy,1)
      call cpu_time(time_start)
      ncount_old =  inow
      dx         =  all_dxzy(ii_d,1)
      dz         =  all_dxzy(ii_d,2)
      dy         =  all_dxzy(ii_d,3)
      if (max(abs(dx),abs(dz),abs(dy)).ne.1) then 
        write(6,*) 'max(abs(dx),abs(dz),abs(dy)).ne.1',dx,dz,dy;flush(6)
        error stop 'max(abs(dx),abs(dz),abs(dy)).ne.1'
      end if 
      do     i=0,nx-1 
        do   j=0,nz-1 
          do k=0,ny-1
            if (is_bottom) then ; is_fluid = (hxz(i,j) < y(k))
            else                ; is_fluid = (hxz(i,j) > y(k)) ; end if
            if (is_fluid) then 
              i2 = modulo(i+dx, nx)
              j2 = modulo(j+dz, nz)
              k2 =        k+dy
              if ((0<=k2).and.(k2<ny)) then 
              if (.not.is_ok(i2,j2,k2)) then 
                if      ((i+dx)<0)   then ; sx =  nx
                else if ((i+dx)>=nx) then ; sx = -nx
                else                      ; sx =  0  ; end if

                if      ((j+dz)<0)   then ; sz =  nz
                else if ((j+dz)>=nz) then ; sz = -nz
                else                      ; sz =  0  ; end if

                all_ijk_t (inow,0:2) =  (/   i + dx + sx ,   j + dz + sz ,   k + dy  /)
                all_ijk_p1(inow,0:2) =  (/   i      + sx ,   j      + sz ,   k       /)
                all_ijk_p2(inow,0:2) =  (/   i - dx + sx ,   j - dz + sz ,   k - dy  /)
                xzy_t (0:2)          =  (/ x(i)+ dx*mag_x, z(j)+ dz*mag_z, y(k + dy) /)
                xzy_p1(0:2)          =  (/ x(i)          , z(j)          , y(k)      /)
                xzy_p2(0:2)          =  (/ x(i)- dx*mag_x, z(j)- dz*mag_z, y(k - dy) /)
                ht                   =  hxz(i2,j2)
                h2                   =  hxz(modulo(i-dx,nx), modulo(j-dz,nz))
                call apply_Fadlun(ca, cb, cw, xzy_t, xzy_p1, xzy_p2, ht, h2, obj, is_bottom)
                all_ca(inow) =  ca
                all_cb(inow) =  cb
                all_cw(inow) =  cw*bc_w
                it           =  all_ijk_t(inow,0)
                jt           =  all_ijk_t(inow,1)
                if (is_bottom) then ; slab(it,jt) = max(all_ijk_t(inow,2), slab(it,jt))
                else                ; slab(it,jt) = min(all_ijk_t(inow,2), slab(it,jt)) ; end if
                inow = inow + 1
                is_ok(i2,j2,k2)  =  .true.
                if (inow>Lmax) then 
                  write(6,*) 'ERROR: inow>Lmax',inow,Lmax,nx,nz,ny,i,j,k,i2,j2,k2,hxz(i,j),hxz(i2,j2),is_bottom;flush(6)
                  error stop 'ERROR: inow>Lmax'
                end if
              end if
              end if
            end if 
          end do 
        end do 
      end do 
      call cpu_time(time_finish)
      write(6,'(A4,A2,A2,L2,3I4,I12,F12.6)') '    ',tag,' :',is_bottom,dx,dz,dy,(inow - ncount_old), &
      (time_finish-time_start); flush(6)
    end do 
  end subroutine

  ! utilities (for main program): is_greater, mergesort_rows
  recursive subroutine mergesort_rows(A,i,j,buffer)
    implicit none
    integer :: i,j,L,mid,i0,j0,k,k2,pos,buffer(0:,:),A(0:,:)
    L = j-i+1
    if (L>1) then
      if (L==2) then
        if (is_greater(A(i,:),A(j,:))) then
          buffer(0,:) = A(i,:)
          A(i,:)      = A(j,:)
          A(j,:)      = buffer(0,:)
        end if
      else
        mid = (i+j)/2
        call mergesort_rows(A,     i, mid, buffer)
        call mergesort_rows(A, mid+1,   j, buffer)
        i0 = i
        j0 = mid+1
        k  = 0
        do while (((i0<=mid).or.(j0<=j)).and.(k<=L))
          if ((j0>j).or.((i0<=mid).and.(.not.is_greater(A(i0,:),A(j0,:))))) then
            buffer(k,:) = A(i0,:)
            i0          = i0 + 1
          else
            buffer(k,:) = A(j0,:)
            j0          = j0 + 1
          end if
          k = k + 1
        end do
        pos = i
        do k2=0,k-1
          A(pos,:) = buffer(k2,:)
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
    ! reched this line ---> all equal ---> is_greater=false
  end function
end module