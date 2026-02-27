program main
  use Mod_Height_Function
  
  implicit none

  type(type_surf_obj)  ::  surf_obj
  integer              ::  i,n
  real*8               :: y_bot, y_top, max_error
  real*8, allocatable  :: xzy_bt(:,:)

  call init_surf_obj(surf_obj)

  open(19,file='./fortran_surf_files/verif_xzy_bt.dat')
    read(19,*) n
    allocate(xzy_bt(0:n-1,0:3))
    do i=0,n-1
     read(19,*) xzy_bt(i,:)
    end do 
  close(19)

  max_error = 0

  do i=0,n-1
    y_bot      =  get_h(xzy_bt(i,0), xzy_bt(i,1), surf_obj, .true.)
    y_top      =  get_h(xzy_bt(i,0), xzy_bt(i,1), surf_obj, .false.)
    max_error  =  max(max_error, &
                      abs(y_bot-xzy_bt(i,2)), &
                      abs(y_top-xzy_bt(i,3)))
    if (max_error>1e-10) then 
      write(6,*) 'ERROR: Mismatch in verification xzy_bt', i, xzy_bt(i,:), y_bot, y_top ; flush(6)
      error stop 'ERROR: Mismatch in verification xzy_bt'
    end if 
  end do 

  write(6,*) 'SUCCESS: Finished verification xzy_bt', max_error ;flush(6)
end program

! module load 2023
! cd /scratch-shared/diezsanhuesa/rafael/Dimples_asymrotated_onlydimple_0.1_corners_rand_321_vrun/input
! gfortran -fdefault-real-8 rough_surf_define.f90 main_verif_rough_surf_define.f90; ./a.out