
# ------------------------------------------------------------------------
#                  Basic Library
# ------------------------------------------------------------------------

from  os.path            import  join                                       as  pjoin
from  os                 import  listdir                                    as  os_listdir
from  os                 import  system                                     as  os_system
from  time               import  time
from  os.path            import  abspath, dirname, basename, isfile, isdir
from  sys                import  argv
from  string             import  ascii_lowercase, ascii_uppercase

try:
    _folder_
except:
    _folder_  =  dirname(abspath(__file__))

# ------------------------------------------------------------------------
#                  Basic Functions
# ------------------------------------------------------------------------

format_dt             =  lambda x :  "%02d:%02d:%02d [hh:mm:ss]" % (x//3600, (x%3600)//60, x%60)
lmap                  =  lambda f,x:  list   (map   (f,x))
lfilter               =  lambda f,x:  list   (filter(f,x))
listdir_full          =  lambda x  :  sorted (pjoin(x,y) for y in os_listdir(x))
listdir_full_files    =  lambda x  :  lfilter(isfile, listdir_full(x))
listdir_full_folders  =  lambda x  :  lfilter(isdir , listdir_full(x))

def writer(fname, A):
    with open(fname, 'w') as f:
        for x in A:
            f.write(x.rstrip('\n') + '\n')

def reader(fname):
    with open(fname,'r') as f:
        return [x.strip('\n') for x in f]

# ------------------------------------------------------------------------
#                  Fortran code
# ------------------------------------------------------------------------

def get_fort_code(params):
    result = """
module fav_utils
  contains
  subroutine read_arr(arr,fname)
    integer           ::  nx, nz, ny, i, j, k
    character(len=*)  ::  fname
    real*8            ::  arr(0:,0:,0:)
    open(19,file=fname, action='read',status='old',form='unformatted',access='stream')
      read(19) nx, nz, ny
      do      k=0,ny-1
        do    j=0,nz-1
          do  i=0,nx-1
            read(19) arr(i,j,k)
          end do
        end do
      end do
    close(19)
  end subroutine

  subroutine write_arr(fname,arr,nx,nz,ny)
    integer           ::  nx, nz, ny, i, j, k
    character(len=*)  ::  fname
    real*8            ::  arr(0:,0:,0:)
    open(19,file=fname, action='write', status='replace', form='unformatted',access='stream')
      write(19) nx, nz, ny
      do      k=0,ny-1
        do    j=0,nz-1
          do  i=0,nx-1
            write(19) arr(i,j,k)
          end do
        end do
      end do
    close(19)
  end subroutine

    subroutine write_dadx_fav_rab(fname,temp_arr,Ra,Rb,R,Rab,dx,nx,nz,ny)
      integer           ::  nx, nz, ny, i, j, k
      character(len=*)  ::  fname
      real*8            ::  Ra(0:,0:,0:), Rb(0:,0:,0:), R(0:,0:,0:), Rab(0:,0:,0:), dx, &
                            fav_a_xp, fav_a_xm, fav_a, fav_b, fav_rab, temp_arr(0:,0:,0:)
      do      k=0,ny-1
        do    j=0,nz-1
          do  i=0,nx-1
            fav_a_xp         =  Ra(modulo(i+1,nx),j,k)/R(modulo(i+1,nx),j,k)
            fav_a_xm         =  Ra(modulo(i-1,nx),j,k)/R(modulo(i-1,nx),j,k)
            fav_a            =  Ra(i,j,k)/R(i,j,k)
            fav_b            =  Rb(i,j,k)/R(i,j,k)
            fav_rab          =  Rab(i,j,k) - R(i,j,k)*fav_a*fav_b
            temp_arr(i,j,k)  =  -(fav_a_xp - fav_a_xm)/(2*dx)*fav_rab
          end do
        end do
      end do
      call write_arr(fname,temp_arr,nx,nz,ny)
    end subroutine

    subroutine write_dadz_fav_rab(fname,temp_arr,Ra,Rb,R,Rab,dz,nx,nz,ny)
      integer           ::  nx, nz, ny, i, j, k
      character(len=*)  ::  fname
      real*8            ::  Ra(0:,0:,0:), Rb(0:,0:,0:), R(0:,0:,0:), Rab(0:,0:,0:), dz, &
                            fav_a_zp, fav_a_zm, fav_a, fav_b, fav_rab, temp_arr(0:,0:,0:)
      do      k=0,ny-1
        do    j=0,nz-1
          do  i=0,nx-1
            fav_a_zp         =  Ra(i,modulo(j+1,nz),k)/R(i,modulo(j+1,nz),k)
            fav_a_zm         =  Ra(i,modulo(j-1,nz),k)/R(i,modulo(j-1,nz),k)
            fav_a            =  Ra(i,j,k)/R(i,j,k)
            fav_b            =  Rb(i,j,k)/R(i,j,k)
            fav_rab          =  Rab(i,j,k) - R(i,j,k)*fav_a*fav_b
            temp_arr(i,j,k)  =  -(fav_a_zp - fav_a_zm)/(2*dz)*fav_rab
          end do
        end do
      end do
      call write_arr(fname,temp_arr,nx,nz,ny)
    end subroutine

    subroutine write_dady_fav_rab(fname,temp_arr,Ra,Rb,R,Rab,yp_ext,nx,nz,ny)
      integer           ::  nx, nz, ny, i, j, k
      character(len=*)  ::  fname
      real*8            ::  Ra(0:,0:,0:), Rb(0:,0:,0:), R(0:,0:,0:), Rab(0:,0:,0:), &
                            fav_a_yp, fav_a_ym, fav_a, fav_b, fav_rab, temp_arr(0:,0:,0:), yp_ext(0:)
      do      k=0,ny-1
        do    j=0,nz-1
          do  i=0,nx-1
            if (k==0) then      ; fav_a_ym = 0
            else                ; fav_a_ym  = Ra(i,j,k-1) / R(i,j,k-1)
            end if
            if (k==(ny-1)) then ; fav_a_yp  = 0
            else                ; fav_a_yp  = Ra(i,j,k+1) / R(i,j,k+1)
            end if
            fav_a            =  Ra(i,j,k)/R(i,j,k)
            fav_b            =  Rb(i,j,k)/R(i,j,k)
            fav_rab          =  Rab(i,j,k) - R(i,j,k)*fav_a*fav_b
            temp_arr(i,j,k)  =  -(fav_a_yp - fav_a_ym)/(yp_ext(k+2)-yp_ext(k))*fav_rab
          end do
        end do
      end do
      call write_arr(fname,temp_arr,nx,nz,ny)
    end subroutine
end module fav_utils

program main
    use fav_utils
    implicit none
    integer              ::  nx, nz, ny, i, j, k
    real*8, allocatable  ::  R(:,:,:), &
                             Ra(:,:,:), &
                             Rb(:,:,:), &
                             a(:,:,:), &
                             b(:,:,:), &
                             ab(:,:,:), &
                             Rab(:,:,:), temp_arr(:,:,:), yu(:), yp_ext(:),dx,dz
   nx = [nx]
   nz = [nz]
   ny = [ny]
   allocate(R(0:nx-1,0:nz-1,0:ny-1))
   allocate(Ra       ,mold=R)
   allocate(Rb       ,mold=R)
   allocate(a        ,mold=R)
   allocate(b        ,mold=R)
   allocate(ab       ,mold=R)
   allocate(Rab      ,mold=R)
   allocate(temp_arr ,mold=R)
   call read_arr(R   ,'[file_R]')
   call read_arr(Ra  ,'[file_Ra]')
   call read_arr(Rb  ,'[file_Rb]')
   call read_arr(a   ,'[file_a]')
   call read_arr(b   ,'[file_b]')
   call read_arr(ab  ,'[file_ab]')
   call read_arr(Rab ,'[file_Rab]')

  block
    real*8 :: fav_a, fav_b, fav_rab, fav_ab, fav_rab_Rm, fav_rab_Rf

    if ([bool_write_fav_a]) then
      do      k=0,ny-1
        do    j=0,nz-1
          do  i=0,nx-1
            fav_a           = Ra(i,j,k)/R(i,j,k)
            temp_arr(i,j,k) = fav_a
          end do
        end do
      end do
      call write_arr('[file_fav_a]',temp_arr,nx,nz,ny)
    end if

    if ([bool_write_fav_b]) then
      do      k=0,ny-1
        do    j=0,nz-1
          do  i=0,nx-1
            fav_b           = Rb(i,j,k)/R(i,j,k)
            temp_arr(i,j,k) = fav_b
          end do
        end do
      end do
      call write_arr('[file_fav_b]',temp_arr,nx,nz,ny)
    end if

    do      k=0,ny-1
      do    j=0,nz-1
        do  i=0,nx-1
          fav_a           =  Ra(i,j,k)/R(i,j,k)
          fav_b           =  Rb(i,j,k)/R(i,j,k)
          fav_rab         =  Rab(i,j,k) - R(i,j,k)*fav_a*fav_b
          temp_arr(i,j,k) =  fav_rab
        end do
      end do
    end do
    call write_arr('[file_fav_rab]',temp_arr,nx,nz,ny)

    do      k=0,ny-1
      do    j=0,nz-1
        do  i=0,nx-1
          fav_a            =  Ra(i,j,k)/R(i,j,k)
          fav_b            =  Rb(i,j,k)/R(i,j,k)
          fav_rab_Rm       =  R(i,j,k) * (ab(i,j,k) - fav_a*b(i,j,k) - fav_b*a(i,j,k) + fav_a*fav_b)
          temp_arr(i,j,k)  =  fav_rab_Rm
        end do
      end do
    end do
    call write_arr('[file_fav_rab_Rm]',temp_arr,nx,nz,ny)

    do      k=0,ny-1
      do    j=0,nz-1
        do  i=0,nx-1
          fav_a            =  Ra(i,j,k)/R(i,j,k)
          fav_b            =  Rb(i,j,k)/R(i,j,k)
          fav_rab          =  Rab(i,j,k) - R(i,j,k)*fav_a*fav_b
          fav_rab_Rm       =  R(i,j,k) * (ab(i,j,k) - fav_a*b(i,j,k) - fav_b*a(i,j,k) + fav_a*fav_b)
          fav_rab_Rf       =  fav_rab - fav_rab_Rm
          temp_arr(i,j,k)  =  fav_rab_Rf
        end do
      end do
    end do
    call write_arr('[file_fav_rab_Rf]',temp_arr,nx,nz,ny)
  end block
    dx = [dx]
    dz = [dz]
    allocate(yu(0:ny), yp_ext(0:ny+1))
    open(19,file='[file_yu]')
      do k=0,ny
        read(19,*) yu(k)
      end do
    close(19)
    yp_ext(0   ) = yu(0 )
    yp_ext(ny+1) = yu(ny)
    do k=0,ny-1
      yp_ext(k+1) = (yu(k)+yu(k+1))/2.
    end do

    if ([bool_dadx_fav_rab]) then ; call write_dadx_fav_rab(&\n'[file_dadx_fav_rab]',&\n temp_arr,Ra,Rb,R,Rab,dx,nx,nz,ny) ; end if
    if ([bool_dbdx_fav_rab]) then ; call write_dadx_fav_rab(&\n'[file_dbdx_fav_rab]',&\n temp_arr,Rb,Ra,R,Rab,dx,nx,nz,ny) ; end if

    if ([bool_dadz_fav_rab]) then ; call write_dadz_fav_rab(&\n'[file_dadz_fav_rab]',&\n temp_arr,Ra,Rb,R,Rab,dz,nx,nz,ny) ; end if
    if ([bool_dbdz_fav_rab]) then ; call write_dadz_fav_rab(&\n'[file_dbdz_fav_rab]',&\n temp_arr,Rb,Ra,R,Rab,dz,nx,nz,ny) ; end if

    if ([bool_dady_fav_rab]) then ; call write_dady_fav_rab(&\n'[file_dady_fav_rab]',&\n temp_arr,Ra,Rb,R,Rab,yp_ext,nx,nz,ny) ; end if
    if ([bool_dbdy_fav_rab]) then ; call write_dady_fav_rab(&\n'[file_dbdy_fav_rab]',&\n temp_arr,Rb,Ra,R,Rab,yp_ext,nx,nz,ny) ; end if

  deallocate(R, Ra, Rb, Rab, a, b, ab, temp_arr, yu, yp_ext)
end program
"""
    assert params['[bool_write_fav_a]'] in ['.true.', '.false.']
    assert params['[bool_write_fav_b]'] in ['.true.', '.false.']
    for key in ['[file_R]',
               '[file_Ra]',
               '[file_Rb]',
               '[file_Rab]',
               '[file_a]',
               '[file_b]',
               '[file_ab]',
               '[bool_write_fav_a]',
               '[bool_write_fav_b]',
               '[file_fav_a]'      ,
               '[file_fav_b]'      ,
               '[file_fav_rab]'    ,
               '[file_fav_ab]'     ,
               '[file_fav_rab_Rm]' ,
               '[file_fav_rab_Rf]' ,
               '[file_fav_ab_Rm]'  ,
               '[file_fav_ab_Rf]'  ,
               '[file_yu]'         ,
               '[dx]'              ,
               '[dz]'              ,
               '[nx]'              ,
               '[nz]'              ,
               '[ny]'              ,
               '[bool_dadx_fav_rab]', '[file_dadx_fav_rab]',
               '[bool_dbdx_fav_rab]', '[file_dbdx_fav_rab]',
               '[bool_dadz_fav_rab]', '[file_dadz_fav_rab]',
               '[bool_dbdz_fav_rab]', '[file_dbdz_fav_rab]',
               '[bool_dady_fav_rab]', '[file_dady_fav_rab]',
               '[bool_dbdy_fav_rab]', '[file_dbdy_fav_rab]']:
        result = result.replace(key,str(params[key]))
    assert not '[' in result
    assert not ']' in result
    return result.split('\n')

# ------------------------------------------------------------------------
#                  Bulk calculations (Fortran)
# ------------------------------------------------------------------------

def main_runner_fav_avg():
    runtime_folder =  dirname(_folder_)
    fname_fav_avg  =  pjoin(_folder_, 'mk_temp_fav_avg.f90')
    all_iter_avg   =  sorted(set([int(basename(x).split('iter_')[1].split('.')[0]) for x in listdir_full_files(pjoin(runtime_folder,'avg')) if x.endswith('.dat')]))
    full_tags      =  sorted(set(basename(x).split('array_arr_3d_')[1].split('_avg_iter_')[0] for x in listdir_full_files(pjoin(runtime_folder,'avg'))))

    memo_written = set()
    def check_write(fname):
        result = '.true.' if (not (basename(fname) in memo_written)) else '.false.'
        memo_written.add(basename(fname))
        return result

    def get_fname_avg(iter_avg, tags):
        result = f"../avg/array_arr_3d_{''.join(sorted(tags))}_avg_iter_{iter_avg:09d}.dat"
        assert isfile(pjoin(_folder_,result))
        return result

    params        =  lfilter(None, [x.split('!')[0].split('#')[0].strip(' \n') for x in reader(pjoin(_folder_,"../../input/params.dat"))])
    params        =  {key.strip(' '): eval(val) for row in params for key,val in [row.split('=')]}
    if not ('dx' in params): params['dx']  =  params['Lx']/params['nx']
    if not ('dz' in params): params['dz']  =  params['Lz']/params['nz']
    as_ftr_bool = lambda x: '.true.' if x else '.false.'
    memo_pTK = set()
    for iter_avg in all_iter_avg:
        for tag in full_tags:
            if len(tag) == 3 and ('R' in tag):
                i_R = [i for i,c in enumerate(tag) if c=='R'][0]
                a,b = [tag[(i_R+d)%3] for d in [1,2]]
                file_fav_a      = f"../avg/array_arr_3d_fav_{a}_avg_iter_{iter_avg:09d}.dat"
                file_fav_b      = f"../avg/array_arr_3d_fav_{b}_avg_iter_{iter_avg:09d}.dat"
                file_fav_rab    = f"../avg/array_arr_3d_fav_R{a}{b}_avg_iter_{iter_avg:09d}.dat"
                file_fav_ab     = f"../avg/array_arr_3d_fav_{a}{b}_avg_iter_{iter_avg:09d}.dat"
                file_fav_rab_Rm = f"../avg/array_arr_3d_fRm_R{a}{b}_avg_iter_{iter_avg:09d}.dat"
                file_fav_rab_Rf = f"../avg/array_arr_3d_fRf_R{a}{b}_avg_iter_{iter_avg:09d}.dat"
                file_fav_ab_Rm  = f"../avg/array_arr_3d_fRm_{a}{b}_avg_iter_{iter_avg:09d}.dat"
                file_fav_ab_Rf  = f"../avg/array_arr_3d_fRf_{a}{b}_avg_iter_{iter_avg:09d}.dat"
                if not all((c in 'UVW') for c in [a,b]): 
                    memo_pTK.add((a,b)) ; memo_pTK.add((b,a))
                need_pTK_ab = not ((a,b) in memo_pTK) ; memo_pTK.add((a,b))
                need_pTK_ba = not ((b,a) in memo_pTK) ; memo_pTK.add((b,a))
                code_fort = get_fort_code({'[file_R]'           : get_fname_avg(iter_avg,['R'    ]) ,
                                           '[file_Ra]'          : get_fname_avg(iter_avg,['R',a  ]) ,
                                           '[file_Rb]'          : get_fname_avg(iter_avg,['R',b  ]) ,
                                           '[file_Rab]'         : get_fname_avg(iter_avg,['R',a,b]) ,
                                           '[file_a]'           : get_fname_avg(iter_avg,[a  ]) ,
                                           '[file_b]'           : get_fname_avg(iter_avg,[b  ]) ,
                                           '[file_ab]'          : get_fname_avg(iter_avg,[a,b]) ,
                                           '[bool_write_fav_a]' : check_write(file_fav_a)       ,
                                           '[bool_write_fav_b]' : check_write(file_fav_b)       ,
                                           '[file_fav_a]'       : file_fav_a                    ,
                                           '[file_fav_b]'       : file_fav_b                    ,
                                           '[file_fav_rab]'     : file_fav_rab                  ,
                                           '[file_fav_ab]'      : file_fav_ab                   ,
                                           '[file_fav_rab_Rm]'  : file_fav_rab_Rm               ,
                                           '[file_fav_rab_Rf]'  : file_fav_rab_Rf               ,
                                           '[file_fav_ab_Rm]'   : file_fav_ab_Rm                ,
                                           '[file_fav_ab_Rf]'   : file_fav_ab_Rf                ,
                                           '[file_yu]'          : "../../input/yu.txt"          ,
                                           '[dx]'               : str(params['dx'])             ,
                                           '[dz]'               : str(params['dz'])             ,
                                           '[nx]'               : str(params['nx'])             ,
                                           '[nz]'               : str(params['nz'])             ,
                                           '[ny]'               : str(params['ny'])             ,
               '[bool_dadx_fav_rab]': as_ftr_bool((b=='U') and need_pTK_ab) , '[file_dadx_fav_rab]': f"../avg/array_arr_3d_pTK_{a}{b}_avg_iter_{iter_avg:09d}.dat",
               '[bool_dbdx_fav_rab]': as_ftr_bool((a=='U') and need_pTK_ba) , '[file_dbdx_fav_rab]': f"../avg/array_arr_3d_pTK_{b}{a}_avg_iter_{iter_avg:09d}.dat",
               '[bool_dadz_fav_rab]': as_ftr_bool((b=='W') and need_pTK_ab) , '[file_dadz_fav_rab]': f"../avg/array_arr_3d_pTK_{a}{b}_avg_iter_{iter_avg:09d}.dat",
               '[bool_dbdz_fav_rab]': as_ftr_bool((a=='W') and need_pTK_ba) , '[file_dbdz_fav_rab]': f"../avg/array_arr_3d_pTK_{b}{a}_avg_iter_{iter_avg:09d}.dat",
               '[bool_dady_fav_rab]': as_ftr_bool((b=='V') and need_pTK_ab) , '[file_dady_fav_rab]': f"../avg/array_arr_3d_pTK_{a}{b}_avg_iter_{iter_avg:09d}.dat",
               '[bool_dbdy_fav_rab]': as_ftr_bool((a=='V') and need_pTK_ba) , '[file_dbdy_fav_rab]': f"../avg/array_arr_3d_pTK_{b}{a}_avg_iter_{iter_avg:09d}.dat",
                                           })
                writer(fname_fav_avg,code_fort)
                commands=f"cd {_folder_};gfortran {basename(fname_fav_avg)} -o a.out;./a.out"
                print(f"Running Favre-averaging: {tag}")
                t0 = time()
                os_system(commands)
                print(f"Done! (Favre-averaging: {tag}) (elapsed_time: {format_dt(time()-t0)})")
    commands=f"cd {_folder_};rm -f a.out;rm -f {basename(fname_fav_avg)}"
    os_system(commands)

# ------------------------------------------------------------------------
#                  Direct Runner
# ------------------------------------------------------------------------

if __name__ == "__main__":
    t0 = time()
    print('----> Begin: main_runner_fav_avg()')
    main_runner_fav_avg()
    print(f'----> End: main_runner_fav_avg() (elapsed_time: {format_dt(time()-t0)})')


