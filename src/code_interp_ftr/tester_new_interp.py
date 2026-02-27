
# ------------------------------------------------------------------------
#                  Generic Libraries
# ------------------------------------------------------------------------

import  sympy                                                       as  sp
import  numpy                                                       as  np
from    os.path  import  join                                       as  pjoin
from    os.path  import  basename, dirname, isfile, isdir, abspath
from    os       import  listdir                                    as  os_listdir
import  random
import  struct
from    os       import  system                                     as  os_system
random.seed(0)
np.random.seed(0)

# ------------------------------------------------------------------------
#                  Basic Functions
# ------------------------------------------------------------------------

_folder_              =  dirname(abspath(__file__))

lmap                  =  lambda f,x: list(map(f,x))
lfilter               =  lambda f,x: list(filter(f,x))
listdir_full          =  lambda   x: sorted([pjoin(x,y) for y in os_listdir(x)])
listdir_full_files    =  lambda   x: lfilter(isfile, listdir_full(x))
listdir_full_folders  =  lambda   x: lfilter(isdir, listdir_full(x))
rand_between          =  lambda a,b: a + (b-a)*float(np.random.rand())

def writer(fname, A):
    if type(A) == str: A = A.split('\n')
    os_system(f'mkdir -p "{dirname(fname)}"')
    with open(fname,'w') as f:
        for x in A:
            f.write(x.rstrip('\n')+'\n')

def reader(fname):
    with open(fname,'r') as f:
        return [x.rstrip('\n') for x in f]

# ------------------------------------------------------------------------
#                  Utilities
# ------------------------------------------------------------------------

def reshape_arrs(result,g):
    if not (g is None):
        if g=='x': return [arr.reshape(-1, 1, 1) for arr in result]
        if g=='z': return [arr.reshape( 1,-1, 1) for arr in result]
        if g=='y': return [arr.reshape( 1, 1,-1) for arr in result]
    return result

def mk_grid_pu(L,n):
    u = np.linspace(0,L,n+1)
    p = u[:-1] + 0.5*np.diff(u)
    return [p, u[:-1]]

def mk_y_pu(L, n, growth=1.1):
    y   =  np.array([0.] + np.cumsum(np.linspace(1,growth,n//2)).tolist())
    y  *=  (L/2)/y.max()
    y   =  y.tolist()
    for _ in range(random.randint(1,3)):
        y = [2*y[0] - y[1]] + y
    y = np.array( y + (L-np.array(y)[:-1][::-1]).tolist() )
    return [y[:-1] + 0.5*np.diff(y), y]

def writerr_3d_arr(fname, A, args):
    shape  =  args['nx'], args['nz'], args['ny']
    os_system(f'mkdir -p "{dirname(fname)}"')
    with open(fname,'wb') as fw:
        [fw.write(struct.pack('i',s)) for s in shape]
        for          k in range(shape[2]):
            for      j in range(shape[1]):
                for  i in range(shape[0]):
                    fw.write(struct.pack('d',float(A[i,j,k])))

def reader_2d_arr_bin(fname):
    with open(fname, "rb") as f:
        n0,n1  =  [ struct.unpack('i',f.read(4))[0] for _ in range(2)]
        A      =  np.zeros((n0,n1), dtype=float)
        for      j in range(n1):
            for  i in range(n0):
                A[i,j]  =  struct.unpack('d',f.read(8))[0]
    return A

# ------------------------------------------------------------------------
#                  Exact Solution
# ------------------------------------------------------------------------

def mk_exact_sol():
    params           =  {f"L{c}":rand_between(0.75,1.25) for c in 'xyz'}
    params.update({f'wall_BC_{key}_{side}':rand_between(0.5,1.5) for key in 'U,V,W,T,P,mu,cond,rho'.split(',') for side in ['top','bot']})
    _2pi             =  2*float(np.pi)
    rand_sin         =  lambda x_adim,allow_phase=True: (rand_between(0.7,1.3)*int(allow_phase) + rand_between(0.5,0.8)*sp.sin((_2pi)*x_adim + rand_between(-2,2)*int(allow_phase)))
    x,y,z            =  sp.symbols('x,y,z',real=True,positive=True,nonzero=True)
    x_adim           =  x/params['Lx'] 
    z_adim           =  z/params['Lz']
    delta            =  params['Ly']/2
    mk_h             =  lambda: delta*(rand_between(0.3,0.35) + rand_between(0.05,0.07)*sp.sin((_2pi)*x_adim + rand_between(-2,2))*
                                                                                        sp.sin((_2pi)*z_adim + rand_between(-2,2)))
    h_bot            =            mk_h()
    h_top            =  2*delta - mk_h()
    theta_adim       =  (y-h_bot)/(h_top - h_bot)
    mk_field         =  lambda: rand_sin(x_adim) * rand_sin(z_adim) * rand_sin(theta_adim, allow_phase=False)

    mk_lin_bc_var    =  lambda key: ((params[f'wall_BC_{key}_top'] - params[f'wall_BC_{key}_bot'])*theta_adim + 
                                      params[f'wall_BC_{key}_bot'])
    result           =  {key: (mk_field()+mk_lin_bc_var(key)) for key in 'UVWTP'}
    result['h_top']  =  h_top
    result['h_bot']  =  h_bot
    result.update({f"{key}_dd{c0}": result[key].diff(c1) for c0,c1 in zip('xyz',[x,y,z]) for key in result})
    return params, {k:f"lambda x,z,y,sin=np.sin,cos=np.cos: ({v})" for k,v in result.items()}

# ------------------------------------------------------------------------
#                  Height Function
# ------------------------------------------------------------------------

def ftn_line_split(A,lim=70,last_cont=False):
    A     =  str(A).split(' ')
    B     =  []
    line  =  ''
    for x in A:
        if (len(line)+len(x)+3)>=lim:
            B.append(line+' &')
            line = ''
        line  = f"{line} {x}".strip(' ')
    B.append((line+' &') if last_cont else line)
    return '\n'.join(B)

def mk_height_function_frtn(d):
    return ("""
module Mod_Height_Function
  implicit none
  private
  public :: init_surf_obj, type_surf_obj, get_h
  type type_surf_obj
    real*8               :: Lx, Lz, Ly
    integer              :: nreps_x, nreps_z
  end type
  contains
  subroutine init_surf_obj(surf_obj)
    implicit none    
    type(type_surf_obj) :: surf_obj
    surf_obj%Lx       =  [HOLDER_Lx]
    surf_obj%Lz       =  [HOLDER_Lz]
    surf_obj%Ly       =  [HOLDER_Ly]
    surf_obj%nreps_x  =  1
    surf_obj%nreps_z  =  1
  end subroutine 
  function get_h(x, z, surf_obj, is_bottom) result(H)
    implicit none    
    type(type_surf_obj)  ::  surf_obj
    logical              ::  is_bottom
    real*8               ::  x,z,H, theta,m,n,a,b
    integer              ::  i
    real*8, parameter    ::  const_2pi =  8.D0*DATAN(1.D0)
    if (is_bottom) then 
      H = ([HOLDER_H_BOT]) + 0*(x+z)
    else 
      H = ([HOLDER_H_TOP]) + 0*(x+z)
    end if 
  end function
end module
""".replace('[HOLDER_Lx]'     , str(d['Lx']))\
   .replace('[HOLDER_Lz]'     , str(d['Lz']))\
   .replace('[HOLDER_Ly]'     , str(d['Ly']))\
   .replace('[HOLDER_H_BOT]', ftn_line_split(str(d['H_bot'])))\
   .replace('[HOLDER_H_TOP]', ftn_line_split(str(d['H_top']))))

def mk_interp_prop(key,params):
    x0,x1 = params['wall_BC_T_bot'], params['wall_BC_T_top']
    y0,y1 = params[f'wall_BC_{key}_bot'], params[f'wall_BC_{key}_top']
    m     = (y1-y0)/(x1-x0)
    return f"((T_ext_ijk-{x0})*&\n ({m}) +&\n ({y0}))"

# ------------------------------------------------------------------------
#                  Exact Solution
# ------------------------------------------------------------------------

if __name__ == '__main__':
    os_system(f"rm -rf {pjoin(_folder_,'trials')}")

    n_trials = 1
    for i_trial in range(n_trials):
        params, exact_sol  =  mk_exact_sol()
        all_fname_copy     =  [x for x in listdir_full_files(_folder_) if x!=abspath(__file__)]
        n_xzy_base         =  [random.randint(32,36) for _ in range(3)]
        for K_grid in [2,4]:
            trial_subfolder = pjoin(_folder_,'trials', f'trial_{i_trial:09d}',f'refine_{K_grid}')
            args    =  {f"n{c}": (K_grid*n) for c,n in zip('xyz',n_xzy_base)}
            coords  =  {}
            (coords['xp'],coords['xu']),(coords['zp'],coords['zu'])  =  [reshape_arrs(mk_grid_pu(params[f'L{c}'], args[f'n{c}']),c) for c in 'xz']
            coords['yp'],coords['yu']                                =   reshape_arrs(mk_y_pu(params['Ly'], args['ny'], growth=1.1),'y')
            args['ny']                                               =  len(coords['yp'].reshape(-1))
            iter_avg                                                 =  random.randint(1,13)
            def quick_write_avg():
                for tag, mm_xzy in [['U', 'upp'],
                                    ['W', 'pup'],
                                    ['V', 'ppu'],
                                    ['T', 'ppp'],
                                    ['P', 'ppp']]:
                    writerr_3d_arr(pjoin(trial_subfolder,'runtime','avg',f"array_arr_3d_{tag}_avg_iter_{iter_avg:09d}.dat"),
                                  eval(exact_sol[tag])(coords[f"x{mm_xzy[0]}"],coords[f"z{mm_xzy[1]}"],coords[f"y{mm_xzy[2]}"]), args)
            quick_write_avg()

            writer(pjoin(trial_subfolder,'input','rough_surf_define.f90'),
                   mk_height_function_frtn(dict(Lx    = str(params['Lx']) ,
                                                Lz    = str(params['Lz']) ,
                                                Ly    = str(params['Ly']) ,
                                                H_bot = exact_sol['h_bot'].split(':')[-1],
                                                H_top = exact_sol['h_top'].split(':')[-1])))
            writer(pjoin(trial_subfolder,'input','yu.txt'),
                   lmap(str,coords['yu'].reshape(-1).tolist()))
            writer(pjoin(trial_subfolder,'input','params.dat'),
                   [f"nx = {args['nx']}",
                    f"nz = {args['nz']}",
                    f"ny = {args['ny']}",
                    f"Lx = {params['Lx']}",
                    f"Lz = {params['Lz']}",
                    f"Ly = {params['Ly']}",
                    f"use_rough_surf = 1",
                    f"use_rough_surf_quady = 0",
                    *[f"wall_BC_{tag}_{side} = {params[f'wall_BC_{tag}_{side}']}" for tag in 'UVWT' for side in ['top','bot']]])
            [writer(pjoin(trial_subfolder,'runtime','code_interp_ftr',basename(fname)), reader(fname)) for fname in all_fname_copy]
            writer(pjoin(trial_subfolder,'input','fortran_surf_files','empty.txt'), [''])
            writer(pjoin(trial_subfolder,'runtime','thermophysical_properties','expressions_rho_mu_cond_buoyancy_cp_div_jacobian.f90'),
                   [ "T_ext_ijk=T_ext(i,j,k)",
                    f"rho_ijk = {mk_interp_prop('rho',params)}",
                    f"mu_ijk = {mk_interp_prop('mu',params)}",
                    f"cond_ijk = {mk_interp_prop('cond',params)}",
                     "rho_arr(i,j,k)=rho_ijk", "mu_arr(i,j,k)=mu_ijk", "cond_arr(i,j,k)=cond_ijk"])
            os_system(f"cd {pjoin(trial_subfolder,'runtime','code_interp_ftr')} ; bash commands.txt")
            
            def run_verif_surf_hxz():
                max_error = float('-inf')
                for fname in listdir_full_files(pjoin(trial_subfolder,'runtime','code_interp_ftr','data_output','surfs_hxz')):
                    side   = basename(fname).split('_')[-1].split('.')[0]
                    mm     = basename(fname).split('_')[-2]
                    h_ref  = eval(exact_sol[f'h_{side}'])(coords[f"x{mm[0]}"],coords[f"z{mm[1]}"],None)[:,:,0]
                    h_read = reader_2d_arr_bin(fname)
                    max_error = max(max_error, np.fabs(h_ref - h_read).max())
                print(f"run_verif_surf_hxz: {max_error = }")
                assert max_error < 1e-10
            run_verif_surf_hxz()

            def run_verif_interp_2d():
                max_error = float('-inf')
                for fname in listdir_full_files(pjoin(trial_subfolder,'runtime','code_interp_ftr','data_output','interp_2d')):
                    ftype     =  basename(fname).split('_')[-2]
                    side      =  basename(fname).split('_')[-3]
                    tag       =  basename(fname).split('_')[-4]
                    key       =  tag + ('' if ftype == 'val' else f"_{ftype}")
                    mm        =  dict(U = 'up', W = 'pu', V = 'pp', T = 'pp', P = 'pp')[tag]
                    yy        =  eval(exact_sol[f'h_{side}'])(coords[f"x{mm[0]}"],coords[f"z{mm[1]}"],None)
                    arr_ref   =  eval(exact_sol[key])(coords[f"x{mm[0]}"],coords[f"z{mm[1]}"],yy)[:,:,0]
                    arr_read  =  reader_2d_arr_bin(fname)
                    max_error =  max(max_error, np.fabs(arr_ref - arr_read).max())
                    print(basename(fname), np.fabs(arr_ref - arr_read).max())
                print(f"run_verif_interp_2d: {max_error = }")
                # assert max_error < 1e-10
            run_verif_interp_2d()


