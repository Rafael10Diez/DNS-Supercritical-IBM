# ------------------------------------------------------------------------
#                  Notes
# ------------------------------------------------------------------------
#
# Motivation:
#    - Initializing a DNS simulation from scratch (flat plate channel flow or U=0) can be expensive.
#    - The development of flow detachment near the rough surfaces can require very small time steps, and
#      long flow-through times.
#
#    - To avoid these issues, this script does the following:
#      1) Read the final state of a previous DNS case (closely related).
#      2) Re-interpolate the solution to the current grid (fitted within the rough surface boundaries).
#      3) Using previous rough surface data as input is useful, even if the initial data used a different rough surface, or
#         Reynolds/Prandtl number.
#      4) Ideally, the initial data should fulfill the periodic buboundaryondary conditions. However, this has not been found to be mandatory, because
#         the DNS solver quickly reaches equilibrum in either case.
#
# Please note that:
#    1) This script can still be RAM memory intensive, even if it uses less memory than previous versions.
#
#    2) It is important to use at least quadratic interpolation for the DNS data.
#       Otherwise, the DNS data can be non-physical (for linear interpolation), because mu*laplacian(U)=0.
#       This implies that the diffusivity in the Navier-Stokes equations would be (initially) missing, 
#       and sometimes it can lead to divergence.
# 
# ------------------------------------------------------------------------
#                  Generic Libraries
# ------------------------------------------------------------------------

from    sys         import  path                                       as  sys_path
from    os          import  listdir                                    as  os_listdir
from    os.path     import  join                                       as  pjoin
from    os          import  system                                     as  os_system 
import  numpy                                                          as  np
from    os.path     import  basename, dirname, abspath, isdir, isfile
import  zipfile
from    time        import  time
from    sys         import  argv
import  subprocess

# ------------------------------------------------------------------------
#                  Basic Definitions
# ------------------------------------------------------------------------

_folder_            =  dirname(abspath(__file__))
lmap                =  lambda f,x:  list(map(f,x))
lfilter             =  lambda f,x:  list(filter(f,x))
deepdirname         =  lambda x,n:  x if n<1 else deepdirname(dirname(x), n-1)

format_dt           =  lambda x :  "%02d:%02d:%02d [hh:mm:ss]" % (x//3600, (x%3600)//60, x%60)

listdir_full        =  lambda   x:  sorted(pjoin(x,y) for y in os_listdir(x))
listdir_full_files  =  lambda   x:  lfilter(isfile, listdir_full(x))

tolist              =  lambda x: x.tolist() if hasattr(x,'tolist') else x

def pop1(A):
    assert type(A) == list
    x, = A 
    return x

def reader(fname):
    with open(fname, 'r') as f:
        return [x.rstrip('\n') for x in f]

def writer(fname, A):
    assert type(A) == list
    with open(fname, 'w') as f:
        for x in A: f.write(x.strip('\n')+'\n')

def import_fort_dat(fname):
    result = {}
    for line in lfilter(None, [x.split('#')[0].split('!')[0].strip(' ') for x in reader(fname)]):
        assert line.count('=') == 1
        key,val = [x.strip(' ') for x in line.split('=')]
        if val in ['.true.','.false.']: val = val.strip('.').title()
        result[key] = eval(val)
    return result

def mk_pu(L,n):
    u = np.linspace(0,L,n+1)
    return u[:-1] + 0.5*np.diff(u), u[:-1]

def load_params(fname):
    return  eval(' '.join(x.split('#')[0] for x in reader(fname)))

def prod(A):
    y = 1 + (0*A[0])
    for x in A:
        y *= x
    return y

# ------------------------------------------------------------------------
#                  Run Commands in Parallel
# ------------------------------------------------------------------------

def run_commands(all_commands, nproc):
    assert nproc >= 1
    running       =  []
    all_commands  =  list(all_commands)[::-1]
    while all_commands or running:
        # run pending command, if processes are available
        if all_commands and (len(running) < nproc):  
            running.append(subprocess.Popen(all_commands.pop(),shell=True))
        else:
            for i in reversed(range(len(running))):
                # remove process, if finished
                if not (running[i].poll() is None):           
                    print(f"Finished running: (running[{i}].args = {running[i].args})")
                    running.pop(i)
    print('Done! (@run_commands)')

# ------------------------------------------------------------------------
#                  Basic Data
# ------------------------------------------------------------------------

seen_codes = set()
class Basic_Data:
    def __init__(self, input_folder):
        
        self.input_folder =  input_folder
        self.param        =  import_fort_dat(pjoin(input_folder, 'params.dat'))
        self.sheight      =  self.import_height_function(input_folder)

        self.xp, self.xu  =  mk_pu(self.param['Lx'], self.param['nx'])
        self.zp, self.zu  =  mk_pu(self.param['Lz'], self.param['nz'])

        self.yu           =  np.array(lmap(float,lfilter(None,reader(pjoin(input_folder, 'yu.txt')))))
        self.yp           =  self.yu[:-1] + 0.5*np.diff(self.yu)
        self.yu           =  self.yu[:-1]

        assert  len(self.xp)  == len(self.xu) == self.param['nx']
        assert  len(self.zp)  == len(self.zu) == self.param['nz']
        assert  len(self.yu)  == len(self.yp) == self.param['ny']

        self.xu  =  self.xu.reshape(-1)
        self.xp  =  self.xp.reshape(-1)
        self.zu  =  self.zu.reshape(-1)
        self.zp  =  self.zp.reshape(-1)
        self.yu  =  self.yu.reshape(-1)
        self.yp  =  self.yp.reshape(-1)



        assert all(abs(self.param[f'L{c}']-getattr(self.sheight,f'L{c}'))<1e-10 for c in 'xyz')

        self.nreps      = {c: getattr(self.sheight,f"nreps_{c}") for c in 'xz'}
        self.nreps['y'] = 1

        for side in ['bottom', 'top']:
            for cx,cz in ['up','pu','pp']:
                x       =  getattr(self, f'x{cx}').reshape(-1, 1)
                z       =  getattr(self, f'z{cz}').reshape( 1,-1)
                setattr(self, f'h_{side}_{cx}{cz}', self.sheight.get_h(x, z, side, lib=np))

    @staticmethod
    def import_height_function(input_folder):
        # why this?
        #     -  because two different files with the same name cannot be imported.
        #     -  the second import (Height_Function) will be ignored, and replaced by the first (Height_Function).
        orig_fname  =  pjoin(input_folder,'rough_surf_define.py')
        new_fname   =  lambda i: pjoin(input_folder,f'copy_{i}_rough_surf_define.py')
        ii          =  0
        while ii in seen_codes: ii += 1
        seen_codes.add(ii)
        assert ii in [0,1]
        os_system(f'cp -f "{orig_fname}" "{new_fname(ii)}"')
        A     =  reader(new_fname(ii))
        jj,   =  [k for k in range(len(A)) if 'class Height_Function:'==A[k][:22]]
        A[jj] =  f'class Height_Function_copy_{ii}:'
        A     =  [x.replace('Height_Function(', f'Height_Function_copy_{ii}(') for x in A]
        A     =  ['# ********************** File Automatically Generated by Python (please ignore or erase this file) **********************\n',*A]
        writer(new_fname(ii), A)
        sys_path.append(input_folder)
        if ii == 0:
            from  copy_0_rough_surf_define  import  Height_Function_copy_0
            sheight = Height_Function_copy_0()
        if ii == 1:
            from  copy_1_rough_surf_define  import  Height_Function_copy_1
            sheight = Height_Function_copy_1()
        sys_path.pop()
        return sheight
    
    def write_ginfo(self, fname, tag):
        mm         =  dict(U = 'upp',
                           W = 'pup', 
                           V = 'ppu',
                           T = 'ppp')[tag]
        x,z,y      =  [getattr(self,f"{c}{m}").reshape(-1) for c,m in zip('xzy',mm)]
        nx,nz,ny   =  len(x), len(z), len(y)
        A          =  [f"{nx} {nz} {ny} ! nx nz ny"]
        A_y        =  lmap(str,y.tolist())
        A_y[0]    += ' ! y(0:ny-1)'
        A.extend(A_y)
        h_bot, h_top  =  [getattr(self, f'h_{side}_{mm[0]}{mm[1]}') for side in ['bottom','top']]
        assert h_bot.shape == h_top.shape == (nx,nz), [h_bot.shape, h_top.shape, (nx,nz)]
        for     i in range(nx):
            for j in range(nz):
                A.append(f"{i} {j} {h_bot[i,j]} {h_top[i,j]}")
                if i==j==0:
                    A[-1] += ' ! i j h_bot(i,j) h_top(i,j)'
        writer(fname, A)
    
    def norm(self):
        self.xu /= (self.sheight.Lx/self.sheight.nreps_x)
        self.xp /= (self.sheight.Lx/self.sheight.nreps_x)
        self.zu /= (self.sheight.Lz/self.sheight.nreps_z)
        self.zp /= (self.sheight.Lz/self.sheight.nreps_z)
        self.yu /=  self.sheight.Ly
        self.yp /=  self.sheight.Ly
        for side in ['bottom', 'top']:
            for cx,cz in ['up','pu','pp']:
                setattr(self, f'h_{side}_{cx}{cz}', 
                getattr(self, f'h_{side}_{cx}{cz}') / self.sheight.Ly)

# ------------------------------------------------------------------------
#                  Smooth Wall Data
# ------------------------------------------------------------------------

class Load_Smooth_Wall_Data(Basic_Data):
    def __init__(self, input_folder):
        
        super().__init__(input_folder)
        self._smooth_input_folder  =  input_folder
        locate_file                =  lambda folder, tag: pop1(lfilter(lambda x: tag in basename(x), 
                                                             listdir_full_files(folder)))
        self.fname_array_ijk       =  {tag: locate_file(self._smooth_input_folder, f'array_ini_{tag}.') for tag in 'UVWT'}

# ------------------------------------------------------------------------
#                  Adjusted Target Coordinates
# ------------------------------------------------------------------------

class Adjustable_Target_Data(Basic_Data):
    def __init__(self, input_folder, first_quadrant = False):
        
        super().__init__(input_folder)

        if first_quadrant:
            for c in 'xz':
                L = self.param[f'L{c}']/self.nreps[c]
                for t in 'up':
                    vname = f"{c}{t}"
                    setattr(self, vname, getattr(self,vname)%L)

# ------------------------------------------------------------------------
#                  reinterp_field
# ------------------------------------------------------------------------

extrap_1d = lambda x: [2*x[0] - x[1]] + tolist(x) + [2*x[-1] - x[-2]]

def get_per_inds(xt, xs, nper = 10):
    # build indexes and coefficients to interpolate from "xs" (source) to "xt" (target)
    # periodic boundary conditions are considered
    xt, xs   = np.array(tolist(xt)), np.array(tolist(xs))
    Ls_orig  =  len(xs)
    for _ in range(nper): xs  =  extrap_1d(xs)
    inds, coeffs  =  [], []
    for i in range(len(xt)):
        all_j = np.argmin(np.fabs(xs-xt[i]))
        all_j = np.array([all_j-1, all_j, all_j+1])
        assert (min(all_j)>=0) and (max(all_j)<len(xs)),(all_j,len(xs), xt, xs)
        all_x  = np.array([xs[j]-xt[i] for j in all_j])
        all_x /= np.fabs(np.diff(all_x)).max()
        coeffs.append(np.linalg.pinv(np.array([all_x**i for i in range(3)]).transpose())[0,:])
        inds  .append((all_j-nper)%Ls_orig)
    return inds, coeffs # [n,3]

def mk_coeffs_interp(smooth, target, tag, all_fname):
    mm        =  dict(U = 'upp',
                      W = 'pup', 
                      V = 'ppu',
                      T = 'ppp')[tag]
    for ind,c in enumerate('xz'):
        fname        =  all_fname[c]
        t            =  getattr(target,f"{c}{mm[ind]}").reshape(-1)
        s            =  getattr(smooth,f"{c}{mm[ind]}").reshape(-1)
        inds,coeffs  =  get_per_inds(t, s)
        n            =  len(t)
        assert len(inds) == len(coeffs) == n
        A            = [str(n)]
        for j in range(n):
            A.append(' '.join(str(val) for i,c in zip(inds[j],coeffs[j]) for val in [i,c]))
        writer(fname, A)

# ------------------------------------------------------------------------
#                  Main Program
# ------------------------------------------------------------------------

def main_reinterpolator(ref_smooth_subfolder):
    fort_compiler   = 'gfortran'
    t00             =  time()
    first_quadrant  =  True
    smooth          =  Load_Smooth_Wall_Data(ref_smooth_subfolder)
    target          =  Adjustable_Target_Data(dirname(abspath(__file__)), first_quadrant = first_quadrant)
    smooth.norm()
    target.norm()
    all_commands = []
    for tag in 'UVWT':
        all_commands.append([])
        add_command = lambda x: all_commands[-1].append(x)
        work_dir        =  pjoin(_folder_,f'mk_ini_{tag}')
        os_system(f"rm -rf {work_dir}")
        os_system(f"mkdir  {work_dir}")

        # copy & compile fortran
        fname_src_fortran   =  pjoin(work_dir, 'mk_ini_data_sparser_bin.f90')
        fname_exec_fortran  =  pjoin(work_dir, 'a.out')

        os_system(f"cp -f {pjoin(_folder_,basename(fname_src_fortran))} {fname_src_fortran}")
        #assert isfile(fname_src_fortran)
        os_system(f"{fort_compiler} {fname_src_fortran} -o {fname_exec_fortran}")
        #assert isfile(fname_exec_fortran)

        src_tag_fname = pjoin(work_dir, f'data_{tag}', f'source_array_{tag}.bin')
        os_system(f"mkdir -p {dirname(src_tag_fname)}")
        add_command(f"cp -f {smooth.fname_array_ijk[tag]} {src_tag_fname}")
        #assert isfile(src_tag_fname)

        smooth.write_ginfo(pjoin(work_dir, f'data_{tag}',f'gs_info_{tag}.dat'), tag)
        target.write_ginfo(pjoin(work_dir, f'data_{tag}',f'gt_info_{tag}.dat'), tag)
        mk_coeffs_interp(smooth, target, tag, dict(x = pjoin(work_dir, f'data_{tag}',f'coeffs_{tag}_ic_x.dat'),
                                                   z = pjoin(work_dir, f'data_{tag}',f'coeffs_{tag}_jc_z.dat')))
        add_command(f"cd {dirname(fname_exec_fortran)};{fname_exec_fortran} {tag}")
        target_fname = pjoin(work_dir, f'data_{tag}', f'target_array_{tag}.bin')
        #assert isfile(target_fname), target_fname
        add_command(f"mv {target_fname} {pjoin(_folder_,basename(smooth.fname_array_ijk[tag]))}")

        add_command(f"rm -rf {work_dir}")
        print(f'Done! (script: {__file__}) (total elapsed time: {format_dt(time()-t00)})')
        all_commands[-1] = ';'.join(x.strip(' ') for x in all_commands[-1])
    run_commands(all_commands, 4)

# ------------------------------------------------------------------------
#                  Direct Runner
# ------------------------------------------------------------------------

if __name__ == '__main__':
    main_reinterpolator(argv[1])