# Purpose
#    Full tester for "mk_ini_data_sparser.py"
#    A sub-folder with the tests is created & erased (at the end)

# ------------------------------------------------------------------------
#                  Generic Libraries
# ------------------------------------------------------------------------

from    os.path  import  join                                       as  pjoin
from    os       import  system                                     as  os_system 
import  sympy                                                       as  sp
import  numpy                                                       as  np
from    os.path  import  basename, dirname, abspath
import  random
import  subprocess
import  struct

random   .seed(0)
np.random.seed(0)

# ------------------------------------------------------------------------
#                  Utilities
# ------------------------------------------------------------------------

_folder_   =  dirname(abspath(__file__))
lmap       =  lambda f,x: list(map(f,x))
lfilter    =  lambda f,x: list(filter(f,x))

def mk_pu(a,b,n,is_y=False):
    u = np.linspace(a,b,n+1)
    if is_y:
        dy  =  np.average(np.diff(u))
        for i in range(1,len(u)-1): u[i] += 0.1*dy*(1 if i%2 else -1)
    p = u[:-1] + 0.5*np.diff(u)
    if not is_y: u = u[:-1]
    return p,u

def writer(fname, A):
    assert type(A) == list
    os_system(f'mkdir -p "{dirname(fname)}"')
    with open(fname,'w') as f:
        for x in A: f.write(x.rstrip('\n') + '\n')

def reader(fname):
    with open(fname,'r') as f:
        return [x.strip('\n') for x in f]

deepdirname  =  lambda x,n: deepdirname(dirname(x),n-1) if n>0 else x

def pop1(A):
    x, = list(A)
    return x

def writer_bin(fname, A, shape):
    os_system(f'mkdir -p "{dirname(fname)}"')
    with open(fname,'wb') as fw:
        [fw.write(struct.pack('i',s)) for s in shape]
        for          k in range(shape[2]):
            for      j in range(shape[1]):
                for  i in range(shape[0]):
                    fw.write(struct.pack('d',float(A[i,j,k])))


# ------------------------------------------------------------------------
#                  Analytical Solution
# ------------------------------------------------------------------------

def manufactured_solution():
    x_adim,z_adim,y_adim =  sp.symbols('x_adim,z_adim,y_adim',real=True)
    _2pi                 =  float(2*np.pi)
    rand_sign            =  lambda                    :  2  *(random.randint(0,1) - 0.5)
    rand_num             =  lambda a=0.1              :      (random.random()*(1-a) + a)*rand_sign() # return a number in the ranges [a,1]*{1,-1}
    f                    =  lambda x,phase,all_K=[1,2]:  sum(rand_num()*sp.sin(_2pi*K*x + (rand_num()*phase)) for K in all_K)/float(len(all_K))
    U,V,W,T              =  [(1.5*f(x_adim,2)*f(z_adim,2)*f(y_adim,0)) for _ in range(4)]
    h1                   =      f(x_adim,2)*f(z_adim,2)*0.6
    h2                   =      f(x_adim,2)*f(z_adim,2)*0.6
    return {tag: str(f) for tag,f in [['U'  , U   ],
                                      ['V'  , V   ],
                                      ['W'  , W   ],
                                      ['T'  , T   ],
                                      ['h1' , h1  ],
                                      ['h2' , h2  ]]}

class CFD_Folder:
    def __init__(self, h_tag, nreps_x, nreps_z, exact_sol, trial_folder, nn, write_fields):
        self.Lx, self.Lz, self.Ly   =  [(2 + 0.5*random.random()) for _ in range(3)]
        self.nreps_x, self.nreps_z  =  nreps_x, nreps_z
        self.exact_sol              =  exact_sol
        self.h_tag                  =  h_tag
        self.pad_y                  =  0.4 + 0.1*random.random()
        self.xp, self.xu            =  mk_pu( 0         , self.Lx           , nn[0])
        self.zp, self.zu            =  mk_pu( 0         , self.Lz           , nn[1])
        self.yp, self.yu            =  mk_pu(-self.pad_y, self.Ly+self.pad_y, nn[2], is_y=True)
        self.h_min, self.h_max      =  float('inf'), float('-inf')
        self.nx, self.nz, self.ny   =  nn
        self.trial_folder           =  trial_folder
        self.case_folder            =  pjoin(self.trial_folder, self.h_tag)

        self.U                      =  self.get_term('U')
        self.W                      =  self.get_term('W')
        self.V                      =  self.get_term('V')
        self.T                      =  self.get_term('T')

        self.write_folder(write_fields)

    def get_term(self, tag, is_fluid = False):
        mm                =  dict(U='upp', W='pup', V='ppu', T='ppp')[tag]
        x,z,y             =  [getattr(self, f'{c}{m}') for c,m in zip('xzy',mm)]
        x                 =  x.reshape(-1, 1, 1)
        z                 =  z.reshape( 1,-1, 1)
        y                 =  y[:self.ny].reshape( 1, 1,-1)
        x_adim            =  x/(self.Lx/self.nreps_x)
        z_adim            =  z/(self.Lz/self.nreps_z)
        caller            =  'lambda x_adim,z_adim,y_adim,sin=np.sin:'
        h_bot             =  eval(caller+self.exact_sol[self.h_tag])(x_adim,z_adim,np.nan)
        h_top             =  self.Ly - h_bot
        self.h_min        =  min(self.h_min, h_bot.min())
        self.h_max        =  max(self.h_max, h_bot.max())
        y_adim            =  (y-h_bot)/(h_top-h_bot)
        field             =  eval(caller+self.exact_sol[tag])(x_adim,z_adim,y_adim)
        field[y_adim<=0]  =  0.
        field[y_adim>=1]  =  0.
        if is_fluid:
            field = (0.1<=y_adim)*(y_adim<=0.9)
        return field
    
    def write_folder(self, write_fields):
        
        if write_fields:
            [writer_bin(pjoin(self.case_folder, f'array_ini_{tag}.dat'), getattr(self,tag), [self.nx, self.nz, self.ny]) for tag in 'UVWT']
        writer(pjoin(self.case_folder, 'yu.txt'              ), lmap(str,self.yu.tolist()))
        writer(pjoin(self.case_folder, 'params.dat'),[f"Lx  =  {self.Lx}",
                                                      f"Lz  =  {self.Lz}",
                                                      f"Ly  =  {self.Ly}",
                                                      f"nx  =  {self.nx}",
                                                      f"nz  =  {self.nz}",    
                                                      f"ny  =  {self.ny}"])
        writer(pjoin(self.case_folder, 'rough_surf_define.py'),
'''import numpy as np
class Height_Function:
    def __init__(self):
        self.Lx = [Lx]
        self.Ly = [Ly]
        self.Lz = [Lz]
        self.nreps_x = [nreps_x]
        self.nreps_z = [nreps_z]
    def get_h(self,x,z,side,lib=np):
        assert (type(side)==str) and (side in ['bottom', 'top'])
        x_adim            =  x/(self.Lx/self.nreps_x)
        z_adim            =  z/(self.Lz/self.nreps_z)
        sin               =  lib.sin
        result            =  [expr] + 0*(x+z)
        if side == 'top': result = self.Ly - result
        return result
'''.replace('[Lx]',str(self.Lx))\
   .replace('[Ly]',str(self.Ly))\
   .replace('[Lz]',str(self.Lz))\
   .replace('[nreps_x]', str(self.nreps_x))\
   .replace('[nreps_z]', str(self.nreps_z))\
   .replace('[expr]'   , str(self.exact_sol[self.h_tag])).split('\n'))

    def check_fields(self, prev, header=False):
        new                         =  {tag: self.import_array_bin(pjoin(self.case_folder, f'array_ini_{tag}.dat'))  for tag in 'UVWT'}
        if prev is None:      prev  =  {tag: np.nan for tag in new.keys()}
        get_err =  lambda a,b,tag: np.fabs((a-b)[self.get_term(tag,is_fluid=True)]).max()
        errors  =  {tag: get_err(getattr(self,tag),new[tag],tag) for tag in 'UVWT'}
        ratios  =  {tag: prev[tag]/errors[tag] for tag in errors}
        fmt     =  lambda tag: f"{tag}: {errors[tag]:.6e} | {ratios[tag]:7.3f}"
        if header:
            fmt_s = lambda tag: f'{tag}: error_{tag}      | ratio_{tag}'
            print(f"{'nx':5s} {'nz':5s} {'ny':5s} {', '.join(lmap(fmt_s,errors.keys()))}")
        print(f"{self.nx:5d} {self.nz:5d} {self.ny:5d} {', '.join(lmap(fmt,errors.keys()))}")
        return errors
    @staticmethod
    def import_array_bin(fname):
        with open(fname, "rb") as f:
            n0,n1,n2  =  [ pop1(struct.unpack('i',f.read(4))) for _ in range(3)]
            A         =  np.zeros((n0,n1,n2), dtype=float)
            for          k in range(n2):
                for      j in range(n1):
                    for  i in range(n0):
                        A[i,j,k]  =  pop1(struct.unpack('d',f.read(8)))
        return A 

# ------------------------------------------------------------------------
#                  Main Runner
# ------------------------------------------------------------------------

if __name__ == '__main__':
    exact_sol  =  manufactured_solution()
    os_system(f"rm -rf {pjoin(_folder_, 'trials_mkini')}")
    errors       =  None
    mk_ini_file  =  'mk_ini_data_sparser_bin.py'
    files_copy   =  [mk_ini_file, 'mk_ini_data_sparser_bin.f90']
    for nn in [16, 32, 64]:
        K               =  2*nn//16
        trial_folder    =  pjoin(_folder_, 'trials_mkini', f'trials_nn_{nn}')
        nreps1, nreps2  =  [[random.randint(1,2), random.randint(1,2), 1] for _ in range(2)]
        folder1         =  CFD_Folder('h1', *nreps1[:2], exact_sol, trial_folder, [random.randint(nn-K,nn+K)*(nreps1[i]) for i in range(3)], write_fields=True )
        folder2         =  CFD_Folder('h2', *nreps2[:2], exact_sol, trial_folder, [random.randint(nn-K,nn+K)*(nreps2[i]) for i in range(3)], write_fields=False)
        [os_system(f"cp {pjoin(_folder_,xx)} {pjoin(folder2.case_folder,xx)}") for xx in files_copy]
        os_system(f"python3 {pjoin(folder2.case_folder,mk_ini_file)} {folder1.case_folder}")
        #subprocess.Popen(f"python3 {pjoin(folder2.case_folder,mk_ini_file)} {folder1.case_folder}",shell=True, stdout=subprocess.PIPE).wait()
        errors = folder2.check_fields(errors, errors is None)
    print(f"folder1  {folder1.h_min:.5e} {folder1.h_max:.5e} {folder1.case_folder}")
    print(f"folder2  {folder2.h_min:.5e} {folder2.h_max:.5e} {folder2.case_folder}")
    os_system(f"rm -rf {pjoin(_folder_, 'trials_mkini')}")