# ------------------------------------------------------------------------
#                  Generic Libraries
# ------------------------------------------------------------------------

from    os.path  import  abspath, dirname, basename
from    os.path  import  join                                       as  pjoin
from    os       import  system                                     as  os_system
import  numpy                                                       as  np
from    math     import  lcm                                        as  math_lcm
from    sys      import  path                                       as  sys_path
import  random
import  struct

SEED = 0
random.seed(SEED)
np.random.seed(SEED)

_folder_ = dirname(abspath(__file__))

# ------------------------------------------------------------------------
#                  Custom Libraries
# ------------------------------------------------------------------------

sys_path.append(pjoin(dirname(dirname(_folder_)),'src','pycode'))
from  calc_bands  import  mk_P_abc
sys_path.pop()

sys_path.append(_folder_)
from  py_poi_solver  import  Solver_Poisson_FFT_Thomas_3d_rfftx2, Solver_Poisson_FFT_Thomas_3d
sys_path.pop()

# ------------------------------------------------------------------------
#                  Basic Functions
# ------------------------------------------------------------------------

lmap       =  lambda   f,x:  list(map(f,x))
rand_val   =  lambda a,b: float(a + (b-a)*np.random.rand())
rand_sign  =  lambda    : 1 if np.random.rand()>0.5 else -1

def rand_multiple(n):
    A = []
    for i in range(1,n+1):
        if (n%i)==0:
            A.append(i)
    return A[random.randint(0,len(A)-1)]

def midpoints(a,b,n,also_yu=False):
    xu  =  np.linspace(a,b,n+1)
    xp  =  xu[:-1] + 0.5*np.diff(xu)
    if also_yu:
        return xp,xu
    else:
        return xp

def writer(fname,A):
    with open(fname,'w') as f:
        for x in A:
            f.write(str(x).strip('\n') + '\n')

def writer_bin_arr(fname, A):
    n0,n1,n2  =  A.shape[0], A.shape[1], A.shape[2]
    with open(fname,'wb') as fw:
        [fw.write(struct.pack('i',s)) for s in [n0,n1,n2]]
        for          k in range(n2):
            for      j in range(n1):
                for  i in range(n0):
                    fw.write(struct.pack('d',float(A[i,j,k])))

def writer_bin_rows(fname, A):
    conv  =  {int:'i', float:'d'}
    with open(fname,'wb') as fw:
        fw.write(struct.pack('i',len(A)))
        for row in A:
            [fw.write(struct.pack(conv[type(val)],val)) for val in row]

def ceil_div(a,b):
    c = a//b
    while a>(b*c): c += 1
    return c*b

def ceil_multiple(A,B):
    if type(B) == int: B = [B]
    B = math_lcm(*B)
    return ceil_div(A,B)

def mk_rand_ibm(shape, n, fname = None):
    used, used_g  =  set(), set()
    result        =  []
    pfix          =  lambda x: tuple([(i%s) for i,s in zip(x,shape)])
    while len(result) < n:
        ijk_g  =  tuple([random.randint(0,s-1) for s in shape])
        d_ijk  = [ 0,0,0]
        while sum(map(abs,d_ijk)) == 0:
            d_ijk = tuple([random.randint(-1,1) for _ in shape])
        ijk_1, ijk_2 = [pfix((np.array(ijk_g) + i*np.array(d_ijk)).tolist()) for i in [1,2]]
        all_ijk = ijk_g, ijk_1, ijk_2
        if ijk_g in used: continue
        if any(ijk in used_g for ijk in [ijk_1, ijk_2]): continue
        _ = [used.add(p) for p in all_ijk]
        used_g.add(ijk_g)
        result.append([*ijk_g, rand_val(-2,2), *ijk_1, rand_val(-2,2), *ijk_2, rand_val(-2,2)])
    if fname: writer_bin_rows(fname, result)
    return result

# ------------------------------------------------------------------------
#                  Direct Runner
# ------------------------------------------------------------------------

if __name__ == '__main__':
    nproc_max           =    8
    n_trials            =   10
    base_trial_folder   =  pjoin(_folder_, 'trials')
    all_mpi_divs_zy     =  [(i,j) for i in range(1,nproc_max+1) for j in range(1,nproc_max+1) if ((i*j)<=nproc_max)]
    os_system(f"rm -rf   {base_trial_folder}")
    os_system(f"mkdir -p {base_trial_folder}")
    for itrial in range(n_trials):
        mpi_divs_z, mpi_divs_y  =  all_mpi_divs_zy[random.randint(0,len(all_mpi_divs_zy)-1)]
        nproc                   =  mpi_divs_z*mpi_divs_y
        vmax_B  =  float('inf')
        while np.isnan(vmax_B) or (vmax_B>1e3):
            nx, nz, ny   =  shape  =  [random.randint(max(mpi_divs_y,mpi_divs_z),max(10,5*max(mpi_divs_y,mpi_divs_z))) for _ in range(2)] + [random.randint(2*mpi_divs_y,max(10,5*max(mpi_divs_y,mpi_divs_z)))]
            trial_folder           =  pjoin(base_trial_folder, f'trial_{itrial:04d}')
            dx, dz                 =  [float(rand_val(0.1,0.2)) for _ in range(2)]
            A_input                =  5*(np.random.rand(*shape)-0.5)
            yu                     =  np.array([0] + np.cumsum([rand_val(0.1,0.2) for _ in range(ny)]).tolist()) + rand_val(-0.2,0.2)
            assert len(yu) == (ny+1)
            B_output  =  Solver_Poisson_FFT_Thomas_3d_rfftx2(shape, dx+0, dz+0, yu.copy())(A_input.copy())
            vmax_B    =  np.fabs(B_output).max()
        vmax_check = np.fabs(B_output - Solver_Poisson_FFT_Thomas_3d(shape, dx+0, dz+0, yu.copy())(A_input.copy())).max()
        assert vmax_check < 1e-10
        os_system(f"mkdir -p {pjoin(trial_folder, 'geom_data')}")
        assert (not np.isnan(B_output).any())
        yp             =  yu[:-1] + 0.5*np.diff(yu)
        yp_ext         =  lmap(float,[yu[0]] + yp.tolist() + [yu[-1]])
        P_bands_abc    =  [' '.join(map(str,row)) for row in mk_P_abc(yp_ext,is_dirich_first=False)]
        P_bands_abc_00 =  [' '.join(map(str,row)) for row in mk_P_abc(yp_ext,is_dirich_first=True )]
        n_ibm          =  random.randint(1,np.prod(shape)//5)
        nhalo, nskip   =  [random.randint(0,1) for _ in range(2)]
        A              =  [str(nproc), ' '.join(map(str,[nx, nz, ny, mpi_divs_z, mpi_divs_y, dx, dz, n_ibm, nhalo, nskip]))]
        A.extend(P_bands_abc)
        A.extend(P_bands_abc_00)
        assert A_input.shape == B_output.shape == tuple(shape)
        writer(pjoin(trial_folder, 'info.dat' ), lmap(str,A))
        writer_bin_arr(pjoin(trial_folder, 'A_input.dat' ), A_input )
        writer_bin_arr(pjoin(trial_folder, 'B_output.dat'), B_output)
        assert n_ibm >= 1
        ibm_coeffs  =  mk_rand_ibm(shape, n_ibm, fname = pjoin(trial_folder, 'geom_data', 'ibm_coeffs_A.dat'))
        A_ibm       =  A_input.copy()
        for ig,jg,kg, c1, i1,j1,k1, c2, i2,j2,k2, cw in ibm_coeffs:
            A_ibm[ig,jg,kg] = c1*A_input[i1,j1,k1] + c2*A_input[i2,j2,k2] + cw
        writer_bin_arr(pjoin(trial_folder, 'A_ibm.dat' ), A_ibm )

