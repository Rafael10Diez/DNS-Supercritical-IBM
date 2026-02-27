
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
#                  Bulk calculations (Fortran)
# ------------------------------------------------------------------------

def get_vars(A):
    assert not (' ' in A)
    result = []
    valid  = set(list(ascii_lowercase + ascii_uppercase + '_')+lmap(str,range(10)))
    for j,c in enumerate(A):
        if c=='(':
            i = j-1
            while (i>=0) and (A[i] in valid): i -= 1
            i += 1
            if (j-i)>0:
                result.append(A[i:j])
    return sorted(set(result))

def get_ijk(A):
    assert not (' ' in A)
    result    = []
    start,end = [[i for i,c in enumerate(A) if c==p] for p in '()']
    for k,c in enumerate(A):
        if c==',':
            i = [i for i in start if i<k  ][-1]
            j = [j for j in end   if   k<j][0]
            result.append(A[i:j+1])
    return sorted(set(result))

def get_tag_type(tag,runtime_folder):
    if tag[:4] in ['fav_', 'fRf_', 'fRm_', 'fdx_', 'fdz_', 'fdy_', 'pTK_']: return 'T'
    A       =  lfilter(None,[x.split('!')[0].split('#')[0].strip(' ') 
                             for x in reader(pjoin(runtime_folder,'metaprogramming','avg_variables.txt'))])
    A       =  {key.strip(' '):val.strip(' ') for key,val in lmap(lambda x: x.split('=>'), A)}
    expr    =  A[f"arr_3d_{tag}_avg"].replace(' ','')
    assert not any(('0' in x) for x in get_ijk(expr))
    if (get_ijk(expr)==['(i,j,k)']) and (get_vars(expr) in [[c] for c in 'UVW']): return get_vars(expr)[0]
    return 'T'

def main_runner_prof_1d():
    runtime_folder             =  dirname(_folder_)
    opts                       =  '-O3 -cpp ' # -g -fbounds-check
    if len(argv)>=2: opts += ' '.join(argv[1:])
    fname_balance_avg_1d_prof  =  pjoin(_folder_, 'balance_avg_1d_prof.f90')
    print(f'INFO: Compilation flags: [{opts}] (only_prof1d.py)')
    assert isfile(fname_balance_avg_1d_prof)

    all_iter_avg  =  sorted(set([int(basename(x).split('iter_')[1].split('.')[0]) for x in listdir_full_files(pjoin(runtime_folder,'avg')) if x.endswith('.dat')]))
    full_tags     =  ','.join(sorted(set(basename(x).split('array_arr_3d_')[1].split('_avg_iter_')[0] for x in listdir_full_files(pjoin(runtime_folder,'avg')))))
    A_txt = []
    n_txt = 0
    for iter_avg in all_iter_avg:
        all_tags = 'U,V,W,T,P'
        if '-D_READ_RMS' in opts: all_tags = full_tags
        for tag in all_tags.split(','):
            n_txt += 1
            A_txt.append(f"{get_tag_type(tag,runtime_folder)} {tag+'_'*(8-len(tag))}")
            A_txt.append(f"../avg/array_arr_3d_{tag}_avg_iter_{iter_avg:09d}.dat")
            A_txt.append(f"../prof_1d_from_avg3d/prof_1d_avg_{tag}_iter_{iter_avg:09d}.dat") 
    A_txt = [str(n_txt), *A_txt]
    writer(pjoin(_folder_,'all_interp_avg.dat'), A_txt)
    
    def run_prof_1d():
        t0 = time()
        print('----> Begin: Running Bulk calculations (Fortran)')
        os_system(f"mkdir -p {pjoin(runtime_folder,'prof_1d_from_avg3d')}")
        os_system(f'cd "{_folder_}"; gfortran {opts} {fname_balance_avg_1d_prof} -o a.out;./a.out')
        print(f'----> End: Running Bulk calculations (Fortran) (elapsed_time: {format_dt(time()-t0)})')
    run_prof_1d()

    os_system(f"rm -f {pjoin(_folder_,'all_interp_avg.dat')}")
    os_system(f"cd {_folder_} ; rm -f *.mod")
    os_system(f"cd {_folder_} ; rm -f a.out")
    
# ------------------------------------------------------------------------
#                  Direct Runner
# ------------------------------------------------------------------------

if __name__ == "__main__":
    t0 = time()
    print('----> Begin: main_runner_prof_1d()')
    main_runner_prof_1d()
    print(f'----> End: main_runner_prof_1d() (elapsed_time: {format_dt(time()-t0)})')


