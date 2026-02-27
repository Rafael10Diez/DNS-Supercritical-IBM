# ------------------------------------------------------------------------
#                  Generic Libraries
# ------------------------------------------------------------------------

from    os.path  import  abspath, basename, dirname, isfile, isdir
from    os.path  import  join                                       as  pjoin
from    os       import  listdir                                    as  os_listdir
from    sys      import  argv

# ------------------------------------------------------------------------
#                  Basic Functions
# ------------------------------------------------------------------------

_folder_             = dirname(abspath(__file__))
lmap                 = lambda f,x: list(map(f,x))
lfilter              = lambda f,x: list(filter(f,x))
listdir_full         = lambda x: sorted([pjoin(x,y) for y in os_listdir(x)])
listdir_full_files   = lambda x: lfilter(isfile, listdir_full(x))
listdir_full_folders = lambda x: lfilter(isdir , listdir_full(x))

def listdir_full_files_deep(x):
    def dfs(x):
        for y in listdir_full_files(x):
            yield y
        for y in listdir_full_folders(x):
           yield from dfs(y)
    return list(dfs(x))

def reader(fname):
  with open(fname,'r') as f:
    return [x.strip('\n') for x in f]

def writer(fname, A):
  assert type(A) == list
  with open(fname,'w') as f:
    for x in A:
        f.write(str(x).strip('\n')+'\n')

import_json  =  lambda fname:  eval(' '.join([x.split('#')[0] for x in reader(fname)]))
lmap         =  lambda   f,x:  list(map(f,x))
lfilter      =  lambda   f,x:  list(filter(f,x))

def dfs_includer(A,root_folder):
    get_i = lambda A: [i for i,x in enumerate(A) if ('include'==x.strip(' ')[:7] and (not 'fftw3.f03' in x))]
    while get_i(A):
        i     = get_i(A)[0]
        fname = A[i].split('include')[1].strip(' "'+"'")
        B     = reader(pjoin(root_folder, fname))
        if dirname(fname):
           folder_B  = pjoin(root_folder,dirname(fname))
        else:
           folder_B  = root_folder
        A     = [*A[:i],
                 f'! ---------------- Begin include: {fname} ---------------- ',
                 '!'+A[i],
                 *dfs_includer(B,folder_B),
                 f'! ---------------- End include: {fname} ---------------- ',
                 *A[i+1:]]
    return A

# ------------------------------------------------------------------------
#                  Meta-programming (schedule variables to average)
# ------------------------------------------------------------------------

def get_pairs_var_expr(fname):
    A = lfilter(None, [x.strip(' \n') for x in reader(fname)])
    return {a.strip(' '):b.strip(' ') for x in A for a,b in [x.split('=>')]}

def mk_prof_1d_commands(all_files):
    pairs  =  get_pairs_var_expr(pjoin(dirname(_folder_),'metaprogramming','prof_1d_variables.txt'))
    holders = {}
    holders['!__[PYTHON_HOLDER_PROF1D_DECLARE]'] = '\n'.join(f" real(dp), allocatable :: {x}(:)" for x,y in pairs.items())

    holders['!__[PYTHON_HOLDER_PROF1D_SET]'] = '\n'.join('\n'.join([f" allocate({x}(0:ny_loc-1))",
                                                                    f" {x} = 0",
                                                                    f"!$acc enter data copyin({x})",
                                                                    f"!$acc wait"]) for x,y in pairs.items())
    
    holders['!__[PYTHON_HOLDER_PROF1D_APPLY_0]'] = '\n'.join('\n'.join([f"block",
                                                                        f"  integer :: i,j,k",
                                                                        f"  !$acc parallel loop  collapse(3) default(present)",
                                                                        f"  do     k = 0, ny_loc-1",
                                                                        f"    do   j = 0, nz_loc-1",
                                                                        f"      do i = 0, nx_glob-1 ",
                                                                        f"        temp_arr_n_xzy(i,j,k) = ({y})",
                                                                        f"      end do",
                                                                        f"    end do",
                                                                        f"  end do",
                                                                        f"!$acc wait",
                                                                        f"  call batched_binary_reduction(temp_arr_n_xzy, prof_1d_temp, nx_glob*nz_loc, ny_loc, .true., nx_glob*nz_glob,&",
                                                                        f"                                tr_red1d_fwd,tr_red1d_bwd,mpi_divs_z,ny_loc_00k,tr_red1d_B,tr_red1d_B_1d,work_tr)",
                                                                        f"  !$acc parallel loop  collapse(1) default(present)",
                                                                        f"  do k = 0, ny_loc-1",
                                                                        f"    {x}(k) = {x}(k) + prof_1d_temp(k)",
                                                                        f"  end do",
                                                                        f"end block"]) for x,y in pairs.items())

    holders['!__[PYTHON_HOLDER_PROF1D_COMPBULK_0]']               =  '\n'.join(f" {x}(k) = {x}(k)*temp_real" for x,y in pairs.items())
    holders['!__[PYTHON_HOLDER_PROF1D_COMPBULK_VARS]']            =  '\n'.join(f" {x}, &" for x,y in pairs.items())
    holders['!__[PYTHON_HOLDER_PROF1D_COMPBULK_FUNC_DECL_VARS]']  =  '\n'.join(f"real(dp) :: {x}(0:)" for x,y in pairs.items())
    holders['!__[PYTHON_HOLDER_PROF1D_COMPBULK_1]']               =  '\n'.join('\n'.join([f"  block",
                                                                                          f"    character(len={(len(x)+37)}) :: my_file",
                                                                                          f"    write(my_file,'(A{(len(x)+24)},I0.9,A4)') './prof_1d/prof_1d_{x}_iter_',iters_full,'.dat'",
                                                                                          f"    call write_prof_1d({x}, &",
                                                                                          f"                       my_file, mpi_divs_y, mpi_pos_z, mpi_pos_y, obj_ranks_0zy, prof_1d_temp, ny_loc)",
                                                                                          f"  end block"]) for x,y in pairs.items())
    holders['!__[PYTHON_HOLDER_PROF1D_COMPBULK_2]'] = '\n'.join(f" {x}(k) = 0." for x,y in pairs.items())
    for fname in all_files:
        A_main = '\n'.join(reader(fname))
        for key,s in holders.items():
            # assert (type(key) == type(s) == str) and (A_main.count(key) >= 1), [A_main.count(key),key,fname,s,A_main]
            A_main = A_main.replace(key,s)
        writer(fname, [A_main])

def mk_avg_commands(all_files):
    pairs   =  get_pairs_var_expr(pjoin(dirname(_folder_),'metaprogramming','avg_variables.txt'))
    holders =  {}
    holders['!__[PYTHON_HOLDER_AVG_DECLARE]'] = '\n'.join([f" real(dp), allocatable :: {x}(:,:,:)"          for x,y in pairs.items()]+
                                                          [f" real(dp), allocatable :: binary_{x}(:,:,:,:)" for x,y in pairs.items()])
    
    holders['!__[PYTHON_HOLDER_AVG_SET]'] = '\n'.join('\n'.join([f" allocate({x}(0:nx_glob-1, 0:nz_loc-1, 0:ny_loc-1))",
                                                                 f" {x} = 0",
                                                                 f"!$acc enter data copyin({x})",
                                                                 f"!$acc wait",
                                                                 f"if (avg_bin_nbins>0) then",
                                                                 f"  allocate(binary_{x}(0:nx_glob-1, 0:nz_loc-1, 0:ny_loc-1, 0:avg_bin_nbins-1))",
                                                                 f"  binary_{x} = 0",
                                                                  "end if"]) for x,y in pairs.items())
    holders['!__[PYTHON_HOLDER_AVG_APPLY_0]'] = '\n'.join(f"{x}(i,j,k) = {x}(i,j,k) + ({y})*avg_inv_nraw" for x,y in pairs.items())
    holders['!__[PYTHON_HOLDER_AVG_APPLY_1]'] = '\n'.join('\n'.join([f"!$acc update host({x})",
                                                                     f"!$acc wait",
                                                                     f" block",
                                                                     f"   character(len={(len(x)+31)})  ::  my_file",
                                                                     f"     associate( A => {x})",
                                                                     f"       write(my_file,'(A{(len(x)+18)},I0.9,A4)') './avg/array_{x}_iter_',iters_full,'.dat'",
                                                                     f"       call cfd_write_arr(my_file, {x}, 0, tr_io_fwd, slab_00k, temp_arr_n_xzy, work_tr, &",
                                                                     f"                          irank_mpi, nproc_mpi, ny_loc, nz_loc, ny_loc_00k, nx_glob, nz_glob, ny_glob)",
                                                                      "     end associate",
                                                                      " end block"]) for x,y in pairs.items())
    holders['!__[PYTHON_HOLDER_AVG_APPLY_2]']         =  '\n'.join(f" {x}(i,j,k) = 0." for x,y in pairs.items())
    holders['!__[PYTHON_HOLDER_BINAVG_0]']            =  '\n'.join(f"!$acc update host({x})\nbinary_{x}(:,:,:, ibin)  =  {x}" for x,y in pairs.items())
    holders['!__[PYTHON_HOLDER_BINAVG_1]']            =  '\n'.join(f"!$acc update host({x})\n{x} =  0.5*({x} + binary_{x}(:,:,:, ibin))\n!$acc update device({x})" for x,y in pairs.items())
    holders['!__[PYTHON_HOLDER_AVG_DECL_FUNC_VARS]']  =  '\n'.join([f" real(dp) :: {x}(0:,0:,0:)"          for x,y in pairs.items()]+
                                                                   [f" real(dp) :: binary_{x}(0:,0:,0:,0:)" for x,y in pairs.items()])
    holders['!__[PYTHON_HOLDER_AVG_VARS]']            =  '\n'.join([f" {x} , &"        for x,y in pairs.items()]+
                                                                   [f" binary_{x}, &"  for x,y in pairs.items()])
    for fname in all_files:
        A_main = '\n'.join(reader(fname))
        for key,s in holders.items():
            # assert (type(key) == type(s) == str) and (A_main.count(key) >= 1), [A_main.count(key),key,fname,s,A_main]
            A_main = A_main.replace(key,s)
        writer(fname, [A_main])

# ------------------------------------------------------------------------
#                  Direct Runner
# ------------------------------------------------------------------------

if __name__ == '__main__':
    root_folder  =  dirname(_folder_)
    all_files    =  lfilter(lambda x: basename(x).endswith('.f90'), (listdir_full_files(root_folder                                    ) + 
                                                                     listdir_full_files(pjoin(root_folder, 'ibm_setup'                )) + 
                                                                     listdir_full_files(pjoin(root_folder, 'poisson_solver_mgpu_dtdma')) ))
    for fname in all_files:
        A            =  reader(fname)
        writer(pjoin(dirname(fname), f'old_{basename(fname)}'), A)
        writer(fname, dfs_includer(A,root_folder))
    if not (len(argv) >= 2 and argv[1]=='avoid_metaprogramming'):
        mk_avg_commands(all_files)
        mk_prof_1d_commands(all_files)
    