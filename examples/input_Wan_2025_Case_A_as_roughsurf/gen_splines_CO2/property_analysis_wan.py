# ------------------------------------------------------------------------
#                  Generic Libraries
# ------------------------------------------------------------------------

import  numpy                                                                 as  np
import  matplotlib.pyplot                                                     as  plt
from    os.path            import  join                                       as  pjoin
from    sys                import  path                                       as  sys_path
from    os                 import  system                                     as  os_system
from    os                 import  listdir                                    as  os_listdir
from    os.path            import  dirname, basename, isfile, isdir, abspath
from    copy               import  deepcopy
from    scipy.interpolate  import  interp1d

# ------------------------------------------------------------------------
#                  Custom Library
# ------------------------------------------------------------------------

_folder_  =  dirname(abspath(__file__))

sys_path.append(_folder_)
from scan_data_wan import imported_functions_wan
sys_path.pop()

# ------------------------------------------------------------------------
#                  Basic Definition
# ------------------------------------------------------------------------

tolist                =  lambda   x:  x.tolist() if hasattr(x,'tolist') else x
lmap                  =  lambda f,x:  list(map(f,x))
lfilter               =  lambda f,x:  list(filter(f,x))

listdir_full          =  lambda x: [pjoin(x,y) for y in os_listdir(x)]
listdir_full_files    =  lambda x: lfilter(isfile, listdir_full(x))
listdir_full_folders  =  lambda x: lfilter(isdir , listdir_full(x))

logspace              =  lambda a,b,n: 10**np.linspace(np.log10(a), np.log10(b), n)

def reader(fname):
    with open(fname,'r') as f:
        return [x.strip('\n') for x in f]

def writer(fname,A):
    if not isdir(dirname(fname)): os_system(f'mkdir -p "{dirname(fname)}"')
    with open(fname,'w') as f:
        for x in A:
            f.write(x+'\n')

def ftn_line_split(A,lim=100,last_cont=False):
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

# ------------------------------------------------------------------------
#                 Make Fortran Code
# ------------------------------------------------------------------------

def mk_fortran_code(f,arr_name):
    import sympy as sp
    to_dbl = lambda s: str(s).replace('E+','D+').replace('E-','D-').replace('e+','D+').replace('e-','D-')
    def decision_tree(spl_c, spl_x, arr_name, lvl):
        pad                =  '  '*lvl
        if spl_c.shape[1] == 1:
            assert len(spl_x) == 2
            xx          =  f"(T_ext_ijk - {float(spl_x[0]):.16e})"
            cc          =  [float(spl_c[i,0]) for i in range(4)]
            return [f"{pad}{arr_name} = ({cc[3]:.16e}) + ({cc[2]:.16e})*({xx}) + ({cc[1]:.16e})*(({xx})**2) + ({cc[0]:.16e})*(({xx})**3)"]
        else:
            mid = len(spl_x)//2
            return [f'{pad}if (T_ext_ijk<=({(spl_x[mid]):.16e})) then ',
                    *decision_tree(spl_c[:,:mid], spl_x[:mid+1], arr_name, lvl+1),
                    f'{pad}else'  ,
                    *decision_tree(spl_c[:,mid:], spl_x[mid:], arr_name, lvl+1),
                    f'{pad}end if']
    return  lmap(to_dbl, [f'! {arr_name}'] + [ftn_line_split(x) for x in decision_tree(f.c, f.x, arr_name, 0)])

# ------------------------------------------------------------------------
#                  Fortran with Python Indentation
# ------------------------------------------------------------------------

def fort2py_indenter(A0):
    A    =  deepcopy(A0)
    lvl  =  0
    for i in range(len(A)):
        A[i]     =  A[i].strip(' ')
        is_if    =  A[i][:2] == 'if'
        is_else  =  A[i][:4] == 'else'
        is_end   =  A[i][:3] == 'end'
        assert (is_if + is_else + is_end) <= 1
        if is_else: lvl -= 1
        if is_end : lvl -= 1
        A[i]  = (' '*(lvl*4)) + A[i]
        if is_if  : lvl += 1
        if is_else: lvl += 1
    return A

# ------------------------------------------------------------------------
#                 Test Fortran Code
# ------------------------------------------------------------------------

def mk_verif_plots(A0, T0_dns, T1_dns, curves_nist, curves_svg, do_plot = True):
    curves_nist, curves_svg = lmap(deepcopy, [curves_nist, curves_svg])
    A = deepcopy(A0)
    def conv_to_python(A):
        x = '\n'.join(fort2py_indenter(A))
        x = x.replace('_dp'    , ''                )\
             .replace('then'   , ':'               )\
             .replace('else'   , 'else:'           )\
             .replace('end if' , ''                )\
             .replace('D+'     , 'e+'              )\
             .replace('&'      , r" \ ".strip(' ') )\
             .replace('D-'     , 'e-'              )\
             .replace('(i,j,k)', '[i,j,k]'         )\
             .replace('!'      , '#'               )
        x = x.split('\n')
        return [f"    {val}" for val in x if val.strip(' ')]
    success = False
    for ii in range(10):
        fname = pjoin(_folder_,'output',f'kernel_ijk_test_{ii}.py')
        if not isfile(fname):
            writer(fname,[f'def kernel_ijk_{ii}(i,j,k,T_ext,rho_arr,mu_arr,cond_arr,cp_arr,div_jacobian_arr,enth_arr,buoyancy_arr):',
                          *conv_to_python(A)])
            sys_path.append(pjoin(_folder_,'output'))
            if ii == 0:
                from kernel_ijk_test_0 import kernel_ijk_0
                my_kernel                  =  kernel_ijk_0
            elif ii==1:
                from kernel_ijk_test_1 import kernel_ijk_1
                my_kernel                  =  kernel_ijk_1
            elif ii==2:
                from kernel_ijk_test_2 import kernel_ijk_2
                my_kernel                  =  kernel_ijk_2
            elif ii==3:
                from kernel_ijk_test_3 import kernel_ijk_3
                my_kernel                  =  kernel_ijk_3
            elif ii==4:
                from kernel_ijk_test_4 import kernel_ijk_4
                my_kernel                  =  kernel_ijk_4
            if ii == 5:
                from kernel_ijk_test_5 import kernel_ijk_5
                my_kernel                  =  kernel_ijk_5
            elif ii==6:
                from kernel_ijk_test_6 import kernel_ijk_6
                my_kernel                  =  kernel_ijk_6
            elif ii==7:
                from kernel_ijk_test_7 import kernel_ijk_7
                my_kernel                  =  kernel_ijk_7
            elif ii==8:
                from kernel_ijk_test_8 import kernel_ijk_8
                my_kernel                  =  kernel_ijk_8
            elif ii==9:
                from kernel_ijk_test_9 import kernel_ijk_9
                my_kernel                  =  kernel_ijk_9
            sys_path.pop()
            success = True
            break
    assert success
    n                 =  10000
    ni,nj,nk          =  n,1,1
    empty             =  lambda : np.zeros((ni,nj,nk), dtype=float)*np.nan
    enth_arr          =  empty()
    rho_arr           =  empty()
    mu_arr            =  empty()
    cond_arr          =  empty()
    cp_arr            =  empty()
    div_jacobian_arr  =  empty()
    buoyancy_arr      =  empty()
    T_ext             =  np.linspace(T0_dns, T1_dns, n).reshape(*empty().shape)
    _ = [my_kernel(i,j,k,T_ext,rho_arr,mu_arr,cond_arr,cp_arr,div_jacobian_arr,enth_arr,buoyancy_arr) 
         for i in range(ni) 
         for j in range(nj)
         for k in range(nk)]
    fields  =  dict( enth          =  enth_arr        .reshape(-1),
                     rho           =  rho_arr         .reshape(-1),
                     mu            =  mu_arr          .reshape(-1),
                     cond          =  cond_arr        .reshape(-1),
                     cp            =  cp_arr          .reshape(-1),
                     div_jacobian  =  div_jacobian_arr.reshape(-1),
                     T_ext         =  T_ext           .reshape(-1))
    if do_plot:
        def add_curves_data(svg_also=False):
            if svg_also:
                T_svg, rho_svg         =  curves_svg['density'].transpose()
                curves_svg['drho_dT']  =  np.array([T_svg[:-1]+0.5*np.diff(T_svg),
                                                    np.diff(rho_svg)/np.diff(T_svg)]).transpose()
            assert 'drho_dT' in curves_nist
            reinterp1d  =  lambda x,y,x0: interp1d(x,y,kind='cubic', fill_value="extrapolate")(x0)
            for C in ([curves_nist,curves_svg] if svg_also else [curves_nist]):
                T, rho        =  C['density'].transpose()
                C['rho']      =  np.array([T,rho]).transpose()
                cp            =  reinterp1d(*C['cp'].transpose()     , T)
                drho_dT       =  reinterp1d(*C['drho_dT'].transpose(), T)
                C['div_jacobian'] = np.array([T,-1/((rho**2)*cp)*drho_dT]).transpose()
        add_curves_data()
        for key in fields.keys():
            if key != 'T_ext':
                fig    =  plt.figure()
                plt.plot(fields['T_ext'], fields[key],'b', label = 'kernel')
                plt.plot(*list(curves_nist[key].transpose()),'r:',label = 'curves_nist')
                #if key !='div_jacobian':
                #    plt.plot(*list(curves_svg[key].transpose()),'g:',label = 'curves_svg')
                title = f'im_kernel_{key}_arr'
                plt.title(title)
                plt.legend()
                plt.savefig(pjoin(_folder_,'output',f'{title}.png')) ; plt.close(fig)
                print('Verif', f"{key:20s}", fields['T_ext'][0], fields[key][0])
    print('\n\nTests passed!')
    return fields

# ------------------------------------------------------------------------
#                  Direct Runner
# ------------------------------------------------------------------------

if __name__ == '__main__':
    try:
        has_cleaned
    except:
        has_cleaned = False
    if not has_cleaned:
        for fname in lfilter(lambda x: (basename(x)[:16]=='kernel_ijk_test_') and (basename(x)[-3:]=='.py'), listdir_full_files(_folder_)):
            print(fname)
            os_system(f"rm -f {fname}")
        has_cleaned = True

    all_f               =  imported_functions_wan()
    A_declare, A_write  =  ['T_ext_ijk'], ['T_ext_ijk = T_ext(i,j,k)']

    for key in all_f.keys():
        if 'scipy.interpolate.' in str(type(all_f[key])):
            A_declare.append(f'{key}_ijk')
            A_write.extend(mk_fortran_code(all_f[key], A_declare[-1]))

    A_write = [*A_write,
                "enth_arr(i,j,k)         = enth_ijk",
                "rho_arr(i,j,k)          = rho_ijk",
                "mu_arr(i,j,k)           = mu_ijk",
                "cond_arr(i,j,k)         = cond_ijk",
                "cp_arr(i,j,k)           = cp_ijk",
                "div_jacobian_arr(i,j,k) = -1_dp/((rho_ijk**2)*cp_ijk)*drho_dT_ijk"]
    
    A_write.append('''buoyancy_arr(i,j,k)     = -146.0390181723833_dp*(rho_ijk-511.8479411035_dp) ! Gb_y is intended to be "Richardson number"

! rho_cold =  780.731303553 # 297.8 
! rho_hot  =  242.964578654 # 317.8 
! rho_pb   =  462.727153614
! h        =  0.001
! U_bulk   =  0.4119731598987495
! Ri       =  1
! g        =  Ri/(((rho_hot-rho_cold)/rho_pb)*h/(U_bulk**2))
! assert abs((((rho_hot-rho_cold)/rho_pb)*g*h/(U_bulk**2)) - Ri) < 1e-10''')
    writer(pjoin(_folder_,'recomputed', 'expressions_rho_mu_cond_buoyancy_cp_div_jacobian.f90'), A_write)

    assert A_declare == ['T_ext_ijk', 'mu_ijk', 'cond_ijk', 'rho_ijk', 'cp_ijk', 'drho_dT_ijk', 'enth_ijk'],A_declare
    fields_now  = mk_verif_plots(A_write, all_f['T0_dns'], all_f['T1_dns'], all_f['curves_nist'], all_f['curves_svg'])
    fields_input = mk_verif_plots(reader(pjoin(dirname(_folder_),'expressions_rho_mu_cond_buoyancy_cp_div_jacobian.f90')), all_f['T0_dns'], all_f['T1_dns'], all_f['curves_nist'], all_f['curves_svg'])

    max_diff = {key: np.fabs(fields_now[key]-fields_input[key]).max() for key in fields_input.keys()}
    max_diff = max(max_diff.values())
    assert max_diff < 1e-9, max_diff
    print('Properties validated!')