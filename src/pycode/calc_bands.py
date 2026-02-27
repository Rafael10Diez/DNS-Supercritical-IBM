# ------------------------------------------------------------------------
#                  Generic Libraries
# ------------------------------------------------------------------------

from    os.path  import  abspath, dirname, basename, isfile, isdir
from    os.path  import  join                                       as  pjoin
from    os       import  mkdir                                      as  os_mkdir
from    os       import  listdir                                    as  os_listdir
from    os       import  system                                     as  os_system
from    sys      import  path                                       as  sys_path
from    copy     import  deepcopy
from    sys      import  argv
import  struct
import  subprocess

# ------------------------------------------------------------------------
#                  Basic Functions
# ------------------------------------------------------------------------

_folder_ = dirname(abspath(__file__))

def reader(fname):
  with open(fname,'r') as f:
    return [x.strip('\n') for x in f]

def writer(fname, A):
  os_system(f"mkdir -p {dirname(fname)}")
  assert type(A) == list
  with open(fname,'w') as f:
    for x in A:
        f.write(str(x).strip('\n')+'\n')

linspace  =  lambda a,b,n: [(a+(b-a)*float(i/(n-1.))) for i in range(n)]

def mk_pu(L,n):
    u = linspace(0,L,n+1)
    p = [(u[i]+0.5*(u[i+1]-u[i])) for i in range(len(u)-1)]
    return p, u[:-1]

tolist       =  lambda     x:  x.tolist() if hasattr(x,'tolist') else x
lmap         =  lambda   f,x:  list(map(f,x))
lfilter      =  lambda   f,x:  list(filter(f,x))

def import_fort_dat(fname):
    result = {}
    for line in lfilter(None, [x.split('#')[0].split('!')[0].strip(' ') for x in reader(fname)]):
        assert line.count('=') == 1
        key,val = [x.strip(' ') for x in line.split('=')]
        if val in ['.true.','.false.']: val = val.strip('.').title()
        result[key] = eval(val)
    return result

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
#                  Utilities
# ------------------------------------------------------------------------

fd_coeffs_ref  =  ' '.join("""lambda xs,xc,xn: {('ddx'  , 'normal'  ): [(-xc**2 + xn**2)/(xc**2*xn - xc**2*xs - xc*xn**2 + xc*xs**2 + xn**2*xs - xn*xs**2), (-xn**2 + xs**2)/(xc**2*xn - xc**2*xs - xc*xn**2 + xc*xs**2 + xn**2*xs - xn*xs**2), (xc**2 - xs**2)/(xc**2*xn - xc**2*xs - xc*xn**2 + xc*xs**2 + xn**2*xs - xn*xs**2)],
                                                ('d2dx2', 'normal'  ): [2*(xc - xn)/(xc**2*xn - xc**2*xs - xc*xn**2 + xc*xs**2 + xn**2*xs - xn*xs**2), 2*(xn - xs)/(xc**2*xn - xc**2*xs - xc*xn**2 + xc*xs**2 + xn**2*xs - xn*xs**2), 2*(-xc + xs)/(xc**2*xn - xc**2*xs - xc*xn**2 + xc*xs**2 + xn**2*xs - xn*xs**2)],
                                                ('ddx'  , 'dirich_s'): [(-xc**2 + xn**2)/(xc**3 - 3*xc**2*xs - xc*xn**2 + 2*xc*xn*xs + 2*xc*xs**2 + xn**2*xs - 2*xn*xs**2), (xc**2 - 2*xc*xs - xn**2 + 2*xs**2)/(xc**3 - 3*xc**2*xs - xc*xn**2 + 2*xc*xn*xs + 2*xc*xs**2 + xn**2*xs - 2*xn*xs**2), (2*xc*xs - 2*xs**2)/(xc**3 - 3*xc**2*xs - xc*xn**2 + 2*xc*xn*xs + 2*xc*xs**2 + xn**2*xs - 2*xn*xs**2)],
                                                ('d2dx2', 'dirich_s'): [2*(xc - xn)/(xc**3 - 3*xc**2*xs - xc*xn**2 + 2*xc*xn*xs + 2*xc*xs**2 + xn**2*xs - 2*xn*xs**2), 2*(xn - xs)/(xc**3 - 3*xc**2*xs - xc*xn**2 + 2*xc*xn*xs + 2*xc*xs**2 + xn**2*xs - 2*xn*xs**2), 2*(-xc + xs)/(xc**3 - 3*xc**2*xs - xc*xn**2 + 2*xc*xn*xs + 2*xc*xs**2 + xn**2*xs - 2*xn*xs**2)],
                                                ('ddx'  , 'dirich_n'): [(2*xc*xn - 2*xn**2)/(xc**3 - 3*xc**2*xn + 2*xc*xn**2 + 2*xc*xn*xs - xc*xs**2 - 2*xn**2*xs + xn*xs**2), (xc**2 - 2*xc*xn + 2*xn**2 - xs**2)/(xc**3 - 3*xc**2*xn + 2*xc*xn**2 + 2*xc*xn*xs - xc*xs**2 - 2*xn**2*xs + xn*xs**2), (-xc**2 + xs**2)/(xc**3 - 3*xc**2*xn + 2*xc*xn**2 + 2*xc*xn*xs - xc*xs**2 - 2*xn**2*xs + xn*xs**2)],
                                                ('d2dx2', 'dirich_n'): [2*(-xc + xn)/(xc**3 - 3*xc**2*xn + 2*xc*xn**2 + 2*xc*xn*xs - xc*xs**2 - 2*xn**2*xs + xn*xs**2), 2*(-xn + xs)/(xc**3 - 3*xc**2*xn + 2*xc*xn**2 + 2*xc*xn*xs - xc*xs**2 - 2*xn**2*xs + xn*xs**2), 2*(xc - xs)/(xc**3 - 3*xc**2*xn + 2*xc*xn**2 + 2*xc*xn*xs - xc*xs**2 - 2*xn**2*xs + xn*xs**2)]}""".split('\n'))

def gen_fd_coeffs():
    import  sympy  as  sp
    result                  =  {}
    x,xs,xc,xn,a,b,c,fs,fc,fn  =  sp.symbols('x,xs,xc,xn,a,b,c,fs,fc,fn',real=True)
    for bc_type in ['normal','dirich_s','dirich_n']:
        f    =  a + b*x + c*(x**2)
        eqs  =  [f.subs(x,xs)- fs,
                 f.subs(x,xc)- fc,
                 f.subs(x,xn)- fn]
        #if bc_type == 'neu_s'  : eqs[ 0] = f.diff(x).subs(x,xs)
        #if bc_type == 'neu_n'  : eqs[-1] = f.diff(x).subs(x,xn)
        if bc_type == 'dirich_s': eqs[ 0] = f.subs(x,2*xs-xc) - (2*fs-fc) #  f.subs(x,xs)
        if bc_type == 'dirich_n': eqs[-1] = f.subs(x,2*xn-xc) - (2*fn-fc) #  f.subs(x,xn)
        sol  =  sp.solve(eqs,a,b,c)
        f    =  f.subs(a,sol[a]).subs(b,sol[b]).subs(c,sol[c])
        result['ddx'  , bc_type] = [f.diff(x)        .subs(x,0).diff(f_) for f_ in [fs,fc,fn]]
        result['d2dx2', bc_type] = [f.diff(x).diff(x).subs(x,0).diff(f_) for f_ in [fs,fc,fn]]
    return f'lambda xs,xc,xn: {result}'

# assert gen_fd_coeffs().replace(' ','') == fd_coeffs_ref.replace(' ','')
# print('INFO: Verified expression for fd_coeffs')

def ibm_deepcopy(ibm_data):
    for tag in ibm_data.keys():
        for ind in range(len(ibm_data[tag]['target'])):
            ibm_data[tag]['target'][ind] = lmap(int  , ibm_data[tag]['target'][ind])
            ibm_data[tag]['coeffs'][ind] = lmap(float, ibm_data[tag]['coeffs'][ind])
            ibm_data[tag]['inds'][ind]   = [lmap(int,ijk) for ijk in ibm_data[tag]['inds'][ind]]
    return ibm_data

def ibm_fix(ibm_data, n):
    for tag in ibm_data.keys():
        fixed_lo,fixed_hi = 0,0
        for ind,ijk_t in enumerate(ibm_data[tag]['target']):
            for ax in range(2):
                if ijk_t[ax] < 0:
                    ijk_t[ax]                                      += n[ax]
                    for ijk in ibm_data[tag]['inds'][ind]: ijk[ax] += n[ax]
                    fixed_lo += 1
                if ijk_t[ax] >= n[ax]:
                    ijk_t[ax]                                      -= n[ax]
                    for ijk in ibm_data[tag]['inds'][ind]: ijk[ax] -= n[ax]
                    fixed_hi += 1
                assert       0 <= ijk_t[ax] <   n[ax]
                assert all((-2 <= ijk  [ax] <= (n[ax]+1)) for ijk in ibm_data[tag]['inds'][ind])
        print(f"(tag, fixed_lo, fixed_hi): {(tag, fixed_lo, fixed_hi)}")
    return ibm_data

def ibm_reorder(ibm_data):
    for key in ibm_data:
        assert sorted(ibm_data[key].keys()) == ['coeffs', 'inds', 'target'], sorted(ibm_data[key].keys())
        coeffs, inds, target = [ibm_data[key][kk] for kk in ['coeffs', 'inds', 'target']]
        re_order =  sorted(range(len(target)), key = lambda i: tuple(target[i][::-1]))
        ibm_data[key]['inds']     =  [inds[i]   for i in re_order]
        ibm_data[key]['coeffs']   =  [coeffs[i] for i in re_order]
        ibm_data[key]['target']   =  [target[i] for i in re_order]
    return ibm_data

def ibm_add_extra_copies(ibm_data, nx_unit, nz_unit, cc_x, cc_z):
    for key in ibm_data:
        assert sorted(ibm_data[key].keys()) == ['coeffs', 'inds', 'target'], sorted(ibm_data[key].keys())
        ibm_data[key]['target'] = [ [i+ii*nx_unit, j+jj*nz_unit, k]                          for (i,j,k) in ibm_data[key]['target'] for ii in range(cc_x) for jj in range(cc_z)]
        ibm_data[key]['inds']   = [[[i+ii*nx_unit, j+jj*nz_unit, k] for (i,j,k) in all_ijk]  for all_ijk in ibm_data[key]['inds']   for ii in range(cc_x) for jj in range(cc_z)]
        ibm_data[key]['coeffs'] = [ row                                                      for row     in ibm_data[key]['coeffs'] for ii in range(cc_x) for jj in range(cc_z)]
    return ibm_data

def ibm_writer(ibm_data, geom_folder):
    for key in ibm_data:
        with open(pjoin(geom_folder,f'ibm_coeffs_{key}.dat'),'wb') as fw:
            fw.write(struct.pack('i',len(     ibm_data[key]['target'])))
            for ijk_t, all_ijk, all_cc in zip(ibm_data[key]['target'],
                                              ibm_data[key]['inds'  ],
                                              ibm_data[key]['coeffs']):
                [fw.write(struct.pack('i',s)) for s in ijk_t]
                assert len(all_ijk) == len(all_cc) == 2
                for cc,ijk in zip(all_cc, all_ijk):
                    fw.write(struct.pack('d',cc))
                    [fw.write(struct.pack('i',s)) for s in ijk]

def ibm_track_slabs(ibm_data,nx,nz,side,geom_folder):
    f = dict(bot=max, top=min)[side]
    for  tag   in ibm_data.keys():
        slab  =  [[None for _ in range(nz)] for _ in range(nx)]
        for i,j,k in ibm_data[tag]['target']:
            if slab[i][j] is None: slab[i][j] = k
            slab[i][j] = f(slab[i][j], k)
        assert not any((val is None) for row in slab for val in row), (tag,nx,nz,slab,ibm_data[tag])
        writer_ibm_slab(geom_folder, tag, side, slab)

def writer_ibm_slab(geom_folder, tag, side, slab):
    nx, nz  =  len(slab), len(slab[0])
    with open(pjoin(geom_folder, f'ibm_inds_{tag}_{side}.dat'),'wb') as fw:
        [fw.write(struct.pack('i',s)) for s in [nx, nz]]
        for     j in range(nz):
            for i in range(nx):
                fw.write(struct.pack('i',slab[i][j]))

def mk_P_abc(yp_ext, is_dirich_first=False):
    assert (type(yp_ext) == list) and (type(yp_ext[0]) == float)
    yp_pad     =  lmap(float,yp_ext)
    ny         =  len(yp_ext)-2
    a,b,c      =  [[0. for _ in range(ny)] for _ in range(3)]
    yp_pad[ 0] =  2*yp_ext[ 0] - yp_ext[ 1]
    yp_pad[-1] =  2*yp_ext[-1] - yp_ext[-2]
    for i in range(ny):
       x0   =  yp_pad[i ]
       x1   =  yp_pad[i+1]
       x2   =  yp_pad[i+2]
       dy   =  0.5*( x2 - x0)
       a[i] =  1.0/((x1 - x0)*dy)
       c[i] =  1.0/((x2 - x1)*dy)
       b[i] = -(a[i]    + c[i])
    if is_dirich_first:  b[0] -= a[0]
    else:                b[0] += a[0]
    b[-1] += c[-1]
    a[0]   = 0
    c[-1]  = 0
    return [lmap(float,row) for row in zip(a,b,c)]

# ------------------------------------------------------------------------
#                  Direct Runner
# ------------------------------------------------------------------------

if __name__ == '__main__':
    geom_folder              =  pjoin(dirname(_folder_)         ,'geom_data')
    input_folder             =  pjoin(dirname(dirname(_folder_)),'input')
    fd_coeffs                =  eval(fd_coeffs_ref)
    yu                       =  lmap(float,lfilter(None,[x.strip(' ') for x in reader(pjoin(input_folder,'yu.txt'))]))
    ny                       =  len(yu) - 1
    json                     =  import_fort_dat(pjoin(input_folder,'params.dat'))

    assert sum((key in json) for key in ['nu'  , 'Ret']) == 1
    assert sum((key in json) for key in ['cond', 'Pr' ]) == 1
    if 'Ret' in json: json['nu']   =         1./float(json['Ret'])
    if 'Pr'  in json: json['cond'] = json['nu']/float(json['Pr'])

    if (json['snap_iter_start']<0) or (json['snap_freq']<0):
        print('INFO: disabling snap_iter_start,snap_freq')
        json['snap_freq']       = -1
        json['snap_iter_start'] = -1
    json['avg_iter_start']   =  max(1,json['avg_iter_start'])
    json['snap_iter_start']  =  max(1,json['snap_iter_start'])
    yp             =  [(yu[i]+yu[i+1])/2 for i in range(ny)]
    yp_ext         =  [yu[0]] + yp + [yu[-1]]
    nstep          =  json['avg_freq']*json['avg_nraw']*max(1,2**json['avg_bin_nbins']) + (json['avg_iter_start']-1)
    assert nstep >= 1
    get_extrap_coeffs = lambda xc,xa,xb: [float((xc  - xa)/(xb - xa)*(-1) + 1), float((xc  - xa)/(xb - xa))]
    # y = (xc  - xa)/(xb - xa)*(yb - ya) + ya
    extrap_T_wall_0 , extrap_T_wall_1   =  get_extrap_coeffs(* yp_ext[  :3])
    extrap_T_wall_m1, extrap_T_wall_m2  =  get_extrap_coeffs(*(yp_ext[-3: ][::-1]))
    def check_forcings():
        bulk_target_Ub, bulk_use_target_Ub,    =  float(json['bulk_target_Ub']) , int(json['bulk_use_target_Ub'])
        bulk_target_Reb, bulk_use_target_Reb   =  float(json['bulk_target_Reb']), int(json['bulk_use_target_Reb'])
        force_Sf_x_pressure, force_Sf_x_accel  =  int(json['force_Sf_x_pressure']), int(json['force_Sf_x_accel'])
        if not bulk_use_target_Ub : assert abs(bulk_target_Ub ) < 1e-10
        if not bulk_use_target_Reb: assert abs(bulk_target_Reb) < 1e-10
        assert (bulk_use_target_Ub  + bulk_use_target_Reb) < 2
        assert (force_Sf_x_pressure + force_Sf_x_accel   ) < 2
        if (bulk_use_target_Ub  + bulk_use_target_Reb):
            assert (force_Sf_x_pressure + force_Sf_x_accel)
        if (force_Sf_x_pressure + force_Sf_x_accel   ):
            assert (bulk_use_target_Ub  + bulk_use_target_Reb)
            assert not int(json['use_incompressible'])
    check_forcings()
    A_base_params = [f"{json['nx']} {json['nz']} {json['ny']} ! nx nz ny_total"      ,
                     f"{int(json['use_rough_surf'])} {int(json['use_rough_surf_quady'])} ! use_rough_surf use_rough_surf_quady",
                     f"{nstep} {json['dtmax']} ! nstep dtmax"                        ,
                     f"{json['Lx']/json['nx']} {json['Lz']/json['nz']} ! dx dz"      ,
                     f"{json['nu']} {json['cond']} ! nu cond"                        ,
                     f"{json['Sf_x']} {json['Sf_z']} {json['Sq']} ! Sf_x Sf_z Sq"    ,
                     f"{json['Gb_x']} {json['Gb_z']} {json['Gb_y']} ! Gb_x Gb_z Gb_y",
                     f"{json['avg_iter_start']} {json['avg_nraw']} {json['avg_freq']} {json['avg_bin_nbins']} ! avg_iter_start avg_nraw avg_freq avg_bin_nbins",
                     f"{json['snap_iter_start']} {json['snap_freq']} ! snap_iter_start snap_freq"                                                      ,
                     f"{json['bulk_nprint']} ! bulk_nprint"                                                                                            ,
                     f"{json['filtering_steps']} {json['filtering_first']} {json['filtering_nstop']} ! filtering_steps filtering_first filtering_nstop",
                     f"{json['wall_BC_U_bot']} {json['wall_BC_V_bot']} {json['wall_BC_W_bot']} {json['wall_BC_T_bot']} ! wall_BC_U/V/W/T_bot",
                     f"{json['wall_BC_U_top']} {json['wall_BC_V_top']} {json['wall_BC_W_top']} {json['wall_BC_T_top']} ! wall_BC_U/V/W/T_top",
                     f"{int(json['use_incompressible'])} ! use_incompressible",
                     f"{int(json['use_utau_work'])} ! use_utau_work",
                     f"{extrap_T_wall_0} {extrap_T_wall_1} {extrap_T_wall_m1} {extrap_T_wall_m2} ! extrap_T_wall_0/1/m1/m2",
                     f"{float(json['bulk_target_Ub'])} {int(json['bulk_use_target_Ub'])} {int(json['bulk_use_rho_weighting'])} ! bulk_target_Ub bulk_use_target_Ub bulk_use_rho_weighting",
                     f"{float(json['bulk_target_Reb'])} {float(json['Ly'])} {int(json['bulk_use_target_Reb'])} ! bulk_target_Reb Ly bulk_use_target_Reb",
                     f"{int(json['force_Sf_x_pressure'])} {int(json['force_Sf_x_accel'])} {int(json['force_Sf_x_semilocal'])} ! force_Sf_x_pressure force_Sf_x_accel force_Sf_x_semilocal",
                     ]
    if     int(json['force_Sf_x_semilocal'])     : assert int(json['force_Sf_x_pressure'])
    if     int(json['use_rough_surf_quady']): assert int(json['use_rough_surf'])
    if not int(json['use_rough_surf'])      : assert not int(json['use_rough_surf_quady'])
    assert not any((c.lower() in line.lower()) for line in A_base_params for c in ['True','False'])
    yp_pad  =  [2*yu[0] - yp[0]] + tolist(yp) + [2*yu[-1] - yp[-1]]
    writer(pjoin(geom_folder,f'base_params.txt'),  A_base_params)
    result         =  dict(inv_dy_p     = [1./(yu[i+1]-yu[i])         for i in range(ny)],
                           inv_dy_u     = [1./(yp_ext[i+1]-yp_ext[i]) for i in range(ny+1)],
                           inv_dy_p_pad = [1./(yp_pad[i+1]-yp_pad[i]) for i in range(ny+1)],
                           yp           = tolist(yp))
    result.update({f"R{x}_{y}_abc"        :[]  for x in ['uw','v']  for y in ['ddy','d2dy2']})
    result.update({'Rv_interp_uw_bottom_top'   : []})
    result.update({"dPdy_ddy_lh": []})
    #
    result["P_band_abc"   ] = mk_P_abc(yp_ext,is_dirich_first=False)
    result["P_band_abc_00"] = mk_P_abc(yp_ext,is_dirich_first=True )
    #
    def add_Ruw():
        for i in range(ny):
            coeffs  =  fd_coeffs(yp_ext[i]   - yp_ext[i+1], # xs
                                 0                        , # xc
                                 yp_ext[i+2] - yp_ext[i+1]) # xn
            if   i == 0     : bc_type = 'dirich_s'
            elif i == (ny-1): bc_type = 'dirich_n'
            else            : bc_type = 'normal'
            result[f"Ruw_ddy_abc"  ].append(coeffs['ddx'  , bc_type])
            result[f"Ruw_d2dy2_abc"].append(coeffs['d2dx2', bc_type])
    add_Ruw()
    #
    def add_Ruw_pad():
        result[f"Ruw_pad_ddy_abc"] = []
        result[f"Ruw_pad_d2dy2_abc"] = []
        for i in range(ny):
            coeffs  =  fd_coeffs(yp_pad[i]   - yp_pad[i+1], # xs
                                 0                        , # xc
                                 yp_pad[i+2] - yp_pad[i+1]) # xn
            bc_type = 'normal'
            result[f"Ruw_pad_ddy_abc"].append(coeffs['ddx'  , bc_type])
            result[f"Ruw_pad_d2dy2_abc"].append(coeffs['d2dx2', bc_type])
    add_Ruw_pad()
    #
    def add_Rv():
        for i in range(ny+1):
            if i in [0,ny]:
                 coeffs = {key:[0.,0.,0.] for key in fd_coeffs(-1,0,1)}
            else:
                coeffs  =  fd_coeffs(yu[i-1] - yu[i], # xs
                                     0              , # xc
                                     yu[i+1] - yu[i]) # xn
            bc_type = 'normal'
            result[f"Rv_ddy_abc"  ].append(coeffs['ddx'  , bc_type])
            result[f"Rv_d2dy2_abc"].append(coeffs['d2dx2', bc_type])
    add_Rv()
    #
    def add_Rv_pad():
        result[f"Rv_pad_ddy_abc"] = []
        yu_pad = [2*yu[0] - yu[1]] + tolist(yu) + [2*yu[-1] - yu[-2]]
        for i in range(ny+1):
            if i in [0,ny]:
                 coeffs = {key:[0.,0.,0.] for key in fd_coeffs(-1,0,1)}
            else:
                coeffs  =  fd_coeffs(yu_pad[i]   - yu_pad[i+1], # xs
                                     0                        , # xc
                                     yu_pad[i+2] - yu_pad[i+1]) # xn
            bc_type = 'normal'
            result[f"Rv_pad_ddy_abc"].append(coeffs['ddx'  , bc_type])
    add_Rv_pad()
    #
    def mk_dPdy_ddy():
        for i in range(ny+1):
            if i in [0,ny]:
                result["dPdy_ddy_lh"].append([0.,0.])
            else:
                dy  =  yp[i] - yp[i-1]
                result["dPdy_ddy_lh"].append([-1/dy, 1/dy])
    mk_dPdy_ddy()
    #
    def mk_Rv_interp():
        n_passed = 0
        for i,yyu in enumerate(yu):
            if i in [0,ny]:
                result['Rv_interp_uw_bottom_top'].append([0. , 0.])
                n_passed += 1
            else:
                dy_top     =  yp[i] - yyu
                dy_bottom  =          yyu - yp[i-1]
                result['Rv_interp_uw_bottom_top'].append([ 1 - (dy_bottom/(dy_top + dy_bottom)) ,
                                                           1 - (dy_top   /(dy_top + dy_bottom)) ])
        assert n_passed == 2
    mk_Rv_interp()
    def mk_geom_vars_write():
        for key,A in result.items():
            if type(A[0]) in [float,int]:
                A = [[x] for x in A]
            A = [str(len(A)), *[' '.join(map(str,row)) for row in A]]
            writer(pjoin(geom_folder,f'{key}.dat'), A)
    mk_geom_vars_write()

    [os_system(f"mkdir -p {pjoin(dirname(geom_folder),kk)}") for kk in ['end','avg','prof_1d']]

    writer(pjoin(dirname(_folder_),"signal_runtime.txt"),
    ["0                         ! signal_runtime (see below)",
     "dtmax                     ! new [real(dp):: dtmax                     ] (for signal_runtime == 1)",
     "bulk_nprint               ! new [integer :: bulk_nprint               ] (for signal_runtime == 1)",
     "snap_iter_start snap_freq ! new [integer :: snap_iter_start, snap_freq] (for signal_runtime == 1)",
     ""                                                                                                 ,
     "! signal_runtime:"                                                                                ,
     "!     0    :  continue"                                                                           ,
     "!     1    :  change parameters (dtmax, bulk_nprint, snap_iter_start, snap_freq)"                 ,
     "!     other:  stop"                                                                               ])

    if json['use_rough_surf']:

        assert json['nx'] % json['nreps_x'] == 0
        assert json['nz'] % json['nreps_z'] == 0

        xp, xu      = mk_pu(json['Lx']/json['nreps_x'], json['nx']//json['nreps_x'])
        zp, zu      = mk_pu(json['Lz']/json['nreps_z'], json['nz']//json['nreps_z'])

        assert (len(xp) == len(xu)) and (len(zp) == len(zu))
        #
        assert len(yu)%2
        coords = dict(xp=xp, xu=xu, yu=yu, yp=yp, zp=zp, zu=zu)
        get_mm = dict(U='upp',
                      W='pup',
                      V='ppu',
                      T='ppp')
        all_commands = []
        for tag in 'UVWT':
            ibm_folder = pjoin(_folder_,f'ibm_fortran_{tag}')
            if isdir(ibm_folder): os_system(f"rm -rf {ibm_folder}")
            os_system(f"mkdir -p {ibm_folder}")
            xx,zz,yy  =  [tolist(coords[f'{c}{t}']) for c,t in zip('xzy',get_mm[tag])]
            writer(pjoin(ibm_folder,'ibm_params.txt'),lmap(str,[f"{tag} {len(xx)} {len(zz)} {len(yy)} {json['nreps_x']} {json['nreps_z']} {int(json['use_rough_surf_quady'])}",
                                                                f"{json[f'wall_BC_{tag}_bot']} {json[f'wall_BC_{tag}_top']}",
                                                                *xx,  *zz,  *yy]))
            all_commands.append(';'.join(lfilter(None,lmap(lambda x: x.strip(' '),
                       """
                       cd [pwd]
                       mkdir [pwd]/fortran_surf_files
                       cp -f ../ibm_coeffs_engine.f90 [pwd]
                       cp -f ../main_runner_ibm_coeffs_engine.f90 [pwd]
                       cp -f ../../../input/params.dat [pwd]
                       cp -f ../../../input/rough_surf_define.f90 [pwd]
                       cp -f ../../../input/fortran_surf_files/* [pwd]/fortran_surf_files
                       gfortran -fdefault-real-8 -fdefault-double-8 rough_surf_define.f90 ibm_coeffs_engine.f90  main_runner_ibm_coeffs_engine.f90
                       ./a.out
                       """.replace('[pwd]',ibm_folder).split('\n')))))
        run_commands(all_commands, nproc=4)
