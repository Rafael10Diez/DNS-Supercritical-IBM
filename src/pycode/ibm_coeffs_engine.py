# ------------------------------------------------------------------------
#                  Generic Libraries
# ------------------------------------------------------------------------

import  numpy                    as  np
from    copy   import  deepcopy

# ------------------------------------------------------------------------
#                  Basic Functions
# ------------------------------------------------------------------------

lmap  =  lambda f,x:  list(map(f,x))
lzip  =  lambda  *x:  list(zip(*x))

class Mk_Extended:
    def __init__(self, x, nper = 10):
        
        self._nper  =  nper

        x    =  x.tolist()
        for _ in range(nper):
            x    =  [2*x[0] - x[1]] + x + [2*x[-1] - x[-2]]
        
        self._x_ext = x
    
    def __call__(self, i):
        i  =  i + self._nper
        assert 0 <= i < len(self._x_ext)
        return self._x_ext[i]

# ------------------------------------------------------------------------
#                  Trilinear Fit
# ------------------------------------------------------------------------

def search_mid_h(xzy_a, xzy_b, get_h):
    xzy_a       =  list(map(np.copy, xzy_a))
    xzy_b       =  list(map(np.copy, xzy_b))

    temp_xzy_a  =  deepcopy(xzy_a)
    temp_xzy_b  =  deepcopy(xzy_b)

    assert len(xzy_a) == len(xzy_b) == 3

    # find points where "a" is higher above the surface than "b" (should not happen)
    m           =  (xzy_a[2] - get_h(*xzy_a[:2])) > (xzy_b[2] - get_h(*xzy_b[:2])) 

    # reverse these points: (a,b = b,a)
    for i in range(3):
        xzy_a[i][m]  =  temp_xzy_b[i][m]
        xzy_b[i][m]  =  temp_xzy_a[i][m]
    del i

    dx, dz, dy  =  [(xzy_b[i] - xzy_a[i]) for i in range(3)]
    x0, z0, y0  =  xzy_a

    eval_h  =  lambda mid: get_h(dx*mid + x0,
                                 dz*mid + z0) - (dy*mid + y0) # positive when below surf (closer to "a")
    
    # define edge case, dx=dz=0
    # in such cases, (y0 = h) and (dy = 0) make code exit at once with right value (search is meaningless)
    m       =  (np.fabs(dx) < 1e-10) * (np.fabs(dz) < 1e-10)
    y0[m]   =  (eval_h(0.5) + (dy*0.5 + y0))[m]
    dy[m]   =  0.

    assert np.all(eval_h(1) <= (eval_h(0)+1e-10))
    assert np.all(         (xzy_b[2] - get_h(*xzy_b[:2])) >= -1e-10) # "b" is above surface
    assert np.all(1e-10 >= (xzy_a[2] - get_h(*xzy_a[:2]))          ) # "a" is below surface

    lo, hi  =  np.zeros_like(dx), np.ones_like(dx)
    iters   =  0

    while np.any((hi-lo) > 1e-10):
        
        iters += 1
        assert iters < 10_000, (lo, hi, xzy_a, xzy_b)

        mid   = (lo + hi) / 2
        new_h =  eval_h(mid)

        m     = (new_h-1e-10) <= 0 # above surf, near b=hi
        hi[m] = mid[m]

        m     = (new_h+1e-10) >= 0 # below surf, near a=lo
        lo[m] = mid[m]

    mid  =  (lo + hi) / 2

    return ( dx*mid + x0,
             dz*mid + z0,
             dy*mid + y0) # eval_h is supposed to be zero (or close)

def apply_Fadlun(fix_points, get_h, force_interior = False):
    print('Note: applying Ghost-version of Fadlun IBM scheme!')
    all_args    =   [None for _ in fix_points]
    all_coeffs  =  [[None for _ in range(2)] for _ in fix_points]
    all_inds    =  [[None for _ in range(2)] for _ in fix_points]
    
    all_xn, all_zn, all_yn  =  search_mid_h([[d[f"{c}2"]               for d in fix_points] for c in 'xzy'], 
                                            [[d['fluid_points'][-1][c] for d in fix_points] for c in 'xzy'],
                                            get_h                                                          )
    for ind in range(len(all_inds)):
        d                 =  fix_points[ind]

        ijk_pi  =  [[arr[c+'_ext'] for c in 'ijk'] for arr in d['fluid_points']]
        xzy_pi  =  [[arr[c       ] for c in 'xzy'] for arr in d['fluid_points']]

        xn,zn,yn = xzy_n  =  all_xn[ind], all_zn[ind], all_yn[ind]
        xzy_g             =  [2*xn - d['x2'],
                              2*zn - d['z2'],
                              2*yn - d['y2']]
        
        xzy_p1, xzy_p2  =  xzy_pi[:2]
        ijk_p1, _       =  ijk_pi[:2]

        KK                    =  1/sum((p2-p1)**2 for p1,p2 in zip(xzy_p1, xzy_p2))
        ee                    =  [(p2-p1)*KK for p1,p2 in zip(xzy_p1, xzy_p2)]
        norm                  =  lambda xzy: sum((c-p1)*val for c,p1,val in zip(xzy, xzy_p1, ee))
        contained             =  lambda target, A , B, norm=norm: [(norm(A)-1e-10) <= norm(target) <= (norm(B)+1e-10),
                                                                   1-(norm(target)-norm(A))/(norm(B)-norm(A))]

        contain_n1, alpha_n   =  contained(xzy_g, xzy_n , xzy_p1)
        contain_12, alpha_p1  =  contained(xzy_g, xzy_p1, xzy_p2)

        # assert (contain_n1 + contain_12) >= 1, (alpha_n, alpha_p1, xzy_g, xzy_n , xzy_p1, xzy_p2, side, ind, d)
        if contain_n1 and contain_12:
            assert abs(alpha_n)<1e-10 and abs(alpha_p1-1)<1e-10
            
        all_coeffs[ind] = [0.     for _ in all_coeffs[ind]]
        all_inds  [ind] = [ijk_p1 for _ in all_inds  [ind]]

        if contain_n1 and (not force_interior):
            # use point, if contained here, and not forcing interior
            #     -> forced point is used later in first triplet (ip,ip+1)
            # all_inds  [ind][0]  =    ijk_p1  # done before
            all_coeffs[ind][0]  =  -(1-alpha_n) # negative signs
        else:
            for ip in range(len(ijk_pi)-1):
                contain__, alpha___  =  contained(xzy_g, xzy_pi[ip], xzy_pi[ip+1])
                if contain_n1 and force_interior:
                    contain__ = True # force previous choice
                if contain__:
                    all_inds  [ind][0]  =    ijk_pi[ip]
                    all_coeffs[ind][0]  =  -(    alpha___) # negative signs
                    all_inds  [ind][1]  =          ijk_pi[ip+1]
                    all_coeffs[ind][1]  =  -(1 - alpha___)
                    break
            assert contain__, (xzy_g,xzy_n,xzy_p1,xzy_p2)
        
        all_args[ind] = {}
    return all_inds, all_coeffs, all_args

# ------------------------------------------------------------------------
#                 Make rough_surface
# ------------------------------------------------------------------------

def mk_d_xzy(mult_distance):
    result = [[dx,dz,dy] for dx in [0, -1, 1]
                         for dz in [0, -1, 1]
                         for dy in [0, -1, 1]]
    result.sort(key = lambda val: [sum(map(abs,val)), -abs(val[2]),-abs(val[0]), -abs(val[1]) ])
    assert result[0] == [0,0,0]
    result =  np.array(result[1:])
    result2 = []
    for mm in sorted(mult_distance):
        result2.extend((mm*result).tolist())
    return result2

def make_rough_surface_ghost(info, mult_distance = [1], force_interior = False, use_tags = 'UVWT', verbose=True):
    final_results = {}
    for tag, x, z, y in [['U', info['xu'], info['zp'], info['yp']],
                         ['W', info['xp'], info['zu'], info['yp']],
                         ['V', info['xp'], info['zp'], info['yu']],
                         ['T', info['xp'], info['zp'], info['yp']]]:
        if not tag in use_tags: continue
        x, z, y         =  lmap(lambda arr: np.array(arr).reshape(-1),[x, z, y])
        x_ext           =  Mk_Extended(x)
        z_ext           =  Mk_Extended(z)
        # a) find points to be fixed
        def get_fix_points():
            fix_points      =  []
            assert y.tolist() == sorted(y)
            is_fluid    =  info['is_fluid'](x.reshape(-1,1,1), z.reshape(1,-1,1),y.reshape( 1, 1,-1))
            is_ok       =  np.copy(is_fluid)
            nx, nz, ny  =  is_fluid.shape
            if verbose: print(f"IBM Stats ({tag}): ([nx,nz,ny]: [{nx:4d},{nz:4d},{ny:4d}])")
            for dx,dz,dy in  mk_d_xzy(mult_distance):
                L_old = len(fix_points)
                for          ii in  range(is_fluid.shape[0]):
                    for      jj in  range(is_fluid.shape[1]):
                        for  kk in  range(is_fluid.shape[2]):
                            if is_fluid[ii,jj,kk]:
                                i2,j2,k2  =  ( (ii + dx)%nx ,
                                               (jj + dz)%nz ,
                                                kk + dy     )
                                if 0<= k2 < ny:
                                    if not is_ok[i2,j2,k2]:
                                        mm                =  max(map(abs, [dx,dz,dy]))
                                        dx_u, dz_u, dy_u  =  dx//mm, dz//mm, dy//mm
                                        assert ((dx_u*mm) == dx) and ((dz_u*mm) == dz) and ((dy_u*mm) == dy) 
                                        fix_points.append(dict(i2_ext       = ii + dx                                                         ,
                                                               j2_ext       = jj + dz                                                         ,
                                                               k2_ext       = kk + dy                                                         ,
                                                               x2           = x_ext(ii + dx)                                                  ,
                                                               z2           = z_ext(jj + dz)                                                  ,
                                                               y2           =     y[k2]                                                       ,
                                                               x_start      = (x_ext(ii) + x_ext(ii + dx))/2                                  ,
                                                               z_start      = (z_ext(jj) + z_ext(jj + dz))/2                                  ,
                                                               fluid_points = [dict(i_ext =      (ii-dx_u*m), j_ext  =      (jj-dz_u*m), k_ext =   kk-dy_u*m ,
                                                                                    x     = x_ext(ii-dx_u*m), z      = z_ext(jj-dz_u*m), y     = y[kk-dy_u*m]) for m in range(mm+1)],
                                                              ))
                                        is_ok[i2,j2,k2] = True
                if verbose: print(f"        ([dx,dz,dy]: [{dx:2d},{dz:2d},{dy:2d}]): {(len(fix_points)-L_old):9d}")
            return fix_points
        fix_points = get_fix_points()
        result_inds, result_coeffs, result_args  =  apply_Fadlun(fix_points, info['get_h'], force_interior = force_interior)
        result_target       =  [[d[f'{c}2_ext'] for c in 'ijk'] for d in fix_points]
        final_results[tag]  =  dict(inds    =  result_inds  , 
                                    coeffs  =  result_coeffs, 
                                    args    =  result_args  , 
                                    target  =  result_target)
    return final_results
