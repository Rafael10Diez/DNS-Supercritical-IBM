
# ------------------------------------------------------------------------
#                  Generic Libraries
# ------------------------------------------------------------------------

from    os.path  import  basename
from    os.path  import  join     as  pjoin
import  numpy                     as  np

# ------------------------------------------------------------------------
#                  Basic Functions
# ------------------------------------------------------------------------

lmap    =  lambda f,x: list(map(f,x))
tolist  =  lambda x: x.tolist() if hasattr(x,'tolist') else x

def reader(fname):
    with open(fname, 'r') as f:
        return [x.rstrip('\n') for x in f]

def find_interval(A, val):
    lo_orig, hi_orig  =  0, len(A)-1
    lo     , hi       =  lo_orig, hi_orig
    while abs(hi-lo)>1:
        mid  =  (lo+hi)//2
        if A[mid] <= val: lo = mid
        if A[mid] >= val: hi = mid
    if lo == hi == lo_orig: hi = lo_orig + 1
    if lo == hi == hi_orig: lo = hi_orig - 1
    assert (lo+1) == hi
    return [lo,hi]

# ------------------------------------------------------------------------
#                  Extrapolation Coefficient
# ------------------------------------------------------------------------

def extrap_coeffs(g, p1, p2):
    mag         =  lambda x: np.sqrt((x**2).sum())
    signed_dist = lambda x: (x*n).sum()
    g, p1, p2   =  lmap(np.array, [g, p1, p2])
    n           =  (p2-p1)/mag(p2-p1)
    x__m_x0     =  signed_dist( g-p1) #    x  - x0 :   g - p1
    x1_m_x0     =  signed_dist(p2-p1) #    x1 - x0 :  p2 - p1
    m           =  x__m_x0 / x1_m_x0
    # m = (y1-y0)/(x1-x0)*(x-x0) + y0
    #     (   -1)/(x1-x0)*(x-x0) + 1  : f0
    #     (1    )/(x1-x0)*(x-x0)      : f1
    return  -m + 1, m

# ------------------------------------------------------------------------
#                  Tridiagonal Solver
# ------------------------------------------------------------------------

class Tridiagonal_Solver:
    def __init__(self, a, b, c):
        self.abc = a+0., b+0., c+0.
        assert np.fabs(a[:,:, 0]).max() < 1e-10
        assert np.fabs(c[:,:,-1]).max() < 1e-10
        self.n = n = a.shape[2]
        assert     a.shape[2] ==     b.shape[2] ==     c.shape[2] == n
        assert len(a.shape)   == len(b.shape)   == len(c.shape)   == 3
        empty          =  lambda: np.zeros_like(b)*np.nan

        self.a         =  a
        self.b0        =  b[:,:,0]

        self.cp        =  empty()
        self.div_bacp  =  empty()

        self.cp[:,:,0] = c[:,:,0]/self.b0

        for i in range(1,n):
            self.div_bacp[:,:,i]  =  b[:,:,i] - a[:,:,i]*self.cp[:,:,i-1]
            self.cp[:,:,i]        =  c[:,:,i]/self.div_bacp[:,:,i]

        self.cp[:,:,-1] *= np.nan

    def __call__(self, d, inplace = True):
        if not inplace: d = d.copy()

        # make d'
        d[:,:,0] /= self.b0
        for i in range(1,self.n):
            d[:,:,i] = (d[:,:,i] - self.a[:,:,i]*d[:,:,i-1])/self.div_bacp[:,:,i]

        x = d
        for i in range(self.n-2,-1,-1):
            x[:,:,i] -= self.cp[:,:,i] * x[:,:,i+1]

        return x

# ------------------------------------------------------------------------
#                  Poisson Solver
# ------------------------------------------------------------------------

class Solver_Poisson_FFT_Thomas_3d:
    def __init__(self, shape, dx, dz, yu):
        yu      =  np.array(yu).reshape(-1)
        mk_p    =  lambda u: u[:-1] + 0.5*np.diff(u)
        yp_ext  =  np.array([yu[0]] + mk_p(yu).reshape(-1).tolist() + [yu[-1]])

        assert shape[2] == (len(yp_ext)-2)

        Nx, Nz, _         =  shape
        k_x               =  np.fft.rfftfreq(Nx, d=1/Nx)
        k_z               =  np.fft. fftfreq(Nz, d=1/Nz)

        mk_ax             =  lambda j,max_tot,dd: -4./(dd*dd)*(np.sin(j*np.pi/max_tot)**2)
        a_xz              =  (mk_ax(k_x, Nx, dx).reshape(-1, 1,1) + # a_x + a_z
                              mk_ax(k_z, Nz, dz).reshape( 1,-1,1))

        a,b,c             =  np.array(self.mk_abc(yp_ext)).transpose()
        a                 =  a.reshape(1,1,-1) + np.zeros_like(a_xz)
        b                 =  b.reshape(1,1,-1) + np.zeros_like(a_xz)
        c                 =  c.reshape(1,1,-1) + np.zeros_like(a_xz)

        a_00, b_00, c_00  =  np.array(self.mk_abc(yp_ext,is_dirich_first=True)).transpose()
        a[0,0,:]          =  a_00
        b[0,0,:]          =  b_00
        c[0,0,:]          =  c_00

        b                +=  a_xz
        self.Nx           =  Nx
        self.Nz           =  Nz

        self.get_rfft2   =  lambda rhs : np.fft. rfft2(rhs ,                      axes = (1, 0))
        self.tridiag     =  Tridiagonal_Solver(a,b,c)
        self.get_irfft2  =  lambda spec: np.fft.irfft2(spec, s=(self.Nz,self.Nx), axes = (1, 0))

    def __call__(self, rhs):
        return self.get_irfft2(self.tridiag(self.get_rfft2(rhs)))

    @staticmethod
    def mk_abc(yp_ext, is_dirich_first=False):
        yp_ext     =  np.array(yp_ext)
        yp_pad     =  yp_ext.copy()
        ny         =  len(yp_ext)-2
        a,b,c      =  [np.zeros_like(yp_ext,shape=ny) for _ in range(3)]
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
#                  Poisson Solver Double RFFT
# ------------------------------------------------------------------------

class Solver_Poisson_FFT_Thomas_3d_rfftx2:
    def __init__(self, shape, dx, dz, yu):
        yu      =  np.array(yu).reshape(-1)
        mk_p    =  lambda u: u[:-1] + 0.5*np.diff(u)
        yp_ext  =  np.array([yu[0]] + mk_p(yu).reshape(-1).tolist() + [yu[-1]])

        assert shape[2] == (len(yp_ext)-2)

        Nx, Nz, _         =  shape
        k_x               =  np.repeat(np.fft.rfftfreq(Nx, d=1/Nx),2)
        k_z               =  np.repeat(np.fft.rfftfreq(Nz, d=1/Nz),2)

        mk_ax             =  lambda j,max_tot,dd: -4./(dd*dd)*(np.sin(j*np.pi/max_tot)**2)
        a_xz              =  (mk_ax(k_x, Nx, dx).reshape(-1, 1,1) + # a_x + a_z
                              mk_ax(k_z, Nz, dz).reshape( 1,-1,1))

        a,b,c             =  np.array(self.mk_abc(yp_ext)).transpose()
        a                 =  a.reshape(1,1,-1) + np.zeros_like(a_xz)
        b                 =  b.reshape(1,1,-1) + np.zeros_like(a_xz)
        c                 =  c.reshape(1,1,-1) + np.zeros_like(a_xz)

        a_00, b_00, c_00  =  np.array(self.mk_abc(yp_ext,is_dirich_first=True)).transpose()
        a[:2,:2,:]          =  a_00
        b[:2,:2,:]          =  b_00
        c[:2,:2,:]          =  c_00

        b                +=  a_xz
        self.Nx           =  Nx
        self.Nz           =  Nz

        self.get_rfft2   =  lambda rhs : np.fft. rfft2(rhs ,                      axes = (1, 0))
        self.tridiag     =  Tridiagonal_Solver(a,b,c)
        self.get_irfft2  =  lambda spec: np.fft.irfft2(spec, s=(self.Nz,self.Nx), axes = (1, 0))

    def __call__(self, x):
        def as_rc_rfft(x,inv):
            if inv is None:
                x = np.fft.rfft(x,axis=0)
                y = np.zeros_like(x,shape=[2*x.shape[0],*x.shape[1:]],dtype=float)
                y[ ::2,...] = x.real
                y[1::2,...] = x.imag
            else:
                y = np.fft.irfft(x[::2,...] + 1j*x[1::2,...], axis=0,n=inv)
            return y
        as_2dfft = lambda x,inv=None: as_rc_rfft(as_rc_rfft(x, self.Nx if inv else None).transpose(1,0,2),
                                                               self.Nz if inv else None).transpose(1,0,2)
        return as_2dfft(self.tridiag(as_2dfft(x)),inv=True)

    @staticmethod
    def mk_abc(yp_ext, is_dirich_first=False):
        yp_ext     =  np.array(yp_ext)
        yp_pad     =  yp_ext.copy()
        ny         =  len(yp_ext)-2
        a,b,c      =  [np.zeros_like(yp_ext,shape=ny) for _ in range(3)]
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
#                  Test Poisson Solver
# ------------------------------------------------------------------------

def test_poisson_solver():
    class Test_Solution:
        def __init__(self):
            import sympy as sp
            sin, cos  =  sp.sin, sp.cos
            def fder(f, der, val):
                for _ in range(der):
                    f = f.diff(val)
                return f
            lapla                      =  lambda f: sum(fder(f,2,val) for val in vars.values())
            f1                         =  lambda x:  2*cos(x+cos(x))
            f2                         =  lambda x:  3*cos(x+sin(x))
            f3                         =  lambda x:  5*sin(x-sin(x)*cos(x))
            vars                       =  {key: sp.symbols(key, real=True) for key in 'xyz'}
            self.Lx, self.Ly, self.Lz  =  5.63, 2, 5.63/2
            self.L_xyz                 =  {key: getattr(self,f"L{key}") for key in 'xyz'}
            p = f1(2*sp.pi * vars['x'] / self.Lx) * \
                f2(2*sp.pi * vars['y'] / self.Ly) * \
                f3(2*sp.pi * vars['z'] / self.Lz)
            b = lapla(p)
            def fder(f, der, val):
                for _ in range(der):
                    f = f.diff(val)
                return f
            def check_derivs():
                for f in [p,b]:
                    for key,val in vars.items():
                        for der in [0,1,2]:
                            assert (fder(f,der,val).subs(val, 0              )  -
                                    fder(f,der,val).subs(val, self.L_xyz[key])  ).simplify() == 0
                val = vars['y']
                f   = p
                for der in [1]:
                    assert (fder(f,der,val).subs(val, 0      ).simplify()  ==
                            fder(f,der,val).subs(val, self.Ly).simplify()  == 0), (der,f)
            check_derivs()
            self.as_str  =  dict(p   = str(p)       ,
                                 rhs = str(b))
            self.as_f    =  {key: eval(f"lambda x,y,z, cos = np.cos, sin = np.sin, pi = np.pi: {val}") for key,val in self.as_str.items()}

    def grid_y(L, n):
        #return np.linspace(0, L, n+1)
        assert n%2
        y, dy  = [0], 1
        all_dy = iter(np.linspace(1,10,n//2).tolist())
        for _ in range(n//2):
            y.append(y[-1] + dy)
            dy = next(all_dy)
        assert len(y) == (n//2+1)
        y          =  np.array(y)/max(y)*0.5
        y[-1]      =  0.5
        other      =  np.cumsum(np.array(np.diff(y).tolist()[::-1]), 0)
        other[-1]  =  0.5
        y          =  y.tolist() + (0.5+other).tolist()
        y[-1]      = 1
        assert len(y) == n
        y = np.array(y)
        y = y[:-1] + 0.5*np.diff(y)
        y = [0.] + y.tolist() + [1.]
        return np.array(y)*L

    mk_p     =  lambda u: np.array(u)[:-1] + 0.5*np.diff(u)
    sol      =  Test_Solution()
    Linf     =  lambda x: float(np.fabs(x).max())
    print(f"{sol.as_str=}")
    prev = np.nan
    for N in [2**i for i in range(2,8)]:
        assert N <= 512
        nx = ny = nz =  N
        yu           =  grid_y(sol.Ly, ny+1)
        xp           =  mk_p(np.linspace(0, sol.Lx, nx+1)).reshape(-1, 1, 1)
        zp           =  mk_p(np.linspace(0, sol.Lz, nz+1)).reshape( 1,-1, 1)
        yp           =  mk_p(yu)                          .reshape( 1, 1,-1)

        dx           =  float(np.mean(np.diff(xp.reshape(-1))))
        dz           =  float(np.mean(np.diff(zp.reshape(-1))))

        rhs          =  sol.as_f['rhs'](xp, yp, zp)
        p_ref        =  sol.as_f['p']  (xp, yp, zp)
        p_ref       -=  np.average(p_ref)

        p_fft        =  Solver_Poisson_FFT_Thomas_3d(tuple(rhs.shape), dx, dz, yu)(rhs)
        p_fft       -=  np.average(p_fft)

        p_fft2       =  Solver_Poisson_FFT_Thomas_3d_rfftx2(tuple(rhs.shape), dx, dz, yu)(rhs)
        p_fft2       -=  np.average(p_fft2)
        error_x2_rfft = Linf(p_fft - p_fft2)
        assert error_x2_rfft < 1e-10
        error        =  Linf(p_ref-p_fft)
        ratio        =  prev/error
        prev         =  error
        print(f"N: {N:3d} (ratio: {ratio:7.3f}) (residual: {error:7.3e}) (error_x2_rfft: {error_x2_rfft:7.3e})")

# ------------------------------------------------------------------------
#                  Direct Runner
# ------------------------------------------------------------------------

if __name__ == '__main__':
    test_poisson_solver()