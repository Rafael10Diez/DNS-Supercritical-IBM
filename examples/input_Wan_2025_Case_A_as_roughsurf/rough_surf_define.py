# ------------------------------------------------------------------------
#                  Generic Libraries
# ------------------------------------------------------------------------

import  numpy                              as  np
from    os.path  import  join              as  pjoin
from    os.path  import  dirname, abspath

# ------------------------------------------------------------------------
#                  Basic Functions
# ------------------------------------------------------------------------

_folder_ = dirname(abspath(__file__))

lmap     =  lambda f,x: list(map(f,x))
lfilter  =  lambda f,x: list(filter(f,x))

def reader(fname):
    with open(fname,'r') as f:
        return [x.rstrip('\n') for x in f]

# ------------------------------------------------------------------------
#                  Height Function
# ------------------------------------------------------------------------

class Height_Function:
    def __init__(self,verbose=True):
        self.import_params(verbose=verbose)
        assert all(hasattr(self,key) for key in ['Lx', 'Lz', 'Ly', 'nreps_x', 'nreps_z'])
        self._avg_height           =  0.
    
    def get_h(self, x, z, side = None, lib = np):
        assert (type(side) == str) and (side in ['bottom', 'top'])
        if not hasattr(x,'shape'):
            x, z    =  np.array(x), np.array(z)
        result  =  self._avg_height + lib.zeros_like(x+z)
        if side == 'top':
            result = self.Ly - result
        else:
            assert side == 'bottom'
        return result
    
    def import_params(self,verbose=True):
        fname_params = pjoin(_folder_, 'params.dat')
        A = reader(fname_params)
        A = lfilter(None,[x.split('!')[0].strip(' ') for x in A])
        A = {key.strip(' '):eval(val) for key,val in [x.split('=') for x in A]}
        print('\n---------------- Height Function ----------------')
        print(f'Scanned Parameters: ({fname_params})')
        for key,f in [['Lx'     , float], 
                      ['Lz'     , float], 
                      ['Ly'     , float], 
                      ['nreps_x', int  ], 
                      ['nreps_z', int  ]]:
            val = f(A[key])
            setattr(self, key, val)
            print(f'    {key:8s}: {val}')
        print('')

# ------------------------------------------------------------------------
#                  Direct Runner
# ------------------------------------------------------------------------

if __name__ == '__main__':
    Height_Function()