# ------------------------------------------------------------------------
#                  Generic Libraries
# ------------------------------------------------------------------------

from    os.path  import  abspath, dirname, basename, isfile, isdir
from    os.path  import  join                                       as  pjoin
from    os       import  mkdir                                      as  os_mkdir
from    os       import  listdir                                    as  os_listdir
from    os       import  system                                     as  os_system
import  sympy                                                       as  sp
import  numpy                                                       as  np

try:
    _folder_
except:
    _folder_  =  dirname(abspath(__file__))

# ------------------------------------------------------------------------
#                  Basic Functions
# ------------------------------------------------------------------------

lmap                  =  lambda f,x: list(map   (f,x))
lfilter               =  lambda f,x: list(filter(f,x))

listdir_full          =  lambda x: sorted(pjoin(x,y) for y in os_listdir(x))
listdir_full_files    =  lambda x: lfilter(isfile, listdir_full(x))
listdir_full_folders  =  lambda x: lfilter(isdir , listdir_full(x))

def pop1(A):
    x, = list(A)
    return x

def writer(fname,A):
    with open(fname,'w') as f:
        for x in A:
            f.write(str(x).strip('\n') + '\n')

def reader(fname):
    with open(fname,'r') as f:
        return [x.strip('\n') for x in f]

# ------------------------------------------------------------------------
#                  Direct Runner
# ------------------------------------------------------------------------

if __name__ == '__main__':
    for tag in ['_poi','_ibm']:
        bad_runs          =  []
        max_error_global  =  float('-inf')
        n_trials          = 0
        for trial_subfolder in listdir_full_folders(pjoin(_folder_,'trials')):
            n_trials     +=  1
            fname_error   =  sorted(lfilter(lambda x: basename(x) == 'max_error.dat', listdir_full_files(trial_subfolder)))
            if len(fname_error) == 1:
                data_error       = lfilter(None,[x.strip(' ') for x in reader(fname_error[0]) if ('=' in x)])
                data_error       = {key.strip(' '):eval(val.strip(' ')) for x in data_error for key,val in [x.split('=')]}
                max_error_global = max(max_error_global, data_error[f'max_error{tag}'])
                if data_error[f'max_error{tag}'] > 1e-10:
                    bad_runs.append([trial_subfolder, data_error[f'max_error{tag}']])
            else:
                bad_runs.append([trial_subfolder, float('nan')])
                max_error_global = float('inf')
        if bad_runs:
            line = '-'*60
            print(f"{line} Bad Runs {line}")
            for trial_subfolder, max_error in bad_runs:
                print(f"max_error: {max_error:12.3e}, {trial_subfolder}")
        assert (not bad_runs) and (max_error_global < 1e-10)
        print(f"SUCCESS! (all tests passed) (tag: {tag}) (n_trials: {n_trials:9d}) (max_error_global: {max_error_global:12.3e}) ({pjoin(_folder_,'trials')})")