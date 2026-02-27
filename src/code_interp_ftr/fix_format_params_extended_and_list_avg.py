
# ------------------------------------------------------------------------
#                  Generic Libraries
# ------------------------------------------------------------------------

from  os.path  import  dirname, abspath, isfile, isdir, basename
from  os.path  import  join                                       as  pjoin
from  os       import  listdir                                    as  os_listdir

# ------------------------------------------------------------------------
#                  Basic Functions
# ------------------------------------------------------------------------

def reader(fname):
    with open(fname,'r') as f:
        return [x.rstrip('\n') for x in f]

def writer(fname,A):
    with open(fname,'w') as f:
        for x in A:
            f.write(x+'\n')

lmap                  =  lambda f,x: list(map(f,x))
lfilter               =  lambda f,x: list(filter(f,x))
listdir_full          =  lambda x: sorted([pjoin(x,y) for y in os_listdir(x)])
listdir_full_files    =  lambda x: lfilter(isfile, listdir_full(x))
listdir_full_folders  =  lambda x: lfilter(isdir , listdir_full(x))

# ------------------------------------------------------------------------
#                  Direct Runner
# ------------------------------------------------------------------------

if __name__ == '__main__':
    _folder_    =  dirname(abspath(__file__))

    def fix_fmt_extended_params():
        fname       =  pjoin(_folder_, 'extended_params.dat')
        A           =  list(filter(None, [x.strip(' ') for x in reader(fname)]))
        A           =  [[x.strip(' ') for x in row.split('=')] for row in A]
        L0          =  max(map(len,[x[0] for x in A]))
        fmt         =  eval('lambda x: f"{x:'+str(L0)+'s}"')
        writer(fname, [f"{fmt(x)} = {eval(y)}" for x,y in A])
    fix_fmt_extended_params()
    
    def list_iters_avg():
        get_int    =  lambda x: int(basename(x).split('avg_iter_')[1].split('.')[0].split('_')[0])
        iters_avg  =  sorted(set(lmap(get_int,listdir_full_files(pjoin(dirname(dirname(_folder_)), 'avg')))))
        fname       =  pjoin(_folder_, 'iters_avg.dat')
        writer(fname, lmap(str,iters_avg))
    list_iters_avg()


    