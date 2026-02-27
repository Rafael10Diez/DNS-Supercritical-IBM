
# ------------------------------------------------------------------------
#                  Generic Libraries
# ------------------------------------------------------------------------

from    os         import  system                                     as  os_system
from    os.path    import  join                                       as  pjoin
from    os         import  listdir                                    as  os_listdir
from    os.path    import  abspath, dirname, basename, isfile, isdir
from    string     import  ascii_letters

# ------------------------------------------------------------------------
#                  Basic Functions
# ------------------------------------------------------------------------

_folder_              =  abspath(dirname(__file__))
lmap                  =  lambda f,x: list(map(f,x))
lfilter               =  lambda f,x: list(filter(f,x))
listdir_full          =  lambda   x: sorted([pjoin(x,y) for y in os_listdir(x)])
listdir_full_files    =  lambda   x: lfilter(isfile, listdir_full(x))
listdir_full_folders  =  lambda   x: lfilter(isdir , listdir_full(x))
get_pad               =  lambda x: len(x)-len(x.lstrip(' '))

def listdir_full_files_deep(x):
    def dfs(x):
        for y in listdir_full_files(x):
            yield y
        for y in listdir_full_folders(x):
            yield from dfs(y)
    return list(dfs(x))

def writer(fname,A):
    os_system(f'mkdir -p "{dirname(fname)}"')
    assert (not isfile(fname)) and (basename(dirname(fname)) != 'redu_fixed')
    with open(fname, 'w') as f:
        for x in A:
            f.write(str(x) + '\n')

def reader(fname):
    with open(fname, 'r') as f:
        return [x.rstrip('\n') for x in f]

def removesuffix(A,x):
    if A.endswith(x): return A[:-len(x)]
    return A

def removeprefix(A,x):
    if A[:len(x)] == x: return A[len(x):]
    return A

# ------------------------------------------------------------------------
#                  Utilities
# ------------------------------------------------------------------------

is_reduction_line  =  lambda x: (does_start(x,'!$acc ')  and  ('reduction(' in x))
valid_chars        =  set(['_',*lmap(str,range(10)),*list(ascii_letters)])
is_word            =  lambda x,s,i: (does_start(x,s,i)                                            and
                                     (((i-1)< 0     )      or (not (x[i-1]      in valid_chars))) and
                                     (((i+len(s))>=len(x)) or (not (x[i+len(s)] in valid_chars))) )

def ftn_line_split(A,lim=128,last_cont=False):
    A     =  str(A).split(' ')
    B     =  []
    line  =  ''
    for x in A:
        if (len(line)+len(x)+4)>=lim:
            B.append(' ' + line+' &')
            line = ''
        line  = f"{line} {x}".strip(' ')
    B.append(' ' + ((line+' &') if last_cont else line))
    return '\n'.join(B)

def ftn_text_file_check(A,lim=128):
    A = [x for line in A for x in line.split('\n')]
    for i in range(len(A)):
        if len(A[i])>lim:
            old  = A[i]
            A[i] = ftn_line_split(A[i],lim=lim,last_cont=False)
            if verbose: print('WARNING: breaking down: ', i, len(old),'\n', old,'\n', A[i])
            assert all(len(line)<=lim for line in A[i].split('\n'))
    return A


def does_start(x,s,i=0):
    x = x.lstrip(' ')
    for j in range(len(s)):
        if ((i+j)>=len(x)) or (x[i+j]!=s[j]): return False
    return True

# clean previous modifications
def clean_prev(A):
    i_def_if_pairs  =  lfilter(lambda x: f'!defined({macro_name})' in A[x[0]] ,
                               list(zip(*[[i for i in range(len(A)) if does_start(A[i],c)] for c in ['#if','#end']])))
    for i,j in reversed(i_def_if_pairs):
        B   =  A[i+1:j]
        k,  =  [i for i in range(len(B)) if does_start(B[i],'#else')]
        A   =  A[:i] + B[:k] + A[j+1:]
    return A

def get_reductions(x,A):
    pad   = len(x)
    x     = x.lstrip(' ')
    pad  -= len(x)
    tag   = 'reduction('
    all_i = [i for i in range(len(x)) if does_start(x,tag,i)]
    assert len(all_i)
    for i in reversed(all_i):
        begin,after =  x[:i], x[i:]
        j           =  [j for j in range(len(after)) if after[j]==')'][0]
        snippet     =  after[:j+1]
        after       =  after[j+1:]
        x           =  begin + after
        A.append(snippet)
    return ' '*pad + x

def outer_commas(x):
    i     = 0
    total = 0
    for c in x:
        if c == '(': i += 1
        if c == ')': i -= 1
        if (c == ',') and (i==1): total += 1
    return total

def replace_holders(x, reduction_vars, reduction_ops, reduction_buffers, pad):
    x     = x.lstrip(' ')
    for word,op,buf in zip(reduction_vars, reduction_ops, reduction_buffers):
        all_i  =  [i for i in range(len(x)) if is_word(x,word,i)]
        for i in reversed(all_i):
            assert x.count('=') == 1
            new            =  buf
            before, after  =  x[:i], x[(i+len(word)):]
            if '=' in before:
                if   (op == '+' ): new  =  '0'
                else:
                    new            =  ''
                    before, after  =  before.strip(' '), after.strip(' ')
                    if before[-1]==after[0]==',':
                        before = removesuffix(before, ',')
                    else:
                        assert (before[-1]+after[0]).count(',') == 1
                        before = removesuffix(before, ',')
                        after  = removeprefix(after, ',')
                    if outer_commas((before+after).split('=')[1]) == 0:
                        assert before.count(f'{op}(') == 1
                        before  = before.replace(f'{op}(', '(')
            x = before + new + after
    return ' '*pad + x

def code_fetch(reduction_vars,reduction_ops,scratch_array,npad):
    pad = ' '*npad
    A = []
    for name,op in zip(reduction_vars,reduction_ops):
        mode_op  =  {'+'  :  '0' ,
                      'max':  '1',
                      'min': '-1'}[op]
        A.append(f"{pad}call quick_binary_scalar_reduction({scratch_array}, {name}, {binred_numel}, {mode_op})")
    # A.append(f'{pad}end block')
    return A

# ------------------------------------------------------------------------
#                  Direct Runner
# ------------------------------------------------------------------------

if __name__ == '__main__':
    verbose     =  True
    macro_name  =  '_BIN_REDUCTION'
    all_fname   =  lfilter(lambda x: basename(x).endswith('.f90'), listdir_full_files_deep(dirname(_folder_)))
    scratch_array     = 'temp_arr_n_xzy'
    for fname in all_fname:
        A  =  reader(fname)
        if any(map(is_reduction_line,A)):
            A      =  clean_prev(A)
            #
            i      =  0
            while i<len(A):
                if is_reduction_line(A[i]):
                    all_reductions  =  []
                    j               =  i
                    # go through !$ lines, and find reductions (erasing the declarations as well)
                    while j<len(A) and does_start(A[j],'!$'):
                        j    +=  1
                    # j should be the start of the nested for-loop (or code will crash)
                    assert does_start(A[j], 'do ')
                    max_d, d, k   = None, 0, j
                    # find max depth of for-loop, as well as loop_vars
                    while (d>0) or (max_d is None):
                        d    += (does_start(A[k], 'do '))
                        d    -= (does_start(A[k], 'end do') or does_start(A[k], 'enddo'))
                        max_d = d if (max_d is None) else max(max_d, d)
                        k    += 1
                    # (k-1) is the end of the nested for-loop
                    k_end =  k-1
                    k_ini = i
                    while ((k_ini-1)>=0) and does_start(A[k_ini],'!$'): k_ini-= 1
                    k_ini += 1
                    k_ini = max(0,k_ini)
                    old_code = A[k_ini:k_end+1]
                    assert (d == 0) and (k_end<len(A))
                    def apply_get_reductions():
                        for j in range(k_ini,k_end+1):
                            if is_reduction_line(A[j]):
                                A[j]  =  get_reductions(A[j], all_reductions) 
                    apply_get_reductions()
                    # prepare variables:
                    reduction_ops, reduction_vars  =  [[s.split(':')[j].strip(' ') for line in all_reductions for s in line.split('(')[1].split(')')[0].split(',')] for j in range(2)]
                    assert len(reduction_vars) == len(set(reduction_vars)) == 1
                    binred_numel  =  'nx_glob*nz_loc*ny_loc'
                    loc_address = 'i,j,k'
                    reduction_buffers = [f"{scratch_array}({loc_address})" for _ in reduction_vars]
                    # loop over body replacing reudction variable(s)
                    for k in range(j,k_end+1):
                        if not ((does_start(A[k], 'do ')) or
                                (does_start(A[k], 'end do') or does_start(A[k], 'enddo'))):
                            A[k] = replace_holders(A[k], reduction_vars, reduction_ops, reduction_buffers,get_pad(A[k]))
                    # add binary reduction at the end
                    A[k_end] = '\n'.join([A[k_end],*code_fetch(reduction_vars,reduction_ops,scratch_array,get_pad(A[k_end]))])
                    # first line line of code before "!$ flags"
                    j = k_ini
                    pad = ' '*get_pad(A[j])
                    code_preamble = [  #f'{pad}block',
                                       #f'{pad}  integer :: ' + ', '.join(var_defs.keys()),
                                       f'{pad}call set_val_slice(temp_arr_n_xzy,0,-1,0.d0)',
                                     #*[f"{pad}  {key}= {val}" for key,val in var_defs.items()]
                                     ]
                    # add block at the start
                    A[j] = '\n'.join([*code_preamble,A[j]])
                    A[j] = '\n'.join([f'#if !defined({macro_name})',*old_code,'#else',A[j]])
                    A[k_end] += '\n#endif'
                    # i continues at k_end+1
                    i    = k_end + 1
                else:
                    i += 1
            writer(pjoin(dirname(fname),'redu_fixed',basename(fname)),ftn_text_file_check(A,lim=128))
