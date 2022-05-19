
import itertools
import logging
import os 
import subprocess
import sys
import tempfile


import pandas as pd

class NonZeroReturnException(Exception):
    """
    Thrown when a command has non-zero return code. 
    """

def readlist(filepath):
    '''
    Assumes file is a list of strings, one per line. 
    Ignores lines beginning with a has '#'
    Ignores characters in a line afeter a '#'
    '''

    if filepath is not None:
        logging.info(f'reading file: {filepath}')
        flist = []
        try:
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    if len(line) > 0:
                        idx = line.find('#')
                        if idx == -1:
                            flist.append(line.strip())
                        elif idx > 0:
                            flist.append(line[:idx].strip())
                    else:
                        pass   # empty line
                        
            logging.debug(f'got list with {len(flist)} items.')
            return flist
        except:
            return []
    else:
        logging.info('no file. return [].')
        return []


def writelist(filepath, dlist, mode=0o644):
    logging.info(f"writing list length={len(dlist)} to file='{filepath}'")
    rootpath = os.path.dirname(filepath)
    basename = os.path.basename(filepath)
    try:
        (tfd, tfname) = tempfile.mkstemp(suffix=None,
                                         prefix=f"{basename}.",
                                         dir=f"{rootpath}/",
                                         text=True)
        logging.debug(f"made temp {tfname}")
        with os.fdopen(tfd, 'w') as f:
            nlines = 0
            for item in dlist:
                f.write(f"{item}\n")
                nlines += 1
        os.rename(tfname, filepath)
        os.chmod(filepath, mode)
        logging.info(f"wrote {nlines} to {filepath}")
    except Exception as ex:
        logging.error(traceback.format_exc(None))

    finally:
        pass



def load_df(filepath):
    """
    Convenience method to load DF consistently accross modules. 
    """
    filepath = os.path.expanduser(filepath)
    df = pd.read_csv(filepath, sep='\t',index_col=0, keep_default_na=False, dtype =str, comment="#")
    df.fillna(value='', inplace=True)
    df = df.astype('str', copy=False)
    return df


def add_rowlist_column(rowlist, colval):
    """
    For use during dataframe construction. Adds col to list of rows with specified.
       
    """
    for row in rowlist:
        row.append(colval)
    return rowlist
    

def merge_write_df(newdf, filepath,  mode=0o644):
    """
    Reads existing, merges new, drops duplicates, writes to temp, renames temp. 
    """
    log = logging.getLogger('utils')
    log.debug(f'inbound new df:\n{newdf}')
    filepath = os.path.abspath( os.path.expanduser(filepath) )
    if os.path.isfile(filepath):
        df = load_df(filepath)
        log.debug(f'read df:\n{df}')
        df = df.append(newdf, ignore_index=True)
        df.fillna(value='', inplace=True)
        df = df.astype('str', copy=False)
        log.debug(f'appended df:\n{df}')
    else:
        df = newdf
        df.fillna(value='', inplace=True)
        df = df.astype('str', copy=False)
    logging.debug(f"df length before dropping dupes is {len(df)}")
    df.drop_duplicates(inplace=True, ignore_index=True, keep='first')
    logging.debug(f"df length after dropping dupes is {len(df)}")
    df = df.reset_index(drop=True)
    rootpath = os.path.dirname(filepath)
    basename = os.path.basename(filepath)
    try:
        (tfd, tfname) = tempfile.mkstemp(suffix=None,
                                         prefix=f"{basename}.",
                                         dir=f"{rootpath}/",
                                         text=True)
        logging.debug(f"made temp {tfname}")
        df.to_csv(tfname, sep='\t')
        os.rename(tfname, filepath)
        os.chmod(filepath, mode)
        logging.info(f"wrote df to {filepath}")

    except Exception as ex:
        logging.error(traceback.format_exc(None))
        raise ex


def matrix2table(df, symmetric=True, merge_label=True, combo_char='x'):
    '''
    Takes  A x B  square matrix of permutations and converts to table of AxB with values.
    Symmetric means A-> B is same as B->A and need appear only once. 
    output is alphabetized on 'comb' column. 
    merge_labels means produce 'comb' column with combo_chr between label names.  

        A   B   C
    A   .8  .2  .1
    B  .2   .6  .3
    C  .1   .3  .5

    first  second  comb   val    
    A      A       AxA    .8
    A      B       AxB    .2
    A      C       AxC    .1
    B      B       BxB    .6
    B      C       BxC    .3
    C      C       CxC    .5
    
    '''
    #tdf = df.stack().reset_index()
    #tdf.columns = ['first','second','val']
    #tdf['combo'] = tdf['first'].astype(str) + combo_char + tdf['second'].astype(str)
    #items = list(set(list(tdf['first'])))
    items = list(df.columns)
    items.sort()

    lol = []
    
    for (f,s) in itertools.combinations_with_replacement(items, 2):
        r = [ f'{f}{combo_char}{s}', df.at[f,s] ]
        logging.debug(f'row is {r}')
        lol.append( r )

    tdf = pd.DataFrame(lol, columns=['pair','val'])
    return tdf




def listdiff(list1, list2):
    logging.debug(f"got list1: {list1} list2: {list2}")
    s1 = set(list1)
    s2 = set(list2)
    sd = s1 - s2
    dl = list(sd)
    dl.sort()
    logging.debug(f"diff has length {len(dl)}")
    return dl


def listmerge(list1, list2):
    logging.debug(f"got list1: {list1} list2: {list2}")
    s1 = set(list1)
    s2 = set(list2)
    sd = s1 | s2
    dl = list(sd)
    dl.sort()
    logging.debug(f"merged has length {len(dl)}")
    return dl


def run_command(cmd):
    """
    cmd should be standard list of tokens...  ['cmd','arg1','arg2'] with cmd on shell PATH.
    
    """
    cmdstr = " ".join(cmd)
    logging.info(f"command: {cmdstr} running...")
    start = dt.datetime.now()
    cp = subprocess.run(cmd, 
                    text=True, 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.STDOUT)
    end = dt.datetime.now()
    elapsed =  end - start
    logging.debug(f"ran cmd='{cmdstr}' return={cp.returncode} {elapsed.seconds} seconds.")
    
    if cp.stderr is not None:
        logging.debug(f"got stderr: {cp.stderr}")
    if cp.stdout is not None:
        logging.debug(f"got stdout: {cp.stdout}")   
    if str(cp.returncode) == '0':
        logging.info(f'successfully ran {cmdstr}')
        return(cp.stderr, cp.stdout, cp.returncode)
    else:
        logging.warn(f'non-zero return code for cmd {cmdstr}')
        raise NonZeroReturnException()


def run_command_shell(cmd):
    """
    maybe subprocess.run(" ".join(cmd), shell=True)
    cmd should be standard list of tokens...  ['cmd','arg1','arg2'] with cmd on shell PATH.
    
    """
    cmdstr = " ".join(cmd)
    #logging.info(f"running command: {cmdstr} ")
    start = dt.datetime.now()
    cp = subprocess.run(" ".join(cmd), 
                    shell=True, 
                    stdout=subprocess.PIPE, 
                    stderr=subprocess.STDOUT)
    #cp = subprocess.run(cmd, 
    #                shell=True, 
    #                stdout=subprocess.PIPE, 
    #                stderr=subprocess.STDOUT)
    end = dt.datetime.now()
    elapsed =  end - start
    #logging.debug(f"ran cmd='{cmdstr}' return={cp.returncode} {elapsed.seconds} seconds.")
    
    if cp.stderr is not None:
        #logging.debug(f"got stderr: {cp.stderr}")
        pass
    if cp.stdout is not None:
        #logging.debug(f"got stdout: {cp.stdout}")
        pass
    if str(cp.returncode) == '0':
        #logging.info(f'successfully ran {cmdstr}')
        return(cp.stderr, cp.stdout, cp.returncode)
    else:
        logging.error(f'non-zero return code for cmd {cmdstr}')
        raise NonZeroReturnException(f'For cmd {cmdstr}')
