#!/usr/bin/env python

import logging


def get_default_config():
    cp = ConfigParser()
    cp.read(os.path.expanduser("~/git/cshlwork/etc/phmmer.conf"))
    return cp


def execute_phmmer(config, queryfile, database=None):
    """
    cpus = 8
    eval_threshold = 1.0e-120
    database=~/data/uniprot/uniprot_sprot.fasta
    remove_self_hits = True
    
    *Excludes geneids ( <protein>_<species> ) of sample from phmmer hit results*
    So the inbound sequence files *must* contain correct geneids (or some other
    method must be used to exclude self-hits). 

    """
    logging.debug(f"executing with filename={filename} version={version} ")
  
    outdir = os.path.expanduser(config.get('phmmer','cachedir'))
    filename =os.path.expanduser(filename)
    outpath = os.path.dirname(filename)
    filebase = os.path.splitext(os.path.basename(filename))[0]
    outfile = "%s/%s.phmmer.tbl.txt" % (outdir, filebase)
    logging.debug(f"outfile={outfile}")
    cpus = config.get('phmmer','cpus')
    eval_threshold = config.get('phmmer','eval_threshold')
    database = os.path.expanduser(config.get('phmmer', 'database'))
    logging.debug(f"Using non-current version of uniprot for phmmer database: {database}")
    
    cmd = ["/usr/bin/which","phmmer"]
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    res = p.stdout.readlines()
    logging.debug(f'which phmmer gave {res}')

    cmd = [ 'phmmer',
           '--tblout', outfile , 
           '--noali',
           '--cpu', cpus,
           '-E', eval_threshold,
           filename,
           database 
           ]
    
    logging.debug(f"Running: {cmd}")
    cp = subprocess.run(cmd)
   
    logging.debug(f"Ran cmd='{cmd}' outfile={outfile} returncode={cp.returncode} " )
    logging.debug(f"returning outfile={outfile} exclude+list={exclude_list}")
    
    return (outfile, exclude_list, cidgidmap)




if __name__ == '__main__':
    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.WARN)
    
    parser = argparse.ArgumentParser()
      
    parser.add_argument('-d', '--debug', 
                        action="store_true", 
                        dest='debug', 
                        help='debug logging')

    parser.add_argument('-v', '--verbose', 
                        action="store_true", 
                        dest='verbose', 
                        help='verbose logging')
    
    parser.add_argument('infile', 
                        metavar='infile', 
                        type=str, 
                        help='Uniprot .dat file') 


    
    
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)   
    
    c = get_default_config()
    execute_phmmer(config, )
