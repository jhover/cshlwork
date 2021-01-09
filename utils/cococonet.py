import logging
import h5py
import pandas as pd

def read_network_hdf5(filename):
    """
    Loads data in file to dataframe.
    """
    with h5py.File(filename, 'r') as f:
        logging.debug("reading matrix...")
        matrix = f['agg'][:]
        logging.debug("reading rows. converting to unicode.")
        rows = [ s.decode() for s in  f['row'][:] ]
        logging.debug("reading columns. converting to unicode")
        columns = [ s.decode() for s in  f['col'][:] ]
        logging.debug("making dataframe...")
        df = pd.DataFrame(matrix,  index=rows, columns = columns )
        logging.debug(f"network shape: {df.shape}")
        logging.debug(f"network:\n {df}")        
    return df        