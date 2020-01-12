#!/usr/bin/env python
#
#
#
import argparse 
import logging

import numpy as np
from scipy import sparse

#              A B C D E F        A B C D E F
#            A 0 1 0 0 0 0      A 0 1 1 1 1 1
#            B 0 0 1 1 0 0      B 0 0 1 1 1 1
#            C 0 0 0 0 1 0  --> C 0 0 0 0 1 1 
#            D 0 0 0 0 1 0      D 0 0 0 0 1 1
#            E 0 0 0 0 0 1      E 0 0 0 0 0 1 
#            F 0 0 0 0 0 0      F 0 0 0 0 0 0
#


DATA = [ [ 0, 1, 0, 0, 0, 0 ],
         [ 0, 0, 1, 1, 0, 0 ],
         [ 0, 0, 0, 0, 1, 0 ],
         [ 0, 0, 0, 0, 1, 0 ],
         [ 0, 0, 0, 0, 0, 1 ],
         [ 0, 0, 0, 0, 0, 0 ] ]


def testmatrix():
    gomatrix = np.array(DATA)
    #gomatrix = np.zeros( shape, dtype=np.int8 )
    #gomatrix = np.full( shape, False, dtype=bool    )
    #gomatrix = sparse.csr_matrix( shape, dtype=bool )
    logging.debug(f"filling in parent matrix for all goterms...")
    #for gt in godict.keys():
    #    for parent in godict[gt]['is_a']:
    #            gomatrix[termidx[parent]][termidx[gt]] = 1     
    logging.debug(f"starting matrix: \n{gomatrix}")
    logging.debug("Calculating sparsity...")
    sparsity = 1.0 - np.count_nonzero(gomatrix) / gomatrix.size
    logging.debug(f"sparsity = {sparsity}")
    # print(gomatrix)
    logging.debug("converting to sparse matrix.")
    gomatrix = sparse.lil_matrix(gomatrix, dtype=bool)
    #gomatrix = sparse.lil_matrix(gomatrix, dtype=bool)
    #logging.debug(f" gomatrix= {gomatrix}")
    logging.debug("converging matrix...")
    gomatrix = converge_matrix(gomatrix)
    logging.debug("convert to dense matrix")
    gomatrix = gomatrix.todense()
    logging.debug(f"got converged matrix:\n{gomatrix}")
    logging.debug("Calculating sparsity...")
    sparsity = 1.0 - np.count_nonzero(gomatrix) / gomatrix.size
    logging.debug(f"sparsity = {sparsity}")        
    #for k in dict.keys():
    #    po = dict[k]['part_of']
    #    print(po)
    #logging.debug(f"got dict: {dict}")

def converge_matrix(mat):
    oldval = 0
    newval = np.sum(mat)
    logging.debug(f"initial sum is {newval}")    
    
    logging.debug("multiplying matrix by itself...")
    #mat2 = np.matmul(mat, mat)
    #mat2 = mat.multiply(mat)
    mat2 = mat @ mat
    logging.debug("adding back original matrix...")
    mat2 = mat2 + mat
    #mat2 = np.where(mat2 > 0, 1, 0)
    logging.debug("calculating matrix sum...")
    
    while newval != oldval:
        logging.debug(f"newval {newval} != oldval {oldval}")
        oldval = newval
        #mat2 = np.matmul(mat, mat)
        #mat2 = mat.multiply(mat)
        mat2 = mat @ mat
        mat2 = mat2 + mat
        #mat2 = np.where(mat2 > 0, 1, 0)
        newval = np.sum(mat2)
        logging.debug(f" newval={newval} oldval={oldval}")
    logging.debug(f"done. values converged with newval {newval}")   
    return mat2
    


def step(mat):
    # Multiply the matrices together to get only the
    # relationships that propagate.
    mat2 = numpy.matmul(mat, mat)
    # Add the new relationships to the original
    # state of the matrix
    mat2 = mat + mat2
    # Flatten values > 0 to 1
    mat2 = np.where(mat2 > 0, 1, 0)
    return(mat2)

if __name__ == '__main__':
    FORMAT='%(asctime)s (UTC) [ %(levelname)s ] %(filename)s:%(lineno)d %(name)s.%(funcName)s(): %(message)s'
    logging.basicConfig(format=FORMAT)
    
    parser = argparse.ArgumentParser()
      
    parser.add_argument('-d', '--debug', 
                        action="store_true", 
                        dest='debug', 
                        help='debug logging')

    parser.add_argument('-v', '--verbose', 
                        action="store_true", 
                        dest='verbose', 
                        help='verbose logging')
          
    args= parser.parse_args()
    
    if args.debug:
        logging.getLogger().setLevel(logging.DEBUG)
    if args.verbose:
        logging.getLogger().setLevel(logging.INFO)
    testmatrix()