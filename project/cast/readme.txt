Barseq2 expression data for CAST
---------------------------------------
file: barseq2.rds
to access: 
R 
> require(SingleCellExperiment)
> barseq2 = readRDS('barseq2.rds')
> # sparse expression matrix
> X = SingleCellExperiment::assay(barseq2)  

file: barseq2.hdf5
R 
> require(rhdf5)
> require(Matrix)
> i = as.numeric(rhdf5::h5read(file = 'barseq2.hdf5', name ='matrix/i') )
> j = as.numeric(rhdf5::h5read(file = 'barseq2.hdf5', name ='matrix/j'))
> x = as.numeric(rhdf5::h5read(file = 'barseq2.hdf5', name ='matrix/x'))
> dim = rhdf5::h5read(file = 'barseq2.hdf5', name ='dim')
> genes = rhdf5::h5read(file = 'barseq2.hdf5', name ='dimnames/genes')
> cells = rhdf5::h5read(file = 'barseq2.hdf5', name ='dimnames/cells')
> M = Matrix::sparseMatrix(i = i, j = j, x=x,dims = dim, dimnames = list(genes,cells ) )

python3
>>> import h5py
>>> from scipy import sparse
>>>
>>> with h5py.File('barseq2.hdf5') as f:
>>>     dim = f['dim'][()]
>>>     genes = [g.decode('utf-8') for g in f['dimnames/genes'][()]]
>>>     cells = [c.decode('utf-8') for c in f['dimnames/cells'][()]]
>>>     i = f['matrix/i'][()] - 1
>>>     j = f['matrix/j'][()] - 1  
>>>     x = f['matrix/x'][()]
>>>     M = sparse.csc_matrix((x, (i,j)) , shape =dim )
