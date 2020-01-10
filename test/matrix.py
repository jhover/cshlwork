import numpy






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