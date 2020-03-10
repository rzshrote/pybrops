import numpy

def matrix_is_sorted(mat):
    for i in range(mat.size-1):
        if mat[i] > mat[i+1]:
            return False
    return True
