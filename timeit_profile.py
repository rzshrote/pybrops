import numpy
import timeit

n_markers = 10000
vector1 = numpy.arange(n_markers)
vector2 = numpy.arange(n_markers)
exec1 = 'vector1[:,None] @ vector1[None,:]'
exec2 = """\
vector2.shape = (len(vector2),1)
vector2 @ vector2.T
"""
exec1_list = list()
exec2_list = list()
for i in range(25):
    exec1_list.append(timeit.timeit(exec1, number=100, globals=globals()))
    exec2_list.append(timeit.timeit(exec2, number=100, globals=globals()))
with open("matrix_mult2.txt", "w") as handle:
    for a,b in zip(exec1_list, exec2_list):
        handle.write("%s\t%s\n" % (a,b))
