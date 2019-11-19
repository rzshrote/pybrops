import numpy

# import meiosis simulating engine
import pybropt.sim.meiosis as meiosis

# generate binary marker data
n_lines = 20
n_markers = 100
n_phases = 2
n_selected = 10

def print_bmat(mat):
    for row in range(mat.shape[0]):
        for col in range(mat.shape[1]):
            print(mat[row,col], sep='', end='')
        print('')
    print('')

# seed random number
numpy.random.seed(116019)
# make markers
markers = numpy.random.binomial(1, 0.3, (n_phases,n_lines,n_markers))
# add 4 to the second phase to easily see introgressions
markers[1,:,:] += 4

print("Phase 1")
print_bmat(markers[0])
print("Phase 2")
print_bmat(markers[1])

# generate a linkage map (2 chromosomes; 50 markers each)
d = numpy.empty(n_markers, dtype=numpy.float64)
for b in range(2):
    d[(b*50):(b*50)+50] = numpy.sort(numpy.random.uniform(0.0, 4.0, 50))
d_size = numpy.array([50,50])
sources = numpy.arange(n_lines)

gout = meiosis.meiosis(markers, d, d_size, sources, verbose=False)

print("Gametes")
print_bmat(gout)
