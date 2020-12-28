class GenomicModelFactory:
    def __init__(self, c, **kwargs):
        self.c = c
        self.calibrate = c
    def calibrate(self, **kwargs):
        raise NotImplementedError("method is abstract")

class A:
    def __init__(self, **kwargs):
        self.data = 'A'

class B:
    def __init__(self, **kwargs):
        self.data = 'B'

f = GenomicModelFactory(A)
g = f.calibrate()

type(g)
