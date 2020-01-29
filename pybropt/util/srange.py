# define a range function that included the stop variable too.

def srange(start, stop, step):
    yield from range(start, stop, step)
    yield stop
