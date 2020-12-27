from itertools import chain

def ixfty1(start, stop, step, rst, rsp):
    """
    Factory for creating start, stop index generator functions.
    Start at 'start'; end at 'stop'; advance by 'step'.
    For 'sym' column functions.
    """
    return zip(
        range(start, stop, step),                       # start indices
        chain(range(start+step, stop, step), (stop,))   # stop indices
    )

def ixfty2(start, stop, step, rst, rsp):
    """
    Factory for creating start, stop index generator functions.
    Start at 'start'; end at 'rsp'; advance by 'step'.
    For 'tril' column functions.
    """
    return zip(
        range(start, rsp, step),                        # start indices
        chain(range(start+step, rsp, step), (rsp,))     # stop indices
    )

def ixfty3(start, stop, step, rst, rsp):
    """
    Factory for creating start, stop index generator functions.
    Start at 'rst'; end at 'stop'; advance by 'step'.
    For 'triu' column functions.
    """
    return zip(
        range(rst, stop, step),                         # start indices
        chain(range(rst+step, stop, step), (stop,))     # stop indices
    )
