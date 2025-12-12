import numpy as np

def addLog(x, y):
    """
    Computes log(exp(x) + exp(y)) safely.
    equivalent to numpy.logaddexp(x, y)
    """
    return np.logaddexp(x, y)

def addLog2(x, y):
    """
    Computes log2(2**x + 2**y) safely.
    equivalent to numpy.logaddexp2(x, y)
    """
    return np.logaddexp2(x, y)

def sumLog2(arr):
    """
    Computes log2(sum(2**x)) for an array x.
    """
    # np.logaddexp2.reduce is not available in all versions, but let's check or implement manually
    # Actually logaddexp2 ufunc has a reduce method.
    return np.logaddexp2.reduce(arr)

def log2(x):
    return np.log2(x)

def exp2(x):
    return np.exp2(x)
