import numpy as np

def vangle(a, b, amagn=None, bmagn=None):
    """
    Calculates angle(s) between two vectors in arrays a and b.
    The shape of both is (n,3), n is arbitrary number of 3d vectors.
    """
    sha = a.shape
    shb = b.shape
    if len(sha) == 1: # a and b are just 3d vectors, not arrays of vectors
        if amagn == None: amagn = np.sqrt(dot(a,a))
        if bmagn == None: bmagn = np.sqrt(dot(b,b))
        return np.arccos(np.dot(a,b)/(amagn*bmagn))
    else:
        n = sha[0]
        theta = np.empty(n, dtype=np.double)
    
    for i in xrange(n):
        if amagn == None: amagn = np.sqrt(np.dot(a[i,:],a[i,:]))
        if bmagn == None: bmagn = np.sqrt(np.dot(b[i,:],b[i,:]))
        theta[i] = np.arccos(np.dot(a[i,:],b[i,:])/(amagn*bmagn))

    return theta
