# -*- coding: utf-8 -*-
"""
Created on Fri Nov 25 12:30:49 2016

@author: https://gist.github.com/djbarker/d52e1c2891b4ff2c3e17
"""

def lowbit0(n):
    """
    Returns the index of the least significant bit that is a 0 in the given int.
    """
    m = 1
    j = 1
    while m&(~n)==0:
        m *= 2
        j += 1
    return j

def sobolSeq(poly,init,bits=32):
    """
    Sets up a generator for producing a maximum of 2**bits Sobol quasi-random (sub-
    random) numbers.
      poly - Generating polynomical coeffieients, including 1 (x^0) but
             excluding x^(N+1). Element of {0,1}^N
      init - Initialising values fed to the recurrance relation. Length must
             match poly. Element of {2*k_i+1|i=1...N}
      bits - Accuracy of the sequence (number of bits).
    """

    # check input
    if len(poly)!=len(init):
        raise ValueError('Length of polynomial coeffs must match length of initializing values')

    for p in poly:
        if p!=0 and p!=1:
            raise ValueError('Polynomial coeffs must be zero or one.')

    for m in init:
        if m%2==0:
            raise ValueError('Initial values must be odd.')

    if bits<0:
        raise ValueError('Cannot have negative accuracy')

    # initialize direction numbers
    q = len(init)
    M = init
    V = [ M[i]*2**(bits-i-1) for i in range(q) ]
    A = poly

    for i in range(len(M),bits):
        m = M[-q]^((2**q)*M[-q])
        for j in range(1,q):
            m ^= (2**j)*M[-j]*A[j]
        M.append( m )
        V.append( m*2**(bits-i-1) )

    # generate
    n = 1
    S = 0
    while True:
        S = S ^ V[lowbit0(n)]
        n += 1
        yield float(S)/(2**bits)

# EXAMPLES

if __name__ == '__main__':

    import matplotlib.pyplot as plt
    import numpy as np
    from math import sqrt

    # Visualize a 2D Sobol sequence
    
    s1 = sobolSeq([1,1],[1,1],32)
    s2 = sobolSeq([1,0,1],[1,3,7],32)

    N = 100
    X =  [ next(s1) for _ in range(N) ]
    Y =  [ next(s2) for _ in range(N) ]

    plt.subplot('222')
    plt.scatter(X,Y,marker='.')
    plt.xlim(0,1)
    plt.ylim(0,1)
    plt.gca().set_aspect('equal')

    plt.subplot('221')
    plt.hist(Y,bins=sqrt(N)//2,orientation='horizontal',normed=True)

    plt.subplot('224')
    plt.hist(X,bins=sqrt(N)//2,normed=True)

    plt.show()

    # 2D Numerical integration example

    def f(x,y):
        return -(x**2)*(y**2)

    F_ana = -1./9. # analytical integral over the unit square

    s1 = sobolSeq([1,1],[1,1],32)
    s2 = sobolSeq([1,0,1],[1,3,7],32)

    N = 100
    X =  [ next(s1) for _ in range(N) ]
    Y =  [ next(s2) for _ in range(N) ]

    # pseudo-random sequences
    X2 = np.random.uniform(size=(N,))
    Y2 = np.random.uniform(size=(N,))

    errSobol = []
    errMC    = []
    accSobol = 0
    accMC    = 0

    idxs = []
    ni=0
    for i in range(N):
        accSobol += f(X[i], Y[i])
        accMC    += f(X2[i],Y2[i])

        if i==ni: 
            errSobol.append(abs( accSobol / (i+1) - F_ana ))
            errMC.append(abs( accMC / (i+1) - F_ana ))
            ni += int(np.log10(i+1)+1)
            idxs.append( i )
            
    n = np.array(range(N),dtype=np.float_)+1
    L = 2*(n)**-1
    S = 2*(n)**-0.5
    plt.loglog(idxs,errSobol)
    plt.loglog(idxs,errMC)
    plt.loglog(n,L,'k--')
    plt.loglog(n,S,'k--')
    plt.ylim(min(L),max(L))
    plt.show()