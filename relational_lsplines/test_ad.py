##
##  Unit testing AD IA
##
## Luke McCulloch
## September 25 2015
##
"""
cd Documents/computational_naval_architecture/alt_implementation/
"""
import numpy as np
import unittest#2 as unittest
import copy
from time_analysis import Timer
from interval_analysis import ia 
import automatic_differentiation as ad
from AD import scalarAD 
from AF import fuzzyNumber

##
## DATA...
##






#not possible in the old system
#xv = np.asarray([
#    AD.scalarAD( ia(0.,1.), np.matrix([ia(1.,1.),ia(0.,0.)]), np.matrix([[ia(0.,0.),ia(0.,0.)],
#                                                               [ia(0.,0.),ia(0.,0.)]]) ),
#    AD.scalarAD( ia(3.,5), np.matrix([ia(0.,0.),ia(1.,1.)]), np.matrix([[ia(0.,0.),ia(0.,0.)],
#                                                               [ia(0.,0.),ia(0.,0.)]]) )
#]) 




def verify_adia_w_AFAD(old, new):
    assert(isinstance(old, fuzzyNumber))
    assert(isinstance(old.min, scalarAD))
    assert(isinstance(old.max, scalarAD))
    assert(isinstance(new, ad))
    assert(isinstance(new, ad))
    assert(isinstance(new.value, ia))
    assert(isinstance(new.value, ia))
    
    
    

class AdditionTests(unittest.TestCase):
    N=3
    Imatrix = np.identity(N)
    gx = ad.adObjectMaker.makeGradient(N,0)#np.matrix([1.,0.,0.])
    gy = ad.adObjectMaker.makeGradient(N,1)#np.matrix([0.,1.,0.])
    gz = ad.adObjectMaker.makeGradient(N,2)#np.matrix([0.,0.,1.])
    H =  ad.adObjectMaker.makeHessian(N)#np.matrix( np.zeros((N,N),float) )
    
    xnew  = ad.ad(3.0, gx, H)
    ynew  = ad.ad(10., gy, H)
    znew  = ad.ad(4.0, gz, H)
    
    x  = scalarAD(3.0,   gx, H)
    y  = scalarAD(10.,   gy, H)
    z  = scalarAD(4.0,   gz, H)
    
    cases = ((x+y, xnew+ynew),
             (x+y+z, xnew+ynew+znew))
    def testOne(self):
        for old, new in self.cases:
            self.failUnlessEqual(old, new)



class MultiplictionTests(unittest.TestCase):
    N=3
    Imatrix = np.identity(N)
    gx = np.matrix([1.,0.,0.])
    gy = np.matrix([0.,1.,0.])
    gz = np.matrix([0.,0.,1.])
    H =  np.matrix( np.zeros((N,N),float) )
    
    xnew  = ad.ad(3.0, gx, H)
    ynew  = ad.ad(10., gy, H)
    znew  = ad.ad(4.0, gz, H)
    
    x  = scalarAD(3.0,   gx, H)
    y  = scalarAD(10.,   gy, H)
    z  = scalarAD(4.0,   gz, H)
    
    cases = ((x*y, xnew*ynew),
                 (x*y*z, xnew*ynew*znew))
    def testOne(self):
        for old, new in self.cases:
            self.failUnlessEqual(old, new)



class matrixMultiplicationTests(unittest.TestCase):
    
    adomg = ad.adObjectMaker.makeGradient
    adomh = ad.adObjectMaker.makeHessian
    
    
    N=3
    xsup = 1.
    xinf = 0.
    ysup = 5.
    yinf = 3.
    zsup = 11.
    zinf = 9.
    
    
    gX = adomg(N,0),#np.matrix([ia(1.,1.),ia(0.,0.),ia(0.,0.)])
    gY = adomg(N,1)#np.matrix([ia(0.,0.),ia(1.,1.),ia(0.,0.)])
    gZ = adomg(N,2)#np.matrix([ia(0.,0.),ia(0.,0.),ia(1.,1.)])
    
#    HX = np.matrix([[ia(0.,0.),ia(0.,0.),ia(0.,0.)],
#                     [ia(0.,0.),ia(0.,0.),ia(0.,0.)],
#                     [ia(0.,0.),ia(0.,0.),ia(0.,0.)]])
    
#    ia_null = ia(0.,0.)
#    
#    m1 = np.matrix(np.zeros((N,N),float))
#    HX=[]
#    for i in range(N):
#        HX.append([])
#        for j in range(N):
#            HX[i].append(ia_null*m1[i,j])
#    HX = np.matrix(HX)
#    
#    HY = copy.deepcopy(HX)
#    HZ = copy.deepcopy(HX)
    
    HX = adomh(N)
    HY = adomh(N)
    HZ = adomh(N)
    
    xvnew = np.asarray([
        ad.ad( ia(xinf, xsup), gX, HX ),
        ad.ad( ia(yinf, ysup), gY, HY),
        ad.ad( ia(zinf, zsup), gZ, HZ)]) 
    
    

    from AD import scalarAD
    from AF import fuzzyNumber
    Imatrix = np.identity(N)
    a = 0.
    b = .5
    c = 1.
    mina    = scalarAD(  a,np.matrix(Imatrix[0]),np.matrix(np.zeros((N,N),float)))
    reala   = scalarAD(  b,np.matrix(Imatrix[0]),np.matrix(np.zeros((N,N),float)))
    maxa    = scalarAD(  c,np.matrix(Imatrix[0]),np.matrix(np.zeros((N,N),float)))
    x1 = fuzzyNumber(mina, reala, maxa)
    
    a = 3.
    b = 4.
    c = 5.
    mina    = scalarAD(  a,np.matrix(Imatrix[1]),np.matrix(np.zeros((N,N),float)))
    reala   = scalarAD(  b,np.matrix(Imatrix[1]),np.matrix(np.zeros((N,N),float)))
    maxa    = scalarAD(  c,np.matrix(Imatrix[1]),np.matrix(np.zeros((N,N),float)))
    x2 = fuzzyNumber(mina, reala, maxa)
    
    a = 9.
    b = 10.
    c = 11.
    mina    = scalarAD(  a,np.matrix(Imatrix[2]),np.matrix(np.zeros((N,N),float)))
    reala   = scalarAD(  b,np.matrix(Imatrix[2]),np.matrix(np.zeros((N,N),float)))
    maxa    = scalarAD(  c,np.matrix(Imatrix[2]),np.matrix(np.zeros((N,N),float)))
    x3 = fuzzyNumber(mina, reala, maxa)
    
    X = [x1,x2,x3]
    

    with Timer() as t:
        r = np.dot(X,X)
    print "=> elasped time for dot product: %s s" % t.secs
    
    
    
    xv1 = xvnew+xvnew
    print 'Time analysis'
    with Timer() as t:
        r = np.dot(xvnew,xvnew)
    print "=> elasped time for dot product: %s s" % t.secs
    
    

def main():
    unittest.main()

if __name__ == '__main__':
    main()