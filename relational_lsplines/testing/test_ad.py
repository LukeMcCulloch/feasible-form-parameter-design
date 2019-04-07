##
##  Unit testing AD IA
##
## Luke McCulloch
## September 25 2015
##
"""
cd Documents/computational_naval_architecture/alt_implementation/
"""

import relational_lsplines as rlspline 
#from relational_lsplines import simple_hull_rules_language as srl

#from relational_lsplines.simple_hull_rules_language import HullGeometryGenerator


import numpy as np
import unittest#2 as unittest
import copy
from timer import Timer

ia = rlspline.ia
ad = rlspline.ad
from AD import scalarAD 
from AF import fuzzyNumber

##
## DATA...
##
N=3
Imatrix = np.identity(N)
gx = np.matrix([1.,0.,0.])
gy = np.matrix([0.,1.,0.])
gz = np.matrix([0.,0.,1.])
H =  np.matrix( np.zeros((N,N),float) )

xnew  = ad(3.0, gx, H)
ynew  = ad(10., gy, H)
znew  = ad(4.0, gz, H)

x  = scalarAD(3.0,   gx, H)
y  = scalarAD(10.,   gy, H)
z  = scalarAD(4.0,   gz, H)


xsup = 1.
xinf = 0.
ysup = 5.
yinf = 3.
zsup = 11.
zinf = 9.


gX = np.matrix([ia(1.,1.),ia(0.,0.),ia(0.,0.)])
gY = np.matrix([ia(0.,0.),ia(1.,1.),ia(0.,0.)])
gZ = np.matrix([ia(0.,0.),ia(0.,0.),ia(1.,1.)])

HX = np.matrix([[ia(0.,0.),ia(0.,0.),ia(0.,0.)],
                 [ia(0.,0.),ia(0.,0.),ia(0.,0.)],
                 [ia(0.,0.),ia(0.,0.),ia(0.,0.)]])
                 
HY = copy.deepcopy(HX)
HZ = copy.deepcopy(HX)

xvnew = np.asarray([
    ad( ia(xinf, xsup), gX, HX ),
    ad( ia(yinf, ysup), gY, HY),
    ad( ia(zinf, zsup), gZ, HZ)])  



N=3
Null = np.zeros((N,N),float)
Iarray = np.identity(N)
x1min = scalarAD(1.,np.matrix(Iarray[0]),np.matrix(Null))
x1mid = scalarAD(1.5,np.matrix(Iarray[0]),np.matrix(Null))
x1max = scalarAD(2.,np.matrix(Iarray[0]),np.matrix(Null))
x1 = fuzzyNumber(x1min,x1mid,x1max)

x2min = scalarAD(2.,np.matrix(Iarray[1]),np.matrix(Null))
x2mid = scalarAD(2.5,np.matrix(Iarray[1]),np.matrix(Null))
x2max = scalarAD(3.,np.matrix(Iarray[1]),np.matrix(Null))
x2 = fuzzyNumber(x2min,x2mid,x2max)

x3min = scalarAD(1.1,np.matrix(Iarray[2]),np.matrix(Null))
x3mid = scalarAD(2.1,np.matrix(Iarray[2]),np.matrix(Null))
x3max = scalarAD(3.1,np.matrix(Iarray[2]),np.matrix(Null))
x3 = fuzzyNumber(x3min,x3mid,x3max)


cmin = scalarAD(2.,np.matrix(Null[0]),np.matrix(Null))
cmid = scalarAD(2.5,np.matrix(Null[0]),np.matrix(Null))
cmax = scalarAD(3.,np.matrix(Null[0]),np.matrix(Null))
c = fuzzyNumber(cmin,cmid,cmax)

V = np.asarray([x1,x2,x3])

Imatrix = np.identity(N)
a = 1.
b = 1.5
c = 2.
mina    = scalarAD(  a,np.matrix(Imatrix[0]),np.matrix(np.zeros((2,2),float)))
reala   = scalarAD(  b,np.matrix(Imatrix[0]),np.matrix(np.zeros((2,2),float)))
maxa    = scalarAD(  c,np.matrix(Imatrix[0]),np.matrix(np.zeros((2,2),float)))
x1 = fuzzyNumber(mina, reala, maxa)

a = 0.
b = .5
c = 1.
mina    = scalarAD(  a,np.matrix(Imatrix[1]),np.matrix(np.zeros((2,2),float)))
reala   = scalarAD(  b,np.matrix(Imatrix[1]),np.matrix(np.zeros((2,2),float)))
maxa    = scalarAD(  c,np.matrix(Imatrix[1]),np.matrix(np.zeros((2,2),float)))
x2 = fuzzyNumber(mina, reala, maxa)

X = [x1,x2]

x = copy.deepcopy(X[0])
y = copy.deepcopy(X[1])
x.min.value = 0.
x.real.value = .5
x.max.value = 1.
y.min.value = 3.
y.real.value = 4.
y.max.value = 5.
A = np.asarray([x,y])
with Timer() as t:
    r = np.dot(A,A)
print "=> elasped time for dot product: %s s" % t.secs



xv1 = xvnew+xvnew
print 'Time analysis'
with Timer() as t:
    r = np.dot(xvnew,xvnew)
print "=> elasped time for dot product: %s s" % t.secs

#not possible in the old system
#xv = np.asarray([
#    AD.scalarAD( ia(0.,1.), np.matrix([ia(1.,1.),ia(0.,0.)]), np.matrix([[ia(0.,0.),ia(0.,0.)],
#                                                               [ia(0.,0.),ia(0.,0.)]]) ),
#    AD.scalarAD( ia(3.,5), np.matrix([ia(0.,0.),ia(1.,1.)]), np.matrix([[ia(0.,0.),ia(0.,0.)],
#                                                               [ia(0.,0.),ia(0.,0.)]]) )
#]) 



def verify_adia_w_AFAD(old, new):
    assert(isinstance(old, AF.fuzzyNumber))
    assert(isinstance(old.min, AD.scalarAD))
    assert(isinstance(old.max, AD.scalarAD))
    assert(isinstance(new, ad))
    assert(isinstance(new, ad))
    assert(isinstance(new.value, ia))
    assert(isinstance(new.value, ia))
    
    
    

class AdditionTests(unittest.TestCase):
    cases = ((x+y, xnew+ynew),
             (x+y+z, xnew+ynew+znew)
             )
             
    def testOne(self):
        for old, new in self.cases:
            self.failUnlessEqual(old, new)



class MultiplictionTests(unittest.TestCase):
    cases = ((x*y, xnew*ynew),
             (x*y*z, xnew*ynew*znew)
             )
             
    def testOne(self):
        for old, new in self.cases:
            self.failUnlessEqual(old, new)

#def main():
#    unittest.main()

if __name__ == '__main__':
    main()