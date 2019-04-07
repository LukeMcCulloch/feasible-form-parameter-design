##
## Interval Analysis testing
##
## Luke McCulloch
## Spetember 25 2015
##


import relational_lsplines as rlspline 
#from relational_lsplines import simple_hull_rules_language as srl

#from relational_lsplines.simple_hull_rules_language import HullGeometryGenerator


import numpy as np
import unittest#2 as unittest
#import copy
from timer import Timer

#my interval library
ia = rlspline.ia
#my auto diff library
ad = rlspline.ad

import interval #test against another library


    
    
class TestArithmeticMethods(unittest.TestCase):

    
    def sp_equal(self, mine,their):
        mu,ms = mine[0],mine[1]
        tu,ts = their[0][0], their[0][1]
        self.assertEqual(mu,tu)
        self.assertEqual(ms,ts)
        return 
    
    
    
    def inf_equal(self, mine,their):
        mu,ms = mine[0],mine[1]
        tu,ts = their[0][0], their[0][1]
        self.assertEqual(mu,-ia.inf)
        self.assertEqual(ms,ia.inf)
        self.assertEqual(tu,-interval.inf)
        self.assertEqual(ts,interval.inf)
        return 

    def test_add(self):
        a = ia(1.,2.)
        b = interval.interval([1.,2.])
        self.sp_equal(a+a,b+b)
        
        
        a = ia(-1.,2.)
        b = interval.interval([-1.,2.])
        self.sp_equal(a+a,b+b)
        
        
        a = ia(0.,0.)
        b = interval.interval([0.,0.])
        self.sp_equal(a+a,b+b)
        
        
    def test_sub(self):
        a = ia(1.,2.)
        b = interval.interval([1.,2.])
        self.sp_equal(a-a,b-b)
        
        
        a = ia(-1.,2.)
        b = interval.interval([-1.,2.])
        self.sp_equal(a-a,b-b)
        
        
        a = ia(0.,0.)
        b = interval.interval([0.,0.])
        self.sp_equal(a-a,b-b)
        

    def test_mul(self):
        a = ia(1.,2.)
        b = interval.interval([1.,2.])
        self.sp_equal(a*a,b*b)
        
        
        a = ia(-1.,2.)
        b = interval.interval([-1.,2.])
        self.sp_equal(a*a,b*b)
        
        
        a = ia(0.,0.)
        b = interval.interval([0.,0.])
        self.sp_equal(a*a,b*b)
        
        
    def test_div(self):
        a = ia(1.,2.)
        b = interval.interval([1.,2.])
        self.sp_equal(a/a,b/b)
        
        
        a = ia(1.,2.)
        a1 = ia(-1.,2.)
        
        #really should build an infinity object
        # instead of doing this!:
        atest = a/a1
        self.assertEqual(atest[0],ia(-2e+25, -1.0) )
        self.assertEqual(atest[1],ia(0.5, 2e+25) )
        
        

                  
class TestStringMethods(unittest.TestCase):

    def test_upper(self):
        self.assertEqual('foo'.upper(), 'FOO')

    def test_isupper(self):
        self.assertTrue('FOO'.isupper())
        self.assertFalse('Foo'.isupper())

    def test_split(self):
        s = 'hello world'
        self.assertEqual(s.split(), ['hello', 'world'])
        # check that s.split fails when the separator is not a string
        with self.assertRaises(TypeError):
            s.split(2)
            
            
if __name__ == '__main__':
    
    xia = ia(0.,1.)
    x   = ia( ad(0.,1.,0.), ad(1.,1.,0.)  )
    
    a = xia+xia
    b = x+x
    
    c = x*x
    
    grad1 = np.matrix([1.,0.])
    hess1 = np.matrix([[0.,0.],[0.,0.]])
    x1sup = ad(0., grad1, hess1)
    x1inf = ad(1., grad1, hess1)
    
    grad2 = np.matrix([0.,1.])
    hess2 = np.matrix([[0.,0.],[0.,0.]])
    x2sup = ad(3., grad2, hess2)
    x2inf = ad(5., grad2, hess2)
    
    xv = np.asarray([ia(x1sup,x1inf), ia(x2sup,x2inf)])
    
    xv1 = xv+xv
    print 'Time analysis'
    #with Timer() as t:
    #    r = np.dot(xv,xv)
    #print "=> elasped time for dot product: %s s" % t.secs
    
    nu = ia(0.,0.)
    
    
    print '\nJaulin Inversion'
    print 'Initial Box'
    x = ia(.5,1.)
    y = ia(2.,3.)
    
    print 'x = {}'.format(x)
    print 'y = {}'.format(y)
    print 'Function := y=exp(x)'
    cy = y & x.exp()
    cx = x & y.log()
    print 'minimal, idempotent, contraction:'
    print 'c([x]) = {}'.format(cx)
    print 'c([y]) = {}'.format(cy)
    
    unittest.main()