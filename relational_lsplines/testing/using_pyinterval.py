# -*- coding: utf-8 -*-
"""
Created on Fri Nov 27 09:45:26 2015

@author: lukemcculloch
"""
import numpy as np
import interval as intpy
from interval_arithmetic import ia
#import intpy.interval as intpy
#import intpy.imath as intmath

interval = intpy.interval
intmath = intpy.imath


class IATest(object):
    def __init__(self, clist, alist):
        self.clist = clist
        self.mylist = alist
        self.fail_dict = {}
        return
        
    def __call_(self):
        return
        
    def iseq(self):
        return
    
    def equal(self, a, b):
        extr = a.extrema
        lex = len(extr)
        if lex == 2:
            o1 = a.extrema[0][0]
            o2 = a.extrema[1][1]
            print '    ',o1, b.inf
            print '    ',o2, b.sup
        if lex == 1:
            print '    scalar result in their stuff'
            o1 = a.extrema[0][0]
            o2 = a.extrema[0][0]
            print '    ',o1, b.inf
            print '    ',o2, b.sup
        return o1 == b.inf, o2 == b.sup
        
    def test_add(self):
        print ''
        print 'ADDITION'
        self.fail_dict['addition'] = []
        for o,m in zip(self.clist, self.mylist):
            print ''
            eqtpl = self.equal(o+o,m+m)
            if not (eqtpl[0] and eqtpl[1]):
                print 'Failed!!!!'
                print 'FAILED addition failed on {}+{} = {}+{}'.format(o,o,m,m)
                print '--------------------------------------------------------'
                self.fail_dict['addition'].append([o,m])
            else:
                print '    passed  {}+{} = {}+{}'.format(o,o,m,m)
        return
        
    def test_subtract(self):
        print ''
        print 'SUBTRACTION'
        self.fail_dict['subtraction'] = []
        for o,m in zip(self.clist, self.mylist):
            print ''
            eqtpl = self.equal(o-o,m-m)
            if not (eqtpl[0] and eqtpl[1]):
                print 'Failed!!!!'
                print 'FAILED subtraction failed on {}-{} = {}-{}'.format(o,o,m,m)
                print '--------------------------------------------------------'
                self.fail_dict['subtraction'].append([o,m])
            else:
                print '    passed  {}-{} = {}-{}'.format(o,o,m,m)
        return
        

    def test_multiply(self):
        print ''
        print 'MULTIPLICATION'
        self.fail_dict['multiplication'] = []
        for o,m in zip(self.clist, self.mylist):
            print ''
            eqtpl = self.equal(o*o,m*m)
            if not (eqtpl[0] and eqtpl[1]):
                print 'Failed!!!!'
                print 'FAILED multiplication failed on {}*{} = {}*{}'.format(o,o,m,m)
                print '--------------------------------------------------------'
                self.fail_dict['multiplication'].append([o,m])
            else:
                print '    passed  {}*{} = {}*{}'.format(o,o,m,m)
        return
        
    def test_division(self):
        print ''
        print 'DIVISION'
        self.fail_dict['division'] = []
        for o,m in zip(self.clist, self.mylist):
            print ''
            eqtpl = self.equal(o/o,m/m)
            if not (eqtpl[0] and eqtpl[1]):
                print 'Failed!!!!'
                print 'FAILED division failed on {}/{} = {}/{}'.format(o,o,m,m)
                print '--------------------------------------------------------'
                self.fail_dict['division'].append([o,m])
            else:
                print '    passed  {}/{} = {}/{}'.format(o,o,m,m)
        return
        
    def test_pow(self):
        print ''
        print 'power'
        self.fail_dict['power'] = []
        powers = [0,1,2,3,4,5,6]
        for p in powers:
            for o,m in zip(self.clist, self.mylist):
                print ''
                eqtpl = self.equal(o**p,m**p)
                if not (eqtpl[0] and eqtpl[1]):
                    print 'Failed!!!!'
                    print 'FAILED Power failed on {}**{} = {}**{}'.format(o,p,m,p)
                    print '--------------------------------------------------------'
                    self.fail_dict['power'].append([o,m,p])
                else:
                    print '    passed  {}**{} = {}**{}'.format(o,p,m,p)
        return
        
    def test_arctan(self):
        print ''
        print 'arctan'
        self.fail_dict['arctan'] = []
        for o,m in zip(self.clist, self.mylist):
            print ''
            eqtpl = self.equal(intmath.atan(o),m.arctan2(1.))
            if not (eqtpl[0] and eqtpl[1]):
                print 'Failed!!!!'
                print 'FAILED arctan failed on intmath.atan({}/{}) = arctan2({}/{})'.format(o,1,m,1)
                print '--------------------------------------------------------'
                self.fail_dict['arctan'].append([o,m])
            else:
                print '    passed  intmath.atan({}/{}) = arctan2({}/{})'.format(o,1,m,1)
        return
        
    
    def test_sqrt(self):
        print ''
        print 'sqrt'
        self.fail_dict['sqrt'] = []
        for o,m in zip(self.clist, self.mylist):
            print ''
            eqtpl = self.equal(  intmath.sqrt(o**2) , (m**2).sqrt() )
            if not (eqtpl[0] and eqtpl[1]):
                print 'Failed!!!!'
                print 'FAILED sqrt failed on intmath.atan({}/{}) = arctan2({}/{})'.format(o,1,m,1)
                print '--------------------------------------------------------'
                self.fail_dict['sqrt'].append([o,m])
            else:
                print '    passed  intmath.sqrt({}) = sqrt({})'.format(o,m)
        return
        
    
    def test_dot_product(self):
        print ''
        print 'dot product'
        self.fail_dict['dot'] = []
        
        theirs = 0.
        for el in self.clist:
            theirs += el*el
        q=np.asarray(self.mylist)
        k = np.dot(q,q)
        eqtpl = self.equal(theirs,k)
        if not (eqtpl[0] and eqtpl[1]):
            print 'Failed!!!!'
            print 'FAILED dot product failed on:'
            for a,b in zip(self.clist, self.mylist):
                a,b
            print '--------------------------------------------------------'
            self.fail_dict['dot'].append([self.clist,self.mylist])
        else:
            print '    passed  dot product'
        return
        



if __name__ == '__main__':
    a = interval([1.,2.])
    b = interval([1.0, 4.0])
    c = interval([0.0, 4.0])
    d = interval([-4.0, -1.0])
    e = interval([-4.0, 4.0])
    g = interval([-1.,4.])
    
    clist = [a,b,c,d,e]
    
    a1 = ia(1.,2.)
    a2 = ia(1.,4.)
    a3 = ia(0.,4.)
    a4 = ia(-4.,-1.)
    a5 = ia(-4.,4.)
    
    mylist = [a1,a2,a3,a4,a5]
    
    TC = IATest(clist, mylist)
    TC.test_add()
    TC.test_subtract()
    TC.test_multiply()
    TC.test_division()
    TC.test_pow()
    TC.test_sqrt()
    TC.test_arctan()
    TC.test_dot_product()
    
    self = TC