##
## Interval Analysis
##
## Luke McCulloch
## Spetember 25 2015
##
import numpy as np
import operator as op #ipow for relational power equations
import copy
import matplotlib.pyplot as plt
from time_analysis import Timer

#get rid of later:
import inspect #http://stackoverflow.com/questions/2654113/python-how-to-get-the-callers-method-name-in-the-called-method
  

import fpu  

try:
    import mpmath  as mp #for outward rounding:
    has_mp = True
except:
    has_mp = False
    
    
runtime_testing = False

def df(a):
    #def dff(a):
    #    return lambda: a
    return lambda: a #dff(a)
    
class ia(object):
    infinity = 1.e25
    def __init__(self, inf, sup, has_mp=True, isempty = False, name='unknown_var'):
        self.inf                    = inf
        self.sup                    = sup
        self.isempty                = isempty
        self.has_mp                 = has_mp
        self.infinity               = 1.e25
        self.name                   = name
        #try:
        self.inf = fpu.down(df(self.inf))#
        self.sup = fpu.up(df(self.sup))#
        #self.inf = float(mp.mpf(mp.iv.mpf(str(self.inf)).a))
        #self.sup = float(mp.mpf(mp.iv.mpf(str(self.sup)).b))
#        except:
#            print 'warning, using old style IA with AD'
#            print 'No validation provided for derivatives in his implementation!'
#            self.inf = float(mp.mpf(mp.iv.mpf(str(self.inf.value)).a))
#            self.sup = float(mp.mpf(mp.iv.mpf(str(self.sup.value)).b))
        return
        
    #@classmethod
    #def infinty(self):
    #    return 1.e25
    
    #def __cmp__(self, other):
    #    #if hasattr(other, 'inf'):
    #    return self.inf.__cmp__(other.inf)
        
    
    def __call__(self ):
        a = self.inf
        c = self.sup 
        return a,c

    def __repr__(self):
        return "ia({}, {})".format(self.inf, self.sup)
        
    def __str__(self):
        return "ia({},{})".format(self.inf, self.sup)
    
    
    def roundedreturn(self,p):
        return ia(round(self.inf,p),round(self.sup,p))


    #def __deepcopy__(self, memo):
    #    return ia(copy.deepcopy(self))
    
    def __neg__(self):
        return ia(-self.sup,-self.inf)
            
        
    def __add__(self, other):
        if self.inf>self.sup:
            print 'error, interval is inverted'+'__add__'
        try:
            return ia(self.inf + other.inf, self.sup + other.sup)
        except:
            return ia(self.inf + other, self.sup + other)

    def __radd__(self, other):
        return self.__add__(other)
        #        if self.inf>self.sup:
        #            print 'error, interval is inverted'+'__radd__'
        #        try:
        #            return ia(self.inf + other.inf, self.sup + other.sup)
        #        except:
        #            return ia(self.inf + other, self.sup + other)
            
    def __sub__(self, other):
        if self.inf>self.sup:
            print 'error, interval is inverted '+'__sub__'
        try:
            return ia(self.inf - other.sup, self.sup - other.inf)
        except:
            return ia(self.inf - other, self.sup - other)
    
    def __rsub__(self, other):
        return self.__sub__(other)
        #        if self.inf>self.sup:
        #            print 'error, interval is inverted'+'__rsub__'
        #        return (-self) + other
        
    
    def __mul__(self, other):
        """Moore, page 23
        """
        if self.inf>self.sup:
            print 'error, interval is inverted '+'__mul__'
            print 'self = {}'.format(self.name)
            print 'self could be value, grad, or hessian!'
            print 'exception:', 1./0.
        try:
            S = (self.inf*other.inf, 
                 self.inf*other.sup, 
                 self.sup*other.inf,
                 self.sup*other.sup)
            return ia(min(S),max(S))
        except:
            S = (self.inf*other, 
                 self.sup*other)
            return ia(min(S),max(S))
        return
            
                  
    def Conditional__mul__(self, other):
        """ Hansen page 8, Moore page 24, Hansen page 
        """
        try:
            a = self.inf
            b = self.sup
            c = other.inf
            d = other.sup
            
            if  a >= 0. and c >= 0.:
                return ia(a*c,b*d)
            elif a >= 0. and c<0.<d:
                return ia(b*c,b*d)
            elif a >= 0. and d <= 0.:
                return ia(b*c,a*d)
            elif a<0.<b and c>=0.:
                return ia(a*d,b*d) #corrected
            elif a<0.<b and d<=0.:
                return ia(b*c,a*c) # I say ad  hansen says ac for sup #corrected# ia(b*d,a*d)- I agree with Eldon
            elif b<=0. and c>=0.:
                return ia(a*d,b*c)
            elif b<=0. and c<0.<d:
                return ia(a*d,a*c)
            elif b<=0. and d<=0.:
                return ia(b*d,a*c)
            elif a<0.<b and c<0.<d:
                return ia(min(b*c,a*d,a*c,b*d),max(b*c,a*d,a*c,b*d))
                #return ia(min(b*c,a*d),max(a*c,b*d))
            else:
                print 'error in multiplication'
                c1 = a*c
                c2 = a*d
                c3 = b*c
                c4 = b*d
                return ia(min(c1,c2,c3,c4),max(c1,c2,c3,c4))
                
            
        except:
            a = self.inf
            b = self.sup
            c = other
            d = other
            
            if  a >= 0. and c >= 0.:
                return ia(a*c,b*d)
            elif a >= 0. and c<0.<d:
                return ia(b*c,b*d)
            elif a >= 0. and d <= 0.:
                return ia(b*c,a*d)
            elif a<0.<b and c>=0.:
                return ia(a*d,b*d) #corrected
            elif a<0.<b and d<=0.:
                return ia(b*c,a*c) # I say ad  hansen says ac for sup #corrected# ia(b*d,a*d)- I agree with Eldon
            elif b<=0. and c>=0.:
                return ia(a*d,b*c)
            elif b<=0. and c<0.<d:
                return ia(a*d,a*c)
            elif b<=0. and d<=0.:
                return ia(b*d,a*c)
            elif a<0.<b and c<0.<d:
                return ia(min(b*c,a*d,a*c,b*d),max(b*c,a*d,a*c,b*d))
                #return ia(min(b*c,a*d),max(a*c,b*d))
            else:
                print 'error in multiplication'
                c1 = a*c
                c2 = a*d
                c3 = b*c
                c4 = b*d
                return ia(min(c1,c2,c3,c4),max(c1,c2,c3,c4))
            
    def My__mul__(self, other):
        """ Moore page 24 (I could use but didn't)
        """
        try:
            a = self.inf*other.inf
            b = self.inf*other.sup
            c = self.sup*other.inf
            d = self.sup*other.sup
            v_min = min(a,b,c,d)
            v_max = max(a,b,c,d)
            return ia(v_min, v_max)
        except:
            a = self.inf*other
            b = self.sup*other
            v_min = min(a,b)
            v_max = max(a,b)
            return ia(v_min, v_max)
    
    def Mat__mul__(self, other, invert = False):
        """ Moore page 24 (I could use but didn't)
        """
        if invert == False:
            try:
                a = np.dot(self.inf,other.inf)
                b = np.dot(self.inf,other.sup)
                c = np.dot(self.sup,other.inf)
                d = np.dot(self.sup,other.sup)
                v_min = min(a,b,c,d)
                v_max = max(a,b,c,d)
                return ia(v_min, v_max)
            except:
                a = np.dot(self.inf,other)
                b = np.dot(self.sup,other)
                v_min = min(a,b)
                v_max = max(a,b)
                return ia(v_min, v_max)
        else:
            return self.__div__(other, invert = False)
    
    def __rmul__(self, other, invert = False):
        return self.__mul__(other)
        #        if self.inf>self.sup:
        #            print 'error, interval is inverted '+'__rmul__'
        #        if invert == False:
        #            return self.__mul__(other)
        #        else:
        #            return self.__div__(other, invert = False)
        #    
    
    
    def NOinvert_one_oldest(self):
        if self.inf>self.sup:
            print 'error, interval is inverted '+'invert_one'
        d1 = copy.deepcopy(self) 
        if self.contains(0.):
            d1.sup = self.infinity
            if self.inf == 0.:
                d1.inf = 1./self.sup
            elif self.sup==0:
                d1.inf = 1./self.inf
            else:
                d1.inf = min(1./self.sup,1./self.inf)
            ans = d1 
            ans.is_extended = True
            return ans
        elif (0.<self.inf or self.sup<0.):
            a=1./self.inf
            b=1./self.sup
            return ia(b,a)
        else:
            print 'error in inversion'
            return
            
    def NOinvert_one_old(self):
        d1 = copy.deepcopy(self) 
        if self.contains(0.):
            print 'other contains 0 - should have caught in __div__ yet we are in invert_one'
            print 'ia self = {}'.format(self)
            d1.sup = self.infinity
            if self.inf == 0.:
                #d1.inf = 1./self.sup
                d1.inf = min(0.,1./self.sup,1./self.inf)
            elif self.sup==0:
                #d1.inf = 1./self.inf
                d1.inf = min(0.,1./self.sup,1./self.inf)
            else:
                d1.inf = min(0.,1./self.sup,1./self.inf)
            ans = d1 
            ans.is_extended = True
            return ans
        elif (0.<self.inf or self.sup<0.):
            a=1./self.inf
            b=1./self.sup
            return ia(b,a) #is this part captured in you new version?
        else:
            print 'error in inversion'
            return
    """
                    if self.contains(0.):
                        ans = copy.deepcopy(self)
                        ans.inf = -ia.infinity
                        ans.sup = ia.infinity
                        return ia(-ia.infinity,ia.infinity)
    #"""
    def invert_one(self):
        if self.contains(0.):
            print 'other contains 0'
            c = self.inf
            d = self.sup
            # if 1. < 0.:
            """
            if d == 0:
                return ia(1./c,ia.infinity)
            elif c == 0:
                return ia(-ia.infinity,1./d)
            else:
                assert(c<0.<d):
                print 'NEED TO SPLIT THE INTERVAL (2.4.1 Hansen page 10, union for 1 on top)'
                return ia(-ia.infinity,ia.infinity)
            #"""
            # if 1. > 0.:
            if d == 0.:
                if c != 0.:
                    return ia(-ia.infinity,1./c)
                else:
                    return ia(-ia.infinity,ia.infinity)
            elif c == 0.:
                return ia(1./d, ia.infinity)
            else:
                assert(c<0.<d)
                #minm = min(0., 1./c, 1./d)
                minm = min(1./c, 1./d)
                print 'NEED TO SPLIT THE INTERVAL (2.4.1 Hansen page 10, 2nd union, but 1 on top!)'
                #curframe = inspect.currentframe()
                #calframe = inspect.getouterframes(curframe, 2)
                #print 'caller name:', calframe[1][3]
                #return ia(-ia.infinity,1./c),ia(1./d,ia.infinity)
                #return ia(-ia.infinity,ia.infinity) #standard union
                return ia(1./d,ia.infinity) ##tlm experiment
                #return ia(minm,ia.infinity)
        else:
            assert(0.<self.inf or self.sup<0.)
            a=1./self.inf
            b=1./self.sup
            return ia(b,a)
            #return (1./self)
    
    def invert_split(self):
        print 'split inversion possible'
        if self.contains(0.):
            print 'other contains 0'
            c = self.inf
            d = self.sup
            # if 1. < 0.: (why was this not a thing?)
            """
            if d == 0:
                return ia(1./c,ia.infinity)
            elif c == 0:
                return ia(-ia.infinity,1./d)
            else:
                assert(c<0.<d):
                print 'NEED TO SPLIT THE INTERVAL (2.4.1 Hansen page 10, union for 1 on top)'
                return ia(-ia.infinity,ia.infinity)
            #"""
            # if 1. > 0.:
            if d == 0.:
                return [ia(-ia.infinity,1./c)]
            elif c == 0.:
                return [ia(1./d, ia.infinity)]
            else:
                assert(c<0.<d)
                minm = min(0., 1./c, 1./d)
                print 'splitting! (2.4.1 Hansen page 10, 2nd union, but 1 on top!)'
                #new_interval       = ia(1./d,ia.infinity)
                #new_interval.split = True
                #new_interval.second = ia(-ia.infinity,1./c)
                #return [ia(-ia.infinity,1./c),ia(1./d,ia.infinity)]
                return [ia(1./d,ia.infinity),ia(-ia.infinity,1./c)]
                #return ia(-ia.infinity,ia.infinity)

        else:
            assert(0.<self.inf or self.sup<0.)
            a=1./self.inf
            b=1./self.sup
            return [ia(b,a)]
            #return (1./self)
    
    def NO__div__old(self, other, invert = False):
        if self.inf>self.sup:
            print 'error, interval is inverted '+'__div__'
        if invert == False:
            if isinstance(other, ia):
                if other.contains(0.):
                    #a = self.inf
                    #b = self.sup
                    #c = other.inf
                    #d = other.sup
                    #if
                    return ia(-ia.infinity,ia.infinity)
                
                else:
                    try:
                        return self *other.invert_one()
                    except:
                        return self*(1./other)
            else:
                return self*(1./other)
        else:
            return self.__mul__(other)
    
    def __div__(self, other, invert = False):
        if invert == False:
            if isinstance(other, ia):
                if other.contains(0.):
                    if self.contains(0.):
                        return ia(-self.infinity,self.infinity)
                    else:
                        #print 'other splits! = {}'.format(other)
                        #print'split from simple __div__, frame is:'
                        #curframe = inspect.currentframe()
                        #calframe = inspect.getouterframes(curframe, 2)
                        #print 'caller name:', calframe[1][3]
                        #print 'do inver_one:'
                        return self *other.invert_one()#_old()
                    """
                    print 'other contains 0'
                    a = self.inf
                    b = self.sup
                    c = other.inf
                    d = other.sup
                    
                    if (a < 0. < b): #self.contains(0.):
                        print '0 is in the numertor and denominator'
                        return ia(-ia.infinity,ia.infinity)
                    if b <= 0.:
                        if d == 0:
                            return ia(b/c,ia.infinity)
                        elif c == 0:
                            return ia(-ia.infinity,b/d)
                        else:
                            assert(c<0.<d)
                            minm = min(0.,a/c,a/d,b/c,b/d)
                            print 'NEED TO SPLIT THE INTERVAL (2.4.1 Hansen page 10, 1st union)'
                            #return ia(-ia.infinity,b/d), ia(b/c,ia.infinity)
                            #return ia(-ia.infinity,ia.infinity)
                            return ia(minm,ia.infinity)
                    elif a >= 0.:
                        if d == 0.:
                            return ia(-ia.infinity,a/c)
                        elif c == 0.:
                            return ia(a/d, ia.infinity)
                        else:
                            assert(c<0.<d)
                            minm = min(0.,a/c,a/d,b/c,b/d)
                            print 'NEED TO SPLIT THE INTERVAL (2.4.1 Hansen page 10, 2nd union)'
                            #return ia(-ia.infinity,a/c),ia(a/d,ia.infinity)
                            #return ia(-ia.infinity,ia.infinity)
                            return ia(minm,ia.infinity)
                    else:
                        print '-----------------------------'
                        print 'Should Never Get Here! because self contains 0!'
                        print 'self = {}'.format(self)
                        print 'other = {}'.format(other)
                        print '-----------------------------'
                        #return self *other.invert_one()
                    #"""
                else: #other doesn't contain 0
                    #print 'other is simple = {}'.format(other)
                    assert(not other.contains(0.))
                    try:
                        #print'split from __div__, frame is:'
                        #curframe = inspect.currentframe()
                        #calframe = inspect.getouterframes(curframe, 2)
                        #print 'caller name:', calframe[1][3]
                        #print 'do inver_one:'
                        return self *other.invert_one()#_old()
                    except:
                        print 'something is wierd'
                        print 'other is almost simple = {}'.format(other)
                        return self*(1./other)
                        
            else: #other is not interval
                return self*(1./other)
                
        else: #invert is true -> an experimental thing..
            return self.__mul__(other)
            
#            try:
#                return self *other.invert_one()
#            except:
#                return self*(1./other)
#        else:
#            return self.__mul__(other)
    def __rdiv__(self, other, invert = False):
        #return self.__div__(other)
        if self.inf>self.sup:
            #print 'other = {}'.format(other)
            #print 'self = {}'.format(self)
            #return self.invert_one()*other#other/self#
            print 'error, interval is inverted '+'__rdiv__'
        if invert == False:
            return self.invert_one()*other
        else:
            return self.invert_one()*other
    
    def __truediv__(self, other, invert = False):
        if isinstance(other, ia):
            if other.contains(0.):
                if self.contains(0.):
                    return ia(-self.infinity,self.infinity)
                else:
                    #print 'other splits! = {}'.format(other)
                    #print'split from simple __div__, frame is:'
                    #curframe = inspect.currentframe()
                    #calframe = inspect.getouterframes(curframe, 2)
                    #print 'caller name:', calframe[1][3]
                    #print 'do inver_one:'
                    return self *other.invert_one()#_old()
                """
                print 'other contains 0'
                a = self.inf
                b = self.sup
                c = other.inf
                d = other.sup
                
                if (a < 0. < b): #self.contains(0.):
                    print '0 is in the numertor and denominator'
                    return ia(-ia.infinity,ia.infinity)
                if b <= 0.:
                    if d == 0:
                        return ia(b/c,ia.infinity)
                    elif c == 0:
                        return ia(-ia.infinity,b/d)
                    else:
                        assert(c<0.<d)
                        minm = min(0.,a/c,a/d,b/c,b/d)
                        print 'NEED TO SPLIT THE INTERVAL (2.4.1 Hansen page 10, 1st union)'
                        #return ia(-ia.infinity,b/d), ia(b/c,ia.infinity)
                        #return ia(-ia.infinity,ia.infinity)
                        return ia(minm,ia.infinity)
                elif a >= 0.:
                    if d == 0.:
                        return ia(-ia.infinity,a/c)
                    elif c == 0.:
                        return ia(a/d, ia.infinity)
                    else:
                        assert(c<0.<d)
                        minm = min(0.,a/c,a/d,b/c,b/d)
                        print 'NEED TO SPLIT THE INTERVAL (2.4.1 Hansen page 10, 2nd union)'
                        #return ia(-ia.infinity,a/c),ia(a/d,ia.infinity)
                        #return ia(-ia.infinity,ia.infinity)
                        return ia(minm,ia.infinity)
                else:
                    print '-----------------------------'
                    print 'Should Never Get Here! because self contains 0!'
                    print 'self = {}'.format(self)
                    print 'other = {}'.format(other)
                    print '-----------------------------'
                    #return self *other.invert_one()
                #"""
            else: #other doesn't contain 0
                #print 'other is simple = {}'.format(other)
                assert(not other.contains(0.))
                try:
                    #print'split from __div__, frame is:'
                    #curframe = inspect.currentframe()
                    #calframe = inspect.getouterframes(curframe, 2)
                    #print 'caller name:', calframe[1][3]
                    #print 'do inver_one:'
                    return self *other.invert_one()#_old()
                except:
                    print 'something is wierd'
                    print 'other is almost simple = {}'.format(other)
                    return self*(1./other)
                    
        else: #other is not interval
            return self*(1./other)
        
    def __rtruediv__(self, other, invert = False):
        if self.inf>self.sup:
            print 'error, interval is inverted '+'__rdiv__'
        if invert == False:
            return self.invert_one()*other
        else:
            return self.invert_one()*other
            
    def div_split(self, other, invert = False):
        assert(isinstance(self, ia))
        if invert == False:
            if other.contains(0.):
                if self.contains(0.):
                    return [ia(-self.infinity,self.infinity)]
                else:
                    return [self *el for el in other.invert_split()]
            else: #other doesn't contain 0
                #try:
                assert(not other.contains(0.))
                #return self *other.invert_one()#_old()
                return [self*el for el in other.invert_split()]
                #except:
                #    print 'something is wierd'
                #    print 'other is almost simple = {}'.format(other)
                #    return [self*(1./other)]
                        
                
        else: #invert is true -> an experimental thing..
            return [self.__mul__(other)]
            
    
    def isodd(self, n):
        return n%2==1
    def iseven(self, n):
        return n%2==0
    
    def __pow__(self, n, invert = False):
        """D,L,M Verified Real Number Calculation
        A library for Interval Arithmetic
        page 3
        """
        assert(isinstance(n, int))
        #if not isinstance(n, int):
        #    return 
        if self.inf>self.sup:
            print 'error, interval is inverted '+'__pow__'
        if n<0:
            invs = 1./self
            assert(n>=0.)
            return invs.__pow__(abs(n))
        if n==0:
            return ia(1.,1.)
        elif self.inf>=0. or ( self.isodd(n) ):
            return ia(self.inf**n, self.sup**n)
        elif self.sup<=0. and self.iseven(n):
            return ia(self.sup**n, self.inf**n)
        else:
            return ia(0., max(self.inf**n, self.sup**n) )
            
            
        
        
            
    def A1__pow__(self, n, invert = False):
        """Hansen page 9
        """
        if invert == False:
            #if n==2:
            #    return self.pow_of_2()
            if n==0:
                return ia(1.,1.) #crucial!
            elif self.inf>=0. or ( self.contains(0.) and self.isodd(n) ):
                return ia(self.inf**n, self.sup**n)
            elif self.sup<=0.:
                return ia(self.sup**n, self.inf**n)
            elif self.contains(0.) and (self.iseven(n) and n >0):
                return ia(0., max(self.inf**n, self.sup**n) )
            #elif n==.5:
            #    return self.sqrt()
            #elif 0.<=self.inf:
            #    return ia(self.inf**n, self.sup**n)
        else:
            n = 1./n
            return self.__pow__(n,invert = False)
        
    def Moore__pow__(self, n, invert = False):
        """Moore page 60
        """
        if invert == False:
            if n==2:
                return self.pow_of_2()
            elif n==0:
                return ia(1.,1.)
            #        try:
            #            assert(n>=0.)
            #        except:
            #            print 'pow issue: n>=2?, n = {}'.format(n)
            # consider pow n=.3, immediately, the gradient has a negvative pow!
            if self.inf>0. or self.isodd(n):
                return ia(self.inf**n, self.sup**n)
            elif self.sup<0. and self.iseven(n):
                return ia(self.inf**n, self.sup**n)
            elif self.contains(0.) and self.iseven(n):
                return ia(0., abs(self)**n)
            elif n==.5:
                return self.sqrt()
            elif 0.<=self.inf:
                return ia(self.inf**n, self.sup**n)
        else:
            n = 1./n
            return self.__pow__(n,invert = False)

    def pow_of_2(self):
        """Moore page 49
        """
        if 0.<= self.inf:
            return ia(self.inf**2, self.sup**2)
        elif self.sup <=0.:
            return ia(self.sup**2, self.inf**2)
        elif self.inf < 0. < self.sup:#self.contains(0.):
            return ia(0., max(self.inf**2, self.sup**2))
        
    def sqrt(self):
        if self.inf>self.sup:
            print 'error, interval is inverted '+'sqrt'
        assert(0.<=self.inf)
        a = np.sqrt(self.inf)
        b = np.sqrt(self.sup)
        return ia(min(a,b), max(a,b))
        #??
        #if extended interval arithmetic, then:
        #?
        #c = -b
        #d = -a
        #return ia(min(a,b,c,d), max(a,b,c,d))
        
    def exp(self, invert = False):
        if self.inf>self.sup:
            print 'error, interval is inverted '+'exp'
        if invert==False:
            a = np.exp(self.inf)
            b = np.exp(self.sup)
            return ia(a, b)
        else:
            return self.log()
        
    def log(self, invert = False):
        if self.inf>self.sup:
            print 'error, interval is inverted '+'log'
        assert(self.inf > 0.)
        if invert==False:
            a = np.log(self.inf)
            b = np.log(self.sup)
            return ia(a, b)
        else:
            return self.exp()
        
    def sin(self):
        """
            TODO: check interval sign before use!
        """
        total_range = self.sup - self.inf
        the_min     = np.sin(self.inf)
        the_max     = np.sin(self.sup)
        
        the_min2     = np.sin(self.inf)
        the_max2     = np.sin(self.sup)
        
        # evaluate these if's on x, not the sin(x)! :
        if total_range >= 2.*np.pi:
            the_max   = 1.0
            the_min   = -1.0
        else:
            if ((self.inf%(2.*np.pi) < np.pi/2.) and (total_range >= (np.pi/2.-self.inf%(2.*np.pi)))) :
                the_max2   = 1.0
                
            if ((self.inf%(2.*np.pi) < 3.*np.pi/2.) and (total_range >= (3.*np.pi/2.-self.inf%(2.*np.pi))) ):
                the_min2   = -1.0
                
        v_min = min(the_min, the_max, the_min2, the_max2)
        v_max = max(the_min, the_max, the_min2, the_max2)
        return ia( v_min,v_max)
    
    def cos(self):
        """
            TODO: recheck ia cosine
        """
        total_range = self.sup - self.inf
        the_min     = np.cos(self.inf)
        the_max     = np.cos(self.sup)
        the_min2    = np.cos(self.inf)
        the_max2    = np.cos(self.sup)
        # evaluate these if's on x, not the sin(x)! :
        if total_range >= 2.*np.pi:
            the_max   =  1.0
            the_min   = -1.0
        else:
            if (   ( (self.inf+np.pi/2.)%(2.*np.pi) < np.pi/2. ) and 
                ( total_range >= (np.pi/2.-(self.inf+np.pi/2.)%(2.*np.pi)) )  ):
                the_max2   = 1.0
                #flag = True
                
            if (   ( (self.inf+np.pi/2.)%(2.*np.pi) < 3.*np.pi/2. ) and 
                ( total_range >= (3.*np.pi/2.-(self.inf+np.pi/2.)%(2.*np.pi)) )  ):
                the_min2   = -1.0

        v_min = min(the_min, the_max, the_min2, the_max2)
        v_max = max(the_min, the_max, the_min2, the_max2)
        return ia( v_min,v_max)
    
    def tan(self):
        if self.inf>self.sup:
            print 'error, interval is inverted '+'tan'
        total_range = self.sup - self.inf
        assert(total_range <= np.pi/2.)
        the_min     = self.inf.tan()
        the_max     = self.sup.tan()  
        v_min = min(the_min, the_max)
        v_max = max(the_min, the_max)
        return ia(v_min, v_max)
    
    def atan(self):
        if self.inf>self.sup:
            print 'error, interval is inverted '+'atan'
        the_min     = self.inf.atan()
        the_max     = self.sup.atan()
        return ia(min(the_min, the_max),
                            max(the_min, the_max))
        
    def arctan2(self, other):
        if self.inf>self.sup:
            print 'error, interval is inverted '+'arctan2'
        try:
            the_min1     = np.arctan2(self.inf, other.inf)
            the_min3     = np.arctan2(self.inf, other.sup)
            the_max1     = np.arctan2(self.sup, other.inf)
            the_max3     = np.arctan2(self.sup, other.sup)
            v_min = min(the_min1, the_max1,
                        the_min3, the_max3)
            v_max = max(the_min1, the_max1,
                        the_min3, the_max3)
        except:
            the_min1     = np.arctan2(self.inf, other)
            the_max1     = np.arctan2(self.sup, other)
            v_min = min(the_min1, the_max1)
            v_max = max(the_min1, the_max1)
        return ia(v_min, v_max)
    
    
    def le_or(self, other): 
        """the less than or equal to part of __or__
        """
        if self.inf>self.sup:
            print 'error, interval is inverted '+'__or__'
        return ia(min(self.inf, other.inf), min(self.sup, other.sup))
        
    def ge_or(self, other): 
        """the greater than or equal to part of __or__
        """
        if self.inf>self.sup:
            print 'error, interval is inverted '+'__or__'
        return ia(max(self.inf, other.inf), max(self.sup, other.sup))
        
    def __or__(self, other): 
        """actually hull, not union
        """
        if self.inf>self.sup:
            print 'error, interval is inverted '+'__or__'
        return ia(min(self.inf, other.inf), max(self.sup, other.sup))
    
    def __and__(self, other): 
        """intersection
        """        
        if self.inf>self.sup:
            print 'error, interval is inverted '+'__and__'
        if self.sup<other.inf or self.inf>other.sup:
            return ia(0.,0., isempty=True)
        elif self.sup==other.inf :
            return ia(self.sup,self.sup)
        elif self.inf==other.sup:
            return ia(self.inf,self.inf)
        else:
            return ia(max(self.inf, other.inf), min(self.sup, other.sup))

    def hull(self, other):
        return self.__or__(other)
    
    def __eq__(self, other):
        if isinstance(other, ia):
            if self.inf>self.sup:
                print 'ERROR:__eq__, inf>sup ia({},{})'.format(self.inf,self.sup)
                #print 'error, interval is inverted '+'__eq__'
                #print 'ia({},{})'.format(self.inf,self.sup)
            if (self.inf == other.inf) & (self.sup == other.sup):
                return True
            else:
                return False
        elif isinstance(other, float):
            return False
            
#    def __req__(self, other):
#        if self.inf>self.sup:
#            print 'error, interval is inverted '+'__eq__'
#        if (self.inf == other.inf) & (self.sup == other.sup):
#            return True
#        else:
#            return False
#    def eq(self,other):
#        return 
        
    def __lt__(self, other):
        if self.inf>self.sup:
            print 'error, interval is inverted '+'__lt__'
        try:
            return self.sup < other.inf
        except:
            return self.sup < other
    
    def __le__(self, other):
        if self.inf>self.sup:
            print 'error, interval is inverted '+'__lt__'
        try:
            return self.sup <= other.inf
        except:
            return self.sup <= other
            
    def __gt__(self, other):
        if self.inf>self.sup:
            print 'error, interval is inverted '+'__gt__'
        try:
            return self.inf > other.sup
        except:
            return self.inf > other
    
    def __ge__(self, other):
        if self.inf>self.sup:
            print 'error, interval is inverted '+'__gt__'
        try:
            return self.inf >= other.sup
        except:
            return self.inf >= other
            
    def __abs__(self):
        if self.inf>self.sup:
            print 'error, interval is inverted '+'__abs__'
        return max(abs(self.inf), abs(self.sup))
        
    def width(self):
        if self.inf>self.sup:
            print 'error, interval is inverted '+'width'
        return self.sup - self.inf
    
    def midpoint(self):
        if self.inf>self.sup:
            print 'error, interval is inverted '+'midpoint'
        return (self.inf + self.sup)*.5
    
    def getpoint(self, pt):
        if self.inf>self.sup:
            print 'error, interval is inverted '+'getpoint'
        assert (0.<=pt<=1.)
        p1 = 1. - pt
        p2 = pt
        return (self.inf*p1 + self.sup*p2)
        
    def dist(self,other): #page 63
        if self.inf>self.sup:
            print 'error, interval is inverted '+'dist'
        return  max(abs(self.inf - other.inf), abs(self.sup-other.sup))
        
    def __contains__(self, other):
        """Overide the intrinsic
        'in' operator
        """
        sinf = fpu.down(df(self.inf))
        ssup = fpu.up(df(self.sup))
        if self.inf>self.sup:
            print 'error, interval is inverted '+'__contains__'
        if isinstance(other, ia):
            #if other.inf>=self.inf and other.sup<=self.sup:return True:
            if other.inf>=sinf and other.sup<=ssup:return True
        elif isinstance(other, ad):
            if self.inf <= other.value <=self.sup:return True
        elif isinstance(other, float):
            if self.inf <= other <=self.sup:return True
        else:
            return False
        return
        
    def contains(self, number):
        """Check if a scalar 
        is in the interval
        """
        sinf = fpu.down(df(self.inf))
        ssup = fpu.up(df(self.sup))
        if self.inf>self.sup:
            print 'error, interval is inverted '+'contains'
        result = None
        #if (self.inf <= number) and (number<=self.sup):
        if (sinf <= number) and (number<=ssup):
            result = True
        elif self.inf == number:
            result = True
        elif self.sup == number:
            result = True
        else:
            result = False
        return result
        
                  
from automatic_differentiation import ad
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
    with Timer() as t:
        r = np.dot(xv,xv)
    print "=> elasped time for dot product: %s s" % t.secs
    
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
    
    if runtime_testing:
        import interval as I
        from interval import interval as pyinterval
        
        pya = pyinterval([0.,1.])
        pyinf = pyinterval([2.,I.inf])
        pyb = pyinterval([1.,2.])
        pyc = pyinterval([0.,1.])
        #pyd = pyinterval([0.,1.])
        
        a = ia(0.,1.)
        binf = ia(2.,ia.infinity)
        b = ia(1.,2.)
        c = ia(0.,1.)
        d = ia(1.,2.)
        
        
        
        print '[0,1]*[2, inf] = {}'.format(a*b)
        
        def rump1(a,b):
            """not meant to show a tight answer!
            see
            http://raim2013.lip6.fr/theme/PDF/Intervalles/[Goldsztejn]%20Interval%20Analysis%20-%20Nonlinear%20Computer-Assisted%20Proofs.pdf
            """
            return ((333.75 - a**2)*(b**6)
                    + (a**2) * (11.*(a**2)*(b**2)-121.*(b**4)-2.) 
                    + 5.5*(b**8) + a/(2.*b))
        
        def rump2(a,b):
            a1 = (333.75 - a**2)*(b**6)
            a2 = (a**2) * (11.*(a**2)*(b**2)-121.*(b**4)-2.) 
            a3 = 5.5*(b**8) + a/(2.*b)
            ans = a1+a2+a3
            return ans
        

        a = ia(77617.9999999999999,77617.0000000000000001)
        b = ia(33096.9999999999999 , 33096.0000000000000001)
        print 'rump func = {}'.format(rump1(a,b))
        
        #a = 77617.
        #b = 33096.
        print 'rump func = {}'.format(rump1(77617.,33096.))

        a = ia(77617.,77617.)
        b = ia(33096.,33096.)
        a1 = (333.75 - a**2)*(b**6)
        a2 = (a**2) * (11.*(a**2)*(b**2)-121.*(b**4)-2.) 
        a3 = 5.5*(b**8) + a/(2.*b)
        print a1
        print a2
        print a3
        print a1+a2+a3
        print 'rump func = {}'.format(rump2(a,b))
        
        a = 77617.
        b = 33096.
        print 'rump func = {}'.format(rump2(a,b))
        
        
        pa = pyinterval([77617.,77617.])
        pb = pyinterval([33096.,33096.])
        a1 = (333.75 - pa**2)*(pb**6)
        a2 = (pa**2) * (11.*(pa**2)*(pb**2)-121.*(pb**4)-2.) 
        a3 = 5.5*(pb**8) + pa/(2.*pb)
        print a1
        print a2
        print a3
        print a1+a2+a3
        rump1(pa,pb)
        rump2(pa,pb)
        
        a = ia(1./3,1./3)
        aa = ia(1.,1.)
        bb = ia(3.,3.)
        cc = aa/bb
        
        
    cb =  ia(0.95,0.95)
    cmidshp =  ia(0.959595959615,0.989583333337)
    cp =  ia(0.96,0.989999999983)
    
    cp = ia(.96,.96)
    cmidshp = cmidshp & (cb/cp)
    print cb & (cmidshp*cp)
    
    inf = float(mp.mpf(mp.iv.mpf(str(cmidshp.inf)).a))
    sup = float(mp.mpf(mp.iv.mpf(str(cmidshp.sup)).b))
    
    
    
    def hansen_p18_v1(x,y):
        return (x-y)/(x+y)
    def hansen_p18_v2(x,y):
        return 1.-(2./(1+(x/y)))
    
    def show_rounding():
        a = ia(1.,1.)
        b = ia(3.,3.)
        c = a/b
        print 'a = ',a
        print 'b = ',b
        print 'c = a/b'
        print 'c = ',c
        print 'c.inf = {:18.17f} '.format(c.inf)
        print 'c.sup = {:18.17f} '.format(c.sup)
        return