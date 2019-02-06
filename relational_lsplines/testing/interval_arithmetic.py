##
## Interval Analysis
##
## Luke McCulloch
## Spetember 25 2015
##
import numpy as np
import copy
import matplotlib.pyplot as plt
#from time_analysis import Timer

try:
    import mpmath  as mp #for outward rounding:
    has_mp = True
except:
    has_mp = False

class ia(object):
    def __init__(self, inf, sup, has_mp=True, isempty = False, name='unknown_var'):
        self.inf                    = inf
        self.sup                    = sup
        self.isempty                = isempty
        self.has_mp                 = has_mp
        self.infinity               = 1.e25
        self.name                   = name
        try:
            self.inf = float(mp.mpf(mp.iv.mpf(str(self.inf)).a))
            self.sup = float(mp.mpf(mp.iv.mpf(str(self.sup)).b))
        except:
            print 'warning, using old style IA with AD'
            print 'No validation provided for derivatives in his implementation!'
            self.inf = float(mp.mpf(mp.iv.mpf(str(self.inf.value)).a))
            self.sup = float(mp.mpf(mp.iv.mpf(str(self.sup.value)).b))
        return
    
    def __call__(self ):
        a = self.inf
        c = self.sup 
        return a,c

    def __repr__(self):
        return "ia[{}, {}]".format(self.inf, self.sup)
        
    def __str__(self):
        return "ia[{},{}]".format(self.inf, self.sup)
    
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
        if self.inf>self.sup:
            print 'error, interval is inverted'+'__radd__'
        try:
            return ia(self.inf + other.inf, self.sup + other.sup)
        except:
            return ia(self.inf + other, self.sup + other)
            
    def __sub__(self, other):
        if self.inf>self.sup:
            print 'error, interval is inverted '+'__sub__'
        try:
            return ia(self.inf - other.sup, self.sup - other.inf)
        except:
            return ia(self.inf - other, self.sup - other)
    
    def __rsub__(self, other):
        if self.inf>self.sup:
            print 'error, interval is inverted'+'__rsub__'
        return (-self) + other
        
    
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
        if self.inf>self.sup:
            print 'error, interval is inverted '+'__rmul__'
        if invert == False:
            return self.__mul__(other)
        else:
            return self.__div__(other, invert = False)
    
    
    
    def invert_one(self):
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
        
    def __div__(self, other, invert = False):
        if self.inf>self.sup:
            print 'error, interval is inverted '+'__div__'
        if invert == False:
            try:
                return self *other.invert_one()
            except:
                return self*(1./other)
        else:
            return self.__mul__(other)
    
    def __rdiv__(self, other, invert = False):
        if self.inf>self.sup:
            print 'error, interval is inverted '+'__rdiv__'
        if invert == False:
            return self.invert_one()*other
        else:
            return self.invert_one()*other
        
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
        total_range = self.sup.value - self.inf.value
        the_min     = self.inf.sin()
        the_max     = self.sup.sin()
        
        the_min2     = self.inf.sin()
        the_max2     = self.sup.sin()
        
        # evaluate these if's on x, not the sin(x)! :
        if total_range >= 2.*np.pi:
            the_max.value   = 1.0
            the_min.value   = -1.0
        else:
            if ((self.inf.value%(2.*np.pi) < np.pi/2.) and (total_range >= (np.pi/2.-self.inf.value%(2.*np.pi)))) :
                the_max2.value   = 1.0
                
            if ((self.inf.value%(2.*np.pi) < 3.*np.pi/2.) and (total_range >= (3.*np.pi/2.-self.inf.value%(2.*np.pi))) ):
                the_min2.value   = -1.0
                
        v_min = min(the_min, the_max, the_min2, the_max2)
        v_max = max(the_min, the_max, the_min2, the_max2)
        return ia( v_min,v_max)
    
    def cos(self):
        """
            TODO: recheck ia cosine
        """
        total_range = self.sup.value - self.inf.value
        the_min     = self.inf.cos()
        the_max     = self.sup.cos()
        the_min2    = self.inf.cos()
        the_max2    = self.sup.cos()
        # evaluate these if's on x, not the sin(x)! :
        if total_range >= 2.*np.pi:
            the_max.value   =  1.0
            the_min.value   = -1.0
        else:
            if (   ( (self.inf.value+np.pi/2.)%(2.*np.pi) < np.pi/2. ) and 
                ( total_range >= (np.pi/2.-(self.inf.value+np.pi/2.)%(2.*np.pi)) )  ):
                the_max2.value   = 1.0
                #flag = True
                
            if (   ( (self.inf.value+np.pi/2.)%(2.*np.pi) < 3.*np.pi/2. ) and 
                ( total_range >= (3.*np.pi/2.-(self.inf.value+np.pi/2.)%(2.*np.pi)) )  ):
                the_min2.value   = -1.0

        v_min = min(the_min, the_max, the_min2, the_max2)
        v_max = max(the_min, the_max, the_min2, the_max2)
        return ia( v_min,v_max)
    
    def tan(self):
        if self.inf>self.sup:
            print 'error, interval is inverted '+'tan'
        total_range = self.sup.value - self.inf.value
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
        else:
            return ia(max(self.inf, other.inf), min(self.sup, other.sup))

    def hull(self, other):
        return self.__or__(other)
    
    def __eq__(self, other):
        if self.inf>self.sup:
            print 'error, interval is inverted '+'__eq__'
        if (self.inf == other.inf) & (self.sup == other.sup):
            return True
        else:
            return False
        
    def __lt__(self, other):
        if self.inf>self.sup:
            print 'error, interval is inverted '+'__lt__'
        try:
            return self.sup < other.inf
        except:
            return self.sup < other
            
    def __gt__(self, other):
        if self.inf>self.sup:
            print 'error, interval is inverted '+'__gt__'
        try:
            return self.inf > other.sup
        except:
            return self.inf > other
            
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
        if self.inf>self.sup:
            print 'error, interval is inverted '+'__contains__'
        if isinstance(other, ia):
            if other.inf>=self.inf and other.sup<=self.sup:return True
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
        if self.inf>self.sup:
            print 'error, interval is inverted '+'contains'
        result = None
        if (self.inf <= number) and (number<=self.sup):
            result = True
        elif self.inf == number:
            result = True
        elif self.sup == number:
            result = True
        else:
            result = False
        return result
        
                  
#from automatic_differentiation import ad
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