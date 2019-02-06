#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 23 18:33:46 2017

@author: luke

building this reverse mode AD tool
is how I learned to 
write the tiny DSL for miniKanren

sources:

from xKanren docs:
    
    For the specifics of how to do it,    
    Bruce Christianson's paper: 'automatic hessians by reverse accumulation'
    helped me decide how to structure the computational graph building
    via returned/stored stuff as you see in PStates.
    
    Some nice tutorials on reverse mode AD helped me more quickly
    check 'this kind of thinking worked for that'.  
    Building confidence to move to other use cases.  ;)  
        e.g. 
        https://justindomke.wordpress.com/2009/03/24/a-simple-explanation-of-reverse-mode-automatic-differentiation/
        https://rufflewind.com/2016-12-30/reverse-mode-automatic-differentiation
        
        better?:
            https://stats.stackexchange.com/questions/224140/step-by-step-example-of-reverse-mode-automatic-differentiation
            http://www.columbia.edu/~ahd2125/post/2015/12/5/




some facinating(!) parallels across many fields:
    
    http://math.ucr.edu/home/baez/networks/networks_3.html
    (start here:
        http://math.ucr.edu/home/baez/networks/networks_1.html)
        
    seems networks are pretty 'well connected'... I couldn't resist

        
"""
import numpy as np
import copy
import inspect


class adObjectMaker(object):
    #ia_null = ia(0.,0.)
    #Im = np.matrix(np.identity((N),float))
    
    @staticmethod
    def makeGradient(N,i):
        """function to
        create a gradient of interval components
        """
        ia_id = ia(1.,1.)
        ia_null = ia(0.,0.)
        if i>= 0:
            Im = np.identity((N),float)#np.matrix(np.identity((N),float))
        elif i==-1:#interval contraint -> gradient is null
            Im = np.zeros((N,N),float)
        GX=[]
        for j in range(N):
            if j==i:
                GX.append(ia_id*Im[i,j])
            else:
                GX.append(ia_null*Im[i,j])
        GX = np.matrix(GX)
        return GX
        
    @staticmethod
    def makeHessian(N):
        """function to
        create a Hessian of interval components
        """
        ia_null = ia(0.,0.)
        m1 = np.zeros((N,N),float)#np.matrix(np.zeros((N,N),float))
        HX=[]
        for i in range(N):
            HX.append([])
            for j in range(N):
                HX[i].append(ia_null*m1[i,j])
        HX = np.matrix(HX)#np.asarray(HX)#
        return HX
      
    @staticmethod
    def scalarGradient(N,i):
        """function to
        create a gradient of real components
        """
        if i>= 0:
            Im = np.identity((N),float)#np.matrix(np.identity((N),float))
        elif i==-1:#interval contraint -> gradient is null
            Im = np.zeros((N,N),float)
        GX=[]
        for j in range(N):
            if j==i:
                GX.append(Im[i,j])
            else:
                GX.append(Im[i,j])
        GX = np.matrix(GX)
        return GX
        
    @staticmethod
    def scalarHessian(N):
        """function to
        create a Hessian of real components
        """
        m1 = np.zeros((N,N),float)#np.matrix(np.zeros((N,N),float))
        HX=[]
        for i in range(N):
            HX.append([])
            for j in range(N):
                HX[i].append(m1[i,j])
        HX = np.matrix(HX)#np.asarray(HX)#
        return HX
        


class Node(object):
    def __init__(self, op, value, args):
        self.op = op
        self.value = value
        self.args = args
        
    #def __str__(self):
    #    return str(self.value.value)
    
    #def __repr__(self):
    #    return self.value



class Graph(object):

    def add_node(self, node):
        return 


class Rad(object):
    count = 0
    def __init__(self, 
                 value,  
                 args=[],
                 N=None,
                 i=None,
                 grad=None,
                 op=None, dop=None, 
                 name=None):
        """ """
        self.value = value
        if op is None:
            self.op = Rad.__init__
            self.dop = Rad.__init__
            self.args = []
        else:
            self.op = op
            self.dop = dop
            self.args = args
        if name is None: 
            Rad.count += 1
            self.name = 'v'+str(Rad.count)
        else:
            self.name = name
        if grad is None:
            if N is not None:
                self.grad = np.matrix(np.zeros(N)).T
                self.N = N
                if i is not None:
                    self.grad[i] = 1.0
            else:
                N=1
                self.grad = np.matrix(np.zeros(N)).T
                self.N = N
        else:
            self.N = np.size(grad)
            self.grad = grad
            
    
    def __str__(self, level=0):
        #return "{}".format(self.value)
        ret = "   "*level+repr(self.name)+""+self.op.__doc__+"\n"
        for el in self.args:
            ret += el.__str__(level+1)
        return ret
    
    def __repr__(self):
        return 'Rad({})'.format(str(self.value))
    
    #@staticmethod
    def __add__(self, other, grad=False):
        """+"""
        if not grad:
            if isinstance(other, Rad):
                ret = self.value + other.value
                return Rad(ret, 
                           N = self.N,
                           op = Rad.__add__, 
                           dop = Rad.Dadd,
                           args=[self, other])
            else:
                ret = self.value + other
                return Rad(ret, 
                           op = Rad.__add__, 
                           dop = Rad.Dadd,
                           args=[self, 0.])
                
        else:
            return self.Dadd(other)
    
    
    def __radd__(self, other, grad=False):
        """+"""
        return self.__add__(other, grad)
    
    
    def Dadd(self,other,which):
        """d+"""
        #return 1. #type error here
        return Rad(1., name = '-', N=3, i=0) #fix..
        # how far can you push this?  -unknown
    
    
        
    #@staticmethod
    def __mul__(self, other, grad=False):
        """*"""
        if not grad:
            if isinstance(other, Rad):
                ret = self.value * other.value
                return Rad( ret, 
                           N = self.N,
                           op = Rad.__mul__, 
                           dop = Rad.Dmul,
                           args=[self, other])
            else:
                ret = self.value * other
                return Rad(ret, 
                           N = self.N,
                           op = Rad.__mul__, 
                           dop = Rad.Dmul,
                           args=[self, 0.])
        else:
            #if isinstance(other, Rad):
            return self.Dmul(other) 
            #else:
            #    return self.Dmul()
    def __rmul__(self, other, grad=False):
        """*"""
        return self.__mul__(self,other)
        #        if not grad:
        #            if isinstance(other, Rad):
        #                ret = self.value * other.value
        #                return Rad( ret, 
        #                           N = self.N,
        #                           op = Rad.__mul__, 
        #                           dop = Rad.Dmul,
        #                           args=[self, other])
        #            else:
        #                ret = self.value * other
        #                return Rad(ret, 
        #                           N = self.N,
        #                           op = Rad.__mul__, 
        #                           dop = Rad.Dmul,
        #                           args=[self, 0.])
        #        else:
        #            #if isinstance(other, Rad):
        #            return self.Dmul(other) 
        #            #else:
        #            #    return self.Dmul()
        #        
    
    def Dmul(self, other, which):
        """D*"""
        if which == 0:
            return other
        else:
            return self
        
        
        
        
class BackPropagator(object):
    """On a given level
    The derivative of the levels func
    with respect to the levels primatives
    is always such that a 
    derivative of the levels primatives
    w/r/t themselves is always 1.
    https://stats.stackexchange.com/questions/224140/step-by-step-example-of-reverse-mode-automatic-differentiation
    
    The general rule is to 
    sum over all possible paths 
    from one node to the other, 
    multiplying the derivatives 
    on each edge of the path together. 
    http://colah.github.io/posts/2015-08-Backprop/#fnref1
    """
    def __call__(self, fwdnode, val=0.):
        return self.backprop(fwdnode, val)
        
    @staticmethod
    def backprop(fwdnode, val=0., grad=np.matrix(np.asarray([0.,0.,0.])).T):
        dop = fwdnode.dop
        args = fwdnode.args
        if len(args)>0:
            #val = val + op(args[0],args[1],grad=True)
            val0 = val + dop(args[0],args[1],0)
            val1 = val + dop(args[0],args[1],1)
            #acc = 0.
            #for el in args:
            acc0 = BackPropagator.backprop(args[0], val)
            acc1 = BackPropagator.backprop(args[1], val)
            val0 = val0 * acc0
            val1 = val1 * acc1
            return val0+val1
        else:
            val = fwdnode.grad
            return val
        
        

         
from interval_arithmetic import ia   
if __name__ == "__main__":
    
    bp = BackPropagator()
    
    a = Rad(3., name = 'a', N=3, i=0)
    b = Rad(4., name = 'b', N=3, i=1)
    c = Rad(5., name = 'c', N=3, i=2)
    d = a+b
    e = d+c
    f = (d*e)+a
    f.name = 'f'
    e.name = 'e'
    d.name = 'd'
    
#    
    g2 = a*a
    g2.name = 'g2'
    r2 = bp(g2)  #[6,0,0]
    assert(np.any(r2.value.T == [6.,0.,0.]))
    
    g3 = a+b+c
    r3 = bp(g3) #[1,1,1]
    assert(np.any(r3.value.T == [1.,1.,1.]))
    
    g4 = a*b*c
    r4 = bp(g4) 
    assert(np.any(r4.value.T == [20.,15.,12.])),'Multiplication check r4'
    
    g1 = a*a + a*b + a*c
    g1.name = 'g1'

    r1 = bp(g1)  
    assert(np.any(r1.value.T == [15.,3.,3.]))
    
    g2 = a*a
    g2.name = 'g2'
    r2 = bp(g2)  
    assert(np.any(r2.value.T == [6.,0.,0.]))
    
    g3 = a+b+c
    r3 = bp(g3) 
    assert(np.any(r3.value.T == [1.,1.,1.]))
    #note you have a type issue here!
    # fixed quickly Feb 23, 2018 (8 months after the fact???)
    # anyway - if this ever needs to be right, come back check it!
    # (algorithmically)
    
    g4 = a*b*c
    r4 = bp(g4) 
    assert(np.any(r4.value.T == [20.,15.,12.]))
    
    print 'so for instance,'
    print 'when we have g4 = a*b*c'
    print 'g4 = {}'.format(g4.value)
    print ''
    print 'and computational graph looks like:'
    print 'g4 = {}'.format(g4)
    print 'The gradient of g4 is r4 = bakprop(g4)'
    print ' r4 = '
    print '{}'.format(r4.value)
    
    #g5 = 1.*a #okay don't push it -- Feb 2018, this puppy is limited.
    #assert(one day this will be fixed)