#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue May 23 18:33:46 2017

@author: luke
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
    def __init__(self, value, grad=None, 
                 op=None, args=None,name=None):
        """ """
        if op is None and args is None:
            self.value = value
            self.grad = 1.
            self.op = Rad.__init__
            self.args = []
        else:
            self.value = value
            self.grad = 1.
            self.op = op
            self.args = args
        if name is None: 
            Rad.count += 1
            self.name = 'v'+str(Rad.count)
        else:
            self.name = name
    
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
                           op = Rad.__add__, 
                           args=[self, other])
            else:
                ret = self.value + other
                return Rad(ret, 
                           op = Rad.__add__, 
                           args=[self, 0.])
                
        else:
            return self.Dadd(other)
    
    
    def __radd__(self, other, grad=False):
        """+"""
        if not grad:
            ret = self.value + other
            return Rad(ret, 
                       op = Rad.__add__, 
                       args=[self, 0.])
                
        else:
            return self.Dadd(other)
    
    
    def Dadd(self, other):
        """d+"""
        if isinstance(other, Rad):
            ret = self.grad + other.grad
            return Rad(ret, 
                       op = Rad.Dadd, 
                       args=[self, other])
        else:
            ret = self.grad + 0.
            return Rad(ret, 
                       op = Rad.Dadd, 
                       args=[self, 0.])
            
        
    #@staticmethod
    def __mul__(self, other, grad=False):
        """*"""
        if not grad:
            if isinstance(other, Rad):
                ret = self.value * other.value
                return Rad( ret, 
                           op = Rad.__mul__, 
                           args=[self, other])
            else:
                ret = self.value * other
                return Rad(ret, 
                           op = Rad.__mul__, 
                           args=[self, 0.])
        else:
            return self.Dmul(other)
        
    
    def Dmul(self, other):
        """D*"""
        if isinstance(other, Rad):
            ret = self.value*other.grad + self.grad*other.value
            return Rad( ret, 
                       op = Rad.Dmul, 
                       args=[self, other])
        else:
            ret = self.value*0.0 + self.grad*other
            return Rad( ret, 
                       op = Rad.Dmul, 
                       args=[self, 0.])
            
       
        
class BackPropagator(object):
    """On a given level
    The derivative of the levels func
    with respect to the levels primatives
    is always such that a 
    derivative of the levels primatives
    w/r/t themselves is always 1.
    https://stats.stackexchange.com/questions/224140/step-by-step-example-of-reverse-mode-automatic-differentiation
    """
    def __call__(self, fwdnode, val=0.):
        return self.backprop(fwdnode, val)
        
    @staticmethod
    def backprop(fwdnode, val=0.):
        if isinstance(fwdnode, Rad):
            #print 'fwdnode = ',fwdnode.name
            op = fwdnode.op
            args = fwdnode.args
            if len(args)>0:
                #print args
                #print args[0].name, args[1].name
                st = op(args[0],args[1],grad=True)
                val1 = st#op(args[1],grad=True)
                #print 'op i =',st.value, st.name
                st1=0.
                for el in args:
                    #print 'backproping el ',el.name
                    #if isinstance(BackPropagator.backprop(el,val1), Rad):
                        #print el.name, el.value
                    st1 += BackPropagator.backprop(el,val1)
                    
                val = val1 + st1
                #val = st1
                #val = val1
            else:
                val = fwdnode.grad #+ val
                #st = op(args[0],args[1],grad=True)
                #val += st#op(args[1],grad=True)
            return val
        else:
            return 0.
        
        
        
def print_graph(node):
    #print node.name
    if len(node.args)>0:
        for el in node.args:
            print_graph(el)
        print 'nodes = [{},{}]'.format(node.args[0].name,
                        node.args[1].name)
    return

def pg(node):
    """
        TODO: use dict implementation of graph??
    """
    #print node.name
    if len(node.args)>0:
        #for el in node.args:
        picture.append( node.name)
        picture.append( 'nodes = [{},{}]'.format(
                pg(node.args[0]),
                pg(node.args[1])
            ))
        for el in node.args:
            return pg(el)#node.args[0])
    else:    
        picture = [node.name]
        return picture
        
         
from interval_arithmetic import ia   
if __name__ == "__main__":
    
    bp = BackPropagator()
    
    a = Rad(3., name = 'a')
    b = Rad(4., name = 'b')
    c = Rad(5., name = 'c')
    d = a+b
    e = d+c
    f = (d*e)+a
    f.name = 'f'
    e.name = 'e'
    d.name = 'd'
    
    a.grad = np.matrix(np.asarray([1.,0.,0.]).T)
    b.grad = np.matrix(np.asarray([0.,1.,0.]).T)
    c.grad = np.matrix(np.asarray([0.,0.,1.]).T)
    
    g1 = a*a + a*b + a*c
    g1.name = 'g1'
    r1 = bp(g1)  #[15.,3.,3.]
    
    g2 = a*a
    g2.name = 'g2'
    r2 = bp(g2)  #[6,0,0]
    
    g3 = a+b+c
    r3 = bp(g3) #[1,1,1]
    
    g4 = a*b*c
    r4 = bp(g4) #[20.,15.,12.]
    
    fwdnode = g1
    val = 0.
    op = fwdnode.op
    args = fwdnode.args
    if len(args)>0:
        print args[0].name, args[1].name
        st = op(args[0],args[1],grad=True)
        val += st#op(args[1],grad=True)
        print 'op i =',st.value,st.name
    #    for el in args:
    #        st = BackPropagator.backprop(el,val)
    #        val += st#BackPropagator.backprop(el,val)
    #        print 'el i = ',st.value, st.name
    #        #BackPropagator.backprop(el,val)