##
## Automatic Differentiation
##
## Luke McCulloch
## Spetember 25 2015
##
import numpy as np
import matplotlib.pyplot as plt
import sys
import copy

from time_analysis  import Timer

"""  Import Note: 

from interval_analysis import ia 

        is done below the class."""
#"""
import os
try:
    os.environ['PATH'] = os.environ['PATH'] + ':/usr/texbin'
except:
    pass
#"""

class Canvas(object):
    def __init__(self):
        self.fig            = plt.figure(1)
        self.ax             = self.fig.add_subplot(111)
        self.fig.add_subplot(self.ax)
        return
        
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
        
    
    
        
        

class ad(object):
    """
    TODO: document that this uses that idiom
    where a class instance's methods return
    instances of that class again.
    
    Note:  
        of_scalars=True     : means you have a vector of real numbers
        of_scalars=False    : means you have a vector of interval numbers
        
    N   = the total number of numbers in your design vector, or 
            the number of interacting scalars in your automatic differentiation space
            
    dim = the position of 'this number' within the vector of N numbers
    """
    def __init__(self, value, grad=None, hess=None,
                 name='unknown_var',
                 of_scalars=False,
                 N=None, dim=None):
        self.value  = value
        self.name   = name ## TODO: can we introspect pythons variable names?
        if N is None:
            assert(grad is not None and hess is not None)
            self.grad   = grad
            self.hess   = hess
            if not isinstance(grad, float):
                self.n      = grad.size
        if N is not None:
            self.n = N
            assert(dim is not None)
            if of_scalars == False:
                self.grad = adObjectMaker.makeGradient(N,dim)
                self.hess = adObjectMaker.makeHessian(N)
            elif of_scalars == True:
                self.grad = adObjectMaker.scalarGradient(N,dim)
                self.hess = adObjectMaker.scalarHessian(N)
        return
        
    def __call__(self):
        return self.value, self.grad, self.hess
        
    def __str__(self):
        return 'ad({})'.format(self.value)
    
    def __repr__(self):
        return 'ad({})'.format(self.value)
    
    
    
    def roundedreturn(self,p):
        try:
            return round(self.value,p)
        except:
            return self.value.roundedreturn(p)
            


    def eval_monoticity(self):
        if isinstance(self.grad[0,0], ia):
            decreasing = [el[0,0].sup<0. for el in self.grad.T]
            increasing = [el[0,0].inf>0. for el in self.grad.T]
            flat = [ t1==False and t2==False for t1, t2 in zip(decreasing, increasing)]
            return decreasing, increasing, flat
        else:
            print 'This gradient has no interval properties'
            return None
    
    def __neg__(self):
        return ad(-1.*self.value, -1.*self.grad, -1.*self.hess)

        
    def __add__(self,other):
        try:
            return ad(self.value + other.value  , self.grad + other.grad, self.hess + other.hess)
        except:
            try:
                return ad(self.value + other  , self.grad, self.hess)
            except:
                print 'addition failed - check types'
                print self, other
                sys.stdout.flush()
                
    def __radd__(self,other):
        #return ad(self.value + other  , self.grad, self.hess)
        return self.__add__(other)
    
    
#        try:
#            return ad(self.value + other.value  , self.grad + other.grad, self.hess + other.hess)
#        except:
#            try:
#                return ad(self.value + other  , self.grad, self.hess)
#            except:
#                print 'addition failed - check types'
#                print self, other
#                sys.stdout.flush()
            

    def __sub__(self,other):
        try:
            return ad(self.value - other.value  , self.grad - other.grad, self.hess - other.hess)
        except:
            try:
                return ad(self.value - other  , self.grad, self.hess)
            except:
                print ''
                print '\nSubtraction failed - check types'
                print 'data: ', self, other
                sys.stdout.flush()

    def __rsub__(self,other):
        return self.__sub__(other)
        #return ad(other-self.value  , -self.grad, -self.hess)
    
    
    def __mul__(self, other):
        try:
            return ad(self.value*other.value,
                     np.dot(self.grad,other.value) + np.dot(other.grad,self.value),
                     np.dot(self.hess,other.value) + np.dot(np.transpose(self.grad),other.grad) + \
                     np.dot(np.transpose(other.grad),self.grad) + np.dot(other.hess,self.value))
        except:
            return ad(self.value*other, self.grad*other, self.hess*other)
            
    def __rmul__(self, other):
        return self.__mul__(other)
        
    def not_quite__div__(self, other):
        try: 
            return ad((self.value/other.value),
                     (np.dot(other.value,self.grad) - np.dot(self.value,other.grad)) /(other.value**2),
                     (np.dot((other.value**2) , ((np.dot(self.hess,other.value) + 
                     np.dot(self.grad.T,other.grad) ) - \
                     (np.outer( other.grad.T,self.grad) + np.dot(other.hess,self.value)))   ) -
                     (np.dot(  (np.dot(other.value,self.grad) - np.dot(self.value,other.grad) ).T  , \
                     2.*np.dot(other.grad,other.value) )  )) / (other.value**4)
                     )
        except:
            return ad((self.value/other)  ,
                     (other*self.grad)/(other**2),
                     ((other**2)*other*self.hess)/other**4)
                     
    def __div__(self, other):
        try:
            return ad((self.value/other.value)  ,
                     (np.dot(other.value,self.grad) - np.dot(self.value,other.grad))/(other.value**2),
                     (np.dot( (other.value**2) , ( (np.dot(self.hess,other.value) + 
                     np.dot( self.grad.T,other.grad) ) - \
                     (np.dot(other.grad.T,self.grad) + np.dot(other.hess,self.value)))  ) - \
                     ( (np.dot(other.value,self.grad) - np.dot(self.value,other.grad)).T * \
                     2.*np.dot(other.grad,other.value)) ) /(other.value**4)
                     )
        except:
            return ad((self.value/other)  ,
                     (other*self.grad)/(other**2),
                     ((other**2)*other*self.hess)/other**4) 
                     
                         
    def SAVEthis__div__(self, other):
        try:
            return ad((self.value/other.value)  ,
                     (np.dot(other.value,self.grad) - np.dot(self.value,other.grad))/(other.value**2),
                     (np.dot( (other.value**2) , ( (np.dot(self.hess,other.value) + 
                     np.outer( self.grad.T,other.grad) ) - \
                     (np.dot(other.grad.T,self.grad) + np.dot(other.hess,self.value)))  ) - \
                     ( (np.dot(other.value,self.grad) - np.dot(self.value,other.grad)).T * \
                     2.*np.dot(other.grad,other.value)) ) /(other.value**4)
                     )
        except:
            return ad((self.value/other)  ,
                     (other*self.grad)/(other**2),
                     ((other**2)*other*self.hess)/other**4) 
                     
    def __rdiv__(self, other):
        return self.__div__(other)
        #return ad((other/self.value)  ,
        #             (-other*self.grad)/(other**2),
        #             (-(other**2)*other*self.hess)/other**4) 
                     
    def sin(self):
        return ad(np.sin(self.value)  , 
                        np.cos(self.value)*self.grad, 
                        -np.sin(self.value)*np.transpose(self.grad)*self.grad + np.cos(self.value)*self.hess )

    def cos(self):
        return ad(np.cos(self.value)  , 
                        -np.sin(self.value)*self.grad , 
                        -np.cos(self.value)*np.transpose(self.grad)*self.grad - np.sin(self.value)*self.hess)
    #"""

    def tan(self):
        return ad(np.tan(self.value)  , self.grad/((np.cos(self.value))**2),
                     ( ((np.cos(self.value))**2)*self.hess - (np.transpose(self.grad)*(2.*np.cos(self.value))*(-np.sin(self.value))*self.grad) )/
                     ((np.cos(self.value))**4)
                     )


    def asin(self): #hessian not really checked!
        return ad(np.arcsin(self.value) ,
                        self.grad/np.sqrt(1.-(self.value*self.value)),
                        ( (np.sqrt(1.-(self.value*self.value)))*(self.hess) - \
                          (np.transpose(self.grad))*((-2.*self.value*self.grad)*(((1.-(self.value*self.value))**(-1./2.))*0.5)) /
                          ((np.sqrt(1.-(self.value*self.value)))**2) )
                        )

    def acos(self): #hessian not really checked!
        return ad(np.arccos(self.value) ,
                        -self.grad/np.sqrt(1.-(self.value*self.value)),
                        ( (np.sqrt(1.-(self.value*self.value)))*(-self.hess) - \
                          (-np.transpose(self.grad))*((-2.*self.value*self.grad)*(((1.-(self.value*self.value))**(-1./2.))*0.5)) /
                          ((np.sqrt(1.-(self.value*self.value)))**2) )
                        )

    def atan(self):
        return ad(np.arctan(self.value) ,
                     self.grad/(1.+self.value*self.value),
                     ( ((1.+self.value*self.value)*self.hess)-(np.transpose(self.grad)*(self.grad*self.value*2.)) ) /(((1.+self.value*self.value))**2)
                     )
    
    def arctan2(self, other): #quandrant aware
        store = self/other
        return ad(np.arctan2(self.value,other.value) ,
                     store.grad/(store.value*store.value +1.),
                     ( \
                     np.dot((store.value*store.value + 1.),store.hess) - \
                     np.dot(np.transpose(store.grad),np.dot(store.grad,store.value*2.))  
                     ) \
                     /(((1.+store.value*store.value))**2)
                     )

    def exp(self, invert = False):
        #if invert==False:
        return ad(np.exp(self.value)  , self.grad*np.exp(self.value), self.hess*np.exp(self.value) + np.transpose(self.grad)*self.grad*np.exp(self.value))
        #else:
        #    return self.log()

    def log(self, invert = False):
        #if invert==False:
        return ad(np.log(self.value)  , self.grad/self.value,
                    (self.value*self.hess - np.transpose(self.grad)*self.grad)/(self.value*self.value))
        #else:
        #    return self.exp()
    def sqrt(self):
        dummy = copy.deepcopy(self)
        if dummy.value < 0.:
            print 'warning, non positive square root encountered in scalarAD'
            print 'TODO: make sure check if less than 0. works for interval and float'
            if isinstance(dummy.value, float):
                dummy.value = 0.
            elif isinstance(dummy.value, ia):
                dummy.value.inf = 0.
                if dummy.value.sup <0.:
                    dummy.value.sup = 0.
        return ad(np.sqrt(dummy.value) , 0.5*dummy.grad/np.sqrt(dummy.value),
                     ( np.dot(np.sqrt(dummy.value),0.5*dummy.hess) - \
                     (0.5*np.transpose(dummy.grad)*0.5*dummy.grad/np.sqrt(dummy.value)))/dummy.value
                     )
                     
    def __abs__(self):
        return ad(abs(self.value),self.grad,self.hess)
    
    def __pow__(self,n):
        assert (type(n) == int or type(n) == float); #"n should be scalar valued."
        return ad(self.value**n,
                     np.dot( (self.value**(n-1))*n, self.grad ),
                     np.dot( (self.grad.T), np.dot( (self.value**(n-2))*(n*n-n) , self.grad ) ) +\
                     np.dot( (self.value**(n-1))*n, self.hess ) #+ \
                     
                     )
        
    def __eq__(self,other):
        try:
            return self.value == other.value
        except:
            return self.value == other

    def __lt__(self, other):
        try:
            return self.value < other.value
        except:
            return self.value < other
            
    def __gt__(self, other):
        try:
            return self.value > other.value
        except:
            return self.value > other
        
    def __le__(self, other):
        try:
            return self.value <= other.value
        except:
            return self.value <= other
    
    def __ge__(self, other):
        try:
            return self.value >= other.value
        except:
            return self.value >= other    
    
    def print_components(self):
        print '\nComponents of {}:'.format(self.name)
        print '{}.value ='.format(self.name)
        print '           ',self.value
        print '{}.grad  ='.format(self.name)
        for i in range(self.grad.size):
            print '           ', self.grad[0,i]
        print '{}.hess  ='.format(self.name)
        for i in range(self.grad.size):
            print '           ', self.hess[i]
        print ''
        return
    
    @staticmethod
    def get_which_dimension(x):
        not_found = -1
        for i in range(np.size(x.grad)):
            if x.grad[0,i].contains(1.):
                return i
        return not_found
            
            
    @staticmethod
    def plot(x, y, canvas = None, get = True, label_x='x', label_y = 'y'):
        """TODO:
        Nice plotting of inerval boxes and their ranges
        a la Jaulin's contractor videos
        This should be in ad and not ia 
        since ad naturally specifies the structure of the space
        """
        if canvas is None:
            canvas = Canvas()
            canvas.max_ = max(1.1*x.value.sup, 1.1*y.value.sup)
            canvas.ax.arrow(0., 0., canvas.max_, 0., head_width=0.05, head_length=0.1, fc='k', ec='k')
            canvas.ax.arrow(0., 0., 0., canvas.max_, head_width=0.05, head_length=0.1, fc='k', ec='k')
            name = str('x')
            plt.annotate(name, xy = (canvas.max_,-.1*canvas.max_), xycoords='data')
            name = str('y')
            plt.annotate(name, xy = (-.1*canvas.max_,canvas.max_), xycoords='data')
            
        #x_dim = ad.get_which_dimension(x)
        #y_dim = ad.get_which_dimension(y)
        #if x_dim < y_dim:
        if isinstance(x,ad):
            canvas.ax.plot(x.value.inf, 0., marker = 'o', color = 'black')
            canvas.ax.plot(x.value.sup, 0., marker = 'o', color = 'black')
            canvas.ax.plot(0., y.value.inf, marker = 'o', color = 'black')
            canvas.ax.plot(0., y.value.sup, marker = 'o', color = 'black')
            canvas.ax.plot([x.value.inf,
                            x.value.sup,
                            x.value.sup,
                            x.value.inf,
                            x.value.inf], 
                           [y.value.inf,
                            y.value.inf,
                            y.value.sup,
                            y.value.sup,
                            y.value.inf], marker = 'o', color = 'black')
        if isinstance(x,ia):
            canvas.ax.plot(x.inf, 0., marker = 'o', color = 'red')
            #name = label_x+'.inf'
#            canvas.ax.annotate(name, xy=(x.inf, 0.),  xycoords='data',
#                xytext=(60, -30), textcoords='offset points',
#                arrowprops=dict(arrowstyle="->")
#                )
            
            #name = label_x+'.sup'
            canvas.ax.plot(x.sup, 0., marker = 'o', color = 'red')
#            canvas.ax.annotate(name, xy=(x.sup, 0.),  xycoords='data',
#                xytext=(40, -20), textcoords='offset points',
#                arrowprops=dict(arrowstyle="->")
#                )
            canvas.ax.plot(0., y.inf, marker = 'o', color = 'red')
            canvas.ax.plot(0., y.sup, marker = 'o', color = 'red')
            canvas.ax.plot([x.inf,
                            x.sup,
                            x.sup,
                            x.inf,
                            x.inf], 
                           [y.inf,
                            y.inf,
                            y.sup,
                            y.sup,
                            y.inf], marker = 'o', color = 'red')
            
        plt.xlim(-1.,canvas.max_+1.)
        plt.ylim(-1.,canvas.max_+1)
        if get==True:
            return canvas
        else:
            return
  
class VectorAD(object):
    """A class for arrays of ad points
    
    -not for arrays nor matrices of ad.grad noror ad.hess 
    of a scalar objective function.
    """
    
    def __init__(self, X):
        assert(isinstance(X, list))
        self.X = np.asarray(X)
        self.n = len(self.X)
        return
    
    @staticmethod
    def generate_mid_vector(X, loc=.5):
        """given a vector of interval components
        return a scalar component vector
        """
        #        mv=[]
        #        for i in range(len(x)):
        #            mv.append(x[i].getpoint(loc))
        #        mv=np.asarray(mv)
        #        return mv
        ## new way:
        xm=copy.deepcopy(X)
        for i in range(len(X)):
            xm[i].value = X[i].value.getpoint(loc)
        return xm
    
    def generate_gradient_from_system():
        """This would be needed were we not using
        Lagrangians.  Had we instead been solving 
        unrelated systems of equations
        
        This might be needed if doing vector
        potentials though
        """
        print 'Not implemented'
        return 'You want a Lagrangian, not this.'

    
    def generate_Hessian_from_system():
        """This would be needed were we not using
        Lagrangians.  Had we instead been solving 
        unrelated systems of equations
        
        This might be needed if doing vector
        potentials though
        """
        print 'Not implemented'
        return 'You want a Lagrangian, not this.'
    
    @staticmethod
    def generate_numpy_repr(X):
        """Take an interval valued matrix
        e.g. a Hessian.
        Return the matrix contraining
        same information (in the same spot)
        """
        n = X.n
        X.igrad = []
        X.ihess = []
        for i in range(n):
            X.igrad.append( [ X.grad[0,i].inf, X.grad[0,i].sup] )
            X.ihess.append([])
            for j in range(n):
                X.ihess[i].append( [ X.hess[i,j].inf, X.hess[i,j].sup ] )
        X.igrad = np.asarray(X.igrad)
        X.ihess = np.asarray(X.ihess)
        return X

 
     

class IntervalAnalysis(object):
    """A class to 
    do interval linear algebra
    """
    def __init__(self, mask = None):
        self.mask = mask
        self.n = len(mask)
        self.check_box = []
        self.has_mp = True
        return
    
    @staticmethod
    def mid_vector(x, loc = .5):
        mv=[]
        for i in range(len(x)):
            mv.append(x[i].value.getpoint(loc))
        mv=np.asarray(mv)
        return mv
        
    @staticmethod
    def generate_mid_gradient(G, size, loc = .5):
        """generate scalar component
        mid vectors 
        from an interval valued gradient
        
        input: 
            F.grad (interval component vector)
        output:
            mid(F.grad) (scalar component vector)
        """
        MG = []
        for i in range(size):
            MG.append(G[0,i].getpoint(loc))
        MG = np.asarray(MG)
        return MG
        
    @staticmethod
    def generate_thin_gradient(G, size):
        """generate interval component
        mid gradient array
        from an interval valued gradient matrix
        
        input: 
            F.grad (interval component matrix vector)
        output:
            mid(F.grad) (interval component array vector)
        """
        MG = []
        for i in range(size):
            MG.append(G[0,i])
        MG = np.asarray(MG)
        return MG
    
    @staticmethod
    def generate_mid_hessian(H, loc = .5):
        """generate mid matrices from hessians of Lagrangians
        
        input: 
            F.hess (interval component matrix)
        output:
            mid(F.hess) (scalar component matrix)
        """
        MXm = []
        for i in range(len(H)):
            MXm.append([])
            for j in range(len(H)):
                MXm[i].append(H[i,j].getpoint(loc))
        MXm = np.asarray(MXm)
        return MXm
        
    
    @staticmethod
    def guass_sidel(A,X,B, nmax=1, tol = 1.e-6):
        """
            A is a square interval matrix
            B is a column vector RHS
            X is a column vector of unknowns
        
        See e.g. Snyder 1992, Gen Model for CAD, page 92
        if ever Z ( X_n, then f has a unique zero in X_n.
        """
        cnvg  = 0
        count = 0
        midA  = IntervalAnalysis.generate_mid_hessian(A)
        try:
            Y     = np.linalg.inv(midA)    
            G     = np.dot(Y,A)
            C     = np.dot(Y,B)
        except:
            print 'no precondition possible'
            G = A
            C = B
        
        size = len(X)
        Xk = copy.deepcopy(X)
        Xold = copy.deepcopy(X)
        while((cnvg < size) and count < nmax):
            Xold = Xk
            for i in range(G.shape[0]):
                s1 = np.dot(G[i, :i], Xk[:i])
                s2 = np.dot(G[i, i + 1:], Xold[i + 1:])
                if (s1 is not None) and (s2 is not None):
                    #Z = ((G[i, i])**(-1))*(C[i] - s1 - s2)
                    Z = (C[i] - s1 - s2)/(G[i, i])
                elif (s1 is not None):
                    s2 = 0.
                    #Z = ((G[i, i])**-1)*(C[i] - s1)
                    Z = (C[i] - s1)/(G[i, i])
                elif (s2 is not None):
                    s1 = 0.
                    #Z = ((G[i, i])**-1)*(C[i] - s2)
                    Z = (C[i] - s2)/(G[i, i])
                else:
                    assert((s1 is not None)or(s2 is not None))

                test_update = Xold[i] & Z       
                #for i, el in enumerate(test_update):
                #if abs(test_update.width().value) > 0.:  #we are ditching perfect answers??
                if test_update.isempty == False:
                    Xk[i] = test_update
                else:
                    Xk[i].isempty = True
                if (abs(Xk[i].width()) < tol): #possible fail for true 0 width intervals? check and Xk[i].isempty = False??
                    cnvg += 1
                
                count += 1
    
        return Xk
        
    @staticmethod
    def split_guass_seidel(A,X,B,
                           nmax=1, tol = 1.e-6):
        """
            
            A       = this_hess
            X       = Vint.T
            B       = grad.T
            nmax    = 1
            tol     = gstol 
            i=0
            
            A is a square interval matrix
            B is a column vector RHS
            X is a column vector of unknowns
            dlist :: list containing the design space of boxes
        
        See e.g. Snyder 1992, Gen Model for CAD, page 92
        if ever Z ( X_n, then f has a unique zero in X_n.
        
        updated to spell G.S. correctly (><)
        
        TODO: capture all split domains just once
        steps:
            -get all the split inervals
            -permute them to create new boxes
        """
        #split_alreay = False
        assert(nmax==1)
        gslist = [] #newly created intervals
        #
        cnvg  = 0
        count = 0
        midA  = IntervalAnalysis.generate_mid_hessian(A)
        try:
            Y     = np.linalg.inv(midA)    
            G     = np.dot(Y,A)
            C     = np.dot(Y,B)
        except:
            print 'no precondition possible'
            G = A
            C = B
        #
        size = len(X)
        Xk = copy.deepcopy(X)
        Xold = copy.deepcopy(X)
        while((cnvg < size) and count < nmax):
            Xold = Xk
            for i in range(G.shape[0]):
                s1 = np.dot(G[i, :i], Xk[:i])
                s2 = np.dot(G[i, i + 1:], Xold[i + 1:])
                assert((s1 is not None) or (s2 is not None))
                if (s1 is not None) and (s2 is not None):
                    #if split_alreay:
                    #    Z = (C[i] - s1 - s2)/(G[i, i])
                    #else:
                    Z = (C[i] - s1 - s2).div_split(G[i, i])
                elif (s1 is not None):
                    s2 = 0.
                    #if split_alreay:
                    #   Z = (C[i] - s1)/(G[i, i])
                    #else:
                    Z = (C[i] - s1).div_split(G[i, i])
                else:# (s2 is not None):
                    s1 = 0.
                    #if split_alreay:
                    #    Z = (C[i] - s2)/(G[i, i])
                    #else:
                    Z = (C[i] - s2).div_split(G[i, i])

                #if not split_alreay:
                if len(Z)>1:
                    #split_alreay = True
                    assert(len(Z)==2)
                    #for i in range(len(Z[1:])):    
                    test_update = Xold[i] & Z[1] 
                    if test_update.isempty == False:
                        Ck = copy.deepcopy(Xk)
                        Ck[i] = test_update #Very Important! lists by refernce
                        #bx = box(Xk)
                        #bx.f = eval_system(bx)
                        gslist.append(Ck)

                test_update = Xold[i] & Z[0]
                #else:
                #    test_update = Xold[i] & Z
                #for i, el in enumerate(test_update):
                #if abs(test_update.width().value) > 0.:  #we are ditching perfect answers??
                if test_update.isempty == False:
                    Xk[i] = test_update
                else:
                    Xk[i].isempty = True
                if (abs(Xk[i].width()) < tol): #possible fail for true 0 width intervals? check and Xk[i].isempty = False??
                    cnvg += 1
                
                count += 1
    
        return Xk, gslist
      
    @staticmethod
    def split_guass_seidel2(A,X,B,
                           nmax=1, tol = 1.e-6):
        """
            
            A       = this_hess
            X       = Vint.T
            B       = grad.T
            nmax    = 1
            tol     = gstol 
            i=0
            
            A is a square interval matrix
            B is a column vector RHS
            X is a column vector of unknowns
            dlist :: list containing the design space of boxes
        
        See e.g. Snyder 1992, Gen Model for CAD, page 92
        if ever Z ( X_n, then f has a unique zero in X_n.
        
        updated to spell G.S. correctly (><)
        
        TODO: capture all split domains just once
        steps:
            -get all the split inervals
            -permute them to create new boxes
            
        idea :: division always returns a list
        """
        #split_alreay = False
        assert(nmax==1)
        #gslist = [] #newly created intervals
        #
        # better way : permute! 2^n
        isize = A.shape[0] #works b/c 1 row per dimension
        ilist = [[]for i in range(isize)]
        #
        cnvg  = 0
        count = 0
        midA  = IntervalAnalysis.generate_mid_hessian(A)
        try:
            Y     = np.linalg.inv(midA)    
            G     = np.dot(Y,A)
            C     = np.dot(Y,B)
        except:
            print 'no precondition possible'
            G = A
            C = B
        #
        size = len(X)
        Xk = copy.deepcopy(X)
        Xold = copy.deepcopy(X)
        while((cnvg < size) and count < nmax):
            Xold = Xk
            for i in range(G.shape[0]):
                s1 = np.dot(G[i, :i], Xk[:i])
                s2 = np.dot(G[i, i + 1:], Xold[i + 1:])
                assert((s1 is not None) or (s2 is not None))
                if (s1 is not None) and (s2 is not None):
                    #if split_alreay:
                    #    Z = (C[i] - s1 - s2)/(G[i, i])
                    #else:
                    Z = (C[i] - s1 - s2).div_split(G[i, i])
                elif (s1 is not None):
                    s2 = 0.
                    #if split_alreay:
                    #   Z = (C[i] - s1)/(G[i, i])
                    #else:
                    Z = (C[i] - s1).div_split(G[i, i])
                else:# (s2 is not None):
                    s1 = 0.
                    #if split_alreay:
                    #    Z = (C[i] - s2)/(G[i, i])
                    #else:
                    Z = (C[i] - s2).div_split(G[i, i])
                    
                for el in Z:
                    ilist[i].append( Xold[i] & el)
                
            results = []
            for i in range(isize-1):
                for el1 in ilist[i]:
                    for j in range(i+1,isize):
                        for el2 in ilist[j]:
                            #print i, el1, j, el2
                            Ck = copy.deepcopy(Xk)
                            Ck[i] = el1
                            Ck[j] = el2
                            results.append(Ck)
                #if (abs(Xk[i].width()) < tol): #possible fail for true 0 width intervals? check and Xk[i].isempty = False??
                #    cnvg += 1
                
            count += 1
    
        return results
      
    @staticmethod
    def some_convergence_criteria(xiold,xi,tol):
        return xi.value.width() < xiold.value.width()
        
    @staticmethod
    def compute_scalar_newton_basic(func,c,x, 
                                    itmax = 1000, 
                                    tol = 1.e-20):
        """The scalar valued interval Newton's method
        """
        xi = copy.deepcopy(x)
        i=0
        not_converged = True
        while(i<itmax and not_converged):
            xm = copy.deepcopy(xi)
            xm.value = xi.value.getpoint(.5)
            f = func(xm,c)
            F = func(xi,c)
            nx = xm - f.value/F.grad[0,0]
            
            xiold = copy.deepcopy(xi)
            xi.value = xi.value & nx.value
            
            not_converged = IntervalAnalysis.some_convergence_criteria(xiold,xi,tol)
            
            i+=1
        return xi, i
     
    
    @staticmethod
    def vector_nx_convergence_criteria(xiold,Xi,tol):
        converged = True
        for i in range(len(xiold)):
            if Xi[i].value.width() < xiold[i].value.width():
                converged = False
        return converged
    
    @staticmethod
    def check_width(X, tol):
        xw = []
        for i in range(len(X)):
            xw.append(X[i].value.width() <= tol)
        if False in xw:
            return False
        else:
            return True
        
    @staticmethod
    def compute_interval_newton_basic(func,X, itmax = 25, tol = 1.e-5, Bspline_style = False):
        
        if Bspline_style:
            x = X[0]
            y = X[1]
            m = len(x)
            x_update = []  
            for el in x:
                x_update.append(el)
            for el in y:
                x_update.append(el)
            Xi = np.asarray(x_update)
        else:
            Xi = copy.deepcopy(X)
        
        n = len(X)
        xiold = [copy.deepcopy(Xi)]
        count=0
        converged = False
        #f_upper_bound = 1.e20
        gstol = tol
        while(count<itmax and not converged):
            
            xm = VectorAD.generate_mid_vector(Xi, loc = .5)

            if Bspline_style:
                f = func( [ xm[0:m], xm[m:2*m] ] )
                F = func([ Xi[0:m], Xi[m:2*m] ])
            else:
                f = func(xm)
                F = func(Xi)
            
            #f_upper_bound = min(f_upper_bound,f.value)
            #flb_check = (F.value.inf > f_upper_bound) #if F.value.inf > f_upper_bound, eliminate the box
            
            #precondition hessian:
            this_hess = F.hess
            this_hess = np.asarray(this_hess)
            #Hm = IntervalAnalysis.generate_mid_hessian(this_hess)
            #iHm = np.linalg.inv(Hm)
            #B = iHm
            #this_hess = np.dot(B,this_hess)
            
            #precondition gradietn:
            grad = IntervalAnalysis.generate_thin_gradient(f.grad, size=f.grad.size)
            #r = np.dot(B,grad)
            #grad = r
            
            Vint = Xi - xm
            for i in range(n):
                Vint[i] = Vint[i].value
            nx = IntervalAnalysis.guass_sidel(this_hess,
                                               Vint.T,
                                               grad.T, 
                                               nmax=1,
                                               tol=gstol )
            nx = xm - nx
            for i in range(n):
                test = Xi[i].value & nx[i].value
                if test in Xi[i].value: 
                    Xi[i].value = test
                    Xi[i].possible_zero = True    
                else:
                    Xi[i].possible_zero = False
                    Xi[i].isempty       = True
            converged = IntervalAnalysis.vector_nx_convergence_criteria(xiold[-1],Xi,tol)
            small_enough = IntervalAnalysis.check_width(Xi, tol)
            if converged or small_enough:
                xiold.append( copy.deepcopy(Xi) )
                #print 'starting with'
                #print xiold[0]
                #print 'now'
                #print Xi
                Xi = copy.deepcopy(xiold[-2])
            else:
                xiold.append( copy.deepcopy(Xi) )
            count += 1
            print count

            
        if Bspline_style:
            X = [[],[]]
            for i in range(m):
                X[0].append(Xi[i])
                X[1].append(Xi[m+i])
            Xi = X
        return Xi, count, F, converged, small_enough
        
    
    @staticmethod
    def compute_interval_newton_basic_split(this,func,
                                            X,
                                            itmax = 1, 
                                            tol = 1.e-5, 
                                            Bspline_style = False):
        
        """
            func    = self.f
            X       = copy.deepcopy(X())
            itmax   = 1
            tol     = self.tol
            Bspline_style = True
        """
        assert(itmax == 1)
        ilist = [] #list of newly generated intervals
        if Bspline_style:
            x = X[0]
            y = X[1]
            m = len(x)
            #print 'm = ', m
            x_update = []  
            for el in x:
                x_update.append(el)
            for el in y:
                x_update.append(el)
            Xi = np.asarray(x_update)
        else:
            Xi = copy.deepcopy(X)
        
        n = len(X)
        #print 'n = ',n
        xiold = [copy.deepcopy(Xi)]
        count=0
        converged = False
        f_upper_bound = 1.e20
        gstol = tol
        CCi = copy.deepcopy(Xi)
        #while(count<itmax and not converged):
            
        xm = VectorAD.generate_mid_vector(Xi, loc = .5)

        if Bspline_style:
            f = func( [ xm[0:m], xm[m:2*m] ] )
            F = func([ Xi[0:m], Xi[m:2*m] ])
        else:
            f = func(xm)
            F = func(Xi)
        
        f_upper_bound = min(f_upper_bound,f.value)
        #flb_check = (F.value.inf > f_upper_bound) #if F.value.inf > f_upper_bound, eliminate the box
        
        #precondition hessian:
        this_hess = F.hess
        this_hess = np.asarray(this_hess)
        #Hm = IntervalAnalysis.generate_mid_hessian(this_hess)
        #iHm = np.linalg.inv(Hm)
        #B = iHm
        #this_hess = np.dot(B,this_hess)
        
        #precondition gradietn:
        grad = IntervalAnalysis.generate_thin_gradient(f.grad, size=f.grad.size)
        #r = np.dot(B,grad)
        #grad = r
        
        Vint = Xi - xm
        for i in range(n):
            Vint[i] = Vint[i].value
            
        nx, gslist = IntervalAnalysis.split_guass_seidel(this_hess,
                                                         Vint.T,
                                                         grad.T, 
                                                         nmax=1,
                                                         tol=gstol )
        print 'nx = ', nx
        print 'glist=',gslist
        nx = xm - nx
        for i in range(n):
            test = Xi[i].value & nx[i].value
            if test in Xi[i].value: 
                Xi[i].value = test
                Xi[i].possible_zero = True    
            else:
                Xi[i].possible_zero = False
                Xi[i].isempty       = True
        converged = IntervalAnalysis.vector_nx_convergence_criteria(xiold[-1],Xi,tol)
        small_enough = IntervalAnalysis.check_width(Xi, tol)
        if converged or small_enough:
            xiold.append( copy.deepcopy(Xi) )
            #print 'starting with'
            #print xiold[0]
            #print 'now'
            #print Xi
            #Xi = copy.deepcopy(xiold[-2])
        else:
            xiold.append( copy.deepcopy(Xi) )
        count += 1
        print count
        
        
        print '\n gslist'
        print gslist, '\n'
        for nx in gslist:
            Ci = copy.deepcopy(CCi)
            nx = xm - nx
            this.ss = [copy.deepcopy(CCi), nx, xm,[]]
            for i in range(n):
                #Ci = copy.deepcopy(CCi) #TLM original fix
                """
                CCi= self.ss[0]
                Ci = copy.deepcopy(CCi)
                nx=self.ss[1]
                xm = self.ss[2]
                #nx = xm-nx
                i=0
                """
                print i
                print Ci[i].value
                print nx[i].value
                test = nx[i].value & Ci[i].value
                this.ss[3].append((i,test))
                """
                    Key aspect ofthe algorithm:
                      only take boxes where something changed
                      once!
                      
                      here denoted by:
                      
                      and Ci[i].value not in Xi[i].value
                """
                if test in Ci[i].value and Ci[i].value not in Xi[i].value: 
                    print 'insert to ilist'
                    Ci[i].value = test
                    Ci[i].possible_zero = True  
                    if Bspline_style:
                        C = [[],[]]
                        for i in range(m):
                            C[0].append(Ci[i])
                            C[1].append(Ci[m+i])
            
                    ilist.append(C)
                    break #hack to prevent multicovers of an interval box
                    
                #else: #drop it
                #    Ci[i].possible_zero = False
                #    Ci[i].isempty       = True
                #    this.boudary.append()
        
        if Bspline_style:
            X = [[],[]]
            for i in range(m):
                X[0].append(Xi[i])
                X[1].append(Xi[m+i])
            Xi = X
        return Xi, count, F, converged, small_enough, ilist


    @staticmethod
    def compute_interval_newton_basic_split2(this,func,
                                            X,
                                            itmax = 1, 
                                            tol = 1.e-5, 
                                            Bspline_style = False):
        
        """
            func    = self.f
            X       = copy.deepcopy(X())
            itmax   = 1
            tol     = self.tol
            Bspline_style = True
        """
        assert(itmax == 1)
        ilist = [] #list of newly generated intervals
        if Bspline_style:
            x = X[0]
            y = X[1]
            m = len(x)
            #print 'm = ', m
            x_update = []  
            for el in x:
                x_update.append(el)
            for el in y:
                x_update.append(el)
            Xi = np.asarray(x_update)
        else:
            Xi = copy.deepcopy(X)
        
        n = len(X)
        #print 'n = ',n
        xiold = [copy.deepcopy(Xi)]
        count=0
        converged = False
        f_upper_bound = 1.e20
        gstol = tol
        CCi = copy.deepcopy(Xi)
        #while(count<itmax and not converged):
            
        xm = VectorAD.generate_mid_vector(Xi, loc = .5)

        if Bspline_style:
            f = func( [ xm[0:m], xm[m:2*m] ] )
            F = func([ Xi[0:m], Xi[m:2*m] ])
        else:
            f = func(xm)
            F = func(Xi)
        
        f_upper_bound = min(f_upper_bound,f.value)
        #flb_check = (F.value.inf > f_upper_bound) #if F.value.inf > f_upper_bound, eliminate the box
        
        #precondition hessian:
        this_hess = F.hess
        this_hess = np.asarray(this_hess)
        #Hm = IntervalAnalysis.generate_mid_hessian(this_hess)
        #iHm = np.linalg.inv(Hm)
        #B = iHm
        #this_hess = np.dot(B,this_hess)
        
        #precondition gradietn:
        grad = IntervalAnalysis.generate_thin_gradient(f.grad, size=f.grad.size)
        #r = np.dot(B,grad)
        #grad = r
        
        Vint = Xi - xm
        for i in range(n):
            Vint[i] = Vint[i].value
            
        
        nx_list = IntervalAnalysis.split_guass_seidel2(this_hess,
                                                         Vint.T,
                                                         grad.T, 
                                                         nmax=1,
                                                         tol=gstol)
        #print '\n nx_list'
        #print nx_list, '\n'
        #nx = copy.deepcopy(nx_list[0])
        for nx in nx_list:
            Ci = copy.deepcopy(CCi)
            nx = xm - nx
            this.ss = [copy.deepcopy(CCi), nx, xm,[]]
            for i in range(n):
                #Ci = copy.deepcopy(CCi) #TLM original fix
                """
                Ci = copy.deepcopy(CCi)
                """
                #print i
                #print Ci[i].value
                #print nx[i].value
                test = nx[i].value & Ci[i].value
                this.ss[3].append((i,test))
                if test in Ci[i].value:
                    #print 'insert to ilist'
                    Ci[i].value = test
                    Ci[i].possible_zero = True  
            if Bspline_style:
                C = [[],[]]
                for i in range(m):
                    C[0].append(Ci[i])
                    C[1].append(Ci[m+i])
    
            ilist.append(C)
            #break #hack to prevent multicovers of an interval box
                    
                #else: #drop it
                #    Ci[i].possible_zero = False
                #    Ci[i].isempty       = True
                #    this.boudary.append()
        
        small_enough = False
        return ilist, count, F, converged, small_enough
        
    
    @staticmethod
    def compute_scalar_newton(func,
                              c,
                              x, 
                              vertex_tpl, 
                              itmax = 1000, 
                              tol = 1.e-20):
        """The scalar valued interval Newton's method
        """
        xi                      = copy.deepcopy(x)
        i                       = 0
        not_converged           = True
        xyz     = vertex_tpl[0]     # which dimension are we in [0,1,2]?
        index   = vertex_tpl[1]   # which vertex are we attempting to contract?
        n       = len(x[0])
        gn      = xyz*n + index      #pick the correct element of the gradient
        while(i<itmax and not_converged):
            midx                    = xi[xyz][index].value.getpoint(.5)
            xm                      = copy.deepcopy(xi)
            xm[xyz][index].value    = midx
            f                       = func(xm,c)
            F                       = func(xi,c)
            
            nx                      = xm - f.value/F.grad[0,gn]
            
            xiold                   = copy.deepcopy(xi)
            
            xi[xyz][index].value    = xi[xyz][index].value & nx.value
            
            not_converged = IntervalAnalysis.some_convergence_criteria(xiold[xyz][index],xi[xyz][index],tol)
            
            i+=1
        return xi, i
        
    @staticmethod
    def montonic_compute_scalar_newton(func,
                                       c,
                                       x, 
                                       vertex_tpl, 
                                       itmax = 1000, 
                                       tol = 1.e-15):
        """The scalar valued interval Newton's method
        TODO: make this exploit monotonicity
        
        
            monotonicity checking box consistency
            (monotonic interval newton iteration on a single DOF)
        """
        xi = copy.deepcopy(x)
        i=0
        not_converged = True
        xyz = vertex_tpl[0]     # which dimension are we in [0,1,2]?
        index = vertex_tpl[1]   # which vertex are we attempting to contract?
        n = len(x[0])
        gn = xyz*n + index      #pick the correct element of the gradient
        while(i<itmax and not_converged):
            midx                    = xi[xyz][index].value.getpoint(.5)
            xm                      = copy.deepcopy(xi)
            xm[xyz][index].value    = midx
            f                       = func(xm,c)
            F                       = func(xi,c)
            
            nx = xm - f.value/F.grad[0,gn]
            
            xiold = copy.deepcopy(xi)
            
            xi[xyz][index].value = xi[xyz][index].value & nx.value
            
            not_converged = IntervalAnalysis.some_convergence_criteria(xiold[xyz][index],xi[xyz][index],tol)
            
            i+=1
        return xi, i
        
        
    def compute_interval_newton(self, F, f, L, vertices, thin_vertices, 
                                thin_L, loc=.5, fjc=False):
        ## compute_interval_newton2 in the '77' code
        x = vertices[0]
        y = vertices[1]
        thinx = thin_vertices[0]
        thiny = thin_vertices[1]
        #loc_=loc
        
        x_update = []
        thin_xupdate = []
        for el, tel in zip(x,thinx):
            x_update.append(el)
            thin_xupdate.append(tel)
        for el,tel in zip(y,thiny):
            x_update.append(el)
            thin_xupdate.append(tel)
        for el,tel in zip(L.equality,thin_L.equality):
            Lv = L.equality[el].interval_Lagrange
            tLv = thin_L.equality[tel].interval_Lagrange
            x_update.append(Lv)
            thin_xupdate.append(tLv)
        if fjc:
            x_update.append(L.fj.interval_Lagrange)
            thin_xupdate.append(thin_L.fj.interval_Lagrange)
        
            
            
        X           = np.asarray(x_update)
        thin_X      = np.asarray(thin_xupdate)
        
        
        #xm          = VectorAD.generate_mid_vector(X, loc=loc_) #correctly returns vector of ad(scalars)
        grad        = IntervalAnalysis.generate_thin_gradient(f.grad, size=f.grad.size)
        this_hess   = F.hess
        this_hess   = np.asarray(this_hess)
        
        usehessian  = this_hess[self.mask,:]
        usehessian  = usehessian[:,self.mask]
        Vint        = X - thin_X #xm
        for i in range(len(Vint)):
            Vint[i] = Vint[i].value
        
            
        nx = IntervalAnalysis.guass_sidel(usehessian,
                                           (Vint[self.mask]).T,
                                           (grad[self.mask]).T, 
                                           nmax=1,
                                           tol = 1.e-6)
        nx = thin_X[self.mask] - nx
        test_updateN    = self.vector_and(X[self.mask], nx)
        return test_updateN
        
    def compute_interval_newton3(self, F, f, L, vertices, thin_vertices, 
                                 thin_L, loc=.5):
        x = vertices[0]
        y = vertices[1]
        thinx = thin_vertices[0]
        thiny = thin_vertices[1]
        
        x_update = []
        thin_xupdate = []
        for el, tel in zip(x,thinx):
            x_update.append(el)
            thin_xupdate.append(tel)
        for el,tel in zip(y,thiny):
            x_update.append(el)
            thin_xupdate.append(tel)
        for el,tel in zip(L.equality,thin_L.equality):
            Lv = L.equality[el].interval_Lagrange
            tLv = thin_L.equality[tel].interval_Lagrange
            x_update.append(Lv)
            thin_xupdate.append(tLv)
            
        X           = np.asarray(x_update)
        thin_X      = np.asarray(thin_xupdate)
        
        
        grad        = IntervalAnalysis.generate_thin_gradient(f.grad, size=f.grad.size)
        this_hess   = F.hess
        this_hess   = np.asarray(this_hess)
        
        usehessian  = this_hess[self.mask,:]
        usehessian  = usehessian[:,self.mask]
        Vint        = X - thin_X #xm
        for i in range(len(Vint)):
            Vint[i] = Vint[i].value
        
            
        A   = usehessian
        Xgs = (Vint[self.mask]).T
        B   = (grad[self.mask]).T
        nmax        = 1
        tol         = 1.e-6
        cnvg        = 0
        count       = 0
        midA = IntervalAnalysis.generate_mid_hessian(A)
        Y = np.linalg.inv(midA)    
        M = np.dot(Y,A)
        r = np.dot(Y,B)
        #size = len(X)
        #Xk = copy.deepcopy(X)
        #Xold = copy.deepcopy(X)
        for i in range(len(M)):
            s1 = np.dot(M[i, :i], Xgs[:i]) #Xnew
            s2 = np.dot(M[i, i + 1:], Xgs[i + 1:]) #Xold
            if (s1 is not None) and (s2 is not None):
                N = thin_X[self.mask][i] + (r[i] - s1 - s2)/(M[i, i])
                #N = (r[i] - s1 - s2)/(M[i, i])
            elif (s1 is not None):
                s2 = 0.
                N = thin_X[self.mask][i] + (r[i] - s1)/(M[i, i])
                #N = (r[i] - s1)/(M[i, i])
            elif (s2 is not None):
                s1 = 0.
                N = thin_X[self.mask][i] + (r[i] - s2)/(M[i, i])
                #N = (r[i] - s2)/(M[i, i])
                
            #d1 = X[self.mask][i].value - \
            #        thin_X[self.mask][i].value
            #test_update = N & d1
            test_update = N.value & X[self.mask][i].value
            if test_update.isempty == False:
                test_update.prove_zero  = test_update in X[self.mask][i].value
                test_update.prove_stuck = test_update.width() == X[self.mask][i].value.width()
                X[self.mask][i].value = test_update
                X[self.mask][i].value.prove_nosol = test_update.isempty
            else:
                X[self.mask][i].value.prove_stuck = test_update.width() == X[self.mask][i].value.width()
                X[self.mask][i].value.isempty = True
                X[self.mask][i].value.prove_nosol = test_update.isempty
            
        return X    
        
    def compute_inewton_gs_intrinsic(self, F, f, L, vertices):
        x = vertices[0]
        y = vertices[1]
        
        x_update = []
        for el in x:
            x_update.append(el)
        for el in y:
            x_update.append(el)
        for el in L.equality:
            Lv = L.equality[el].interval_Lagrange
            x_update.append(Lv)
        X        = np.asarray(x_update)
        
        
        xm          = VectorAD.generate_mid_vector(X, loc=.5) #correctly returns vector of ad(scalars)
        grad        = IntervalAnalysis.generate_thin_gradient(f.grad, size=f.grad.size)
        this_hess   = F.hess
        this_hess   = np.asarray(this_hess)
        
        usehessian  = this_hess[self.mask,:]
        usehessian  = usehessian[:,self.mask]
        Vint        = X - xm
        for i in range(len(Vint)):
            Vint[i] = Vint[i].value
        
            

        A   = usehessian
        Xgs = (Vint[self.mask]).T
        B   = (grad[self.mask]).T
        nmax        = 1
        tol         = 1.e-6
        cnvg        = 0
        count       = 0
        midA = IntervalAnalysis.generate_mid_hessian(A)
        Y = np.linalg.inv(midA)    
        M = np.dot(Y,A)
        r = np.dot(Y,B)
        #size = len(X)
        #Xk = copy.deepcopy(X)
        #Xold = copy.deepcopy(X)
        for i in range(len(M)):
            s1 = np.dot(M[i, :i], Xgs[:i]) #Xnew
            s2 = np.dot(M[i, i + 1:], Xgs[i + 1:]) #Xold
            if (s1 !=None) and (s2!=None):
                N = xm[self.mask][i] + (r[i] - s1 - s2)/(M[i, i])
            elif (s1 !=None):
                s2 = 0.
                N = xm[self.mask][i] + (r[i] - s1)/(M[i, i])
            elif (s2 !=None):
                s1 = 0.
                N = xm[self.mask][i] + (r[i] - s2)/(M[i, i])
                
            test_update = N.value & X[self.mask][i].value
            if test_update.isempty == False:
                test_update.prove_zero  = test_update in X[self.mask][i].value
                test_update.prove_stuck = test_update.width() == X[self.mask][i].value.width()
                X[self.mask][i].value = test_update
                X[self.mask][i].value.prove_nosol = test_update.isempty
            else:
                X[self.mask][i].value.prove_stuck = test_update.width() == X[self.mask][i].value.width()
                X[self.mask][i].value.isempty = True
                X[self.mask][i].value.prove_nosol = test_update.isempty
            
        return X
        
    def vector_and(self,xi, delta):
        return_values = copy.deepcopy(xi)
        for i in range(len(xi)):

            #element of xi
            a = xi[i].value.inf
            c = xi[i].value.sup
            x_simple = ia(a,c)
            
            #element of delta
            d = delta[i].value.inf
            f = delta[i].value.sup
            dz_simple = ia(d,f)
            
            new_values = x_simple & dz_simple
            prove_zero = dz_simple in x_simple
            prove_stuck = new_values.width() == x_simple.width()
            
            return_values[i].value.inf   = new_values.inf
            return_values[i].value.sup   = new_values.sup
            
            return_values[i].isempty     = new_values.isempty
            return_values[i].prove_zero  = prove_zero
            return_values[i].prove_nosol = new_values.isempty
            return_values[i].prove_stuck = prove_stuck
        return return_values
        
    @staticmethod
    def monotonic_contraction(precomputed_f, vertices, vmap):
        """
            see An Interval Constraint Propogation Algorithm
            Exploiting Monoticity page 4, paragraph 5.
            by A.N.T. (ze French guis, oui oui)
            
            f : some constraint function evalutated on the intervals, 
                    especially a monotonic one
            vmap : maps a linear pt index to the curve vertex list of pt lists
        """
        nv_max = copy.deepcopy(vertices)
        nv_min = copy.deepcopy(vertices)
        
        #nv = copy.deepcopy(vertices)

        decreasing, increasing, flat = precomputed_f.eval_monoticity()
        
        for i in range(len(vmap)):
            ivtx = vmap[i]
            deg_max = ia(vertices[ ivtx[0] ][ ivtx[1] ].value.sup,
                         vertices[ ivtx[0] ][ ivtx[1] ].value.sup)
            deg_min = ia(vertices[ ivtx[0] ][ ivtx[1] ].value.inf,
                         vertices[ ivtx[0] ][ ivtx[1] ].value.inf)
            if decreasing[i]:
                nv_min[ ivtx[0] ][ ivtx[1] ].value = deg_max
                nv_max[ ivtx[0] ][ ivtx[1] ].value = deg_min
            elif increasing[i]:
                nv_min[ ivtx[0] ][ ivtx[1] ].value = deg_min
                nv_max[ ivtx[0] ][ ivtx[1] ].value = deg_max
            else: ## not monotonic keep same as old
                pass
        return  nv_min, nv_max
        
    def invertMid(self, f):
        """
            Function to invert and solve the mid(M)
            of an interval matrix equation Mx=b
        """
        mask = self.mask
        SOSC_mid    = f.hess[:,mask]
        SOSC_mid    = SOSC_mid[mask,:]
        return np.linalg.inv(SOSC_mid)
        
    def get_thin_expansion_point(self, pt, vertices, L):
        xpts = []
        ypts = []
        temp_Lagrange = copy.deepcopy(L)
        
        #self.count_fv   = 0  #free variable index (corresponds to a mask location)
        for i in range(L.curve.n):
            xpt = copy.deepcopy(vertices[0][i])
            thin_xpt = vertices[0][i].value.getpoint(pt)
            xpt.value = ia(thin_xpt,thin_xpt)
            xpts.append(xpt)
            
            ypt = copy.deepcopy(vertices[1][i])
            thin_ypt = vertices[1][i].value.getpoint(pt)
            ypt.value = ia(thin_ypt,thin_ypt)
            ypts.append(ypt)
            #self.count_fv +=2

        for el in temp_Lagrange.equality:
            lpt = copy.deepcopy(L.equality[el].interval_Lagrange)
            thin_fp = L.equality[el].interval_Lagrange.value.getpoint(pt)
            lpt.value = ia(thin_fp,thin_fp)
            temp_Lagrange.equality[el].interval_Lagrange = lpt
            #self.count_fv +=1
        
        thin_vertices = [xpts,ypts]
        return thin_vertices, temp_Lagrange
    


def print_components(x):
        print '\nComponents of {}:'.format(x.name)
        print '{}.value ='.format(x.name)
        print '           ',x.value
        print '{}.grad  ='.format(x.name)
        for i in range(x.grad.size):
            print '           ', x.grad[0,i]
        print '{}.hess  ='.format(x.name)
        for i in range(x.grad.size):
            print '           ', x.hess[i]
        print ''
        return
        
##-----------------------------------------------------------------------------
        
        
##-----------------------------Testing-----------------------------------------
        
        
##-----------------------------------------------------------------------------
def compare_multiplication_with_squaring(x):
    print '\ncompare multiplication with squaring'
    x.print_components()
    b=x*x
    c=x**2
    b.name = 'b'
    c.name = 'c'
    print 'b = {}*{}'.format(x.name,x.name)
    print 'c = {}**2'.format(x.name)
    b.print_components()
    c.print_components()
    print '\nResult:'
    print '             the width of a sqared ia object '
    print '             can be smaller that of the same'
    print '             ia object multiplied with itself'
    print 'Is the width of c.value less than that of b.value? (assumes interval type used for the test)'
    print 'c.width = {} <=? b.width = {}'.format(c.value.width(), b.value.width()    )
    print '\nok? = {}'.format(c.value.width() <= b.value.width())
    return 
    

def check_division(x,y):
    print '\ncheck division, q =  x/y:'
    q = x/y
    x.print_components()
    y.print_components()
    print 'q = x/y...'
    q.name = 'q'
    q.print_components()
    print 'verified via hand calc for x=[0,1], y=[3,5]'
    return
    

def Jaulin_contractor_exp_ln():
    #initialize x and y
    xi = ad( ia(.5,1.), name = 'x', N=2, dim = 0)
    yi = ad( ia(2.,3.), name = 'y', N=2, dim = 1)
    #start plotting it
    plt.rc('text', usetex=True)
    plt.plot([xi.value.inf,xi.value.sup],[-.2,-.2], marker = 'o', label=r'$\text{initial x interval}$', color = 'blue')
    plt.plot([-.2,-.2],[yi.value.inf,yi.value.sup], marker = 'o', label=r'$\text{initial y interval}$', color = 'blue')
    plt.title('Contractor for e')
    plt.xlabel('x')
    plt.xlabel('y')
    cn = ad.plot(xi,yi, label_x = xi.name, label_y = yi.name)

    #contract it
    cy = yi.value & xi.value.exp()
    cx = xi.value & yi.value.log()
    #plot it
    cy.name = 'cy'
    cx.name = 'cx'
    dx = copy.deepcopy(cx)
    dy = copy.deepcopy(cy)
    dx.inf = 0.
    dx.sup = 0.
    dy.inf = 0.
    dy.sup = 0.
    dy.name = ''
    dx.name = ''
    #cn = ad.plot(dx,cy,canvas = cn, label_x = cx.name, label_y = cy.name)
    # step 1 y = y & e(x)
    plt.plot([dx.inf,dx.sup],[cy.inf, cy.sup], marker = 'o', label = r'$y = y \cap e([x^{I}])$', color = 'red')
    
    # step 2 x = x & e(x)^-1
    plt.plot([cx.inf,cx.sup],[dy.inf, dy.sup], marker = 'o', label = r'$x = x \cap \left( \exp \left( x^{I} \right)\right)^{-1}$', color = 'red')    
    #plt.annotate('y = exp(x)', xy=dx.value)
    #cn = ad.plot(cx,dy,canvas = cn, label_x = cx.name, label_y = cy.name)
    cn = ad.plot(cx,cy,canvas = cn, label_x = cx.name, label_y = cy.name)
    plt.legend()
    return


def Newt_Func2(x,c):
    return x*x - c
def Newt_Func3(x,c):
    return x**2 -c#c*(-1) + x**2

def some_convergence_criteria(xi,nx,tol):
    return nx.value.width() < xi.value.width() + tol
def scalar_nm(x,c,func, itmax = 1000, tol = 1.e-20):   
    """The scalar valued interval Newton's method
    """
    xi = copy.deepcopy(x)
    i=0
    not_converged = True
    while(i<itmax and not_converged):
        xm = copy.deepcopy(xi)
        xm.value = xi.value.getpoint(.5)
        f = func(xm,c)
        F = func(xi,c)
        nx = xm - f.value/F.grad[0,0]
        not_converged = some_convergence_criteria(xi,nx,tol)
        xi.value = xi.value & nx.value
        i+=1
    return xi, i
    
    
def Moore_NM_example1(func = Newt_Func3):
    """Moore, page 124 example 8.3
    """
    c = ad( ia(2.0,3.0), name = 'c', N=1, dim = -1)
    x = ad( ia(1.0,2.0), name = 'x', N=1, dim = 0 )
    y1, nit = scalar_nm(x,c,func)
    return

def three_camel(X1,X2):
    return 2.*(X1**2) - 1.05*(X1**4) + (X1**6)/6. - X1*X2 + X2**2
    
def hansen_GO_v1_p117(func = three_camel):
    X1 = ad(ia(0.,1.), name = 'X1', N=2, dim = 0)
    X2 = ad(ia(2.,3.), name = 'X2', N=2, dim = 1)
    F = three_camel(X1,X2)
    print 0. in F.grad[0,0]
    print 0. in F.grad[0,1]
    print 'Because 0 is not in the gradient of F w/r/t X2,'
    print 'There is no stationary point in this box'
    return
    
#
# Hansen Test problems, G.O. version 1, page 141
#
def thcf(X):
    X1 = X[0]
    X2 = X[1]
    return 12.*(X1**2) - 6.3*(X1**4) + (X1**6) + 6.*(X2**2 - X1*X2)
    #return 12.*(X1**2) - 6.3*(X1**4) + (X1**6) + 6.*X2*(X2 - X1)
    
def test_problem1():
    X1 = ad(ia(-12.,10.), name = 'X1', N=2, dim = 0)
    X2 = ad(ia(-11.,10.3), name = 'X2', N=2, dim = 1)
    X = np.asarray([X1,X2])
    func = thcf
    tol = 1.e-20
    Xn, nit, F, converged = IntervalAnalysis.compute_interval_newton_basic(thcf,X,itmax=100)
    
    X1 = ad(ia(-1.,1.), name = 'X1', N=2, dim = 0)
    X2 = ad(ia(-1.,1.), name = 'X2', N=2, dim = 1)
    X = np.asarray([X1,X2])
    Xn, nit, F, converged = IntervalAnalysis.compute_interval_newton_basic(thcf,X)
    return

def levy(x,c):#c is not used but fits the format
    return x**6 - 15.*(x**4) + 27.*(x**2) + 250.

def test_problem2(func = levy):
    """Moore, page 124 example 8.3
    """
    c = ad( ia(0.,0.), name = 'c', N=1, dim = -1)
    x = ad( ia(2.,4), name = 'x', N=1, dim = 0 )
    y1, nit = IntervalAnalysis.compute_scalar_newton_basic(levy, c, x)
    return
    

def beale(x):
    x1 = x[0]
    x2 = x[1]
    return ( (1.5 - x1 + x1*x2)**2 + 
             (2.25 - x1 + x1*(x2**2))**2 + 
             (2.625 - x1 + x1*(x2**3))**2 )
    
def test_problem17():
    X1 = ad(ia(2.8,3.2), name = 'X1', N=2, dim = 0)
    X2 = ad(ia(.3,.7), name = 'X2', N=2, dim = 1)
    X = np.asarray([X1,X2])
    Xn, nit, F, converged = IntervalAnalysis.compute_interval_newton_basic(beale,X)
    return
    
    
def schwefel(x):
    f = 0.
    for i in range(3):
        f += ( (x[1] - x[i]**2)**2 + (x[i] - 1.)**2   )
    return f
    
def test_problem18():
    X1 = ad(ia(.1,5.), name = 'X1', N=3, dim = 0)
    X2 = ad(ia(.1,5.), name = 'X2', N=3, dim = 1)
    X3 = ad(ia(.1,5.), name = 'X3', N=3, dim = 2)
    X = np.asarray([X1,X2,X3])
    Xn, nit, F, converged = IntervalAnalysis.compute_interval_newton_basic(schwefel,X)
    return

def simple(x):
    return x[0]**2 + x[0]*x[1]
def test_simple():
    X1 = ad(ia(-100.,100.), name = 'X1', N=2, dim = 0)
    X2 = ad(ia(-100.,100.), name = 'X2', N=2, dim = 1)
    X = np.asarray([X1,X2])
    Xn, nit, F, converged = IntervalAnalysis.compute_interval_newton_basic(simple,X)
    return

def TLM(x):
    return ((x[0]**3)*x[1])-(x[0]*(x[1]**3))-x[1]
def AF_Matrix_Interal_AD_NM_test(func = TLM):
    print '\n\nMatrix Interval AD Newtons Method:'
    X1 = ad(ia(-.6,-.2), name = 'X1', N=2, dim = 0)
    X2 = ad(ia(0.5,0.9), name = 'X2', N=2, dim = 1)
    V = np.asarray([X1,X2])
    Xn, nit, F, converged = IntervalAnalysis.compute_interval_newton_basic(func,V)
    print 'convergence in {} iterations to:'.format(nit)
    for el in Xn:
        print '{}'.format(el)
        print 'and width:'
        print el.value.width()
    F = func(Xn)
    F.name = 'obj func'
    print 'Examine the objective Function:'
    F.print_components()
    print 
    return Xn, F, TLM


def cons1(X):
    x1 = X[0]
    x2 = X[1]
    x3 = X[2]
    x4 = X[3]
    return x1**3 - x2 + x3**2
    
def cons2(X):
    x1 = X[0]
    x2 = X[1]
    x3 = X[2]
    x4 = X[3]
    return x1**2 - x2 - x4**2

def num_ex_12p6page184(X,C):
    def obj(X):
        return -(X[0]**2 + X[1]**2 + X[2]**2 + X[3]**2)
    X1 = ad(ia(-1.1,-0.7), name = 'X1', N=6, dim = 0)
    X2 = ad(ia(-1.2,1.0),  name = 'X2', N=6, dim = 1)
    X3 = ad(ia( 0.0,2.0),  name = 'X3', N=6, dim = 2)
    X4 = ad(ia( 0.0,1.6),  name = 'X4', N=6, dim = 3)
    L1 = ad(ia( -100.0,100.),name = 'L5', N=6, dim = 4)
    L2 = ad(ia( -100.0,100.),name = 'L6', N=6, dim = 5)
    V = np.asarray([X1,X2,X3,X4,L1,L2])
    
    Xi = copy.deepcopy(V)
    n=len(Xi)
    
    xm = VectorAD.generate_mid_vector(Xi)
    f = obj(xm) - V[4]*cons1(xm) - V[5]*cons2(xm)
    F = obj(Xi) - V[4]*cons1(Xi) - V[5]*cons2(Xi)
    this_hess = F.hess
    this_hess = np.asarray(this_hess)
    grad = IntervalAnalysis.generate_thin_gradient(f.grad, size=f.grad.size)
    Vint = Xi - xm
    for i in range(n):
        Vint[i] = Vint[i].value
    nx = IntervalAnalysis.guass_sidel(this_hess,
                                      Vint.T,
                                      grad.T, 
                                      nmax=10,
                                      tol = 1.e-6)
    nx = xm - nx 
    
    
    return 


from interval_arithmetic import ia   
            
if __name__ == '__main__':
    adomg = adObjectMaker.makeGradient #intrinsic to the class now
    adomh = adObjectMaker.makeHessian #intrinsic to the class now
    
    Moore_NM_example1(Newt_Func2)
    Moore_NM_example1(Newt_Func3)
    Jaulin_contractor_exp_ln()
    #Xn, F, TLM = AF_Matrix_Interal_AD_NM_test()
    
    
    N=3
    xsup = 1.
    xinf = 0.
    ysup = 5.
    yinf = 3.
    zsup = 11.
    zinf = 9.
    
    xv = np.asarray([
        ad( ia(xinf, xsup), adomg(N,0), adomh(N)),
        ad( ia(yinf, ysup), adomg(N,1), adomh(N)),
        ad( ia(zinf, zsup), adomg(N,2), adomh(N))
        ])
        
    N=3
    x = ad( ia(xinf, xsup), name = 'x', N=3, dim=0)
    y = ad( ia(yinf, ysup), name = 'y', N=3, dim=1)
    z = ad( ia(zinf, zsup), name = 'z', N=3, dim=2)
    
    sx = ad( ia(-1., 1.), name = 'sx', N=3, dim=0)
    sx.print_components()

    compare_multiplication_with_squaring(sx)

    xv1 = xv+xv
    
    q = xv[0]*xv[1]
    with Timer() as t:
        r = np.dot(xv,xv)
    print "=> elasped time for dot product: %s s" % t.secs
    print 'new 3 vec r = ',r
    r1 = xv[0]*xv[0] + xv[1]*xv[1] + xv[2]*xv[2]
 
    xv = np.asarray([ad( ia(0.,1.), N=2, dim=0 ),
                     ad( ia(3.,5.), N=2, dim = 1)])
    
    with Timer() as t:
        rn = np.dot(xv,xv)
    print "=> elasped time for dot product: %s s" % t.secs
    print 'new 2 vec rn = ',rn

    xv[0].name = 'x'
    xv[1].name = 'y'
    check_division(xv[0],xv[1])

    self = x
    other = y
    store = self/other
    x.arctan2(y)
    
    
    x1 = ia(-1.,1.)
    x2 = ia(-1.,1.)
    x3 = ia(-1.,1.)
    
    def func(x1,x2,x3):
        """ Hansen example f(x1,x2,x3)
        GO 1st edition page 52
        should return [-3.,3.]
        """
        return (x1/(x2+2.))+(x2/(x3+2.))+(x3/(x1+2))
    
    
    x = ad( ia(-6.824, 6.824), name = 'x', N=2, dim=0)
    y = ad( ia(-6.824, 6.824), name = 'y', N=2, dim=1)
    V = [x,y]
    x = V[0]
    y = V[1]
    F = 10.*y + x**3 + 100.*x
    print F.value
    print ''
    print F.grad
    print ''
    print F.hess
    print ''
    print ''
    
    
    def test_HC_BC_2(V):
        """
            constraint:
            f(x,y) = x^3 + 100x + 10y = 0
            
            
        #"""
        x = V[0]
        y = V[1]
        return 10.*y + x**3 + 100.*x
    
    
    xo = ad( ia(-100., 100.), name = 'x', N=2, dim=0)
    yo = ad( ia(-100., 100.), name = 'y', N=2, dim=1)
    X = [xo,yo]
    g = test_HC_BC_2(X)
    print g.value
    print ''
    print g.grad
    print ''
    print g.hess
    
    
    x = ad(0., name='x', N=2, dim=0)
    y = ad(0., name='x', N=2, dim=0)