import numpy as np
from utilities import vector_AND_
#import AF
from interval_arithmetic import ia

class E1(object):
    """The code below is desinged as 
    a one level/one curve procedure.
    
    But, for a given set of projected/truncated/zeroed
    vertices at this level
    It will return the correct (partial) evaluation
    of the total THB spline E value
    -partial because it is 'only' evaluated 
    at this partial level.
    
    In summary, both
    
    the  Boier-Martin evaluation
    using a dictionary of truncated basis per level
    
    and
    
    project all to the finest level 
    
    should work just fine.
    
    
        self.curve = curve
        self.curve.basis_matrix() 
    """
    def __init__(self, curve):
        self.curve = curve
        self.curve.basis_matrix() 
        self.method_dict = {'equality':self.compute,
                            'max':self.compute_max,
                            'min':self.compute_min,
                            'LS':self.compute_LS}
        
    def __call__(self, *args, **kwargs):
        """
            arg[0] : method
            arg[1] : xpts
            arg[2] : ypts
        """
        method = args[0]
        args = args[1:]
        return self.method_dict[method](*args, **kwargs)
        
    def compute(self, *args, **kwargs):
        vertices = args[0]
        #xpts = vertices[0]
        #ypts = vertices[1]
        #self.curve.E1 = np.dot(np.dot(xpts,self.curve.M1),xpts)+\
        #                np.dot(np.dot(ypts,self.curve.M1),ypts)
        self.curve.E1 = 0.
        for pts in vertices:
            self.curve.E1 = self.curve.E1 + \
                        np.dot(np.dot(pts,self.curve.M1),pts) 
        return self.curve.E1
        
    def compute_LS(self, *args, **kwargs):
        desired_fairness = args[1]
        fairness = self.compute(*args, **kwargs)
        return (fairness-desired_fairness)
        
    def compute_max(self, *args, **kwargs):
        max_fairness = args[1]
        fairness = self.compute(*args, **kwargs)
        return (fairness-max_fairness)*(-1.)
        
    def compute_min(self, *args, **kwargs):
        min_fairness = args[1]
        fairness = self.compute(*args, **kwargs)
        return (fairness - min_fairness)
    
    def contractor(self, args, vertices = None):
        return 'Contract on fairness - what are you crazy?'
        
  
class E2(object):
    def __init__(self, curve):
        self.curve = curve
        self.curve.basis_matrix() 
        self.method_dict = {'equality':self.compute,
                            'max':self.compute_max,
                            'min':self.compute_min,
                            'LS':self.compute_LS}
        
    def __call__(self, *args, **kwargs):
        """
            arg[0] : method
            arg[1] : xpts
            arg[2] : ypts
        """
        method = args[0]
        args = args[1:]
        return self.method_dict[method](*args, **kwargs)
        
    def compute(self, *args, **kwargs):
        vertices = args[0]
        #xpts = vertices[0]
        #ypts = vertices[1]
        #self.curve.E2 = np.dot(np.dot(xpts,self.curve.M2),xpts)+\
        #                np.dot(np.dot(ypts,self.curve.M2),ypts)
        self.curve.E2 = 0.
        for pts in vertices:
            self.curve.E2 = self.curve.E2 + \
                        np.dot(np.dot(pts,self.curve.M2),pts) 
        return self.curve.E2
        
    def compute_LS(self, *args, **kwargs):
        desired_fairness = args[1]
        fairness = self.compute(*args, **kwargs)
        return (fairness-desired_fairness)
        
    def compute_max(self, *args, **kwargs):
        max_fairness = args[1]
        fairness = self.compute(*args, **kwargs)
        return (fairness-max_fairness)*(-1.)
        
    def compute_min(self, *args, **kwargs):
        min_fairness = args[1]
        fairness = self.compute(*args, **kwargs)
        return (fairness - min_fairness)
    
    def contractor(self, args, vertices = None):
        return 'Contract on fairness - what are you crazy?'
        

class E3(object):
    def __init__(self, curve):
        self.curve = curve
        self.curve.basis_matrix() 
        self.method_dict = {'equality':self.compute,
                            'max':self.compute_max,
                            'min':self.compute_min,
                            'LS':self.compute_LS}
        
    def __call__(self, *args, **kwargs):
        """
            arg[0] : method
            arg[1] : xpts
            arg[2] : ypts
        """
        method = args[0]
        args = args[1:]
        return self.method_dict[method](*args, **kwargs)
        
    def compute(self, *args, **kwargs):
        vertices = args[0]
        #xpts = vertices[0]
        #ypts = vertices[1]
        #self.curve.E3 = np.dot(np.dot(xpts,self.curve.M3),xpts)+\
        #                np.dot(np.dot(ypts,self.curve.M3),ypts)
        self.curve.E3 = 0.
        for pts in vertices:
            self.curve.E3 = self.curve.E3 + \
                        np.dot(np.dot(pts,self.curve.M3),pts) 
        return self.curve.E3
        
    def compute_LS(self, *args, **kwargs):
        desired_fairness = args[1]
        fairness = self.compute(*args, **kwargs)
        return (fairness-desired_fairness)
        
    def compute_max(self, *args, **kwargs):
        max_area = args[1]
        fairness = self.compute(*args, **kwargs)
        return (fairness-max_area)*(-1.)
        
    def compute_min(self, *args, **kwargs):
        min_area = args[1]
        fairness = self.compute(*args, **kwargs)
        return (fairness - min_area)
    
    def contractor(self, args, vertices = None):
        return 'Contract on fairness - what are you crazy?'
        
        





class ArcLength(object):
    def __init__(self, curve):
        self.curve = curve
        self.curve.pts_M_pts()
        self.method_dict = {'equality':self.compute,
                            'max':self.compute_max,
                            'min':self.compute_min,
                            'LS':self.compute_LS}
        
    def __call__(self, *args, **kwargs):
        """
            arg[0] : method
            arg[1] : xpts
            arg[2] : ypts
        """
        method = args[0]
        args = args[1:]
        return self.method_dict[method](*args, **kwargs)
        
    def compute(self, *args, **kwargs):
        vertices = args[0]
        #xpts = vertices[0]
        #ypts = vertices[1]
        presum = 0.0
        ssum   = 0.0
        for ts in range(1,len(self.curve.s)-1,1):
            M = self.curve.M[:,:,ts] #depends on curve.pts_M_pts()
            presum = 0.
            for pts in (vertices):
                presum   =  np.dot((np.dot(pts,M)),pts)
                presum   =  presum.sqrt()
            #presum is the scalar valued function 
            #we must integrate over the local parameter value:
            a        = self.curve.s[ts-1]   #min 
            b        = self.curve.s[ts]     #max
            ssum     = presum*(b-a) + ssum   # very simple integration
        self.curve.AL=ssum
        return ssum
    
    def compute_LS(self, *args, **kwargs):
        desired_length = args[1]
        length = self.compute(*args, **kwargs)
        return (length-desired_length)
        
    def compute_max(self, *args, **kwargs):
        max_length = args[1]
        length = self.compute(*args, **kwargs)
        return (length-max_length)*(-1.)
        
    def compute_min(self, *args, **kwargs):
        min_length = args[1]
        length = self.compute(*args, **kwargs)
        return (length - min_length)
    
    def contractor(self, args, vertices = None):
        return 'Contract on arc length? ah, no'
        
        
        
        
class ArcLengthApprox(object):
    def __init__(self, curve):
        self.curve = curve
        self.curve.basis_matrix() 
        self.method_dict = {'equality':self.compute,
                            'max':self.compute_max,
                            'min':self.compute_min,
                            'LS':self.compute_LS}
        
    def __call__(self, *args, **kwargs):
        """
            arg[0] : method
            arg[1] : xpts
            arg[2] : ypts
        """
        method = args[0]
        args = args[1:]
        return self.method_dict[method](*args, **kwargs)
        
    def compute(self, *args, **kwargs):
        vertices = args[0]
        #xpts = vertices[0]
        #ypts = vertices[1]
        cheap_sum = 0.
        for i in range(1,self.curve.n):
            test = 0.
            for pts in vertices:
                test += (pts[i] - pts[i-1])**2
            #lx = (xpts[i] -xpts[i-1])
            #ly = (ypts[i] -ypts[i-1])
            #test = lx**2+ly**2#test = lx*lx+ly*ly
            if isinstance(test, ia):
                if abs(test.value.inf <.0):
                    test.value.inf = .0
            #test = test.sqrt() #removing this allows it to converge
            #more often
            cheap_sum = cheap_sum + test
        self.curve.ALapprox = cheap_sum
        return cheap_sum #cheap_sum.sqrt()
    
    def compute_LS(self, *args, **kwargs):
        desired_length = args[1]
        length = self.compute(*args, **kwargs)
        return (length-desired_length)
        
    def compute_max(self, *args, **kwargs):
        max_length = args[1]
        length = self.compute(*args, **kwargs)
        return (length-max_length)*(-1.)
        
    def compute_min(self, *args, **kwargs):
        min_length = args[1]
        length = self.compute(*args, **kwargs)
        return (length - min_length)
    
    def contractor(self, args, vertices = None):
        return 'Contract on arc length? ah, no'