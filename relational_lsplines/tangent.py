import numpy as np
import matplotlib.pyplot as plt
import copy
from utilities import vector_AND_


class tangent(object):
    """
        uses radians
    """
    def __init__(self, curve, s):
        self.curve  = curve
        self.loc    = s
        self.compute_basis()
        self.method_dict = {'equality':self.compute,
                            'max':self.compute_max,
                            'min':self.compute_min,
                            'LS':self.compute_LS}
    
    def __call__(self, *args, **kwargs):
        """
            arg[0] : method
            arg[1] : vertices
        """
        method = args[0]
        args = args[1:]
        return self.method_dict[method](*args, **kwargs)
    
    def plot(self, 
             constraint, 
             canvas = None, 
             Lspline = None, 
             color = 'blue',
             legendloc = None,
             cn = None,
             font_size = 14):
                 
        if canvas == None:print 'tangent plot has no canvas'
        
        if Lspline == None:
            print 'Error in tangent plot:'
            print 'no Lspline defined, returning to last curve'
            curve = self.curve
            vertices = [curve.xpts,curve.ypts]
        else:
            curve = Lspline.curve
            vertices = [curve.xpts, curve.ypts]
        plt.rc('text', usetex=True)
        plt.rc('font', family='sans-serif')
        angle_value   = self.compute(vertices) #tan_value = Lspline.curve.compute_tangent(0.)
        loc         = self.loc
        v0 = curve.CurvePoint(loc)
        mag = .5
        
        v1 = np.asarray([v0[0] + mag,
                         v0[1]]  )
        v2 = np.asarray([v0[0] + mag,
                         v0[1] + np.tan(angle_value.value)*mag ]) 

        V = np.asarray([v0,v1,v2])
        canvas.plot(V[:,0],V[:,1], color, alpha = .4)
        canvas.plot([V[-1,0],V[0,0]],[V[-1,1],V[0,1]], color, alpha = .4)
        #"""
        if self.loc == 0.:
            name = r'$h_{%d} := $' % cn + \
                   ' ' +r'$\alpha_{B}$'+ \
                   ' - ' +r'${}$'.format(abs(np.rad2deg(constraint.pass_value))) + \
                   ' = ' +r'${}$'.format(np.round(abs(np.rad2deg(angle_value.value - constraint.pass_value)),decimals=2))
        elif self.loc == 1.:
            name = r'$h_{%d} := $' % cn + \
                   ' ' +r'$\alpha_{E}$' + \
                   ' - ' +r'${}$'.format(abs(np.rad2deg(constraint.pass_value))) + \
                   ' = ' +r'${}$'.format(np.round(abs(np.rad2deg(angle_value.value - constraint.pass_value)),decimals=2))
        else:
            name = r'$h_{%d} := $' % cn + \
                   ' '+r'$\alpha_{}$'.format(str(cn)) + \
                   ' - ' +r'${}$'.format(abs(np.rad2deg(constraint.pass_value))) + \
                   ' = ' +r'${}$'.format(np.round(abs(np.rad2deg(angle_value.value - constraint.pass_value)),decimals=2))
        #"""
        if legendloc == None:
            canvas.annotate(name, 
                            xy = (v0 +[.0,.5]),
                            xycoords='data')
        else:
            if loc <=.5:
                xloc = legendloc[2]
            else:
                xloc = legendloc[0]
            """
            canvas.annotate(name, 
                            xy = (v0),
                            xytext = (xloc,legendloc[1]),
                            xycoords='data',
                            fontsize=font_size)
            #"""
            if abs(self.loc) < 0.01:
                tstring = r'$\alpha_{B}$' + \
                   ' = ' +r'${}$'.format(abs(np.rad2deg(constraint.pass_value))) 
            elif abs(self.loc-1.) < 0.01:
                tstring = r'$\alpha_{E}$' + \
                   ' = ' +r'${}$'.format(abs(np.rad2deg(constraint.pass_value))) 
            else:
                tstring = r'$\alpha_{}$'.format(str(cn)) + \
                   ' = ' +r'${}$'.format(abs(np.rad2deg(constraint.pass_value))) 
            canvas.annotate(tstring, 
                            xy = (v0+[.0,.5]),
                            xycoords='data',
                            fontsize=font_size)
            #"""
            
            canvas.plot(v0[0],v0[1], 
                        marker = 'o', 
                        color = 'black', 
                        alpha = .4,
                        label=name )
        return canvas
        
    def compute_basis(self):
        localBasis = np.zeros((self.curve.n,self.curve.n),float)
        span = self.curve.FindSpan(self.loc)
        self.curve.DersBasisFunc(span,self.loc,localBasis[:,span-self.curve.p:span+1])
        self.localBasis = localBasis
        return
    
    def compute(self, *args, **kwargs):
        """
            See The NURBS Book, page 91
            
            1st derivativeof a B-spline Curve :=
                1st deriv basis-funcs DOT Contol_Points
                
            We want the angle of this tangent
            Hence the code below:
        """
        vertices = args[0]
        xpts = vertices[0]
        ypts = vertices[1]
        qx = np.dot(xpts,self.localBasis[1]) 
        qy = np.dot(ypts,self.localBasis[1])
        return (qy).arctan2(qx)
        
    def compute_LS(self, *args, **kwargs):
        desired_tangent_angle = args[1]
        tan_angle = self.compute(*args, **kwargs)
        return (tan_angle - desired_tangent_angle)**2
        
    def compute_max(self, *args, **kwargs):
        maxvalue = args[1]
        tan = self.compute(*args, **kwargs)
        return (tan-maxvalue)*(-1.0)
        
    def compute_min(self, *args, **kwargs):
        minvalue = args[1]
        tan = self.compute(*args, **kwargs)
        return (tan-minvalue)
        
    def monotonic_relation(self, *args, **kwargs):
        """
            how to determine reation between points?
            
            use basis functions?
            
            invert the function?
            
        """
        return
        
    def contractor(self, *args, **kwargs):
        """
            a fwd-bckwd contractor 
            for the tangent angle constraint 
        """
        vertices    = copy.deepcopy(args[0])
        constraint  = copy.deepcopy(args[1])
        nrange = len(vertices[0])
        xpts = []
        ypts = []
        for i in range(nrange):
            xpts.append(vertices[0][i].value)
            ypts.append(vertices[1][i].value)
            
            
        qx = np.dot(xpts,self.localBasis[1]) 
        #qx = np.dot(xpts,self.localBasis[0])
        qy = np.dot(ypts,self.localBasis[1])

        ## split the constraint - need to abstract this!    
        #constraint could be an interval or float.  probably not ADscalar
        tanc = np.tan(constraint) 
        if isinstance(tanc, np.float) :
            lhs = qx*float(tanc)
        else:
            lhs = qx*tanc
            
        total_ans       = []
        useful_indices  = []
        bad_indices     = []
        for i in range(len(ypts)):   
            min_ans = 0.
            for j in range(len(ypts)):
                if j == i:  ## in interval analysis we *must* avoid duplicate computation, if we hope to contract
                    pass   ## so I go about inverting the computation by hand in this loop
                elif j!=i:
                    #min_ans = (lhs - ypts[j]*float(self.localBasis[1,j])) + min_ans #WRONG! left for posterity
                    min_ans = min_ans + ypts[j]*float(self.localBasis[1,j])
            min_ans = lhs - min_ans
            if (abs(float(self.localBasis[1,i])) > 0.0):
                min_ans = min_ans/float(self.localBasis[1,i])
                useful_indices.append(i)
            else: #divide by zero:
                bad_indices.append(i)
            total_ans.append(min_ans)
            
        new_ans    = vector_AND_(ypts, total_ans)
        for i in useful_indices:
            if not new_ans[i].isempty:  #  #abs( new_ans[i].width() ) > 0.: #
                ypts[i] = ypts[i] & new_ans[i]
            else:
                print 'warning, possible constraint violation, tangent'
                
                
        for i in range(nrange):
            vertices[0][i].value = xpts[i]
            vertices[1][i].value = ypts[i]
        return vertices#[xpts, ypts]




class vertical_tangent(object):
    """
        uses radians
    """
    def __init__(self, curve, s):
        self.curve  = curve
        self.loc    = s
        self.compute_basis()
        self.method_dict = {'equality':self.compute,
                            'max':self.compute_max,
                            'min':self.compute_min,
                            'LS':self.compute_LS}
    
    def __call__(self, *args, **kwargs):
        """
            arg[0] : method
            arg[1] : vertices
        """
        method = args[0]
        args = args[1:]
        return self.method_dict[method](*args, **kwargs)
    
    def plot(self, 
             constraint, 
             canvas = None, 
             Lspline = None, 
             color = 'blue',
             legendloc = None,
             cn = None,
             font_size = 14):
                 
        if canvas == None:print 'tangent plot has no canvas'
        
        if Lspline == None:
            print 'Error in tangent plot:'
            print 'no Lspline defined, returning to last curve'
            curve = self.curve
            vertices = [curve.xpts,curve.ypts]
        else:
            curve = Lspline.curve
            vertices = [curve.xpts, curve.ypts]
        plt.rc('text', usetex=True)
        plt.rc('font', family='sans-serif')
        angle_value   = self.compute(vertices)
        loc         = self.loc
        v0 = curve.CurvePoint(loc)
        mag = .5#np.linalg.norm(curve.vertices[-1] - curve.vertices[0])
        #mag = 0.1*mag
        
        v1 = np.asarray([v0[0] + mag,
                         v0[1]]  )
        v2 = np.asarray([v0[0] + mag,
                         v0[1] + np.tan(angle_value.value)*mag ]) 
        
        
        V = np.asarray([v0,v1,v2])
        canvas.plot(V[:,0],V[:,1], color, alpha = .4)
        canvas.plot([V[-1,0],V[0,0]],[V[-1,1],V[0,1]], color, alpha = .4)
        """
        if self.loc == 0.:
            name = r'$h_{%d} := $' % cn + \
                   ' ' +r'$\alpha_{B}$'.format(str(cn)) + \
                   ' - ' +r'${}$'.format(abs(np.rad2deg(constraint.pass_value))) + \
                   ' = ' +r'${}$'.format(np.round(abs(np.rad2deg(angle_value.value - constraint.pass_value)),decimals=2))
        elif self.loc == 1.:
            name = r'$h_{%d} := $' % cn + \
                   ' ' +r'$\alpha_{E}$'.format(str(cn)) + \
                   ' - ' +r'${}$'.format(abs(np.rad2deg(constraint.pass_value))) + \
                   ' = ' +r'${}$'.format(np.round(abs(np.rad2deg(angle_value.value - constraint.pass_value)),decimals=2))
        else:
            name = r'$h_{%d} := $' % cn + \
               ' ' +r'$\alpha_{}$'.format(str(cn)) + \
               ' - ' +r'${}$'.format(abs(np.rad2deg(constraint.pass_value))) + \
               ' = ' +r'${}$'.format(np.round(abs(np.rad2deg(angle_value.value - constraint.pass_value)),decimals=2))
        #"""
        #"""
        if abs(self.loc) < 0.01:
            name = r'$h_{%d} := $' % cn + \
                   ' ' +r'$\alpha_{B}$' + \
                   ' - ' +r'${}$'.format(abs(np.rad2deg(constraint.pass_value))) + \
                   ' = ' +r'${}$'.format(np.round(abs(np.rad2deg(angle_value.value - constraint.pass_value)),decimals=2))
        elif abs(self.loc-1.) < 0.01:
            name = r'$h_{%d} := $' % cn + \
                   ' ' +r'$\alpha_{E}$' + \
                   ' - ' +r'${}$'.format(abs(np.rad2deg(constraint.pass_value))) + \
                   ' = ' +r'${}$'.format(np.round(abs(np.rad2deg(angle_value.value - constraint.pass_value)),decimals=2))
        else:
            name = r'$h_{%d} := $' % cn + \
               ' ' +r'$\alpha_{}$'.format(str(cn)) + \
               ' - ' +r'${}$'.format(abs(np.rad2deg(constraint.pass_value))) + \
               ' = ' +r'${}$'.format(np.round(abs(np.rad2deg(angle_value.value - constraint.pass_value)),decimals=2))
        #"""
        if legendloc == None:
            canvas.annotate(name, 
                            xy = (v0 +[0.,0.5]),
                            xycoords='data')
        else:
            if loc <=.5:
                xloc = legendloc[2]
            else:
                xloc = legendloc[0]
            """
            canvas.annotate(name, 
                            xy = (xloc,legendloc[1]),
                            xycoords='data',
                            fontsize=font_size)
            #"""
            if abs(self.loc) < 0.01:
                tstring = r'$\alpha_{B}$' + \
                   ' = ' +r'${}$'.format(abs(np.rad2deg(constraint.pass_value))) 
            elif abs(self.loc-1.) < 0.01:
                tstring = r'$\alpha_{E}$' + \
                   ' = ' +r'${}$'.format(abs(np.rad2deg(constraint.pass_value))) 
            else:
                tstring = r'$\alpha_{}$'.format(str(cn)) + \
                   ' = ' +r'${}$'.format(abs(np.rad2deg(constraint.pass_value))) 
            canvas.annotate(tstring, 
                            xy = (v0+[.0,0.5]),
                            xycoords='data',
                            fontsize=font_size)
            #"""
            canvas.plot(v0[0],v0[1], 
                        marker = 'o', 
                        color = 'black', 
                        alpha = .4,
                        label=name)
        return canvas
        
    def compute_basis(self):
        localBasis = np.zeros((self.curve.n,self.curve.n),float)
        span = self.curve.FindSpan(self.loc)
        self.curve.DersBasisFunc(span,self.loc,localBasis[:,span-self.curve.p:span+1])
        self.localBasis = localBasis
        return
    
    def compute(self, *args, **kwargs):
        """
            See The NURBS Book, page 91
            
            1st derivativeof a B-spline Curve :=
                1st deriv basis-funcs DOT Contol_Points
                
            We want the angle of this tangent
            Hence the code below:
        """
        vertices = args[0]
        xpts = vertices[0]
        ypts = vertices[1]
        qx = np.dot(xpts,self.localBasis[1]) 
        qy = np.dot(ypts,self.localBasis[1])
        return  (qx/qy).atan()# qx.arctan2(qy)# 
        
    def compute_LS(self, *args, **kwargs):
        desired_tangent = args[1]
        tan = self.compute(*args, **kwargs)
        return (tan - desired_tangent)**2
        
    def compute_max(self, *args, **kwargs):
        maxvalue = args[1]
        tan = self.compute(*args, **kwargs)
        return (tan-maxvalue)*(-1.0)
        
    def compute_min(self, *args, **kwargs):
        minvalue = args[1]
        tan = self.compute(*args, **kwargs)
        return (tan-minvalue)
        
    def monotonic_relation(self, *args, **kwargs):
        """
            how to determine reation between points?
            
            use basis functions?
            
            invert the function?
            
        """
        return
        
    def contractor(self, *args, **kwargs):
        """
            a fwd-bckwd contractor 
            for the tangent angle constraint 
        """
        vertices = copy.deepcopy(args[0])
        constraint  = copy.deepcopy(args[1])
        
        nrange = len(vertices[0])
        xpts = []
        ypts = []
        for i in range(nrange):
            xpts.append(vertices[0][i].value)
            ypts.append(vertices[1][i].value)
        
        
        qx = np.dot(xpts,self.localBasis[1]) 
        #qy = np.dot(ypts,localBasis[1])

        ## split the constraint - need to abstract this!    
        #constraint could be an interval or float.  probably not ADscalar
        tanc = np.tan(constraint) 
        if isinstance(tanc, np.float) :
            lhs = qy*float(tanc)
        else:
            lhs = qy*tanc
            
        total_ans       = []
        useful_indices  = []
        bad_indices     = []
        for i in range(len(xpts)):   
            min_ans = 0.
            for j in range(len(xpts)):
                if j == i:  ## in interval analysis we *must* avoid duplicate computation, if we hope to contract
                    pass   ## so I go about inverting the computation y hand in this loop
                else:
                    #min_ans = (lhs - xpts[j]*float(self.localBasis[1,j])) + min_ans
                    min_ans = min_ans + ypts[j]*float(self.localBasis[1,j])
            min_ans = lhs - min_ans
            if (abs(float(self.localBasis[1,i])) > 0.0):
                min_ans = min_ans/float(self.localBasis[1,i])
                useful_indices.append(i)
            else:
                bad_indices.append(i)
            total_ans.append(min_ans)
            
        new_ans    = vector_AND_(ypts, total_ans)
        for i in useful_indices:
            if not new_ans[i].isempty:  #  #abs( new_ans[i].width() ) > 0.: #
                xpts[i] = xpts[i] & new_ans[i]
            else:
                print 'warning, possible constraint violation, tangent'
                
        for i in range(nrange):
            vertices[0][i].value = xpts[i]
            vertices[1][i].value = ypts[i]
        return vertices#[xpts, ypts]



class tangent_XoverZ(object):
    """
        uses radians
        
        -meant for the peculiar 
        coordinate system of the real hull
    """
    def __init__(self, curve, s):
        self.curve  = curve
        self.loc    = s
        self.compute_basis()
        self.method_dict = {'equality':self.compute,
                            'max':self.compute_max,
                            'min':self.compute_min,
                            'LS':self.compute_LS}
    
    def __call__(self, *args, **kwargs):
        """
            arg[0] : method
            arg[1] : vertices
        """
        method = args[0]
        args = args[1:]
        return self.method_dict[method](*args, **kwargs)
    
        
    def compute_basis(self):
        localBasis = np.zeros((self.curve.n,self.curve.n),float)
        span = self.curve.FindSpan(self.loc)
        self.curve.DersBasisFunc(span,self.loc,localBasis[:,span-self.curve.p:span+1])
        self.localBasis = localBasis
        return
    
    def compute(self, *args, **kwargs):
        """
            See The NURBS Book, page 91
            
            1st derivativeof a B-spline Curve :=
                1st deriv basis-funcs DOT Contol_Points
                
            We want the angle of this tangent
            Hence the code below:
        """
        vertices = args[0]
        xpts = vertices[0]
        #ypts = vertices[1]
        zpts = vertices[2]
        qx = np.dot(xpts,self.localBasis[1]) 
        qz = np.dot(zpts,self.localBasis[1])
        #return  (qz/qx).atan()#
        return  qz.arctan2(qx)# 
        
    
        vertices = args[0]
        xpts = vertices[0]
        #ypts = vertices[1]
        zpts = vertices[2]
        qx = np.dot(xpts,self.localBasis[1]) 
        qz = np.dot(zpts,self.localBasis[1])
        return (qx).arctan2(qz)
        
    def compute_LS(self, *args, **kwargs):
        desired_tangent_angle = args[1]
        tan_angle = self.compute(*args, **kwargs)
        return (tan_angle - desired_tangent_angle)**2
        
    def compute_max(self, *args, **kwargs):
        maxvalue = args[1]
        tan = self.compute(*args, **kwargs)
        return (tan-maxvalue)*(-1.0)
        
    def compute_min(self, *args, **kwargs):
        minvalue = args[1]
        tan = self.compute(*args, **kwargs)
        return (tan-minvalue)
    
    
class vertical_XoverZ_tangent(object):
    """
        uses radians
        
        -meant for the peculiar 
        coordinate system of the real hull
    """
    def __init__(self, curve, s):
        self.curve  = curve
        self.loc    = s
        self.compute_basis()
        self.method_dict = {'equality':self.compute,
                            'max':self.compute_max,
                            'min':self.compute_min,
                            'LS':self.compute_LS}
    
    def __call__(self, *args, **kwargs):
        """
            arg[0] : method
            arg[1] : vertices
        """
        method = args[0]
        args = args[1:]
        return self.method_dict[method](*args, **kwargs)
    
        
    def compute_basis(self):
        localBasis = np.zeros((self.curve.n,self.curve.n),float)
        span = self.curve.FindSpan(self.loc)
        self.curve.DersBasisFunc(span,self.loc,localBasis[:,span-self.curve.p:span+1])
        self.localBasis = localBasis
        return
    
    def compute(self, *args, **kwargs):
        """
            See The NURBS Book, page 91
            
            1st derivativeof a B-spline Curve :=
                1st deriv basis-funcs DOT Contol_Points
                
            We want the angle of this tangent
            Hence the code below:
        """
        vertices = args[0]
        xpts = vertices[0]
        #ypts = vertices[1]
        zpts = vertices[2]
        qx = np.dot(xpts,self.localBasis[1]) 
        qz = np.dot(zpts,self.localBasis[1])
        #return  (qz/qx).atan()#
        return  qz.arctan2(qx)# 
        
    def compute_LS(self, *args, **kwargs):
        desired_tangent = args[1]
        tan = self.compute(*args, **kwargs)
        return (tan - desired_tangent)**2
        
    def compute_max(self, *args, **kwargs):
        maxvalue = args[1]
        tan = self.compute(*args, **kwargs)
        return (tan-maxvalue)*(-1.0)
        
    def compute_min(self, *args, **kwargs):
        minvalue = args[1]
        tan = self.compute(*args, **kwargs)
        return (tan-minvalue)














class harries_tangent(object):
    """
        uses radians
    """
    def __init__(self, curve, s, alpha,index,
                 run_index=0, rise_index=1):
        self.curve  = curve
        self.loc    = s
        self.alpha  = alpha
        if s == 0.:
            assert(index == 0)
        elif s==1:
            assert(index == -1 or index == curve.n-1)
            index -= 1
        self.index      = index
        self.run_index  = run_index
        self.rise_index = rise_index
        #self.compute_basis()
        self.method_dict = {'equality':self.compute,
                            'max':self.compute_max,
                            'min':self.compute_min,
                            'LS':self.compute_LS}
    
    def __call__(self, *args, **kwargs):
        """
            arg[0] : method
            arg[1] : vertices
        """
        method = args[0]
        args = args[1:]
        return self.method_dict[method](*args, **kwargs)
    
    def plot(self, 
             constraint, 
             canvas = None, 
             Lspline = None, 
             color = 'blue',
             legendloc = None,
             cn = None,
             font_size = 14):
                 
        if canvas == None:print 'tangent plot has no canvas'
        
        if Lspline == None:
            print 'Error in tangent plot:'
            print 'no Lspline defined, returning to last curve'
            curve = self.curve
            vertices = [curve.xpts,curve.ypts]
        else:
            curve = Lspline.curve
            vertices = [curve.xpts, curve.ypts]
        plt.rc('text', usetex=True)
        plt.rc('font', family='sans-serif')
        angle_value   = self.compute(vertices) #tan_value = Lspline.curve.compute_tangent(0.)
        loc         = self.loc
        v0 = curve.CurvePoint(loc)
        mag = .5
        
        v1 = np.asarray([v0[0] + mag,
                         v0[1]]  )
        v2 = np.asarray([v0[0] + mag,
                         v0[1] + np.tan(angle_value.value)*mag ]) 

        V = np.asarray([v0,v1,v2])
        canvas.plot(V[:,0],V[:,1], color, alpha = .4)
        canvas.plot([V[-1,0],V[0,0]],[V[-1,1],V[0,1]], color, alpha = .4)
        #"""
        if self.loc == 0.:
            name = r'$h_{%d} := $' % cn + \
                   ' ' +r'$\alpha_{B}$'+ \
                   ' - ' +r'${}$'.format(abs(np.rad2deg(constraint.pass_value))) + \
                   ' = ' +r'${}$'.format(np.round(abs(np.rad2deg(angle_value.value - constraint.pass_value)),decimals=2))
        elif self.loc == 1.:
            name = r'$h_{%d} := $' % cn + \
                   ' ' +r'$\alpha_{E}$' + \
                   ' - ' +r'${}$'.format(abs(np.rad2deg(constraint.pass_value))) + \
                   ' = ' +r'${}$'.format(np.round(abs(np.rad2deg(angle_value.value - constraint.pass_value)),decimals=2))
        else:
            name = r'$h_{%d} := $' % cn + \
                   ' '+r'$\alpha_{}$'.format(str(cn)) + \
                   ' - ' +r'${}$'.format(abs(np.rad2deg(constraint.pass_value))) + \
                   ' = ' +r'${}$'.format(np.round(abs(np.rad2deg(angle_value.value - constraint.pass_value)),decimals=2))
        #"""
        if legendloc == None:
            canvas.annotate(name, 
                            xy = (v0 +[.0,.5]),
                            xycoords='data')
        else:
            if loc <=.5:
                xloc = legendloc[2]
            else:
                xloc = legendloc[0]
            """
            canvas.annotate(name, 
                            xy = (v0),
                            xytext = (xloc,legendloc[1]),
                            xycoords='data',
                            fontsize=font_size)
            #"""
            if abs(self.loc) < 0.01:
                tstring = r'$\alpha_{B}$' + \
                   ' = ' +r'${}$'.format(abs(np.rad2deg(constraint.pass_value))) 
            elif abs(self.loc-1.) < 0.01:
                tstring = r'$\alpha_{E}$' + \
                   ' = ' +r'${}$'.format(abs(np.rad2deg(constraint.pass_value))) 
            else:
                tstring = r'$\alpha_{}$'.format(str(cn)) + \
                   ' = ' +r'${}$'.format(abs(np.rad2deg(constraint.pass_value))) 
            canvas.annotate(tstring, 
                            xy = (v0+[.0,.5]),
                            xycoords='data',
                            fontsize=font_size)
            #"""
            
            canvas.plot(v0[0],v0[1], 
                        marker = 'o', 
                        color = 'black', 
                        alpha = .4,
                        label=name )
        return canvas
        
    def compute_basis(self):
        """not used in Harries' method
        """
        localBasis = np.zeros((self.curve.n,self.curve.n),float)
        span = self.curve.FindSpan(self.loc)
        self.curve.DersBasisFunc(span,
                                 self.loc,
                                 localBasis[:,span-self.curve.p:span+1])
        self.localBasis = localBasis
        return
    
    def compute(self, *args, **kwargs):
        """
            See Harries 1998 dissertation 
            page 71
            and Lothar's 6145 class notes
        """
        vertices = args[0]
        xpts = vertices[self.run_index]
        ypts = vertices[self.rise_index]
        #qx = np.dot(xpts,self.localBasis[0]) 
        #qy = np.dot(ypts,self.localBasis[0])
        
        V0 = [xpts[self.index],
              ypts[self.index]]
        V1 = [xpts[self.index+1],
              ypts[self.index+1]]
        #        dx = (V1[0]-V0[0])/np.cos(self.alpha)
        #        dy = (V1[1]-V0[1])/np.sin(self.alpha)
        #        dalpha = np.sqrt(dx**2+dy**2)
        #        return np.arccos( (V1[0]-V0[0]) / dalpha )
        return np.arctan2( (V1[1]-V0[1]) , (V1[0]-V0[0]) )
        
    def compute_LS(self, *args, **kwargs):
        desired_tangent_angle = args[1]
        tan_angle = self.compute(*args, **kwargs)
        return (tan_angle - desired_tangent_angle)**2
        
    def compute_max(self, *args, **kwargs):
        maxvalue = args[1]
        tan = self.compute(*args, **kwargs)
        return (tan-maxvalue)*(-1.0)
        
    def compute_min(self, *args, **kwargs):
        minvalue = args[1]
        tan = self.compute(*args, **kwargs)
        return (tan-minvalue)
        
    def monotonic_relation(self, *args, **kwargs):
        """
            how to determine reation between points?
            
            use basis functions?
            
            invert the function?
            
        """
        return
        