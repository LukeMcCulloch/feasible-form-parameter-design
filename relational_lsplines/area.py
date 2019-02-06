import numpy as np
import copy
from utilities import vector_AND_


class area(object):
    """memoize the curve's integrated basis product matrix
    used for the area computation.
    """
    def __init__(self, curve):
        self.verbose = False
        self.curve = curve
        self.curve.MomentMatrices()
        self.AM1 = copy.deepcopy(self.curve.AM1)
        self.method_dict = {'equality':self.compute,
                            'max':self.compute_max,
                            'min':self.compute_min,
                            'LS':self.compute_LS}
        
    def __call__(self, *args, **kwargs):
        """
            arg[0] : method
            arg[1] : vertice
        """
        method = args[0]
        args = args[1:]
        return self.method_dict[method](*args, **kwargs)
    
    def plot(self, 
             constraint,
             canvas = None, 
             Lspline = None, 
             color = 'green',
             legendloc = None,
             cn = None,
             font_size = 14):
        if canvas == None:print 'error, area plot has no canvas'
        
        if Lspline == None:
            print 'Error in area plot:'
            print 'no Lspline defined, returning to last curve'
            curve = self.curve
            vertices = [curve.xpts,curve.ypts]
        else:
            curve = Lspline.curve
            vertices = [curve.xpts, curve.ypts]
        xi,xe,yi,ye = curve.extreme_C0()
        ymin = min(yi,ye)
        canvas.fill_between(curve.r[:,0],ymin,curve.r[:,1],
                            facecolor = color, alpha=.1,
                            label='area')
        
        curve.compute_moments()
        try:
            Xc = curve.My/curve.area.value
            Yc = curve.Mx/curve.area.value
        except:
            Xc = curve.My/curve.area
            Yc = curve.Mx/curve.area
        name = r'$h_{} := $'.format(str(cn)) + \
            ' '+ r'$A$'  + \
            ' - ' +r'${}$'.format(constraint.pass_value) +\
            ' = ' +r'${}$'.format(np.round(abs(curve.area.value-constraint.pass_value),decimals=2))
        if legendloc == None:
            canvas.annotate(name, 
                            xy = (Xc-1.5,Yc-1.1),
                            xycoords='data')
        else:
            xloc = legendloc[0]
            """
            canvas.annotate(name, 
                            xy = (xloc,legendloc[1]),
                            xycoords='data',
                            fontsize=font_size)
            #"""
            canvas.annotate(r'$A$'+ ' = '+ \
                            r'${}$'.format(constraint.pass_value),
                            xy = (Xc+.5,Yc+.5),
                            xycoords='data',
                            fontsize=font_size)
            #"""
        canvas.plot(Xc,Yc, 
                    marker = 'o', 
                    color = 'black', 
                    alpha = .4,
                    label=name)
        return canvas
        
    def compute(self, *args, **kwargs):
        vertices = args[0]
        xpts = vertices[0]
        ypts = vertices[1]
        #dotprdct   = np.dot(np.dot(xpts,np.transpose(self.curve.AM1)),ypts)
        dotprdct   = np.dot(ypts,np.dot(self.curve.AM1,xpts))
        smple_bdry  = (xpts[-1]*ypts[-1] - xpts[0]*ypts[0])
        self.curve.area = (dotprdct + smple_bdry)*.5
        if self.verbose:
            print 'self.curve.area = {}'.format(self.curve.area )
        return self.curve.area
        
    def compute_LS(self, *args, **kwargs):
        desired_area = args[1]
        area = self.compute(*args, **kwargs)
        return ((area-desired_area)**2)#.sqrt()
        
    def compute_max(self, *args, **kwargs):
        max_area = args[1]
        area = self.compute(*args, **kwargs)
        return (area-max_area)*(-1.)#((area-max_area)**2).sqrt()*(-1.)
        
    def compute_min(self, *args, **kwargs):
        min_area = args[1]
        area = self.compute(*args, **kwargs)
        return (area - min_area) #((area-min_area)**2).sqrt()#
        
    def contractor(self, *args):
        #curve       = copy.deepcopy(args[0])
        vertices = args[0]
        nrange = len(vertices[0])
        xpts = []
        ypts = []
        for i in range(nrange):
            xpts.append(vertices[0][i].value)
            ypts.append(vertices[1][i].value)
        constraint  = args[1]
        #if vertices == None:
        #    xpts = curve.interval_xpts
        #    ypts = curve.interval_ypts
        #else:
        #    xpts = vertices[0]
        #    ypts = vertices[1]
    
        #Area calculations:
        temp1   = np.dot(ypts,np.dot(self.curve.AM1,xpts))
        store2  = (xpts[-1]*ypts[-1] - xpts[0]*ypts[0])
        self.curve.interval_area = (temp1 + store2)*.5
        
        lhs = (store2 - constraint*2.)*(-1.) #this is = to temp1
        #lhs = lhs & temp1
        

        xlist = np.dot(self.curve.AM1,xpts)
        # expanded area computation
        # area_q = 0.
        # for a,b in zip(ypts,xlist):
        #   area_q += a*b
        total_ans       = []
        useful_indices  = []
        bad_indices     = []
        #
        # assemble the dot product with y
        # "by hand"
        for i in range(len(ypts)):   
            min_ans = lhs
            for j in range(len(ypts)):
                if j==i:
                    pass
                else:
                    min_ans = ypts[j]*xlist[j]*(-1) + min_ans
            if xlist[i].contains(0.) == False:
                min_ans = min_ans/xlist[i]
                useful_indices.append(i)
            else:
                bad_indices.append(i)
            total_ans.append(min_ans)
        #print 'ypts = ', ypts
        #print 'total_ans = ', total_ans
        for i in useful_indices:
            if abs( total_ans[i].width() ) > 0.:
                ypts[i] = ypts[i] & total_ans[i]
        #ypts    = vector_AND_(ypts, total_ans)
    
    
        ylist = np.dot(ypts, self.curve.AM1)
        # expanded area computation
        # area_q = 0.
        # for a,b in zip(xpts,ylist):
        #   area_q += a*b
        total_ans       = []
        useful_indices  = []
        bad_indices     = []
        #
        # assemble the dot product with x
        # "by hand"
        for i in range(len(xpts)):   
            min_ans = lhs
            for j in range(len(xpts)):
                if j==i:
                    pass
                else:
                    min_ans = xpts[j]*ylist[j]*(-1.) + min_ans
            #if (abs(ylist[i]) > 0.0):
            if not ylist[i].contains(0.):
                min_ans = min_ans/ylist[i]
                useful_indices.append(i)
            else:
                bad_indices.append(i)
            total_ans.append(min_ans)
        for i in useful_indices:
            if abs( total_ans[i].width() ) > 0.:
                xpts[i] = xpts[i] & total_ans[i]    
        #xpts    = vector_AND_(xpts, total_ans)
        
        for i in range(nrange):
            vertices[0][i].value = xpts[i]
            vertices[1][i].value = ypts[i]
        return vertices
        
  

class area_to_y_axis(object):
    """memoize the curve's integrated basis product matrix
    used for the area (from curve to the y axis) computation."""
    def __init__(self, curve):
        self.verbose = False
        self.curve = curve
        self.curve.MomentMatrices()
        self.AM1 = copy.deepcopy(self.curve.AM1)
        self.method_dict = {'equality':self.compute,
                            'max':self.compute_max,
                            'min':self.compute_min,
                            'LS':self.compute_LS}
        
    def __call__(self, *args, **kwargs):
        """
            arg[0] : method
            arg[1] : vertice
        """
        method = args[0]
        args = args[1:]
        return self.method_dict[method](*args, **kwargs)
    
    def plot(self, 
             constraint,
             canvas = None, 
             Lspline = None, 
             color = 'green',
             legendloc = None,
             cn = None,
             font_size = 14):
        if canvas == None:print 'error, area plot has no canvas'
        
        if Lspline == None:
            print 'Error in area plot:'
            print 'no Lspline defined, returning to last curve'
            curve = self.curve
            vertices = [curve.xpts,curve.ypts]
        else:
            curve = Lspline.curve
            vertices = [curve.xpts, curve.ypts]
        xi,xe,yi,ye = curve.extreme_C0()
        #ymin = min(yi,ye)
        ymax = max(yi,ye)
        ze = np.ones((curve.nump),float)
        ze *= ymax
        canvas.fill_between(curve.r[:,0],curve.r[:,1],ze,
                            facecolor = color, alpha=.1,
                            label='area')
        
        curve.compute_moments_to_y()
        Xc = curve.My_to_y/curve.area_to_y.value
        Yc = curve.Mx_to_y/curve.area_to_y.value
        name = r'$h_{%} := $' % cn + \
            ' '+ r'$A$'  + \
            ' - ' +r'${}$'.format(abs(constraint.pass_value)) +\
            ' = ' +r'${}$'.format(np.round(abs(curve.area_to_y.value-constraint.pass_value)))
        if legendloc==None:
            canvas.annotate(name, 
                            xy = (Xc-2.,Yc),
                            xycoords='data')
        else:
            
            """
            canvas.annotate(name, 
                            xy = (legendloc[0],legendloc[1]),
                            xycoords='data',
                            fontsize=font_size)
            #"""
            canvas.annotate(r'$h_{%d}$' % cn, 
                                xy = (Xc,Yc+.5),
                                xycoords='data',
                                fontsize=font_size)
        #"""
        canvas.plot(Xc,Yc, 
                    marker = 'o', 
                    color = 'black', 
                    alpha = .4,
                    label=name)
        return canvas
        
    def compute(self, *args, **kwargs):
        vertices = args[0]
        xpts = vertices[0]
        ypts = vertices[1]
        #dotprdct   = np.dot(np.dot(xpts,np.transpose(self.curve.AM1)),ypts)
        #        dotprdct   = np.dot(ypts,np.dot(self.curve.AM1,xpts))
        #        smple_bdry  = (xpts[-1]*ypts[-1] - xpts[0]*ypts[0])
        #        self.curve.area = (dotprdct + smple_bdry)*.5
        
        #temp1   = np.dot(np.dot(xpts,(-self.curve.AM1).T),ypts)
        
#        temp1   = np.dot(ypts,np.dot(self.curve.AM1,xpts))
#        store2  = (xpts[0]*ypts[0] - xpts[-1]*ypts[-1])
#        self.curve.area_to_y = (temp1 + store2)*.5
        
        temp1   = np.dot(ypts,np.dot(self.curve.AM1,xpts))
        store2  = (xpts[0]*ypts[0] - xpts[-1]*ypts[-1])
        self.curve.area_to_y = -(temp1 + store2)*0.5
        
        if self.verbose:
            print 'self.curve.area = {}'.format(self.curve.area_to_y )
        return self.curve.area_to_y
        
    def compute_LS(self, *args, **kwargs):
        desired_area = args[1]
        area = self.compute(*args, **kwargs)
        return ((area-desired_area)**2)#.sqrt()
        
    def compute_max(self, *args, **kwargs):
        max_area = args[1]
        area = self.compute(*args, **kwargs)
        return (area-max_area)*(-1.)#((area-max_area)**2).sqrt()*(-1.)
        
    def compute_min(self, *args, **kwargs):
        min_area = args[1]
        area = self.compute(*args, **kwargs)
        return (area - min_area) #((area-min_area)**2).sqrt()#
        
    def contractor(self, *args):
        vertices = args[0]
        nrange = len(vertices[0])
        xpts = []
        ypts = []
        for i in range(nrange):
            xpts.append(vertices[0][i].value)
            ypts.append(vertices[1][i].value)
        constraint  = args[1]

        #temp1   = np.dot(np.dot(xpts,(-self.curve.AM1).T),ypts)
        #temp1   = np.dot(ypts,np.dot(self.curve.AM1,xpts))
        #store2  = (xpts[0]*ypts[0] - xpts[-1]*ypts[-1])
        #self.curve.area_to_y = (temp1 + store2)*.5
        
        temp1   = -np.dot(ypts,np.dot(self.curve.AM1,xpts))
        store2  = -(xpts[0]*ypts[0] - xpts[-1]*ypts[-1])
        self.curve.area_to_y = (temp1 + store2)*0.5
        
        lhs = (store2 - constraint*2.)*(-1.) #this is = to temp1
        #lhs = lhs & temp1
        
        ## Contract on Y:
        xlist = np.dot(xpts,(-self.curve.AM1).T)
        # expanded area computation
        # area_q = 0.
        # for a,b in zip(ypts,xlist):
        #   area_q += a*b
        total_ans       = []
        useful_indices  = []
        bad_indices     = []
        #
        # assemble the dot product with y
        # "by hand"
        for i in range(len(ypts)):   
            min_ans = lhs
            for j in range(len(ypts)):
                if j==i:
                    pass
                else:
                    min_ans = ypts[j]*xlist[j]*(-1) + min_ans
            if xlist[i].contains(0.) == False:
                min_ans = min_ans/xlist[i]
                useful_indices.append(i)
            else:
                bad_indices.append(i)
            total_ans.append(min_ans)
        #print 'ypts = ', ypts
        #print 'total_ans = ', total_ans
        for i in useful_indices:
            if abs( total_ans[i].width() ) > 0.:
                ypts[i] = ypts[i] & total_ans[i]
        #ypts    = vector_AND_(ypts, total_ans)
                
                
        ylist = np.dot((-self.curve.AM1).T,ypts)
        # expanded area computation
        # area_q = 0.
        # for a,b in zip(xpts,ylist):
        #   area_q += a*b
        total_ans       = []
        useful_indices  = []
        bad_indices     = []
        #
        # assemble the dot product with x
        # "by hand"
        for i in range(len(xpts)):   
            min_ans = lhs
            for j in range(len(xpts)):
                if j==i:
                    pass
                else:
                    min_ans = xpts[j]*ylist[j]*(-1.) + min_ans
            #if (abs(ylist[i]) > 0.0):
            if not ylist[i].contains(0.):
                min_ans = min_ans/ylist[i]
                useful_indices.append(i)
            else:
                bad_indices.append(i)
            total_ans.append(min_ans)
        for i in useful_indices:
            if abs( total_ans[i].width() ) > 0.:
                xpts[i] = xpts[i] & total_ans[i]    
        #xpts    = vector_AND_(xpts, total_ans)
        for i in range(nrange):
            vertices[0][i].value = xpts[i]
            vertices[1][i].value = ypts[i]
        return vertices


class area_to_any_x_axis(object):
    def __init__(self, curve, x_axis_loc = 0.):
        self.verbose        = False
        self.curve          = curve
        self.x_axis_loc     = x_axis_loc
        self.curve.MomentMatrices()
        self.AM1 = copy.deepcopy(self.curve.AM1)
        self.method_dict = {'equality':self.compute,
                            'max':self.compute_max,
                            'min':self.compute_min,
                            'LS':self.compute_LS}
        
    def __call__(self, *args, **kwargs):
        """
            arg[0] : method
            arg[1] : vertice
        """
        method = args[0]
        args = args[1:]
        return self.method_dict[method](*args, **kwargs)
    
    def plot(self, 
             constraint,
             canvas = None, 
             Lspline = None, 
             color = 'green',
             legendloc = None,
             cn = None,
             font_size = 14):
        if canvas == None:print 'error, area plot has no canvas'
        
        if Lspline == None:
            print 'Error in area plot:'
            print 'no Lspline defined, returning to last curve'
            curve = self.curve
            vertices = [curve.xpts,curve.ypts]
        else:
            curve = Lspline.curve
            vertices = [curve.xpts, curve.ypts]
        xi,xe,yi,ye = curve.extreme_C0()
        ymin = min(yi,ye)
        canvas.fill_between(curve.r[:,0],ymin,curve.r[:,1],
                            facecolor = color, alpha=.1,
                            label='area')
        
        curve.compute_moments()
        try:
            Xc = curve.My/curve.area.value
            Yc = curve.Mx/curve.area.value
        except:
            Xc = curve.My/curve.area
            Yc = curve.Mx/curve.area
        name = r'$h_{} := $'.format(cn) + \
            ' '+ r'$A$'  + \
            ' - ' +r'${}$'.format(constraint.pass_value) +\
            ' = ' +r'${}$'.format(np.round(abs(curve.area.value-constraint.pass_value),decimals=2))
        if legendloc == None:
            canvas.annotate(name, 
                            xy = (Xc-1.5,Yc-1.1),
                            xycoords='data',
                            fontsize=font_size)
        else:
            xloc = legendloc[0]
            """
            canvas.annotate(name, 
                            xy = (xloc,legendloc[1]),
                            xycoords='data',
                            fontsize=font_size)
            #"""
            canvas.annotate(r'$A$'+ ' = '+ \
                            r'${}$'.format(constraint.pass_value), 
                            xy = (Xc+.5,Yc+.5),
                            xycoords='data',
                            fontsize=font_size)
            #"""
        canvas.plot(Xc,Yc, 
                    marker = 'o', 
                    color = 'black', 
                    alpha = .4,
                    label=name)
        return canvas
        
    def compute(self, *args, **kwargs):
        vertices = args[0]
        xpts = vertices[0]
        ypts = [self.x_axis_loc - el for el in vertices[1]]
        #dotprdct   = np.dot(np.dot(xpts,np.transpose(self.curve.AM1)),ypts)
        dotprdct   = np.dot(ypts,np.dot(self.curve.AM1,xpts))
        smple_bdry  = (xpts[-1]*ypts[-1] - xpts[0]*ypts[0])
        self.curve.area = (dotprdct + smple_bdry)*.5
        if self.verbose:
            print 'self.curve.area = {}'.format(self.curve.area )
        return self.curve.area
        
    def compute_LS(self, *args, **kwargs):
        desired_area = args[1]
        area = self.compute(*args, **kwargs)
        return ((area-desired_area)**2)#.sqrt()
        
    def compute_max(self, *args, **kwargs):
        max_area = args[1]
        area = self.compute(*args, **kwargs)
        return (area-max_area)*(-1.)#((area-max_area)**2).sqrt()*(-1.)
        
    def compute_min(self, *args, **kwargs):
        min_area = args[1]
        area = self.compute(*args, **kwargs)
        return (area - min_area) #((area-min_area)**2).sqrt()#
        
    def contractor(self, *args):
        #curve       = copy.deepcopy(args[0])
        vertices = args[0]
        nrange = len(vertices[0])
        xpts = []
        ypts = []
        for i in range(nrange):
            xpts.append(vertices[0][i].value)
            ypts.append(vertices[1][i].value)
        constraint  = args[1]
        #if vertices == None:
        #    xpts = curve.interval_xpts
        #    ypts = curve.interval_ypts
        #else:
        #    xpts = vertices[0]
        #    ypts = vertices[1]
    
        #Area calculations:
        temp1   = np.dot(ypts,np.dot(self.curve.AM1,xpts))
        store2  = (xpts[-1]*ypts[-1] - xpts[0]*ypts[0])
        self.curve.interval_area = (temp1 + store2)*.5
        
        lhs = (store2 - constraint*2.)*(-1.) #this is = to temp1
        #lhs = lhs & temp1
        

        xlist = np.dot(self.curve.AM1,xpts)
        # expanded area computation
        # area_q = 0.
        # for a,b in zip(ypts,xlist):
        #   area_q += a*b
        total_ans       = []
        useful_indices  = []
        bad_indices     = []
        #
        # assemble the dot product with y
        # "by hand"
        for i in range(len(ypts)):   
            min_ans = lhs
            for j in range(len(ypts)):
                if j==i:
                    pass
                else:
                    min_ans = ypts[j]*xlist[j]*(-1) + min_ans
            if xlist[i].contains(0.) == False:
                min_ans = min_ans/xlist[i]
                useful_indices.append(i)
            else:
                bad_indices.append(i)
            total_ans.append(min_ans)
        #print 'ypts = ', ypts
        #print 'total_ans = ', total_ans
        for i in useful_indices:
            if abs( total_ans[i].width() ) > 0.:
                ypts[i] = ypts[i] & total_ans[i]
        #ypts    = vector_AND_(ypts, total_ans)
    
    
        ylist = np.dot(ypts, self.curve.AM1)
        # expanded area computation
        # area_q = 0.
        # for a,b in zip(xpts,ylist):
        #   area_q += a*b
        total_ans       = []
        useful_indices  = []
        bad_indices     = []
        #
        # assemble the dot product with x
        # "by hand"
        for i in range(len(xpts)):   
            min_ans = lhs
            for j in range(len(xpts)):
                if j==i:
                    pass
                else:
                    min_ans = xpts[j]*ylist[j]*(-1.) + min_ans
            #if (abs(ylist[i]) > 0.0):
            if not ylist[i].contains(0.):
                min_ans = min_ans/ylist[i]
                useful_indices.append(i)
            else:
                bad_indices.append(i)
            total_ans.append(min_ans)
        for i in useful_indices:
            if abs( total_ans[i].width() ) > 0.:
                xpts[i] = xpts[i] & total_ans[i]    
        #xpts    = vector_AND_(xpts, total_ans)
        
        for i in range(nrange):
            vertices[0][i].value = xpts[i]
            vertices[1][i].value = ypts[i]
        return vertices