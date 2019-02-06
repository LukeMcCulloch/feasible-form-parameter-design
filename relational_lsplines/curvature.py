import numpy as np
import matplotlib.pyplot as plt
import copy
from utilities import vector_AND_
from routines import crossproduct
#from AF import fuzzyNumber
from  interval_arithmetic import ia
#from interval_analysis import monotonic_contraction


class curvature(object):
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
            arg[1] : xpts
            arg[2] : ypts
        """
        method = args[0]
        args = args[1:]
        return self.method_dict[method](*args, **kwargs)
    
    def plot(self, 
             constraint,
             canvas = None, 
             Lspline = None, 
             scale = 4.0, 
             color = 'blue',
             legendloc = None, 
             cn = None,
             docurvature = False,
             font_size = 14):
                 
        if canvas == None:print 'curvature plot has no canvas'
        
        if Lspline == None:
            print 'Error in curvature plot:'
            print 'no Lspline defined, returning to last curve'
            curve = self.curve
            vertices = [curve.xpts,curve.ypts]
        else:
            curve = Lspline.curve
            vertices = [curve.xpts, curve.ypts]
        plt.rc('text', usetex=True)
        plt.rc('font', family='sans-serif')
        ortho_vector = np.asarray([0.,0.,1.])
        if docurvature:
            x = []
            y = []
            for i in range(curve.nump):
                s = curve.s[i]
                C2 = curve.compute_curvature(s)
                curvature = C2
                slope = curve.compute_tangent(s)
                tanpoint = [np.cos(slope),np.sin(slope),0.]
                vector = crossproduct(ortho_vector, tanpoint)
                tanpoint = [curve.r[i][0]+vector[0],
                            curve.r[i][1]+vector[1]]
                curvature_point = [curve.r[i][0]+scale*curvature*vector[0],
                                   curve.r[i][1]+scale*curvature*vector[1]]
                x.append(curvature_point[0])
                y.append(curvature_point[1])
                canvas.plot([curve.r[i][0],curvature_point[0]],
                            [curve.r[i][1],curvature_point[1]], 
                            color = 'grey', alpha = .4)
            #canvas.plot(x,y)
            C2 = curve.compute_curvature(self.loc)
            x = np.asarray(x)
            y = np.asarray(y)
            canvas.plot(x, y, color = 'grey', alpha = .4)
            canvas.fill_between(curve.r[:,0], 
                                y, 
                                curve.r[:,1], 
                                facecolor = color, 
                                alpha=.0,
                                label='curvature')
        else:
            C2 = curve.compute_curvature(self.loc)
            
            name = r'$h_{%d} := $' % cn + \
                   ' '+r'$c_{%d} $' % cn+ \
                   ' - ' +r'${}$'.format(np.rad2deg(constraint.pass_value)) + \
                   ' = ' +r'${}$'.format(np.round(abs(C2 - constraint.pass_value),decimals=2))
        if legendloc == None:
            canvas.annotate(name, 
                            xy = (curve.CurvePoint(constraint.loc) + [0.,1.3] ),
                            xycoords='data')
        else:
            v0 = curve.CurvePoint(constraint.loc)
            if self.loc <=.5:
                xloc = legendloc[2]
            else:
                xloc = legendloc[0]
            """
            canvas.annotate(name, 
                            xy = (xloc,legendloc[1]),
                            xycoords='data',
                            fontsize=font_size)
            """
            cstrg = r'$c_{%d} $' % cn
            canvas.annotate(cstrg + ' = ' +\
                            r'${}$'.format(np.rad2deg(constraint.pass_value)),                    
                            xy = (v0+[.0,1.5]),
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
        vertices = args[0]
        xpts = vertices[0]
        ypts = vertices[1]
        qxdot   = np.dot(xpts,self.localBasis[1])
        qxddot  = np.dot(xpts,self.localBasis[2])
        qydot   = np.dot(ypts,self.localBasis[1])
        qyddot  = np.dot(ypts,self.localBasis[2])
        store = (qxdot*qyddot - qydot*qxddot)
        temp = np.sqrt(qxdot**2 + qydot**2)
        if isinstance(temp, ia):
            if temp.inf<=0:
                temp.inf = 0.
        denom = temp*((temp)**2)#**.5## #problem foud with sqrt
        """
        if isinstance(denom, fuzzyNumber):
            if denom.min<0:
                denom.min.value = 0.
            #if denom.min.value == 0.:
            #    denom.min.value = 1.e-9
            if store.contains(0.) and (store.width()<1.e-4):
                curvature = store
                dummy = copy.deepcopy(curvature)
                dummy.min.value = 0.
                dummy.real.value = 0.
                dummy.max.value = 0.
                curvature = curvature*dummy/denom
            else:
                curvature = store/denom
        else:
            curvature = store/(np.sqrt(qxdot*qxdot + qydot*qydot)**3.)
        #"""
        curvature = store/denom#((np.sqrt(qxdot*qxdot + qydot*qydot))**3.)
        return curvature
        
    def compute_LS(self, *args, **kwargs):
        desired_tangent = args[1]
        curvature = self.compute(*args, **kwargs)
        return (curvature - desired_tangent)**2
        
    def compute_max(self, *args, **kwargs):
        maxvalue = args[1]
        curvature = self.compute(*args, **kwargs)
        return (curvature-maxvalue)*(-1.0)
        
    def compute_min(self, *args, **kwargs):
        maxvalue = args[1]
        curvature = self.compute(*args, **kwargs)
        return (curvature-maxvalue)
        
    def update_allq(self, xpts, ypts):
        """ compute automatic differentiated curvature:
        """
        qxdot   = np.dot(xpts,self.localBasis[1,:])
        qxddot  = np.dot(xpts,self.localBasis[2,:])
        qydot   = np.dot(ypts,self.localBasis[1,:])
        qyddot  = np.dot(ypts,self.localBasis[2,:]) 
        return qxdot,qxddot,qydot,qyddot
            
    def contractor(self, *args, **kwargs):
        """
            a fwd-bckwd contractor 
            for the curvature constraint 
            
            args[0] : interval xpts
            args[1] : interval ypts
            args[2] : constraint
        """
        vertices    = copy.deepcopy(args[0])
        nrange = len(vertices[0])
        xpts = []
        ypts = []
        for i in range(nrange):
            xpts.append(vertices[0][i].value)
            ypts.append(vertices[1][i].value)
        constraint  = copy.deepcopy(args[1])
        
        
        
        
        qxdot,qxddot,qydot,qyddot = self.update_allq(xpts,ypts)
        
        ## the all important computation split (need to abstract this kind of thing)
        ##lhs = (np.sqrt(qxdot*qxdot + qydot*qydot)**3.) *constraint
        lhs = ( ( np.sqrt(qxdot**2 + qydot**2) )**3 )*constraint
        
        #        check2 = qxdot*qyddot
        #        if check2.width() < 1.e-2:
        #            check2.min.value = check2.real.value
        #            check2.max.value = check2.real.value
        #        t1 = (lhs - check2)/qydot
        
        #
        # qyddot
        #
        check2 = qydot*qxddot
        if check2.width() < 1.e-2 and check2.contains(0.):
            check2.inf = 0.
            check2.sup = 0.
        #if qxdot.contains(0.) and abs(qxdot.min.value)>1.e-6:
        #    print 'qxdot = ',qxdot
        #    print 'qxdot not invertable, implement other logic please'
        if abs(float(qxdot.inf))<1.e-6:
            qxdot.inf = 1.e-10
        print 'invert qxdot'
        print 'qxdot = ', qxdot
        
        #t1 = (lhs + qydot*qxddot)/(qxdot)
        t1 = (lhs + check2)/(qxdot)
        
        t1 = t1 & qyddot # go ahead and shrink t1 to qyddot - they are logically equivalent
        total_ans       = []
        useful_indices  = []
        bad_indices     = []
        for i in range(len(ypts)):   
            min_ans = 0.
            for j in range(len(ypts)):
                if j==i:
                    pass
                else:
                    min_ans = (ypts[j]*float(self.localBasis[2,j])) + min_ans
            min_ans = t1 - min_ans
            if (abs(float(self.localBasis[2,i])) > 0.0):
                min_ans = min_ans/float(self.localBasis[2,i])
                useful_indices.append(i)
            else:
                bad_indices.append(i)
            total_ans.append(min_ans)
            
        new_ans    = vector_AND_(ypts, total_ans)
        for i in useful_indices:
            if new_ans[i].isempty == False:  #  abs( new_ans[i].width() ) > 0.:
                ypts[i] = ypts[i] & new_ans[i]
                qxdot,qxddot,qydot,qyddot = self.update_allq(xpts,ypts)
            else:
                print 'warning, possible constraint violation, curvature 1'
                
        ## 
        ## qxdot
        ##
        check2 = qydot*qxddot
        if check2.width() < 1.e-2 and check2.contains(0.):
            check2.inf = 0.
            check2.sup = 0.
        #if qyddot.contains(0.):
        #    print 'qyddot = ',qyddot
        #    print 'qyddot not invertable, implement other logic please'
        
        if qyddot.contains(0.) and qyddot.width()<1.e-6:
            qxdot.inf = 0.#1.e-10
        print 'invert qyddot'
        print 'qyddot = ',qyddot
        fix =  (lhs + check2)*(1./qyddot)#*(qyddot**-1.)
        fix = fix & qxdot # go ahead and shrink fix to qxdot - they are logically equivalent
        total_ans       = []
        useful_indices  = []
        bad_indices     = []
        
        for i in range(len(xpts)): #contract on x[i]
            min_ans = 0.
            for j in range(len(xpts)): # add up all jth pieces of the dot product except i
                if j==i:
                    pass
                else:
                    
                    min_ans = (xpts[j]*float(self.localBasis[1,j] ) ) + min_ans
            min_ans = fix - min_ans
            if (abs(float(self.localBasis[1,i]) ) >0.0 ):
                min_ans = min_ans/float(self.localBasis[1,i])
                useful_indices.append(i)
            else:
                bad_indices.append(i)
            total_ans.append(min_ans)
        
        new_ans = vector_AND_(xpts, total_ans)
        for i in useful_indices:
            if not new_ans[i].isempty:  #  abs( new_ans[i].width() ) > 0.:
                xpts[i] = xpts[i] & new_ans[i]
                qxdot,qxddot,qydot,qyddot = self.update_allq(xpts,ypts)
            else:
                print 'warning, possible constraint violation, curvature 2'
            
                
        ## switch to the other side
                
        ##
        ## contract on qydot
        ##
        check2 = qxdot*qyddot
        if check2.width() < 1.e-2 and check2.contains(0.):
            check2.inf = 0.
            check2.sup = 0.
#        if qxddot.contains(0.):
#            print 'qxddot = ',qxddot
#            print 'qxddot not invertable, implement other logic please'
#            qxddot.min.value = 0.
        if qxddot.contains(0.):
            qxddot.inf = 0.
            
            print 'invert qxddot'
            print 'qxddot = ',qxddot
            t1 = (lhs - check2)/(-qxddot)#*(-qxddot**-1)
            t1 = t1 & qydot
            total_ans       = []
            useful_indices  = []
            bad_indices     = []
            for i in range(len(ypts)):   
                min_ans = 0.
                for j in range(len(ypts)):
                    if j==i:
                        pass
                    else:
                        #print 't1 = ',t1
                        #print 'ypts[{}] = {}'.format(i,ypts[i])
                        #print 'localbasis[{},{}] = {}'.format(1,i,self.localBasis[1,j])
                        min_ans = (ypts[j]*float(self.localBasis[1,j])) + min_ans
                min_ans = t1 - min_ans
                if (abs(float(self.localBasis[1,i])) > 0.0):
                    min_ans = min_ans/float(self.localBasis[1,i])
                    useful_indices.append(i)
                else:
                    bad_indices.append(i)
                total_ans.append(min_ans)
                
            new_ans    = vector_AND_(ypts, total_ans)
            for i in useful_indices:
                if not new_ans[i].isempty:  #  abs( new_ans[i].width() ) > 0.:
                    ypts[i] = ypts[i] & new_ans[i]
                else:
                    print 'warning, possible constraint violation, curvature 3'
                
        ##contract on qxdot
            
        check2 = qxdot*qyddot
        if check2.width() < 1.e-2 and check2.contains(0.):
            check2.inf = 0.
            check2.sup = 0.
        #contract on qxddot
#        if qydot.contains(0.):
#            print 'qydot = ',qxddot
#            print 'qydot not invertable, implement other logic please'
        if qydot.contains(0.):
            qydot.inf = 0.
            print 'invert qydot'
            print 'qydot = ',qydot
            fix = (lhs - qxdot*qyddot)/(-qydot)#*(-qydot**-1)
            fix = fix & qxddot # go ahead and shrink t1 to quddot - they are logically equivalent
            total_ans       = []
            useful_indices  = []
            bad_indices     = []
            for i in range(len(xpts)):
                min_ans = 0.
                for j in range(len(xpts)):
                    if j==i:
                        pass
                    else:
                        min_ans = (xpts[j]*float(self.localBasis[2,j] ) ) + min_ans
                min_ans = fix - min_ans
                if (abs(float(self.localBasis[2,i]) ) >0.0 ):
                    min_ans = min_ans/float(self.localBasis[2,i])
                    useful_indices.append(i)
                else:
                    bad_indices.append(i)
                total_ans.append(min_ans)
            
            new_ans = vector_AND_(xpts, total_ans)
            for i in useful_indices:
                if not new_ans[i].isempty:  #  abs( new_ans[i].width() ) > 0.:
                    xpts[i] = xpts[i] & new_ans[i]
                else:
                    print 'warning, possible constraint violation, curvature 4'
                    
        for i in range(nrange):
            vertices[0][i].value = xpts[i]
            vertices[1][i].value = ypts[i]
        return vertices

    # monotonic_contraction
    def monotonic_contractor(self, *args, **kwargs):
        """
            a fwd-bckwd contractor 
            for the curvature constraint 
            
            args[0] : interval xpts
            args[1] : interval ypts
            args[2] : constraint
        """
        vertices    = copy.deepcopy(args[0])
        nrange = len(vertices[0])
        xpts = []
        ypts = []
        for i in range(nrange):
            xpts.append(vertices[0][i].value)
            xpts.append(vertices[1][i].value)
        constraint  = copy.deepcopy(args[1])
        
        
        # compute automatic differentiated curvature:
        qxdot   = np.dot(xpts,self.localBasis[1,:])
        qxddot  = np.dot(xpts,self.localBasis[2,:])
        qydot   = np.dot(ypts,self.localBasis[1,:])
        qyddot  = np.dot(ypts,self.localBasis[2,:]) 
        #computation of doubledots is expanded below
        
    
        ## the all important computation split (need to abstract this kind of thing)
        ##lhs = ((np.sqrt(qxdot*qxdot + qydot*qydot) )**3. )*constraint
        lhs = (np.sqrt(qxdot**2 + qydot**2)**3.) *constraint
        
        #        check2 = qxdot*qyddot
        #        if check2.width() < 1.e-2:
        #            check2.min.value = check2.real.value
        #            check2.max.value = check2.real.value
        #        t1 = (lhs - check2)/qydot
        
        #
        # qyddot
        #
        check2 = qydot*qxddot
        #if check2.width() < 1.e-2:
        #    check2.min.value = check2.real.value
        #    check2.max.value = check2.real.value
        if qxdot.contains(0.):
            print 'qxdot = ',qxdot
            print 'qxdot not invertable, implement other logic please'
        else:
            print 'invert qxdot'
            print 'qxdot = ', qxdot
            t1 = (lhs + qydot*qxddot)/(qxdot)#*(qxdot**-1.)
            t1 = t1 & qyddot # go ahead and shrink t1 to qyddot - they are logically equivalent
            total_ans       = []
            useful_indices  = []
            bad_indices     = []
            for i in range(len(ypts)):   
                min_ans = 0.
                for j in range(len(ypts)):
                    if j==i:
                        pass
                    else:
                        min_ans = (t1 - ypts[j]*float(self.localBasis[2,j])) + min_ans
                if (abs(float(self.localBasis[2,i])) > 0.0):
                    min_ans = min_ans/float(self.localBasis[2,i])
                    useful_indices.append(i)
                else:
                    bad_indices.append(i)
                total_ans.append(min_ans)
            new_ans    = vector_AND_(ypts, total_ans)
            for i in useful_indices:
                if not new_ans[i].isempty:  #  abs( new_ans[i].width() ) > 0.:
                    ypts[i] = ypts[i] & new_ans[i]
                
        ## 
        ## qxdot
        ##
        
        if qyddot.contains(0.):
            print 'qyddot = ',qyddot
            print 'qyddot not invertable, implement other logic please'
        else:
            print 'invert qyddot'
            print 'qyddot = ',qyddot
            fix =  (lhs + qydot*qxddot)/(qyddot)#*(qyddot**-1.)
            fix = fix & qxdot # go ahead and shrink fix to qxdot - they are logically equivalent
            total_ans       = []
            useful_indices  = []
            bad_indices     = []
            
            for i in range(len(xpts)): #contract on x[i]
                min_ans = 0.
                for j in range(len(xpts)): # add up all jth pieces of the dot product except i
                    if j==i:
                        pass
                    else:
                        
                        min_ans = (fix - xpts[j]*float(self.localBasis[1,j] ) ) + min_ans
                if (abs(float(self.localBasis[1,i]) ) >0.0 ):
                    min_ans = min_ans/float(self.localBasis[1,i])
                    useful_indices.append(i)
                else:
                    bad_indices.append(i)
                total_ans.append(min_ans)
            
            new_ans = vector_AND_(xpts, total_ans)
            for i in useful_indices:
                if not new_ans[i].isempty:  #  abs( new_ans[i].width() ) > 0.:
                    xpts[i] = xpts[i] & new_ans[i]
                
                
        ## switch to the other side
                
        ##
        ## contract on qydot
        ##
        check2 = qxdot*qyddot
        #if check2.width() < 1.e-2:
        #    check2.min.value = check2.real.value
        #    check2.max.value = check2.real.value
        if qxddot.contains(0.):
            print 'qxddot = ',qxddot
            print 'qxddot not invertable, implement other logic please'
        else:
            print 'invert qxddot'
            print 'qxddot = ',qxddot
            t1 = (lhs - qxdot*qyddot)/(-qxddot)#*(-qxddot**-1)
            t1 = t1 & qydot
            total_ans       = []
            useful_indices  = []
            bad_indices     = []
            for i in range(len(ypts)):   
                min_ans = 0.
                for j in range(len(ypts)):
                    if j==i:
                        pass
                    else:
                        #print 't1 = ',t1
                        #print 'ypts[{}] = {}'.format(i,ypts[i])
                        #print 'localbasis[{},{}] = {}'.format(1,i,self.localBasis[1,j])
                        min_ans = (t1 - ypts[j]*float(self.localBasis[1,j])) + min_ans
                if (abs(float(self.localBasis[1,i])) > 0.0):
                    min_ans = min_ans/float(self.localBasis[1,i])
                    useful_indices.append(i)
                else:
                    bad_indices.append(i)
                total_ans.append(min_ans)
                
            new_ans    = vector_AND_(ypts, total_ans)
            for i in useful_indices:
                if not new_ans[i].isempty:  #  abs( new_ans[i].width() ) > 0.:
                    ypts[i] = ypts[i] & new_ans[i]
                
        ##contract on qxdot
            
            
        #contract on qxddot
        if qydot.contains(0.):
            print 'qydot = ',qxddot
            print 'qydot not invertable, implement other logic please'
        else:
            print 'invert qydot'
            print 'qydot = ',qydot
            fix = (lhs - qxdot*qyddot)/(-qydot)#*(-qydot**-1)
            fix = fix & qxddot # go ahead and shrink t1 to quddot - they are logically equivalent
            total_ans       = []
            useful_indices  = []
            bad_indices     = []
            for i in range(len(xpts)):
                min_ans = 0.
                for j in range(len(xpts)):
                    if j==i:
                        pass
                    else:
                        
                        min_ans = (fix - xpts[j]*float(self.localBasis[2,j] ) ) + min_ans
                if (abs(float(self.localBasis[2,i]) ) >0.0 ):
                    min_ans = min_ans/float(self.localBasis[2,i])
                    useful_indices.append(i)
                else:
                    bad_indices.append(i)
                total_ans.append(min_ans)
            
            new_ans = vector_AND_(xpts, total_ans)
            for i in useful_indices:
                if not new_ans[i].isempty:  #  abs( new_ans[i].width() ) > 0.:
                    xpts[i] = xpts[i] & new_ans[i]
                
        
        for i in range(nrange):
            vertices[0][i].value = xpts[i]
            vertices[1][i].value = ypts[i]
   
        return vertices
    
    

class curvature_XoverZ(object):
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
            arg[1] : xpts
            arg[2] : ypts
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
        """3D curvature, we are only interested 
        in longitudinal-transverse flatness..
        i.e. the z-x plane in my coordinate system
        """
        vertices = args[0]
        xpts = vertices[2] # z plays the 'x' part
        ypts = vertices[0] # x plays the 'y' part
        #zpts = vertices[1]
        #********************************************
        # switcharoo:  using z in place of x
        # using x in place of y
        # i.e.
        #
        #  y <- x
        #  x <- z
        #
        qxdot   = np.dot(xpts,self.localBasis[1])
        qxddot  = np.dot(xpts,self.localBasis[2])
        qydot   = np.dot(ypts,self.localBasis[1])
        qyddot  = np.dot(ypts,self.localBasis[2])
        store = (qxdot*qyddot - qydot*qxddot)
        temp = np.sqrt(qxdot**2 + qydot**2)
        if isinstance(temp, ia):
            if temp.inf<=0:
                temp.inf = 0.
        denom = temp*((temp)**2)#**.5## #problem foud with sqrt
        #
        curvature = store/denom#((np.sqrt(qxdot*qxdot + qydot*qydot))**3.)
        return curvature
        
    def compute_LS(self, *args, **kwargs):
        desired_tangent = args[1]
        curvature = self.compute(*args, **kwargs)
        return (curvature - desired_tangent)**2
        
    def compute_max(self, *args, **kwargs):
        maxvalue = args[1]
        curvature = self.compute(*args, **kwargs)
        return (curvature-maxvalue)*(-1.0)
        
    def compute_min(self, *args, **kwargs):
        maxvalue = args[1]
        curvature = self.compute(*args, **kwargs)
        return (curvature-maxvalue)