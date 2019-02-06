"""
    TODO:
        dimension agnostic
        outward rounding
        full suit
"""
import numpy as np
import copy
from  automatic_differentiation import ad
from   interval_arithmetic import ia

try: #outward rounding:
    import mpmath  as mp
    has_mp = True
except:
    has_mp = False

class IntervalAnalysis(object):
    
    def __init__(self, mask = None):
        self.mask = mask
        self.n = len(mask)
        self.check_box = []
        self.has_mp = has_mp
        return
        

    def round_up(self,num):
        return float(mp.mpf(mp.iv.mpf(str(num)).a))
        
    def guass_sidel(self, A,X,B, nmax=1, tol = 1.e-6): #dino gs
        """
            A is a square interval matrix
            B is a column vector RHS
            X is a column vector of unknowns
        """
        cnvg        = 0
        count       = 0
        midA = self.matrix_mid(A)
        Y = np.linalg.inv(midA)    
        G = np.dot(Y,A)
        C = np.dot(Y,B)
        """
            Snyder, 92, Generative Modeling for CAD
            page 148:
            
            Interval Newtons Method
            Remarkably, 
            if at any stage the iteration  we find Z ( X_n,
            then f has a unique zero in X_n...
        """
        size = len(X)
        Xk = copy.deepcopy(X)
        Xold = copy.deepcopy(X)
        while((cnvg < size) and count < nmax):
            Xold = Xk
            for i in range(G.shape[0]):
                s1 = np.dot(G[i, :i], Xk[:i])
                s2 = np.dot(G[i, i + 1:], Xold[i + 1:])
                if (s1 !=None) and (s2!=None):
                    #Z = ((G[i, i])**(-1))*(C[i] - s1 - s2)
                    Z = (C[i] - s1 - s2)/(G[i, i])
                elif (s1 !=None):
                    s2 = 0.
                    #Z = ((G[i, i])**-1)*(C[i] - s1)
                    Z = (C[i] - s1)/(G[i, i])
                elif (s2 !=None):
                    s1 = 0.
                    #Z = ((G[i, i])**-1)*(C[i] - s2)
                    Z = (C[i] - s2)/(G[i, i])
    
                #Xk[i] = Xold[i] & Z
                test_update = Xold[i] & Z       
                #for i, el in enumerate(test_update):
                #if abs(test_update.width().value) > 0.:  #we are ditching perfect answers??
                if test_update.isempty == False:
                    Xk[i] = test_update
                else:
                    Xk[i].isempty = True
                if (abs(Xk[i].width().value) < tol): #possible fail for true 0 width intervals? check and Xk[i].isempty = False??
                    cnvg += 1
                
                try:
                    Xk[i].real = Xk[i].midpoint()
                except:
                    print 'error, could not update midpoint'
                    print 'some other error occured'
                    pass
                count += 1
    
        return Xk
        
    def compute_scalar_newton(self, func, constraint, vertices, vertex_tpl, max_it = 30, tol = 1.e-6):
        """For Constraints Only!
        """
        x_update = copy.deepcopy(vertices)
        not_conv = True
        count = 0
        xyz = vertex_tpl[0]     # which dimension are we in [0,1,2]?
        index = vertex_tpl[1]   # which vertex are we attempting to contract?
        n = len(vertices[0])
        gn = xyz*n + index
        while(not_conv and count < max_it):
                midx                    = x_update[xyz][index].getpoint(.5)
                x_thin                  = copy.deepcopy(vertices)
                x_thin[xyz][index]      = midx
                f_scalar                = func(x_thin, constraint)  #reduces the width in accordance with thick constraint
                f_all                   = func(x_update, constraint)
                
                #f_allmid                = f_all.midpoint()
                #f_scalar                = fuzzyNumber(f_allmid,f_allmid,f_allmid)  #zeros the width! 
                
                this_midx               = midx.value
                this_value              = AF.fuzzyNumber(f_scalar.min.value, 
                                                         f_scalar.real.value, 
                                                         f_scalar.max.value)
                this_gradient           = AF.fuzzyNumber(f_all.min.der[0,gn],
                                                         f_all.real.der[0,gn],
                                                         f_all.max.der[0,gn])
                #nx                      = (this_value*(this_gradient**(-1))*(-1.)) + this_midx #standard
                                                         
                #nx                      = ( (1./this_gradient)*this_value*(-1.) ) + this_midx
                nx                      = ((this_value/this_gradient)*-1) + this_midx
                
                #nx                      = (f_scalar*(this_gradient**(-1))*(-1.)) + midx
                this_update = AF.fuzzyNumber(x_update[xyz][index].min.value, 
                                             x_update[xyz][index].real.value, 
                                             x_update[xyz][index].max.value)
                test_update = nx & this_update     
                if abs(test_update.width()) < tol:
                    not_conv = False
                else:
                    this_update = test_update
                #x_update = nx & x_update
                x_update[xyz][index].min.value  = this_update.min
                x_update[xyz][index].real.value = this_update.real
                x_update[xyz][index].max.value  = this_update.max
                count+=1
        print 'scalar Newtons method complete with {} iterations'.format(count)
        return x_update
        
        
    def montonic_compute_scalar_newton(self, 
                                       func, 
                                       vertices, 
                                       vertex_tpl, 
                                       max_it = 30, 
                                       tol = 1.e-6):
        """
            monotonicity checking box consistency
            (monotonic interval newton iteration on a single DOF)
        """
        x_update = copy.deepcopy(vertices)
        not_conv = True
        count = 0
        xyz = vertex_tpl[0]     # which dimension are we in [0,1,2]?
        index = vertex_tpl[1]   # which vertex are we attempting to contract?
        n = len(vertices[0])
        gn = xyz*n + index
        while(not_conv and count < max_it):
                midx                    = x_update[xyz][index].getpoint(.5)
                x_thin                  = copy.deepcopy(vertices)
                x_thin[xyz][index]      = AF.fuzzyNumber(midx,midx,midx)
                f_scalar                = func(x_thin)  #reduces the width in accordance with thick constraint
                f_all                   = func(x_update)
                
                #f_allmid                = f_all.midpoint()
                #f_scalar                = fuzzyNumber(f_allmid,f_allmid,f_allmid)  #zeros the width! 
                
                this_midx               = midx.value
                this_value              = AF.fuzzyNumber(f_scalar.min.value, 
                                                         f_scalar.real.value, 
                                                         f_scalar.max.value)
                this_gradient           = AF.fuzzyNumber(f_all.min.der[0,gn],
                                                         f_all.real.der[0,gn],
                                                         f_all.max.der[0,gn])
                #nx                      = (this_value*(this_gradient**(-1))*(-1.)) + this_midx
                                                         
                #nx                      = ((1./this_gradient)*this_value*(-1.)) + this_midx
                nx                      = ((this_value/this_gradient)) + this_midx
                
                #nx                      = (f_scalar*(this_gradient**(-1))*(-1.)) + midx
                this_update = AF.fuzzyNumber(x_update[xyz][index].min.value, 
                                             x_update[xyz][index].real.value, 
                                             x_update[xyz][index].max.value)
                if np.isnan(nx.min) or np.isnan(nx.max):
                    test_update = this_update
                    return x_update
                test_update = nx & this_update     
                if abs(test_update.width()) < tol:
                    not_conv = False
                else:
                    this_update = test_update
                #x_update = nx & x_update
                x_update[xyz][index].min.value  = this_update.min
                x_update[xyz][index].real.value = this_update.real
                x_update[xyz][index].max.value  = this_update.max
                count+=1
        print 'scalar Newtons method complete with {} iterations'.format(count)
        return x_update
        
        
    def compute_interval_newton(self, f, thin_f, L, vertices):
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
        x_update        = np.asarray(x_update) 
        
        
        
        midx            = self.mid_vector(x_update)   #vector of scalarAD
        this_hess       = self.generate_hessian(f)    #thick Interval matrix -non AD-   

        #"""
        #key : must pass a thin interval to generate a thin gradient
        #       artificially thinning a thick gradient will not do! 
        # note moreover, and maybe all that above is junk, because a gradient of a thin f, isn't in general thin.
        grad            = self.generate_thin_gradient(f=thin_f, size=self.n)   #thin vector of -non AD- Intervals
        
        Vint            = x_update - midx
        
        usehessian      = this_hess[self.mask,:]
        usehessian      = usehessian[:,self.mask]
        
        nx              = self.guass_sidel(usehessian,np.transpose(Vint[self.mask]),np.transpose(grad[self.mask]), nmax=10, tol = 1.e-6)
        
        nx              = nx*(-1.) + midx[self.mask]           
        test_updateN    = self.vector_and(x_update[self.mask], nx)
        
        return test_updateN
        
        
        
    def interval_krawczyk_inversion(self, f, Lagrangian, vertices):  #|TLM1> #|000> #dino 2
        """
            Invertval Function Inversion for Newton Iteration
        """
        size = self.n

        #f = self.interval_f
        x = vertices[0]
        y = vertices[1]

        identity_matrix = np.identity(size)

        x_update = []
        for el in x:
            x_update.append(el)
        for el in y:
            x_update.append(el)
        for el in Lagrangian.equality:
            l = Lagrangian.equality[el].interval_Lagrange
            x_update.append(l)
        x_update    = np.asarray(x_update)   #vector Interval scalarAD
        midx        = self.mid_vector(x_update)   #vector of scalarAD
        this_hess   = self.generate_hessian(f)    #thick Interval matrix -non AD-    
        
        ##
        #        fy = []
        #        for i in range(size):
        #            fy.append(f.real.der[0,i])  #this isn't necessarily in the middle...
        #        fy = np.asarray(fy)             #vector of -non AD- scalars
        
        #        #fy = generate_thin_gradient(f, size, self.mask, rsize)  #makes a thin Interval -non AD- gradient

        f_mid = f.midpoint()
        fy = []
        for i in range(size):
            fy.append(f_mid.der[0,i])
        fy = np.asarray(fy)
        
        #mJ              = mid_hessian(this_hess, size)
        mJ              = self.generate_preconditioner(f)
        Y               = np.linalg.inv(mJ)
        Yfy             = np.dot(Y,np.transpose(fy))
        p1              = ( np.dot( np.transpose(this_hess),np.transpose(Y)) * (-1.) + identity_matrix )
        p2              = x_update - midx    
        kx              = (np.dot(p1,p2))+ midx - Yfy
        test_updateK    = self.vector_and(x_update[self.mask], kx[self.mask])
    
        return test_updateK
        
        
        
    def mid_vector(self, x):
        mv=[]
        for i in range(len(x)):
            mv.append(x[i].midpoint())
        mv=np.asarray(mv)
        return mv
        
    def matrix_mid(self, A):
        midA = []
        for i in range(len(A)):
            midA.append([])
            for j in range(len(A)):
                midA[i].append(A[i,j].midpoint())
        midA = np.asarray(midA)
        return midA
        
    def generate_hessian(self, f):
        a = f.min.hess
        b = f.real.hess
        c = f.max.hess
        size = len(f.min.hess)
        hessian = []
        for i in range(size):
            hessian.append([])
            for j in range(size):
                gmax = max(a[i,j],b[i,j],c[i,j])
                gmin = min(a[i,j],b[i,j],c[i,j])
                gmid = (gmax + gmin)*0.5
                hessian[i].append(AF.fuzzyNumber(gmin,gmid,gmax))
        hessian = np.asarray(hessian)
        return hessian
        
    def mid_hessian(self, hessian, N=None):
        """
            returns a thin matrix of the
            computed mid of each 
            fuzzy number in the hessian
            
            each element may be of AD type
            
            h.midpoint()
        """
        if N == None:
            N = len(hessian)
        midh = []
        for i in range(N):
            midh.append([])
            for j in range(N):
                midh[i].append(hessian[i,j].midpoint() )
        midh = np.asarray(midh)
        return midh
        
    def generate_gradient(self,f_scalar): #|TLM6>
        """
            scalar Interval of AD objects 
            returns interval vector gradient
        """
        a = f_scalar.min.der
        b = f_scalar.real.der
        c = f_scalar.max.der
        gradient = []
        size = len(f_scalar.min.hess)
        for i in range(size):
            gmax = max(a[0,i],b[0,i],c[0,i])
            gmin = min(a[0,i],b[0,i],c[0,i])
            gmid = (gmin+gmax)*0.5
            gradient.append(AF.fuzzyNumber(gmin,gmid,gmax))
        gradient = np.asarray(gradient)
        return gradient
        
    def hansen_gradient(self,f):  #hansen
        """
            accepts: scalarAD valued objective func
            
            returns: the thin interval gradient
            
            TODO: should account for outward rounding
        """
        hansen_grad = []  #Hansen page 186, page 254
        for i in range(len(f.hess)):
            hansen_grad.append(AF.fuzzyNumber(f.der[0,i],f.der[0,i],f.der[0,i]))
        hansen_grad = np.asarray(hansen_grad)
        return hansen_grad
        
    def generate_thin_gradient(self, f, size=None, mask = None, rsize = None):
        """
            ERROR!
            scalar Interval of AD objects 
            returns an interval gradient vector
            at the midpoint of interval_f
        """

        a = f.min
        b = f.real
        c = f.max
        if size ==None:
            size = len(a.hess)
        agrad = a.der
        bgrad = b.der
        cgrad = c.der
        gradient = []
        for i in range(size):
            ag = agrad[0,i]
            bg = bgrad[0,i]
            cg = cgrad[0,i]
            v_min = min(ag,bg,cg)
            v_max = max (ag,bg,cg)
            v_mid = (v_min+v_max)*0.5
            gradient.append(AF.fuzzyNumber(v_min,v_mid,v_max))
        gradient = np.asarray(gradient)

        return gradient
        
    def generate_preconditioner(self, f_all):
        """
            takes a scalarAD 
                e.g.  F.real.hess 
            returns a scalarized component matrix
        """
        hessian = []
        a = f_all.min.hess
        b = f_all.real.hess
        c = f_all.max.hess
        size = len(f_all.min.hess)
        for i in range(size):
            hessian.append([])
            for j in range(size):
                gmax = max(a[i,j],b[i,j],c[i,j])
                gmin = min(a[i,j],b[i,j],c[i,j])
                gmid = (gmin+gmax)*0.5
                hessian[i].append(gmid)
        hessian = np.asarray(hessian)
        return hessian
        
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
        #FormParam       = L
        xpts = []
        ypts = []
        temp_Lagrange = copy.deepcopy(L)
        #min_temp_Lagrange = []
        #max_temp_Lagrange = []
        
        #self.count_fv   = 0  #free variable index (corresponds to a mask location)
        for i in range(L.curve.n):
            thin_xpt = vertices[0][i].getpoint(pt)
            xpts.append(AF.fuzzyNumber(thin_xpt,thin_xpt,thin_xpt))
        
            thin_ypt = vertices[1][i].getpoint(pt)
            ypts.append(AF.fuzzyNumber(thin_ypt, thin_ypt, thin_ypt))
            #self.count_fv +=2

        for el in temp_Lagrange.equality:
            thin_fp = L.equality[el].interval_Lagrange.getpoint(pt)
            temp_Lagrange.equality[el].interval_Lagrange = AF.fuzzyNumber(thin_fp,thin_fp,thin_fp)
            #self.count_fv +=1
        
        thin_vertices = [xpts,ypts]
        return thin_vertices, temp_Lagrange
        
        
    def vector_and(self,xi, delta):
        return_values = copy.deepcopy(xi)
        for i in range(len(xi)):
            #element of xi
            a = xi[i].min.value
            b = xi[i].real.value
            c = xi[i].max.value
            x_simple = AF.fuzzyNumber(a,b,c)
            
            #element of delta
            d = delta[i].min.value
            e = delta[i].real.value
            f = delta[i].max.value
            dz_simple = AF.fuzzyNumber(d,e,f)
            
            new_values = x_simple & dz_simple
            prove_zero = dz_simple in x_simple
            prove_stuck = new_values.width() == x_simple.width()
            
            return_values[i].min.value   = new_values.min
            return_values[i].max.value   = new_values.max
            return_values[i].real.value  = new_values.midpoint()
            return_values[i].isempty     = new_values.isempty
            return_values[i].prove_zero  = prove_zero
            return_values[i].prove_nosol = new_values.isempty
            return_values[i].prove_stuck = prove_stuck
        return return_values
        
    def monotonic_contraction(self, precomputed_f, vertices, vmap):
        """
            see An Interval Constraint Propogation Algorithm
            Exploiting Monoticity page 4, paragraph 5.
            by A.N.T. 
            
            f : some constraint function evalutated on the intervals, 
                    especially a monotonic one
            vmap : maps a linear pt index to the curve vertex list of pt lists
        """
        nv_max = copy.deepcopy(vertices)
        nv_min = copy.deepcopy(vertices)

        decreasing, increasing, flat = precomputed_f.eval_monoticity()

        for i in range(len(vmap)):
            ivtx = vmap[i]
            deg_max = AF.fuzzyNumber(vertices[ ivtx[0] ][ ivtx[1] ].max,
                                     vertices[ ivtx[0] ][ ivtx[1] ].max,
                                     vertices[ ivtx[0] ][ ivtx[1] ].max)
            deg_min = AF.fuzzyNumber(vertices[ ivtx[0] ][ ivtx[1] ].min,
                                     vertices[ ivtx[0] ][ ivtx[1] ].min,
                                     vertices[ ivtx[0] ][ ivtx[1] ].min)
            if decreasing[i]:
                nv_min[ ivtx[0] ][ ivtx[1] ] = deg_max
                nv_max[ ivtx[0] ][ ivtx[1] ] = deg_min
            elif increasing[i]:
                nv_min[ ivtx[0] ][ ivtx[1] ] = deg_min
                nv_max[ ivtx[0] ][ ivtx[1] ] = deg_max
            else: ## not monotonic keep same as old
                pass
        return     nv_min, nv_max
        
    def check_for_stationary_pt(self, grad_f):
        return np.transpose(np.asarray([grad_f[:,0] <= 0.,0. <= grad_f[:,2]]))
        
    def check_for_convergence(self,grad_f,tol):
        test = np.transpose(np.asarray([-tol <= grad_f[:,0], tol >= grad_f[:,2] ]) )
        if False in test:
            return False
        else:
            return True
        
    def convex_trace_check(self, f):
        x=f.ihess
        ct = []
        for i in range(len(x)):
            ct.append(0<=x[i,i,0])
        return
        