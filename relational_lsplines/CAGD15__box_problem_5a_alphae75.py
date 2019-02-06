# -*- coding: utf-8 -*-
"""
Created on Fri Apr  1 22:36:43 2016

@author: lukemcculloch
"""
import numpy as np
import matplotlib.pyplot as plt
#from mpl_toolkits.axes_grid.axislines import SubplotZero
#import matplotlib.patches as patches
#import matplotlib.path    as path


import curve             as     spline
import knot_operations   as     ko
from   initialValues     import InitializeControlPoints, InitializeControlVertices
import copy
from   ADILS             import IntervalLagrangeSpline, Lagrangian
from   FormParameter     import FormParameterDict
from   initialValues     import interval_bounds, lagrangian_bounds

from minDistance import naive_minimizer

np.set_printoptions(threshold= 200)

class box_constraint(object):
    
    def __init__(self,curve,box):
        self.curve  = curve
        self.box    = box
        self.tol = 1.e-7
        return
    
    def find_box_pts_of_interest(self):
        
        return
        
    def ADILS_optiization(self,curve,yval,xval,x2,y2,
                   tanb=0.,tane=0.,cvb=0.,cve=0.,area=72.):
        interval_data, small = interval_bounds(curve)        
        FPD = FormParameterDict(curve) 
        #---box constraint
        FPD.add_yPointConstraint(kind= 'max', value=4., location=xval[0] )
        FPD.add_xPointConstraint(kind= 'min', value=6., location=xval[0] )
        #FPD.add_yPointConstraint(kind= 'LS', value=4., location=yval[0], weight=-100. )
        #FPD.add_xPointConstraint(kind= 'LS', value=6., location=xval[0], weight=-100. )
        #---
        #---box constraint
        FPD.add_yPointConstraint(kind= 'max', value=8., location=y2[0] )
        FPD.add_xPointConstraint(kind= 'min', value=6., location=x2[0] )
        #FPD.add_yPointConstraint(kind= 'LS', value=8., location=y2[0], weight=-100. )
        #FPD.add_xPointConstraint(kind= 'LS', value=6., location=x2[0], weight=-100. )
        #---
        FPD.add_AreaConstraint(kind='LS',value=area)
        FPD.add_AngleConstraint(kind='equality',location=0.,value=tanb)
        FPD.add_AngleConstraint(kind='equality',location=1.,value=tane)
        FPD.add_CurvatureConstraint(kind='equality',location=0.,value=cvb)
        FPD.add_CurvatureConstraint(kind='equality',location=1.,value=cve)
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize(stop = 50)
        return Lspline
        
    def iteration(self):
        curve = copy.deepcopy(self.curve)
        #curve = Lspline.curve
        y_interest = 4.0
        yindex = 1
        yval = curve.FindPoint(y_interest,yindex)
        
        x_interest = 6.
        x_index = 0
        xval = curve.FindPoint(x_interest,x_index)
        
        y2 = curve.FindPoint(8.,yindex)
        x2 = curve.FindPoint(6.,yindex)
        
        while (abs(xval[0]-yval[0])>self.tol):
            
            
            
            Lspline = self.ADILS_optiization(curve,
                                              xval,
                                              yval,
                                              x2,
                                              y2,
                                              tanb = 0.,
                                              tane = 0.,
                                              cvb = 0.,
                                              cve = 0.,
                                              area = 66.)
            curve = Lspline.curve
            
            y_interest = 4.0
            yindex = 1
            yval = curve.FindPoint(y_interest,yindex)
            
            x_interest = 6.
            x_index = 0
            xval = curve.FindPoint(x_interest,x_index)
            
            y2 = curve.FindPoint(8.,yindex)
            x2 = curve.FindPoint(6.,yindex)
            
        
        return Lspline
        
        
        
class CP_Intersection(object):
    """
        d : curve degree
    """
    def __init__(self, curve1, curve2):
        self.curve1 = curve1
        self.d1 = curve1.k-1
        self.curve2 = curve2
        self.d2 = curve2.k-1
        self.curves = [curve1,curve2]
        return
#    def greville_abscissa(self,i,j):
#        """G.A. : the knot average 
#        to corresponding control poin c_{i}
#            of the
#              - ith knot
#              - jth curve
#        """
#        curve = self.curves[j]
#        d = curve.k-1
#        ga = 0.
#        for t in range(i+1,i+d+1):
#            ga += self.curves.t[t]
#        ga = ga/d
#        return ga
    
    def eval_control_polygon(self, i,j,t):
        """
            i : ith knot
            j : jth curve
            t : greville abcissa
            t E [ t(i)_{ga} , t(i+1)_{ga} ]
        """
        curve = self.curves[j]
        tib = self.greville_abscissa(i,j)
        tibp = self.greville_abscissa(i+1,j)
        ci = curve.vertices[i]
        cip = curve.vertices[i+1]
        g = ((tibp-t)/(tibp-tib))*ci
        g += ((t-tib)/(tibp-tib))*cip
        return g
        
    def intersect_cp_segments(self,i,j):
        """
            i : control vertex num i E {n:nE[1,curve1.n-1}
                ie segment i
            2 : control vertex num j E {n:nE[1,curve2.n-1}
                ie segment j
        """
        b   = self.curve1.vertices
        c   = self.curve2.vertices
        si1 = self.curve1.greville_abscissa(i)
        si2 = self.curve1.greville_abscissa(i-1)
        tj1 = self.curve2.greville_abscissa(j)
        tj2 = self.curve2.greville_abscissa(j-1)
        db  = (b[i]-b[i-1])/(si1-si2)
        dc  = (c[j]-c[j-1])/(tj1-tj2)
        M   = np.asarray([db,-dc]).T
        Mi= np.linalg.inv(M)
        diff = b[i]-c[j]
        N = np.asarray([si1,tj1])
        
        return  N-Mi.dot(diff)
        
    def do_box_constraint(self, box):
        return

class box(object):
    
    def __init__(self, pts):
        assert(len(pts)==4), 'generalized box takes 4 pts'
        self.pts = pts
        self.k = 2
        self.nump = 15
        self.curves = [] 
        self.make_box()
        self.interest_x = self.pts[1][0]
        self.interest_y = self.pts[0][1]
        return
        
    def make_box(self):
        c1 = spline.Bspline(np.asarray([self.pts[0],self.pts[1]]),self.k,self.nump)
        c2 = spline.Bspline(np.asarray([self.pts[1],self.pts[2]]),self.k,self.nump)
        c3 = spline.Bspline(np.asarray([self.pts[2],self.pts[3]]),self.k,self.nump)
        c4 = spline.Bspline(np.asarray([self.pts[3],self.pts[0]]),self.k,self.nump)
        self.curves = [c1,c2,c3,c4]
        return
        
    def pt(self, s):
        if 0.<=s<=.25:
            return self.curves[0].CurvePoint(u=s/.25)
        if .25<=s<=.5:
            return self.curves[0].CurvePoint(u=(s-.25)/.25)
        if .5<=s<=.75:
            return self.curves[0].CurvePoint(u=(s-.5)/.25)
        if .75<=s<=1.:
            return self.curves[0].CurvePoint(u=(s-.75)/.25)
        else:
            print 'error, s = {u : u E [0.,1.]}, but s = {}'.format(s)
        return
    
    def plot(self):
        for cv in self.curves:
            cv.plot()
        return
        
if __name__ == "__main__":
    #curve.FindNearestPoint(point = pt_to_find)
    a=[0.,4.]
    b=[4.,4.]
    c=[8.,8.]
    d=[0.,8.]
    box1 = box([a,b,c,d])
    do_all=False

    test = 'smart_section'
    #test = 'cusp'
    
    if do_all:
        test = 'cusp'
    if test == 'cusp':
        print 'test = {}'.format(test)
        #print 'option = {}'.format(option)
        k=4
        nump=30
        xb = 0.
        yb = 0.
        xe = 12.
        ye = 12.
        alphab = 0.#-5.
        alphae = 0.
        Cab = 0.
        Cae = 0.
        xc = 4.
        yc = 4.
        curve_area = 72.
        slope = 'up'
        
        ae = alphae
        ab = alphab
        
        numcv = 9 #8 #10
        ini_v = InitializeControlVertices(xb,yb,xe,ye,alphae=ae, alphab=ab,
                                          Cab_given=Cab,Cae_given=Cae,
                                           nCV = numcv, slope = 'up')

        
        curve = spline.Bspline(ini_v.vertices, k, nump)  
        curve.knot_insertion(knot=.55)
        curve.knot_insertion(knot=.55)
        curve.knot_insertion(knot=.55)
        curve.knot_insertion(knot=.55)
        curve.knot_insertion(knot=.55)
        curve.knot_insertion(knot=.55)
        curve.knot_insertion(knot=.55)
        curve.knot_insertion(knot=.55)
        curve.knot_insertion(knot=.55)
        c = spline.Bspline(curve.vertices, k, 3*nump)  
        curve.plotcurve_detailed()

    if do_all:
        test = 'smart_section'
    if test == 'smart_section':
        print 'test = {}'.format(test)
        #print 'option = {}'.format(option)
        k=4
        nump=30
        xb = 0.
        yb = 0.
        xe = 10. #12.
        ye = 12. #10.
        alphab = 0.#-5.
        alphae = 75.#0.
        Cab = 0.
        Cae = 0.
        xc = 4.
        yc = 4.
        curve_area = 31. #55.#72.
        slope = 'up'
        
        ae = alphae
        ab = alphab
        
        numcv = 9 #8 #10
        ini_v = InitializeControlVertices(xb,yb,xe,ye,alphae=ae, alphab=ab,
                                          Cab_given=Cab,Cae_given=Cae,
                                           nCV = numcv, slope = 'up')

        
        curve = spline.Bspline(ini_v.vertices, k, nump)  
        #curve.compute_area()
        #curve_area = curve.area
        c = spline.Bspline(ini_v.vertices, k, 3*nump)  
        curve.plotcurve_detailed()
        box1.plot()
        
        
        
        #intersect c.p.'s
        #
        ilist = []
        i=3
        for j in range(1,4+1):
            self = CP_Intersection(curve,box1.curves[j-1])
            p1,p2=self.intersect_cp_segments(i,1)
            ilist.append([p1,p2])
        i=3
        j=1
        #
        # use the 
        #     intersect_cp_segments
        # to find the orrect box segment instead
        #
#        self = CP_Intersection(curve,box1.curves[0])
#        p1,p2=self.intersect_cp_segments(i=3,j=1)
#        
#        c1 = copy.deepcopy(curve)
#        cspan = c1.knot_insertion(knot=p1)
#        b1 = copy.deepcopy(box1.curves[0])
#        bspan = b1.knot_insertion(knot=p2)
#        
#        self = CP_Intersection(c1,box1.curves[0])
#        p1,p2=self.intersect_cp_segments(i=4,j=1)
        
        
        
        #"""
        #my way
        #"""
        y_interest = 4.0
        yindex = 1
        yval0 = curve.FindPoint(y_interest,yindex)
        
        x_interest = 8.
        x_index = 0
        xval0 = curve.FindPoint(x_interest,x_index)
        
        
        
        #print pval[0]
        #print curve.CurvePoint(yval[0])
        #"""
        ccr = copy.deepcopy(curve)
        interval_data, small = interval_bounds(ccr)        
        FPD = FormParameterDict(ccr) 
        #"""
        """
        FPD.add_yPointConstraint(kind= 'max', value=4., location=.4 )
        FPD.add_xPointConstraint(kind= 'min', value=6., location=.45 )
        FPD.add_yPointConstraint(kind= 'LS', value=4., location=.4,weight=-100. )
        FPD.add_xPointConstraint(kind= 'LS', value=6., location=.45,weight=-100. )
        #"""
        FPD.add_yPointConstraint(kind= 'max', value=4., location=yval0[0] )
        FPD.add_xPointConstraint(kind= 'min', value=8., location=xval0[0] )
        #
        #FPD.add_yPointConstraint(kind= 'min', value=8., location=xval0[0] )#new and not needed
        #"""
        FPD.add_AreaConstraint(kind='equality',value=curve_area)
        FPD.add_AngleConstraint(kind='equality',location=0.,value=alphab)
        FPD.add_AngleConstraint(kind='equality',location=1.,value=alphae)
        FPD.add_CurvatureConstraint(kind='equality',location=0.,value=Cab)
        FPD.add_CurvatureConstraint(kind='equality',location=1.,value=Cae)
        #"""
        #"""
        FPD.add_E1(kind='LS', weight = .5)#.05)
        FPD.add_E2(kind='LS', weight = .5)#1.)
        FPD.add_E3(kind='LS', weight = .5)#01)
        
        
        #FPD.add_E1(kind='LS', weight = .5)
        #FPD.add_E2(kind='LS', weight = .5)
        #FPD.add_E3(kind='LS', weight = .5)
        
        FPD.add_ArcLengthApprox(kind='LS', weight = 100.)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(ccr, L, data = interval_data)
        Lspline.optimize(stop = 50)
        #"""
        L0 = copy.deepcopy(Lspline)
        L0.curve.plotcurve_detailed()
        L0.curve.plot_area_to_x_shiftaxis(shift_axis=12.)
        L0.curve.plotCurvature_expensive(scale=3.5,
                                         nump=5,
                                         flip=True)
        box1.plot()
        #L0.curve.plotCurvature_spines()
        #L0.curve.plotCurvature_nospines(alpha = 0.)
        
        
        
        Lspline.curve.pts_M_pts()
        Lspline.curve.compute_arclength()
        
        E1_0 = copy.deepcopy(Lspline.curve.E1)
        E2_0 = copy.deepcopy(Lspline.curve.E2)
        E3_0 = copy.deepcopy(Lspline.curve.E3)
        S_0 = copy.deepcopy(Lspline.curve.AL)
        
        
        
        
        
        
        
        
        
        
        
        if True:
            #my way
            #"""
            x_interest = 4.
            x_index = 0
            xval4 = curve.FindPoint(x_interest,x_index)
            
            y_interest = 4.0
            yindex = 1
            yval4 = curve.FindPoint(y_interest,yindex)
            
            
            x_interest = 6.
            x_index = 0
            xval6 = curve.FindPoint(x_interest,x_index)
            
            y_interest = 6.0
            yindex = 1
            yval6 = curve.FindPoint(y_interest,yindex)
            
            x_interest = 8.
            x_index = 0
            xval8 = curve.FindPoint(x_interest,x_index)
            
            y_interest = 8.0
            yindex = 1
            yval8 = curve.FindPoint(y_interest,yindex)
            print '\n------------------------------------------\n'
            print '(1)  yval8:',yval8
            print '\n------------------------------------------\n'
            
            
            
            ##
            ##  45 area best box tangent
            ##
            #print pval[0]
            #print curve.CurvePoint(yval[0])
            #"""
            c1 = copy.deepcopy(curve)
            interval_data, small = interval_bounds(c1)        
            FPD = FormParameterDict(c1) 
            #"""
            """
            FPD.add_yPointConstraint(kind= 'max', value=4., location=.4 )
            FPD.add_xPointConstraint(kind= 'min', value=6., location=.45 )
            FPD.add_yPointConstraint(kind= 'LS', value=4., location=.4,weight=-100. )
            FPD.add_xPointConstraint(kind= 'LS', value=6., location=.45,weight=-100. )
            #"""
            FPD.add_yPointConstraint(kind= 'max', value=4., location=yval4[0] )
            FPD.add_yPointConstraint(kind= 'min', value=8., location=yval8[0] )
            #FPD.add_xPointConstraint(kind= 'min', value=4., location=yval4[0] )
            FPD.add_xPointConstraint(kind= 'min', 
                                     value=8., 
                                     location=yval8[0],
                                     weight=1.)
            
            
            #FPD.add_xPointConstraint(kind= 'min', value=4., location=xval4[0] )
            ##
            #TLM Sept good one:
            #FPD.add_xPointConstraint(kind= 'LS', 
            #                         value=8., 
            #                         location=yval8[0],
            #                         weight=10.)
            
            
            #FPD.add_xPointConstraint(kind= 'equality', value=8., location=xval8[0] )
            
            #FPD.add_yPointConstraint(kind= 'min', value=8., location=.76)#yval8[0] )
            #FPD.add_xPointConstraint(kind= 'min', value=4., location=yval4[0] )
            #FPD.add_xPointConstraint(kind= 'min', value=8., location=.76)#yval8[0] )
            
            #FPD.add_yPointConstraint(kind= 'max', value=6., location=yval6[0] )
            #FPD.add_xPointConstraint(kind= 'min', value=6., location=yval6[0] )
            
            
            
            #FPD.add_yPointConstraint(kind= 'max', value=4., location=yval1[0] )
            #FPD.add_xPointConstraint(kind= 'min', value=8., location=xval1[0] )
            #"""
            #FPD.add_AreaConstraint(kind='equality',value=48.)
            FPD.add_Area_to_y_Constraint(kind='equality',value=72.)
            FPD.add_AngleConstraint(kind='equality',location=0.,value=alphab)
            FPD.add_AngleConstraint(kind='equality',location=1.,value=alphae)
            FPD.add_CurvatureConstraint(kind='equality',location=0.,value=Cab)
            FPD.add_CurvatureConstraint(kind='equality',location=1.,value=Cae)
            #"""
            #"""
            FPD.add_E1(kind='LS', weight = .5/E1_0.value)
            FPD.add_E2(kind='LS', weight = .5/E2_0.value)
            #FPD.add_E3(kind='LS', weight = .01)#works
            #FPD.add_E3(kind='LS', weight = .5)
            FPD.add_E3(kind='LS', weight = .5/E3_0.value)
            
            #FPD.add_E1(kind='LS', weight = .5)
            #FPD.add_E2(kind='LS', weight = .5)
            #FPD.add_E3(kind='LS', weight = .5)
            
            FPD.add_ArcLengthApprox(kind='LS', weight = .5/S_0)
            #FPD.add_ArcLengthApprox(kind='LS', weight = .5)
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data)#, small, 1.e4)
            Lspline = IntervalLagrangeSpline(c1, L, data = interval_data)
            Lspline.optimize(stop = 50)
            
            L1 = copy.deepcopy(Lspline)
            L1.curve.plotcurve_detailed()
            L1.curve.plotCurvature_expensive(scale=3.5,
                                             nump=5,
                                             flip=True)
            L1.curve.plot_area_to_x_shiftaxis(shift_axis=12.)
            box1.plot()
            """This one just touches the box edge"""
            #L1.curve.plotCurvature_spines()
            #L1.curve.plotCurvature_nospines(alpha = 0.)
            
            #cv = spline.Bspline(L1.curve.vertices,k=4,nump=L1.curve.nump*5)
            #cv.plotCurvature_nospines(scale=1.0, factor=1, alpha=0.)
            #cv.plotCurvature_spines(scale=1.0)
            
        
        
        
        
        
        
        
        
        
        
        
        
        
        if False:
            #my way
            #"""
            x_interest = 4.
            x_index = 0
            xval4 = curve.FindPoint(x_interest,x_index)
            
            y_interest = 4.0
            yindex = 1
            yval4 = curve.FindPoint(y_interest,yindex)
            
            
            x_interest = 6.
            x_index = 0
            xval6 = curve.FindPoint(x_interest,x_index)
            
            y_interest = 6.0
            yindex = 1
            yval6 = curve.FindPoint(y_interest,yindex)
            
            x_interest = 8.
            x_index = 0
            xval8 = curve.FindPoint(x_interest,x_index)
            
            y_interest = 8.0
            yindex = 1
            yval8 = curve.FindPoint(y_interest,yindex)
            print '\n------------------------------------------\n'
            print '(2)  yval8:',yval8
            print '\n------------------------------------------\n'
            
            ##
            ##  45 area best box tangent
            ##
            #print pval[0]
            #print curve.CurvePoint(yval[0])
            #"""
            c1 = copy.deepcopy(curve)
            interval_data, small = interval_bounds(c1)        
            FPD = FormParameterDict(c1) 
            #"""
            """
            FPD.add_yPointConstraint(kind= 'max', value=4., location=.4 )
            FPD.add_xPointConstraint(kind= 'min', value=6., location=.45 )
            FPD.add_yPointConstraint(kind= 'LS', value=4., location=.4,weight=-100. )
            FPD.add_xPointConstraint(kind= 'LS', value=6., location=.45,weight=-100. )
            #"""
            FPD.add_yPointConstraint(kind= 'max', value=4., location=yval4[0] )
            FPD.add_yPointConstraint(kind= 'min', value=8., location=yval8[0] )
            #FPD.add_xPointConstraint(kind= 'min', value=4., location=xval4[0] )
            FPD.add_xPointConstraint(kind= 'min', value=8., location=yval8[0] ) #orginal good one TLM Sept
            
            
            #FPD.add_xPointConstraint(kind= 'min', value=4., location=xval4[0] )
            ##
            #TLM Sept good one:
            #FPD.add_xPointConstraint(kind= 'equality', value=8., location=yval8[0] )
            #
            #
            #FPD.add_xPointConstraint(kind= 'equality', value=8., location=xval8[0] )
            
            #FPD.add_yPointConstraint(kind= 'min', value=8., location=.76)#yval8[0] )
            #FPD.add_xPointConstraint(kind= 'min', value=4., location=yval4[0] )
            #FPD.add_xPointConstraint(kind= 'min', value=8., location=.76)#yval8[0] )
            
            #FPD.add_yPointConstraint(kind= 'max', value=6., location=yval6[0] )
            #FPD.add_xPointConstraint(kind= 'min', value=6., location=yval6[0] )
            
            
            
            #FPD.add_yPointConstraint(kind= 'max', value=4., location=yval1[0] )
            #FPD.add_xPointConstraint(kind= 'min', value=8., location=xval1[0] )
            #"""
            #FPD.add_AreaConstraint(kind='equality',value=48.)
            FPD.add_Area_to_y_Constraint(kind='equality',value=72.)
            FPD.add_AngleConstraint(kind='equality',location=0.,value=alphab)
            FPD.add_AngleConstraint(kind='equality',location=1.,value=alphae)
            FPD.add_CurvatureConstraint(kind='equality',location=0.,value=Cab)
            FPD.add_CurvatureConstraint(kind='equality',location=1.,value=Cae)
            #"""
            #"""
            FPD.add_E1(kind='LS', weight = .5)
            FPD.add_E2(kind='LS', weight = .5)
            FPD.add_E3(kind='LS', weight = .01)
            
            #FPD.add_E1(kind='LS', weight = .5)
            #FPD.add_E2(kind='LS', weight = .5)
            #FPD.add_E3(kind='LS', weight = .5)
            
            FPD.add_ArcLengthApprox(kind='LS', weight = .5)
            #FPD.add_ArcLengthApprox(kind='LS', weight = .5)
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data)#, small, 1.e4)
            Lspline = IntervalLagrangeSpline(c1, L, data = interval_data)
            Lspline.optimize(stop = 50)
            
            L1 = copy.deepcopy(Lspline)
            L1.curve.plotcurve_detailed()
            L1.curve.plotCurvature_expensive(scale=3.5,
                                             nump=5,
                                             flip=True)
            L1.curve.plot_area_to_x_shiftaxis(shift_axis=12.)
            box1.plot()
            """This one just touches the box edge"""
            #L1.curve.plotCurvature_spines()
            #L1.curve.plotCurvature_nospines(alpha = 0.)
            
            #cv = spline.Bspline(L1.curve.vertices,k=4,nump=L1.curve.nump*5)
            #cv.plotCurvature_nospines(scale=1.0, factor=1, alpha=0.)
            #cv.plotCurvature_spines(scale=1.0)
            
        
        
        
        
        
        
        
        
        
        
        if False:
            #my way
            #"""
            x_interest = 4.
            x_index = 0
            xval4 = curve.FindPoint(x_interest,x_index)
            
            y_interest = 4.0
            yindex = 1
            yval4 = curve.FindPoint(y_interest,yindex)
            
            
            x_interest = 6.
            x_index = 0
            xval6 = curve.FindPoint(x_interest,x_index)
            
            y_interest = 6.0
            yindex = 1
            yval6 = curve.FindPoint(y_interest,yindex)
            
            x_interest = 8.
            x_index = 0
            xval8 = curve.FindPoint(x_interest,x_index)
            
            y_interest = 8.0
            yindex = 1
            yval8 = curve.FindPoint(y_interest,yindex)
            print '\n------------------------------------------\n'
            print '(3)  yval8:',yval8
            print '\n------------------------------------------\n'
            
            ##
            ##  45 area best box tangent
            ##
            #print pval[0]
            #print curve.CurvePoint(yval[0])
            #"""
            c1 = copy.deepcopy(curve)
            interval_data, small = interval_bounds(c1)        
            FPD = FormParameterDict(c1) 
            #"""
            """
            FPD.add_yPointConstraint(kind= 'max', value=4., location=.4 )
            FPD.add_xPointConstraint(kind= 'min', value=6., location=.45 )
            FPD.add_yPointConstraint(kind= 'LS', value=4., location=.4,weight=-100. )
            FPD.add_xPointConstraint(kind= 'LS', value=6., location=.45,weight=-100. )
            #"""
            FPD.add_yPointConstraint(kind= 'max', value=4., location=yval4[0] )
            FPD.add_yPointConstraint(kind= 'min', value=8., location=yval8[0] )
            #FPD.add_xPointConstraint(kind= 'min', value=4., location=xval4[0] )
            FPD.add_xPointConstraint(kind= 'min', value=8., location=yval8[0] ) #orginal good one TLM Sept
            
            
            #FPD.add_xPointConstraint(kind= 'min', value=4., location=xval4[0] )
            ##
            #TLM Sept good one:
            #FPD.add_xPointConstraint(kind= 'equality', value=8., location=xval8[0] )
            #
            #
            #FPD.add_xPointConstraint(kind= 'equality', value=8., location=xval8[0] )
            
            #FPD.add_yPointConstraint(kind= 'min', value=8., location=.76)#yval8[0] )
            #FPD.add_xPointConstraint(kind= 'min', value=4., location=yval4[0] )
            #FPD.add_xPointConstraint(kind= 'min', value=8., location=.76)#yval8[0] )
            
            #FPD.add_yPointConstraint(kind= 'max', value=6., location=yval6[0] )
            #FPD.add_xPointConstraint(kind= 'min', value=6., location=yval6[0] )
            
            
            
            #FPD.add_yPointConstraint(kind= 'max', value=4., location=yval1[0] )
            #FPD.add_xPointConstraint(kind= 'min', value=8., location=xval1[0] )
            #"""
            #FPD.add_AreaConstraint(kind='equality',value=48.)
            FPD.add_Area_to_y_Constraint(kind='equality',value=72.)
            FPD.add_AngleConstraint(kind='equality',location=0.,value=alphab)
            FPD.add_AngleConstraint(kind='equality',location=1.,value=alphae)
            FPD.add_CurvatureConstraint(kind='equality',location=0.,value=Cab)
            FPD.add_CurvatureConstraint(kind='equality',location=1.,value=Cae)
            #"""
            #"""
            FPD.add_E1(kind='LS', weight = .5)
            FPD.add_E2(kind='LS', weight = .5)
            FPD.add_E3(kind='LS', weight = .02)
            
            #FPD.add_E1(kind='LS', weight = .5)
            #FPD.add_E2(kind='LS', weight = .5)
            #FPD.add_E3(kind='LS', weight = .5)
            
            FPD.add_ArcLengthApprox(kind='LS', weight = .5)
            #FPD.add_ArcLengthApprox(kind='LS', weight = .5)
            L = Lagrangian(FPD)
            interval_data, small = lagrangian_bounds(L, interval_data)#, small, 1.e4)
            Lspline = IntervalLagrangeSpline(c1, L, data = interval_data)
            Lspline.optimize(stop = 50)
            
            L1 = copy.deepcopy(Lspline)
            L1.curve.plotcurve_detailed()
            L1.curve.plotCurvature_expensive(scale=3.5,
                                             nump=5,
                                             flip=True)
            L1.curve.plot_area_to_x_shiftaxis(shift_axis=12.)
            box1.plot()
            """This one just touches the box edge"""
            #L1.curve.plotCurvature_spines()
            #L1.curve.plotCurvature_nospines(alpha = 0.)
            
            #cv = spline.Bspline(L1.curve.vertices,k=4,nump=L1.curve.nump*5)
            #cv.plotCurvature_nospines(scale=1.0, factor=1, alpha=0.)
            #cv.plotCurvature_spines(scale=1.0)
        
        
        
        
        
        
        
        
        
        
#==============================================================================
#         if True:
#             ##
#             ##  50 attempt
#             ##  add a vertex...
#             #print pval[0]
#             #print curve.CurvePoint(yval[0])
#             #"""
#             c2 = copy.deepcopy(L1.curve)
#             c2.knot_insertion(.15)
#             interval_data, small = interval_bounds(c2)        
#             FPD = FormParameterDict(c2) 
#             
#             #"""
#             """
#             FPD.add_yPointConstraint(kind= 'max', value=4., location=.4 )
#             FPD.add_xPointConstraint(kind= 'min', value=6., location=.45 )
#             FPD.add_yPointConstraint(kind= 'LS', value=4., location=.4,weight=-100. )
#             FPD.add_xPointConstraint(kind= 'LS', value=6., location=.45,weight=-100. )
#             #"""
#             FPD.add_yPointConstraint(kind= 'max', value=4., location=yval4[0] )
#             FPD.add_yPointConstraint(kind= 'max', value=8., location=yval8[0] )
#             FPD.add_xPointConstraint(kind= 'min', value=4., location=yval4[0] )
#             FPD.add_xPointConstraint(kind= 'min', value=8., location=yval8[0] )
#             
#             #FPD.add_yPointConstraint(kind= 'max', value=6., location=yval6[0] )
#             #FPD.add_xPointConstraint(kind= 'min', value=6., location=yval6[0] )
#             
#             
#             
#             #FPD.add_yPointConstraint(kind= 'max', value=4., location=yval1[0] )
#             #FPD.add_xPointConstraint(kind= 'min', value=8., location=xval1[0] )
#             #"""
#             FPD.add_AreaConstraint(kind='equality',value=50.)
#             FPD.add_AngleConstraint(kind='equality',location=0.,value=alphab)
#             FPD.add_AngleConstraint(kind='equality',location=1.,value=alphae)
#             FPD.add_CurvatureConstraint(kind='equality',location=0.,value=Cab)
#             FPD.add_CurvatureConstraint(kind='equality',location=1.,value=Cae)
#             #"""
#             #"""
#             FPD.add_E1(kind='LS', weight = .5)
#             FPD.add_E2(kind='LS', weight = 1.)
#             #FPD.add_E3(kind='LS', weight = .1)
#             FPD.add_ArcLengthApprox(kind='LS', weight = 400.)
#             L = Lagrangian(FPD)
#             interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
#             Lspline = IntervalLagrangeSpline(c2, L, data = interval_data)
#             Lspline.optimize(stop = 50)
#             
#             L2 = copy.deepcopy(Lspline)
#             L2.curve.plotcurve_detailed()
#             L2.curve.plotCurvature_spines()
#             L2.curve.plotCurvature_nospines(alpha = 0.)
#             box1.plot()
#==============================================================================
