##
## 20140228
##
## TLM
## Compute the point on a curve
## which is the minimum distance
## from a point.
#
## Does not use AD:
## but it could - just extend an AD inerface to use ALLDERS from the curve class
##  in its obj.ders and obj.hess
##    but maybe the instatiation isn't worth the time - lots of calls with lots of point
##      new AD objects each time..???
#


def plot_all(plotlist):

    curveSize = 0.
    fig = plt.figure()
    ax = fig.gca(projection='3d')

    ## Test Objects:
    k=4
    vertices = np.array([[0.,0.,0.],
                        [3.,4.,1.],
                        [-4.,0.,1.],
                        [-1.,4.,2.]])
    #curve  = Bspline(vertices, k)
    exlist = [1,1]

    print ' Plot Objects in plotlist: '
    for item in plotlist:
        try:
            if True:#type(item) == type(curve):
                dummy = item.vertices[-1]-item.vertices[0]
                curveSize = max(curveSize,np.linalg.norm(dummy))
                dim = item.dim
                if dim == 3:
                    print 'plotting 3D curve'
                    try:
                        ax.plot(item.ivertices[:,0],item.ivertices[:,1],item.ivertices[:,2],marker = "o",  linestyle = "none", color='k')
                    except:
                        pass
                    ax.plot(item.vertices[:,0],item.vertices[:,1],item.vertices[:,2],marker = "o",  linestyle = "--", color='k')
                    ax.plot(item.r[:,0], item.r[:,1], item.r[:,2], label='parametric curve', color='b')

                if dim == 2:
                    print 'plotting 2D curve'
                    try:
                        ax.plot(item.ivertices[:,0],item.ivertices[:,1],marker = "o",  linestyle = "none", color='k')
                    except:
                        pass
                    ax.plot(item.vertices[:,0],item.vertices[:,1],marker = "o",  linestyle = "--", color='k')
                    ax.plot(item.r[:,0], item.r[:,1],  label='parametric curve', color='b')
                
                
            elif type(item) == type(exlist):
                print 'plotting points'
                ax.plot(item[0],item[1],item[2], marker = "o",  linestyle = "none", color='r')
        except:
            if len(item)==3:
                ax.plot(item[0],item[1],item[2], marker = "o",  linestyle = "none", color='r')
            elif len(item)==2:
                ax.plot(item[0],item[1], marker = "o",  linestyle = "none", color='r')
#        elif type(item) == type(vertices[0]):
#            
#            if len(vertices[0])==3:
#                print 'plotting 3D point'
#                #ax.plot(item[0],item[1],item[2])
#                ax.plot(item,item,item, marker = "o",  linestyle = "none", color='r')
#            elif len(vertices[0])==2:
#                print 'plotting 2D point'
#                #ax.plot(item[0],item[1])
#                ax.plot(item,item,item, marker = "o",  linestyle = "none", color='r')
        else:
            print 'Type unknown'
    plt.show()
    return

def sqr_diff(x1,x2):
    sqdf = (x1-x2)**2
    return sqdf
def diff_sqr(x1,x2):
    return x1**2 - x2**2
def diff(x1,x2):
    return x1-x2

def ObjFunction(curve, s, p):
    """
        Accepts:
            curve   : Bsplineclass curve
            s       : parameter location on the curve
            p       : vector location in R2 or R3
        Returns:
            D an AD object of the min distance from P to the curve
    """
    D = 0.
    X = curve.CurvePoint(s)
    for a,b in zip(X,p):
        D += (a-b)**2
    #D = (X[0]-p[0])**2 + (X[1]-p[1])**2 + (X[2]-p[2])**2
    return D #np.sqrt(D)


def dobj(X,p):
    Dk = (X[0][0]-p[0])**2 + (X[0][1]-p[1])**2 + (X[0][2]-p[2])**2
    dDk = (X[0][0]-p[0])*X[1][0] + (X[0][1]-p[1])*X[1][1] + (X[0][2]-p[2])*X[1][2]
#    dDk = (  X[0][0]-p[0])*X[1][0]/X[0][0]   + \
#            (X[0][1]-p[1])*X[1][1]/X[0][1]   + \
#            (X[0][2]-p[2])*X[1][2]/X[0][2]
    dDk *= 2.
    return .5*((Dk)**(-.5))*dDk
    
def ddobj(X,p):
    Dk = (X[0][0]-p[0])**2 + (X[0][1]-p[1])**2 + (X[0][2]-p[2])**2
    dDk = (X[0][0]-p[0])*X[1][0] + (X[0][1]-p[1])*X[1][1] + (X[0][2]-p[2])*X[1][2]
#    dDk = (  X[0][0]-p[0])*X[1][0]/X[0][0]   + \
#            (X[0][1]-p[1])*X[1][1]/X[0][1]   + \
#            (X[0][2]-p[2])*X[1][2]/X[0][2]
    dDk *= 2.
    
    ddDk =   X[1][0]*X[1][0]+(X[0][0]-p[0])*X[2][0] \
           + X[1][1]*X[1][1]+(X[0][1]-p[1])*X[2][1] \
           + X[1][2]*X[1][2]+(X[0][2]-p[2])*X[2][2]
#    ddDk =   X[1][0]*X[1][0]/X[0][0]/X[0][0] + (X[0][0]-p[0])*X[2][0]/X[0][0] \
#           + X[1][1]*X[1][1]/X[0][1]/X[0][1] + (X[0][1]-p[1])*X[2][1]/X[0][1] \
#           + X[1][2]*X[1][2]/X[0][2]/X[0][2] + (X[0][2]-p[2])*X[2][2]/X[0][2]
    ddDk *= .2
         
    s1 = 0.5*((Dk)**-.5)*ddDk 
    s2 = -.25*((Dk)**(-3./2.))*(dDk**2)
    return s1+s2

def naive_minimizer(curve, point):
    """
        return the poin s 
        on the pre-parameterized curve
        that is min distancefrom point
    """
    check = []
    for i in range(curve.nump):
        check.append(ObjFunction(curve,curve.s[i], point))
    for i in range(len(check)):
        if check[i]==min(check):
            store = i
    min_s = curve.s[store]
    return min_s
  
def quadratic_minimizer(curve, p, slist = [[0.,.15,.3],[.3,.5,.6],[.6,.75,1.]]):    
#def quadratic_minimizer(curve, p, slist = [[0.,0.5,1.]]):#[[0.,.15,.3],[.3,.5,.6],[.6,.75,1.],]):
    """
        Accepts:
            curve   : Bspline class curve
            p       : vector location in R2 or R3
        Returns:
            s       : parameter location on the curve
    """
    truemin = []
    stored = []
    for sv in slist:
        for i in range(5):
            y23 = diff_sqr(sv[1],sv[2])
            y31 = diff_sqr(sv[2],sv[0])
            y12 = diff_sqr(sv[0],sv[1])
            s23 = diff(sv[1],sv[2])
            s31 = diff(sv[2],sv[0])
            s12 = diff(sv[0],sv[1])
            
            d = []
            d.append(ObjFunction(curve,sv[0],p))
            d.append(ObjFunction(curve,sv[1],p))
            d.append(ObjFunction(curve,sv[2],p))
            
            sp = 0.5*(y23*d[0]+y31*d[1]+y12*d[2])/(s23*d[0]+s31*d[1]+s12*d[2])
            #smap = {0:d[0],1:d[1],2:d[2]}
            
            #evaluate the new point:
            if curve.trange[0] <= sp and sp<= curve.trange[-1]:
                
                dsp = ObjFunction(curve,sp,p)
        
                
                for j in range(len(d)):
                    if d[j] == min(d):
                        storemin = j
                    if d[j] == max(d):
                        storemax = j
                
                if d[storemax] > dsp:
                    d[storemax] = dsp
                    sv[storemax] = sp
        
        for i in range(len(d)):
            if d[i] == min(d):
                store = i
        
        truemin.append(sv[store])
        stored.append(min(d))
    #print stored
    #print truemin
    for i in range(len(stored)):
        if stored[i] == min(stored):
            store = i
    
    return truemin[store]

def nm_minimizer(curve, point, si):
    
    #print si, ObjFunction(curve, si, point)
    Allderivs = np.zeros((curve.k, curve.dim))
    X=Allderivs
    for i in range(10):
        curve.CurveDerivsAlg1(si,X)
        #dD = dobj( X, point)
        #ddD = ddobj( X, point)
        store = si
        #si = si - dD/ddD
        e = X[0,:]-point[:]
        if curve.dim==3:
            dD = (e[0]*X[1][0]+e[1]*X[1][1]+e[2]*X[1][2])
            ddD = (X[1][0]**2 + X[1][1]**2 + X[1][2]**2 + \
                        e[0]*X[2][0]+e[1]*X[2][1]+e[2]*X[2][2])
        elif curve.dim==2:
            dD = (e[0]*X[1][0]+e[1]*X[1][1])
            ddD = (X[1][0]**2 + X[1][1]**2 +\
                        e[0]*X[2][0]+e[1]*X[2][1])
                    
        si = si -  dD/ddD
                    

        #print 'dD, ddD, si {}, {}, {}'.format(dD, ddD, si)
        #print dD
        if (si<0.):
            print 'Outside Parameter Range! si = {}'.format(si)
            si = 0.
            break
        elif (si>1.):
            print 'Outside Parameter Range! si = {}'.format(si)
            si = 1.0
            break
        
        #print si, ObjFunction(curve, si, point)
    
    return si

def nearest_param(curve, point, start=None):
    """
        curve : Bspline
        point : xyz
        start : parameter loc, s
        smin : parameter loc of min distance
    """

    if start == None:
        print 'using quadratic starting value'
        start1 = quadratic_minimizer(curve,point)
        smin = nm_minimizer(curve, point, start1)
        start = start1
    else:
        print 'using outside starting value'
        smin = nm_minimizer(curve, point, start)

    return smin, start
    
##---------------------------------------------------------------------
##---------------------------------------------------------------------
##---------------------------------------------------------------------
def main():
    #point = np.array(([1.,2.,4.]),float)
    for z in range(6):
        #point[2] = float(z)/2.
        

        vertices = np.array([[0.,0.,0.],
                               [3.,4.,1.],
                               [3.1,4.,2.],
                               [4.,2.,3.],
                               [3.,3.,4.],
                               [4.,4.,5.]])
        point = vertices[z]+.5

        k=4
        curve1  = Bspline(vertices, k)


        check = []
        time1=time.clock()
        for i in range(curve1.nump):
            check.append(ObjFunction(curve1,curve1.s[i], point))
        time2=time.clock()
        print 'Approx way took ',str(time2-time1),' sec.'
        for i in range(len(check)):
            if check[i]==min(check):
                store = i

        print '\nUsing a Naive Search, the objective function is:'
        print curve1.s[store], ObjFunction(curve1, curve1.s[store], point)


        start = curve1.s[store]
        time1=time.clock()
        smin, start = nearest_param(curve1, point, start)
        #smin, start = nearest_param(curve1, point)
        time2=time.clock()

        print '\nAt initialization, the objective function is:'
        print start, ObjFunction(curve1, start, point)

        print '\nAfter newton iteration, the objective function is:'
        print smin, ObjFunction(curve1, smin, point)
        print '\nOpt way took ',str(time2-time1),' sec.'



        a=curve1.CurvePoint(smin)
        vector = [[a[0],point[0]],[a[1],point[1]],[a[2],point[2]]]
        plot_all([curve1,vector])
        curve1.FormParam = {}
        #curve1.xpts, curve1.ypts, curve1.zpts = constrainedAD3D(curve1, curve1.FormParam)
    return
def main2():
    point = np.array([ 22.,  12.,   2.],float)
    

    vertices = np.array([[0.,0.,0.],
                        [10.,5.,1.],
                        [15.,7.5,1.5],
                        [20.,10.,2.],
                        [30.,15.,3.],
                        [40.,20.,4.]])
    k=4
    nump=5
    curve = Bspline(vertices, k, nump)
    check = []
    time1=time.clock()
    for i in range(curve.nump):
        check.append(ObjFunction(curve,curve.s[i], point))
    time2=time.clock()
    print 'Approx way took ',str(time2-time1),' sec.'
    for i in range(len(check)):
        if check[i]==min(check):
            store = i

    print '\nUsing a Naive Search, the objective function is:'
    print curve.s[store], ObjFunction(curve, curve.s[store], point)


    start = curve.s[store]
    time1=time.clock()
    smin, start = nearest_param(curve, point, start)
    #smin, start = nearest_param(curve, point)
    time2=time.clock()

    print '\nAt initialization, the objective function is:'
    print start, ObjFunction(curve, start, point)

    print '\nAfter newton iteration, the objective function is:'
    print smin, ObjFunction(curve, smin, point)
    print '\nOpt way took ',str(time2-time1),' sec.'

    a=curve.CurvePoint(smin)
    vector = [[a[0],point[0]],[a[1],point[1]],[a[2],point[2]]]
    plot_all([curve,vector])
    curve.FormParam = {}
    
import numpy as np
import matplotlib.pyplot as plt
import time
from curve import *#Bspline, interpolatedBspline
#from AD import constrainedAD3D
#from EqualityFunctions import ADevalx, ADevaly, ADevalz

if __name__ == '__main__':
    option = 'original test'
    #option = 'test2'
    if option == 'original test':
        main()
    if option =='test2':
        main2()
