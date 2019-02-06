# Luke McCulloch
# B-Spline Knot Insertion Algorithm
# August 2012

import numpy as np
import copy

def distance(point1,point2):
    """
        generic distance formula
        takes:
            any dimension as long asthey match
            tuples
            lists
            arrays
    """
    point1 = np.asarray(point1)
    point2 = np.asarray(point2)
    dist = 0.
    for a,b in zip(point1,point2):
        dist += ((a-b)**2)
    return np.sqrt(dist)

def distance2D(point1,point2):

    x = point1[0]-point2[0]
    y = point1[1]-point2[1]
    dist = np.sqrt((x*x) + (y*y))

    return dist

def distance3D(point1,point2):

    x = point1[0]-point2[0]
    y = point1[1]-point2[1]
    z = point1[2]-point2[2]
    dist = np.sqrt((x*x) + (y*y) + (z*z))

    return dist

def Distance4D(point1,point2):
    """Should be changed to some P&T algorithm to compute a NURBS distance"""
    x = point1[0]-point2[0]
    y = point1[1]-point2[1]
    z = point1[2]-point2[2]
    dist = np.sqrt((x*x) + (y*y) + (z*z))

    return dist


def intKnot(t):
    """
        Documentation needed!
    """
    ti = list(t)
    ti.insert(0,t[0]) #inserts t[0] into position 0 in ti
    ti.append(t[-1])
    ti = np.asarray(ti)
    return ti

def exteriorKnotsSet(t):
    knot_list           = list(t)          
    knot_set            = set(t)
    outer_knot_set      = set([min(knot_set),max(knot_set)])
    return outer_knot_set

def interiorKnotsSet(t):
    knot_list   = list(t)          
    knot_set    = set(t)
    outer_knot_set = set([min(knot_set),max(knot_set)])
    interior_knotset = knot_set - outer_knot_set
    return interior_knotset

def multiplicity(t):
    """
        Finds the multiplicity of each knot in the knot vector
            uses array to list
            set of array
            {} dictionary - knot multiplicity pairs
    """
    ti = list(t)
    knotset=set(t)  #each entry in a set is unique
    mult = {}       #Dictionary: For each unique knot, store the mutiplicity
    for knot in knotset:
        mult[knot]=ti.count(knot)

    return mult


def curve_pts_dict(N_span,k):
    """
        Use bezier decomposition code to make a hash table (dict)
        containing all pts in each bezier curve
        also the knots of each bezier curve
        {curve:pts}
        {curve:knots} -> not needed! we find the product knot set
        and the old knot set after forming B-Splines into Bezier form

        function would:

            return the dict "bezier_CurveToPoint"

        require that you pass:
            N_span
            k
    """
    bezier_CurveToPoint={} #our dictionary

    for i in range(N_span):
        span_num=i

        """I should not be creating a new list each time, but I have to clear it"""
        ilist=[]

        for i in range(((span_num+1)*k)-span_num-k,((span_num+1)*k)-span_num,1):
            """Iterate through the points by index, listing the points (via pt index) that are in each curve"""
            ilist.append(i)

        """When the list has every value associated with that span, put them in the dict"""
        bezier_CurveToPoint[span_num] = ilist

    return bezier_CurveToPoint



def new_points(k,r,ta,t,points):
    """
        k       :   curve order
        r       :   is the index location of
                    the new knot
        ta      :   the new knot
        t       :   knot vector
        points  :   control points
    """
    n = len(points)
    points_list = list()

    # Manipulate the Vertices
    for j in range(n+1):
        if j<=r-k+1:
            lambdaj = 1
        elif (r-k+2<=j<=r):
            lambdaj = (ta-t[j])/(t[j+k-1]-t[j])
        elif j>=r+1:
            lambdaj = 0
        if j==0:
            points_list.append(lambdaj*points[j])
        elif j<n:
            points_list.append((1-lambdaj)*points[j-1]+lambdaj*points[j])
        else:
            points_list.append((1-lambdaj)*points[j-1])
    points = np.asarray(points_list)
    return points



def knot_insertion(k,points,t,ta):
    """
        k       :   curve order
        points  :   control vertices
        t       :   knot vector
        ta      :   knot to be inserted
                
        modified on Dec 5 2015. TLM
        
    """
    

    knot_list = list(t)
    
    found = False
    # Find where the new knot will be inserted
    for i in range(len(t)-1):
        if t[i]<ta<=t[i+1]:
            found = True
            r=i+1 
            #print r
            knot_list.insert(r,ta) #inserts ta before r.
            break
    if found == False:
        for i in reversed(range(len(t)-1)):
            if t[i]<=ta<t[i+1]:
                found = True
                r=i+1
                knot_list.insert(r,ta) 
                print 'ir',r
                break
    assert(found==True),'knot insertion failed'
    # L Birk knot vector methods wants index before insertion point so pass r as r-1
    points = new_points(k,r-1,ta,t,points)

    # Overwrite the old knot vector
    t = np.asarray(knot_list)

    
    return t, points
    

def knot_insertion_hierarchical(k,points,t,ta):
    """
        k       :   curve order
        points  :   control vertices
        t       :   knot vector
        ta      :   knot to be inserted
                
        modified on Dec 5 2015. TLM
    """
    knot_list = list(t)
    found = False
    # Find where the new knot will be inserted
    for i in range(len(t)-1):
        if t[i]<ta<=t[i+1]:
            found = True
            r=i+1 
            #print r
            knot_list.insert(r,ta) #inserts ta before r.
            break
    if found == False:
        for i in reversed(range(len(t)-1)):
            if t[i]<=ta<t[i+1]:
                found = True
                r=i+1
                knot_list.insert(r,ta) 
                print 'ineresting:  we are in the back-loop of knot insertion ir = ',r
                break
    assert(found==True),'hierarchical knot insertion failed'
    # L Birk knot vector methods wants index before insertion point so pass r as r-1
    points = new_points(k,r-1,ta,t,points)

    # Overwrite the old knot vector
    t = np.asarray(knot_list)

    
    return t, points, r





        

def knot_removal(curve,index,
                 TOL=.000001):
                 #p,points,
                 #t,ur,
                 #ur_index,
                 #TOL=.000001):
    """
        Lets try to remove a knot and see what the curve looks like.
        1st calculate the points, then remove the knot!

        Inputs:
            p           =  curve degree
            points      =  vertices
            ur          =  knot to be removed
            t           =  knot vector
            store_knot  = ???
            
        Other:
            r       =  last index loc of mult. knot to be removed     
            s       =  multiplicity of ur
            ith     =  the number of times a knot is removed

        mult:   multiplicity of all knots
                    usage: mult[knot] returns the multiplicity of that knot
                    
                    
        
        p = b2.p
        k = b2.k
        
        points = b2.vertices
        
        t = b2.t
        ta = .5
        ur = .5
        ur_index = 4
        
        TOL=.000001
        
    """
    
    p   = curve.p
    points  = curve.vertices
    t = curve.t
    ur = curve.t[index]
    ur_index = index
    
    U   = curve.t
    Pw  = curve.vertices
    
    #temp        = np.zeros((len(points),2))
    temp        =np.zeros((len(points),len(points[0])),float)

    knot_list   = list(t)           #"knot_list.remove(knot)" will remove the knot once.
    knot_set    = set(t)
    outer_knot_set = set([min(knot_set),max(knot_set)])
    interior_knotset = knot_set - outer_knot_set

    if ur in interior_knotset:
        print "interior knot"
    
    
        point_list  = list(points)
        
        mult        = multiplicity(t)   
        s           = mult[ur]          #multiplicity of knot to be removed
        print "mult[ur] = {}".format(s)

        m = len(t)
        #------------------------------------------
        if p==1:
            print "cannot remove a point,"
            print "curve is the control polygon"
            store_knot.remove(ur)
            return t, points

        #------------------------------------------
        # Reverse Loop to Find the knot to be removed
        #for i in range(m-1,0,-1):
        #        if ur_index is None:
        #            r = -1
        #            rfound = False
        #            for i in range(m-p-2,p,-1):
        #                if ur==t[i]:
        #                    rfound = True
        #                    r=i #index of knot to be removed
        #                    #if (len(t)-r)  ?????????????????????????
        #                    break
        #        else:
        r=ur_index
        
        #assert(rfound),'Error, could not find index of knot to remove'
        #k = p+1
        
        i = r-p # knot index - order
        j = r-s # knot index - multiplicity
        print "i = {}, j = {}".format(i,j)

        #temp[0] = points[i-1]
        #temp[j+1] = points[j+1]
        temp[i-1] = points[i-1]
        temp[j+1] = points[j+1]

        print "temp[i-1] = {}, temp[j+1] = {}".format(temp[i-1],temp[j+1])
        
        if (j-i>=0):  #if?
            #for u in interior_knotset:
            
            alphai   = (ur-t[i])/(t[i+p+1]-t[i])
            alphaj   = (ur-t[j])/(t[j+p+1]-t[j])
            temp[i] = (points[i]-(1.0-alphai)*temp[i-1])/alphai
            temp[j] = (points[j]-alphaj*temp[j+1])/(1.0-alphaj)

            print "temp[i] = {}, temp[j] = {}".format(temp[i],temp[j])
                
            i  = i+1
            j  = j-1
            
            
        storei = i-1
        storej = j+1
        print "storei = {}, storej = {}".format(storei, storej)
        #flag=True
        #if ((j-i)<=0):
        #Dedent
        if (distance2D( temp[i-1],temp[j+1] ) <= TOL):
            #flag=True
            print "Knot removal was successful, Save new control points"
    ##            i = r-p
    ##            j = r-s
    ##            while (j-i>0):
    ##                points[i]=temp[i]
    ##                points[j]=temp[j-i+1]
    ##                i = i+1
    ##                j = j-1
            #Update knot list
            knot_list.remove(ur)
            t = np.asarray(knot_list)

            #Update points
            #print "Update Point List, Point_List = {}".format(point_list)
            point_list.insert(storei+1,temp[i])
            #print "Update Point List, Point_List = {}".format(point_list)
            del point_list[storei]
            #print "Update Point List, Point_List = {}".format(point_list)
            del point_list[storej]
            #print "Update Point List, Point_List = {}".format(point_list)

            points=np.asarray(point_list)
                               

        else:
            #flag=False
            print "cannot remove the knot"
            print "temp[i-1] = {}, temp[j+1] = {}".format( temp[i-1],temp[j+1] )
            print "tdistance2D = {}".format(distance2D( temp[i-1],temp[j+1] ))
            #store_knot.remove(ur)
        #end dedent
        return t,points

    elif ur in outer_knot_set:
        print "Error: knot is not on the interior."
        #store_knot.remove(ur)
        return t, points
    else:
        print "knot not found"
        #store_knot.remove(ur)
        return
        

def unsafe_knot_removal(curve,ur,index=None):
    """
        Lets try to remove a knot and see what the curve looks like.
        1st calculate the points, then remove the knot!

        Inputs:
            p           =  curve degree
            points      =  vertices
            ur          =  knot to be removed
            t           =  knot vector
            store_knot  = 
            
        Other:
            r       =  last index loc of mult. knot to be removed     
            s       =  multiplicity of ur
            ith     =  the number of times a knot is removed

        mult:   multiplicity of all knots
                    usage: mult[knot] returns the multiplicity of that knot
        
    """
    p           = curve.p
    points       = curve.vertices
    t           = curve.t
    store_knot  = copy.copy(curve.t)
    #
    temp        = np.zeros((len(points),2))
    #temp        =np.zeros((len(poins),len(points[0])),float)

    knot_list   = list(t)           #"knot_list.remove(knot)" will remove the knot once.
    knot_set    = set(t)
    outer_knot_set = set([min(knot_set),max(knot_set)])
    interior_knotset = knot_set - outer_knot_set

    if ur in interior_knotset:
        print "interior knot"
    
    
        point_list  = list(points)
        
        mult        = multiplicity(t)   
        s           = mult[ur]          #multiplicity of knot to be removed
        print "mult[ur] = {}".format(s)

        m = len(t)
        #------------------------------------------
        if p==1:
            print "remove a point from the cuve,"
            print "curve is the control polygon"
            store_knot.remove(ur)
            return t, points

        #------------------------------------------
        # Reverse Loop to Find the knot to be removed
        #for i in range(m-1,0,-1):
        for i in range(m-p-2,p,-1):
            if ur==t[i]:
                r=i #index of knot to be removed
                #if (len(t)-r)  ?????????????????????????
                break
        
        if index is not None:
            if index is i:
                print 'would have worked anyway'
            else:
                print 'found index = {}'.format(i)
                print 'actual index = {}'.format(index)
            r=index
            
        k = p+1
        
        i = r-p # knot index - order
        j = r-s # knot index - multiplicity
        print "i = {}, j = {}".format(i,j)

        #temp[0] = points[i-1]
        #temp[j+1] = points[j+1]
        temp[i-1] = points[i-1]
        temp[j+1] = points[j+1]

        print "temp[i-1] = {}, temp[j+1] = {}".format(temp[i-1],temp[j+1])
        
        #if (j-i>=0):  #if?
        for u in interior_knotset:
        
            alphai   = (ur-t[i])/(t[i+p+1]-t[i])
            alphaj   = (ur-t[j])/(t[j+p+1]-t[j])
            temp[i] = (points[i]-(1.0-alphai)*temp[i-1])/alphai
            temp[j] = (points[j]-alphaj*temp[j+1])/(1.0-alphaj)

        print "temp[i] = {}, temp[j] = {}".format(temp[i],temp[j])
            
        i  = i+1
        j  = j-1
            
            
        TOL=.000001
        storei = i-1
        storej = j+1
        print "storei = {}, storej = {}".format(storei, storej)
        #flag=True
        #if ((j-i)<=0):
        #Dedent
        #if (distance( temp[i-1],temp[j+1] ) <= TOL):
        #flag=True
        print "Knot removal was successful, Save new control points"
##            i = r-p
##            j = r-s
##            while (j-i>0):
##                points[i]=temp[i]
##                points[j]=temp[j-i+1]
##                i = i+1
##                j = j-1
        #Update knot list
        knot_list.remove(ur)
        t = np.asarray(knot_list)

        #Update points
        #print "Update Point List, Point_List = {}".format(point_list)
        point_list.insert(storei+1,temp[i])
        #print "Update Point List, Point_List = {}".format(point_list)
        del point_list[storei]
        #print "Update Point List, Point_List = {}".format(point_list)
        del point_list[storej]
        #print "Update Point List, Point_List = {}".format(point_list)

        points=np.asarray(point_list)
                               

        #else:
        #    #flag=False
        #    print "cannot remove the knot"
        #    print "temp[i-1] = {}, temp[j+1] = {}".format( temp[i-1],temp[j+1] )
        #    print "tdistance2D = {}".format(distance( temp[i-1],temp[j+1] ))
        #    store_knot.remove(ur)
        #end dedent
        return t,points

    elif ur in outer_knot_set:
        print "Error: knot is not on the interior."
        store_knot.remove(ur)
        return t, points
    else:
        print "knot not found"
        #store_knot.remove(ur)
        return

def RefineKnotVector(n,p,U,Pw,X,r,Ubar,Qw):
    """P&T algorithm page 164 - Jamie has a version
        n = # vertices - 1
        p = order of the curve
        U = old knot vector
        Pw = control points
        X = knots to be inserted into U
        r = dimensions of the curve space
        Ubar = new knot vector
        Qw = new set of control points"""
    return



def knot_removal3D(curve,index):
                  # p,points,t,ur,store_knot,dim):
    """
        Lets try to remove a knot and see what the curve looks like.
        1st calculate the points, then remove the knot!

        Inputs:
            p           =  curve degree
            points      =  vertices
            ur          =  knot to be removed (float)
            t           =  knot vector
            store_knot  = list of knots to be removed (to restore if they fail to remove)
            dim         = dimensions of the curve space (2 for 2D, 3 for 3D, etc)
            
        Other:
            r       =  last index loc of mult. knot to be removed     
            s       =  multiplicity of ur
            ith     =  the number of times a knot is removed

        mult:   multiplicity of all knots
                    usage: mult[knot] returns the multiplicity of that knot
        
    """
    p = curve.p
    points = curve.vertices
    t = curve.t
    ur = curve.t[index]
    store_knot = curve.t[index]
    dim = curve.dim
    
    temp        = np.zeros((len(points),dim))

    knot_list   = list(t)           #"knot_list.remove(knot)" will remove the knot once.
    knot_set    = set(t)
    outer_knot_set = set([min(knot_set),max(knot_set)])
    interior_knotset = knot_set - outer_knot_set

    if ur in interior_knotset:
        #print "interior knot"
    
    
        point_list  = list(points)
        
        mult        = multiplicity(t)   
        s           = mult[ur]          #multiplicity of knot to be removed
        #print "mult[ur] = {}".format(s)

        m = len(t)
        #------------------------------------------
        if p==1:
            print "cannot remove a point,"
            print "curve is the control polygon"
            store_knot.remove(ur)
            return t, points

        #------------------------------------------
        # Reverse Loop to Find the knot to be removed
        #for i in range(m-1,0,-1):
        #        for i in range(m-p-2,p,-1):
        #            if ur==t[i]:
        #                r=i #index of knot to be removed
        #                #if (len(t)-r)  ?????????????????????????
        #                break
            
        r = index
        
        k = p+1
        
        i = r-p # knot index - order
        j = r-s # knot index - multiplicity
        #print "i = {}, j = {}".format(i,j)

        #temp[0] = points[i-1]
        #temp[j+1] = points[j+1]
        #print temp[i-1] , points[i-1]
        temp[i-1] = points[i-1]
        temp[j+1] = points[j+1]

        #print "temp[i-1] = {}, temp[j+1] = {}".format(temp[i-1],temp[j+1])
        
        if (j-i>=0):  #if?
            #for u in interior_knotset:
            
            alphai   = (ur-t[i])/(t[i+p+1]-t[i])
            alphaj   = (ur-t[j])/(t[j+p+1]-t[j])
            temp[i] = (points[i]-(1.0-alphai)*temp[i-1])/alphai
            temp[j] = (points[j]-alphaj*temp[j+1])/(1.0-alphaj)

            #print "temp[i] = {}, temp[j] = {}".format(temp[i],temp[j])
                
            i  = i+1
            j  = j-1
            
            
        TOL=.000001
        storei = i-1
        storej = j+1
        #print "storei = {}, storej = {}".format(storei, storej)
        #flag=True
        #if ((j-i)<=0):
        #Dedent
        if (distance3D( temp[i-1],temp[j+1] ) <= TOL):
            #flag=True
            print "Knot removal was successful, Save new control points"
    ##            i = r-p
    ##            j = r-s
    ##            while (j-i>0):
    ##                points[i]=temp[i]
    ##                points[j]=temp[j-i+1]
    ##                i = i+1
    ##                j = j-1
            #Update knot list
            knot_list.remove(ur)
            t = np.asarray(knot_list)

            #Update points
            #print "Update Point List, Point_List = {}".format(point_list)
            point_list.insert(storei+1,temp[i])
            #print "Update Point List, Point_List = {}".format(point_list)
            del point_list[storei]
            #print "Update Point List, Point_List = {}".format(point_list)
            del point_list[storej]
            #print "Update Point List, Point_List = {}".format(point_list)

            points=np.asarray(point_list)
                               

        else:
            #flag=False
            print "cannot remove the knot {}".format(ur)
            print "temp[i-1] = {}, temp[j+1] = {}".format( temp[i-1],temp[j+1] )
            print "tdistance3D = {}".format(distance3D( temp[i-1],temp[j+1] ))
            print 'tolerance = {}'.format(TOL)
            store_knot.remove(ur)
        #end dedent
        return t,points

    elif ur in outer_knot_set:
        print "Error: knot is not on the interior."
        store_knot.remove(ur)
        return t, points
    else:
        print "knot not found"
        store_knot.remove(ur)
        return
    




def mult_product_knot_set(k1,k2,m1,m2):
    """
        k1  =   order of curve 1
        k2  =   order of curve 2
        m1  =   multiplicity of knots in t1
        m2  =   multiplicity of knots in t2
        m   =   product knot set
    """
    """ Get the multiplicity
        of the product knot set
    """
    if m1>0 and m2>0:
        return max(k1-1+m2,k2-1+m1)
    elif m1==0 and m2>0:
        return
    elif m1>0 and m2==0:
        return
    elif m1==0 and m2==0:
        return 0
    else:
        print 'error: I did not return the product knot set'
        return




#def Distance4D():
#    """P&T algorithm to compute a NURBS distance"""
#    return
#def RemoveCurveKnot(n,p,U,Pw,u,r,s,num,t):
def RemoveCurveKnot(curve,u,r=None,
                    s=1,num=1,
                    TOL = 1.0000001):
    """P&T algorithm page 185 (page 99 in PDF)
        retuires Distance4D page?
            n   = curve.n-1
            p   = curve.p
            U   = curve.t
            Pw  = curve.vertices
            
            u   = knot to be removed
            r   = index(of u in U)
            
            s   = multiplicity of u
            num = number of times to remove u.
            t   = the number of times removable u is.
        """
    if u is None: u = curve.t[r]
    n   = curve.n-1
    p   = curve.p
    U   = curve.t
    Pw  = curve.vertices
    
    pws = list(np.shape(Pw))
    pws[0] = pws[0]-1
    #Pw_final = np.zeros(pws,float)
    Pw_new = np.copy(Pw)


    temp = np.zeros((2*p+1,curve.dim),float)

    m    = n+p+1
    _ord = p+1          #changed code to avoid python intrinsic function
    fout = (2*r-s-p)/2  #1st control point out
    last  = r-s #error, r
    first = r-p
    
    t = range(num)[0]
    for t in range(num): # Eq. 5.28
        off = first-1 #Difference in index between temp and P
        temp[0] = Pw[off]
        temp[last+1-off]=Pw[last+1]
        i=first
        j=last
        ii=1
        jj=last-off
        remflag=0
        while (j-i > t): #compute new control points for one removal step
            alfi = (u-U[i])/(U[i+_ord+t]-U[i])
            alfj = (u-U[j-t])/(U[j+_ord]-U[j-t])
            temp[ii] = (Pw[i]-(1.0-alfi)*temp[ii-1])/alfi
            temp[jj] = (Pw[j]-alfi*temp[jj+1])/(1.0-alfj)
            i=i+1
            ii=ii+1
            j=j-1
            jj=jj-1
        if (j-i<t): #check if knot is removable
            #check = Distance4D(temp[ii-1],temp[jj+1])
            check = distance(temp[ii-1],temp[jj+1])
            #print '1 Distance = {}'.format(check)
            if ( check <= TOL):
                remflag = 1
        else:
            alfi = (u-U[i])/(U[i+_ord+t]-U[i])
            testPoint=alfi*temp[ii+t+1]+(1.0-alfi)*temp[ii-1]
            #check = Distance4D(Pw[i],testPoint)
            check = distance(Pw[i],testPoint)
            #print '2 Distance = {}'.format(check)
            if (check <=TOL):
                remflag = 1
            else:
                #remflag = 1 #hack on this a bit
                #print Pw[i], testPoint
                pass
        if (remflag ==0): #cannot remove any more knots
            #print 'failed'
            pass
            #remflag = 1 #hack on this a bit
            #break #hack on this a bit too
        else: #successful removal. Save new control pts.
            #print 'success, save new'
            i=first
            j=last
            #Pw_new[:i] = Pw[:i]
            #Pw_new[j-1:] = Pw[j:]
            while(j-i>t):
                print i,j
                Pw_new[i] = temp[i-off]
                Pw_new[j] = temp[j-off]
                i=i+1
                j=j-1
        
        first = first-1
        last = last+1 #End of the for loop
    if (t==0):
        return Pw_new
    for k in range(r+1,m):
        U[k-t]=U[k]
    j=fout
    i=j # Pj through Pi will be overwritten
    for k in range(1,t+1):
        if ((k%2)==1):
            i = i+1
        else:
            j=j-1
    for k in range(i+1,n): #shift
        Pw_new[j] = Pw[k]
        j=j+1
    return Pw_new #,temp



def knot_deformation(kv1,kv2,steps):
    """
    Parameters:
    --------------------
        kv1: initial knot vector
        kv2: final knot vector
        steps: number of steps to complete the transformation
        
    Returns:
    --------------------
        A closure which maps from kv1 to kv2 in steps # of steps
        
    Assumes:
    --------------------
        Knot vectors are the same size (easily achievable via knot insertion)
        just insert a knot in the smaller vector
        (assumes you want to parameterize the larger vector and match
        the shape of the curve using the smaller vector)
        
    TODO:
    --------------------
        -not safe
        -use yield?
    """
    #    difference = []
    #    for k1,k2 in zip(kv1,kv2):
    #        difference.append(k1-k2)
    #        
    #    listoftransforms = []
    #    for d in difference:
    #        listoftransforms.append(np.linspace())
        
    listoftransforms = []
    for k1,k2 in zip(kv1,kv2):
        listoftransforms.append(np.linspace(k1,k2,steps,endpoint=True))
        
    vt = np.asarray(listoftransforms)
    
    def knotstep(n,vector=vt):
        return vector[:,n]
    
    return knotstep


if __name__ =='__main__':
    from curve import Bspline, linear_vertices
    k       = 4
    nump    = 50
    start   = (0.,0.)
    stop    = (10.,10.)
    num     = 7
    vertices = linear_vertices(start, stop, num)
    curve = Bspline(vertices, k, nump)
    

    ta = 0.4
    tnew, newpoints = knot_insertion(curve.k,curve.vertices,curve.t,ta)
    c1 = Bspline(newpoints, k, nump)
    c1.t = tnew
    
    tnew, newpoints = knot_insertion(c1.k,c1.vertices,c1.t,ta)
    c2 = Bspline(newpoints, k, nump)
    c2.t = tnew