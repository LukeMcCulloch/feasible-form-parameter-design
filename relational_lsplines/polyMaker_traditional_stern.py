## B-Spline Hull Surface to VTK polygon surface
## A Wigley code translation: FORTRAN to Python
## Modified for different hull definition
## TLM August 28, 2013


import numpy as np
from routines import crossproduct
from panelCentroids import panelCentroids

def polymake(surf1, FSheight, panelShift, nhl, nhg, nhf, Compute_Curvatures,L,B,T):


    #--------------------------------------------------------------------------------------------------------------
    #FSheight = FSH

    filename = 'TLM'

    polyfile = open("BtoPloyBoat.txt", "w")
    realBoat = "TLM_Boat.dat"

    comment = "A B-spline Surface Discretization"+'\n'

    #Points networks:
    tnet    = [] #Transverse Curve Net: (would be output from Bspline optimization program)
    lnet    = [] #Longitudinal Curve Net
    surfnet = [] #surface point net

    zero = 0.0
    small = 0.000001
    #L = 75.
    #B = 5.
    #T = 10. # a bit misleading - but I want + areas...  This is the max rise of the hull from 0.
    C_L = zero

    longiOffset = L/2.

    # number of hull panels and points (small numbers to test program)
    #nhl = 30  # lengthwise 
    #nhg = 8   # girthwise below WL  ex: for Fn 0.3 wavelength hull of length 1, I need # panels
    #nhf = 1   # girthwise above WL

    freeSurfaceWidthFactor = .02
    



    npanels = nhl*(nhg+nhf)
    npoints = (nhl+1)*(nhg+nhf+1)

    npoints_uv = (nhl+1)*(nhg+1) # number of u,v coords to store. A u,v point for every pt on hull
    
    




    Lfs = 4.*L #free surface discretization length


    nfsl = 4*nhl
    nfst = int(0.75*nhl) # may need to be adjusted

    nfspanels = nfsl*nfst
    npoints = npoints + (nfsl+1)*(nfst+1)

    # position of free surface panels above z=0
    # H. Raven suggests: about 0.8 of panel length
    #     zfspanels = 0.8 * L/ nhl
    zfspanels = FSheight#0.65*L / float(nhl)








    ## Allocate memory------------------------------
    panels = np.zeros((4, npanels+nfspanels),int)
    points = np.zeros((3, npoints),float)

    u_v_points = np.zeros((2,npoints),float)    # store the u,v coordinates of each panel point
    curvature = np.zeros((npanels+nfspanels),float)   # easier to write curavture on panel centroids
    mean_curvature = np.zeros((npanels+nfspanels),float) 
    

    ## Compute hull offsets
    ## We go from bow to stern and from keel to WL
    ipt = 0  # point counter

    dz = T / float(nhg)

    for j in range(nhg + nhf +1 ):  #transverse direction
        #u=float(j)/float(nhg + nhf)
        u=float(j)/float(nhg)
        for i in range(nhl+1):      #longitudinal direction
            v=float(i)/float(nhl)
            if j <= nhg:
                
                u_v_points[:,ipt] = u,v  #save the panel points in u,v space
                points[:,ipt]=surf1.evalSurface(u,v) #evaluate the surface and save the panel points in x,y,z space

                
                
                points[2,ipt] -= T #move down by T
            if j > nhg: #we have topped out... in u space
                u=1.
                points[:,ipt]=surf1.evalSurface(u,v)
                u_v_points[:,ipt] = u,v 
                points[2,ipt]= zfspanels/float(nhf) * float(j-nhg) #LB
                #points[2,ipt]= zfspanels    #match free surf?
            
            points[0,ipt] = -points[0,ipt] + longiOffset
            
            if abs(points[1,ipt])<small:
                points[1,ipt]=zero
                
            polyfile.write(str(points[0,ipt])+', '+str(points[1,ipt])+', '+str(points[2,ipt])+'\n')
            ipt +=1



    # number of points on the hull
    nhpt = ipt



    ## Compute hull panels
    ## We go from bow to stern and from keel to WL
    ipan = 0  # panel counter

    for j in range(1,nhg + nhf +1):
        for i in range(1,nhl+1):
            panels[0,ipan] = i + (j-1)*(nhl+1)
            panels[1,ipan] = i + j*(nhl+1)
            panels[2,ipan] = i + 1 + j*(nhl+1)
            panels[3,ipan] = i + 1+ (j-1)*(nhl+1)
            #WRITE(6,'(I4'': '',4(I5,2X))') ipan, panels(:, ipan)
            #polyfile.write(ipan, panels[:, ipan])
            #polyfile.write(str(ipan)+', '+str(panels[0,ipan])+', '+str(panels[1,ipan])+', '+str(panels[2,ipan])++', '+str(panels[3,ipan])'\n')
            polyfile.write(str(panels[0,ipan])+', '+str(panels[1,ipan])+', '+str(panels[2,ipan])+', '+str(panels[3,ipan])+'\n')
            ipan = ipan+1



    ## Compute free surface offsets
    ## We go from bow to stern and from midships to max. width

    dx = Lfs/float(nfsl)
    ##******************************* Improvising an x shift#*************************
    deltax = panelShift #.8*dx
    ##*******************************Improvising an x shift#**************************
    dy = freeSurfaceWidthFactor * L # width of first free surface strip midships
    yfs = zero  #start at midships... ?

    # width of free surface
    Bfs = 0.
    for j in range(nfst):
        Bfs = Bfs + 1.1**j
    Bfs = dy*Bfs

    loc_flag = 0 #initialize location flag - recognize when we are on or off the hull
    ii = 1.0
    for j in range(nfst + 1):  #transverse
        for i in range(nfsl+1):  #longitudinal
            
            points[1,ipt]=0.
            
            #if loc_flag==0: #before the hull
            x = L - float(i)*dx
            if (abs(x) <= longiOffset):# and j==0): 
                loc_flag=1 # we are on the hull
                u=1.
                v=0.5*(1.-x/longiOffset)
                points[:,ipt] = surf1.evalSurface(u,v) #get x,y,z location on the surface
                points[0,ipt] = -points[0,ipt] # reverse the x orientation
                points[0,ipt] += longiOffset #move the x position to reflect shift of origin in computational space
                x=points[0,ipt]
                store=x

            if x < -longiOffset:#we are at the first location behind the hull
                loc_flag=2
                #print loc_flag , x,store
                x = store - ii*dx
                ii=ii+1.

            #print x

            y = (Bfs-points[1,ipt])/Bfs*yfs + points[1,ipt]
            if abs(y)<small:
                y=zero
                #print y
            points[:,ipt]=[x,y,zfspanels]
            
            
            #points[2,ipt] -= T
            polyfile.write(str(points[0,ipt])+', '+str(points[1,ipt])+', '+str(points[2,ipt])+'\n')
            ipt += 1
        ii=1.
        #print '\n\n'
        yfs = yfs + dy*1.1**(j) 







    ## Compute the free surface panels




    for j in range(1,nfst+1):
        for i in range(1,nfsl+1):
            panels[0,ipan] = nhpt + i + (j-1)*(nfsl+1)
            panels[1,ipan] = nhpt + i + j*(nfsl+1)
            panels[2,ipan] = nhpt + i + 1 + j*(nfsl+1)
            panels[3,ipan] = nhpt + i + 1 + (j-1)*(nfsl+1)
            #WRITE(6,'(I4'': '',4(I5,2X))') ipan, panels(:, ipan)
            #polyfile.write(ipan, panels[:, ipan])
            #polyfile.write(str(ipan)+', '+str(panels[0,ipan])+', '+str(panels[1,ipan])+', '+str(panels[2,ipan])++', '+str(panels[3,ipan])'\n')
            polyfile.write(str(panels[0,ipan])+', '+str(panels[1,ipan])+', '+str(panels[2,ipan])+', '+str(panels[3,ipan])+'\n')
            ipan = ipan+1

    polyfile.close()



    #non - dimensionalize to boat length:
    #points=points/L

    
    if Compute_Curvatures==True:
        # Find the centroids of the hull panels in u,v space.
        corners = np.zeros((2,4),float) # temp storage for panel corners
        ip1 = np.asarray([1,2,3,0])
        
        # loop over all panels
        for i in range(npanels):
            for j in range(0,4):
                corners[:,j]=u_v_points[:,panels[j,i]]  #set cornersLocal
            u_c, v_c = panelCentroids(corners)

            # evaluate curvature at the panel points:
            curvature[i] = surf1.eval_K(u_c,v_c)
            mean_curvature[i] =  surf1.eval_meanCurvature(u_c,v_c)
        
    ##*******************************************************************************

    realFile = open(realBoat, "w")

    realFile.write(comment)
    realFile.write('    '+'0'+'    '+'1'+'\n')
    realFile.write('    '+str(npanels)+'    '+str(nfspanels)+'    '+str(npoints)+'\n')
    realFile.write('    '+str(deltax)+'\n')
    for i in range(npanels+nfspanels):
        realFile.write(str(panels[0,i])+'   '+str(panels[1,i])+'    '+str(panels[2,i])+'    '+str(panels[3,i])+'\n')
    for i in range(npoints):
        realFile.write(str(points[0,i])+'   '+str(points[1,i])+'    '+str(points[2,i])+'\n')
    #if Compute_Curvatures==True:
    for i in range(len(curvature)):
        realFile.write(str(curvature[i])+'\n')
    for i in range(len(mean_curvature)):
        realFile.write(str(mean_curvature[i])+'\n')
    realFile.close()
    #-----------------------------------------------------------------------------------------------------



    return points, panels




def poly_prep(L,B,T,surf1,Compute_Curvatures=False,
                 nhl = 32, nhg=16,nhf=1, 
                 hFraction=.8, xFraction=.65):
                     
    Lfs = 4.*L              #free surface discretization length
    nfsl = 4*nhl            # num longitudinal panels
    dx = Lfs/float(nfsl)
    deltax = xFraction*dx
    panelShift  = deltax
    

    FSH         = hFraction*L / float(nhl)
    
    #print '\n\n'
    
    
    #x = Lfs/float(nfsl)

    xFraction = 0.8 #float(raw_input('New Setting --> '))
    deltax = xFraction*dx
    panelShift  = deltax
    
    points, panels = polymake(surf1, FSH, panelShift, nhl, nhg, nhf, Compute_Curvatures, L,B,T)
