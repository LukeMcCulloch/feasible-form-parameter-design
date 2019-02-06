## B-Spline Hull Surface to VTK polygon surface
## A Wigley code translation: FORTRAN to Python
## Modified for different hull definition
## TLM August 28, 2013


import numpy as np
from routines import crossproduct
from panelCentroids import panelCentroids

def polymake(surf1, FSheight, panelShift, nhl, nhg, nhf, Compute_Curvatures,L,B,T):


    hFraction=FSheight*float(nhl)/L
    #FSH         = hFraction*L / float(nhl)

    Lfs = 4.*L
    nfsl = 4*nhl
    dx = Lfs/float(nfsl)
    xFraction=dx/panelShift
    #panelShift= xFraction*dx

    #--------------------------------------------------------------------------------------------------------------
    #FSheight = FSH

    filename = 'TLM'

    polyfile = open("BtoPloyBoat.txt", "w")
    realBoat = "TLM_Boat.dat"

    comment = "A B-spline Surface Discretization   "+'\n'

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
    



    npanels = nhl*(nhg+nhf)         #number of hull panels
    npoints = (nhl+1)*(nhg+nhf+1)   #number of hull points

    npoints_uv = (nhl+1)*(nhg+1) # number of u,v coords to store. A u,v point for every pt on hull
    
    




    Lfs = 4.*L #free surface discretization length


    nfsl = 4*nhl                    # number of free surface panels longitudinaly
    nfst = int(0.75*nhl)            # number of free surface panels transversly (may adjust)

    nfspanels = nfsl*nfst           # number of free surface panels

    ntransoml = 5*nhl/2             # number of transom fs panels longitudinaly -> assuming nfsl=4*nhl
    ntransomt = nhg#+nhf             # the number of transom panels transversly  -> equal to the hull panels or something else..
    
    ntransompanels = ntransoml*ntransomt # total number of transom panels

    npoints = npoints + (nfsl+1)*(nfst+1) +(ntransoml+1)*(ntransomt+1) #total number of points


    
    # position of free surface panels above z=0
    # H. Raven suggests: about 0.8 of panel length
    #     zfspanels = 0.8 * L/ nhl
    zfspanels = FSheight#0.65*L / float(nhl)








    ## Allocate memory------------------------------
    panels = np.zeros((4, npanels+nfspanels+ntransompanels),int)
    points = np.zeros((3, npoints),float)
    ldx = np.zeros((nfsl),float)



    u_v_points = np.zeros((2,npoints),float)    # store the u,v coordinates of each panel point
    curvature = np.zeros((npanels+nfspanels),float)   # easier to write curavture on panel centroids
    mean_curvature = np.zeros((npanels+nfspanels),float) 
    

    ## Compute hull offsets
    ## We go from bow to stern and from keel to WL
    
    ipt = 0  # point counter - saved after the loop

    dz = T / float(nhg) # vertical step = draft/num of hull panels girthwise


    # u runs along the transverse b-splines from keel to beam
    # v runs along the longitudinal b-splines from bow to stern
    for j in range(nhg + nhf +1 ):                          #transverse direction including free surface panel.
        #u=float(j)/float(nhg + nhf)
        u=float(j)/float(nhg)
        for i in range(nhl+1):                              #longitudinal direction
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
            
    storelastYtransom = points[1,ipt-1] #y at the aft outboard end of the hull panels
    print points[1,ipt-1], points[1,ipt]


    # number of points on the hull
    nhpt = ipt # correct since it starts at zero and we've tacked on one at the end of the loops.

    print nhpt

    ## Compute hull panels
    ## We go from bow to stern and from keel to WL
    ipan = 0  # panel counter - saved after the loops

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

    
    yfs = zero  #start at midships... 

    # width of free surface - somewhat arbitrary:
    Bfs = 0.
    for j in range(nfst):
        Bfs = Bfs + 1.1**j
    Bfs = dy*Bfs


    # ipan
    # ipoint
    loc_flag = 0    # initialize location flag - recognize when we are on or off the hull
    ii = 1.0        # dx multiplier for the distance behind the last hull panel -> initialize to one

    fflag = 0
    storezfs = np.copy(zfspanels) #assignments make a deep copy, but this is more explicit
    for j in range(nfst + 1):  #transverse
        loc_flag = 0
        ii       = 1.0
        for i in range(nfsl+1):  #longitudinal
            
            points[1,ipt]=0.
            
            #if loc_flag==0: #before the hull
            x = L - float(i)*dx
            if j==0 and i<nfsl: #we are on the first line of hull panels
                ldx[i]=dx

            if (abs(x) <= longiOffset):# and j==0): 
                loc_flag=1 # we are on the hull
                u=1.
                v=0.5*(1.-x/longiOffset)
                points[:,ipt] = surf1.evalSurface(u,v) #get x,y,z location on the surface
                points[0,ipt] = -points[0,ipt] # reverse the x orientation
                points[0,ipt] += longiOffset #move the x position to reflect shift of origin in computational space
                x=points[0,ipt]
                store=x         #becomes the final x loc on the hull
                storelastY=points[1,ipt]    # becomes the final y loc on the hull
                #print 'y=', y, 'points[1,ipt] =',points[1,ipt] 
                #storelastYpt=ipt# also the minimum transverse location of the fs panel pts.


                """Set the v vertical offset by the local dx step size """
                #zfspanels=hFraction*abs(points[0,ipt]-points[0,ipt-1])

            if x < -longiOffset:#we are at the first location behind the hull
                loc_flag=2
                #print 'x=', x
                #print 'loc_flag = ', loc_flag
                
                #print loc_flag , x,store
                x = store - ii*dx  #store is the last x location on the hull.
                
                ii=ii+1.

            #print x

            #modify panels for a transom stern
            #aft of the ship, make the original nfs panels start at the beam at the stern:
            if loc_flag==2:
                pass    #wow - use old y - it will be perfect all the way down!
            else:
                y = (Bfs-points[1,ipt])/Bfs*yfs + points[1,ipt]
                
            if abs(y)<small:
                y=zero
                #print y

            #this is what we put into the panel file:
            points[:,ipt]=[x,y,zfspanels] #use a trick - y is saved at the last longitudinal location - just use that every time down the longi row!
            zfspanels=storezfs
 
            if j==0:
                if i>0 and i<nfsl:
                    #print i, x,ipt, points[0,ipt-1]
                    ldx[i]=abs(x-points[0,ipt-1])
            
            
            
            #points[2,ipt] -= T
            polyfile.write(str(points[0,ipt])+', '+str(points[1,ipt])+', '+str(points[2,ipt])+'\n')
            
            ipt += 1

        yfs = yfs + dy*1.1**(j) 


    ## Compute the free surface panels


    for j in range(1,nfst+1):
        for i in range(1,nfsl+1):
            panels[0,ipan] = nhpt + i + (j-1)*(nfsl+1)          #|  top right
            panels[1,ipan] = nhpt + i + j*(nfsl+1)              #<  top left
            panels[2,ipan] = nhpt + i + 1 + j*(nfsl+1)          #<_ bottom left
            panels[3,ipan] = nhpt + i + 1 + (j-1)*(nfsl+1)      #_  bottom right
            #WRITE(6,'(I4'': '',4(I5,2X))') ipan, panels(:, ipan)
            #polyfile.write(ipan, panels[:, ipan])
            #polyfile.write(str(ipan)+', '+str(panels[0,ipan])+', '+str(panels[1,ipan])+', '+str(panels[2,ipan])++', '+str(panels[3,ipan])'\n')
            polyfile.write(str(panels[0,ipan])+', '+str(panels[1,ipan])+', '+str(panels[2,ipan])+', '+str(panels[3,ipan])+'\n')
            ipan = ipan+1


    tnpt = ipt  # should be the total number of points computed so far
    print tnpt









    ## Compute transom panels:
            
    ## Compute free surface offsets
    ## We go from bow to stern and from midships to max. width

    #dx = Lfs/float(nfsl)
    ##******************************* Improvising an x shift#*************************
    #deltax = panelShift #.8*dx
    ##*******************************Improvising an x shift#**************************
    #dy = freeSurfaceWidthFactor * L # width of first free surface strip midships



    #ntransoml - number of transom panels longitudinaly
    #ntransomt - number of trnsom panels transversely



    Bfs = B     #storelastY#storelastYtransom    # width is equal to the transom width
    yfs = zero          # start at midships

    # width of free surface
    #Bfs = 0.
    #for j in range(ntransomt):
    #    Bfs = Bfs + 1.1**j
    #Bfs = dy*Bfs


    x=store  # start at the aft extent of the ship
    dy = Bfs/float(ntransomt)
    print 'dy', dy
    #loc_flag = 0 #initialize location flag - recognize when we are on or off the hull
    
    for j in range(ntransomt + 1):  #transverse
        ii = 0.                                 #reset x direction counter for next longitudinal pass
        for i in range(ntransoml+1):  #longitudinal
            
            points[1,ipt]=0.
            
            #if loc_flag==0: #before the hull
            #x = L - float(i)*dx
            #if j==0 and i<nfsl:
            #    ldx[i]=dx

            #if (abs(x) <= longiOffset):# and j==0): 
            #    loc_flag=1 # we are on the hull
            #    u=1.
            #    v=0.5*(1.-x/longiOffset)
            #    points[:,ipt] = surf1.evalSurface(u,v) #get x,y,z location on the surface
            #    points[0,ipt] = -points[0,ipt] # reverse the x orientation
            #    points[0,ipt] += longiOffset #move the x position to reflect shift of origin in computational space
            #    x=points[0,ipt]
            #    store=x
            #    storelastY=y
                

            #if x < -longiOffset:#we are at the first location behind the hull
            #    loc_flag=2
                
                #print loc_flag , x,store
            x = store - ii*dx
                
            

            #print x

            #modify panels for a transom stern
            #aft of the ship, make the original nfs panels start at the beam at the stern:
            #$if loc_flag==2:
            #   yfs = storelastY  """This will set up the fs to insert transom fs panels!!"""
            y = (Bfs-points[1,ipt])/Bfs*yfs + points[1,ipt]
            if abs(y)<small:
                y=zero
                #print y

            #this is what we put into the panel file:
            points[:,ipt]=[x,y,zfspanels]
 
            #if j==0:
            #    if i>0 and i<nfsl:
            #        #print i, x,ipt, points[0,ipt-1]
            #        ldx[i]=abs(x-points[0,ipt-1])
            
            
            
            #points[2,ipt] -= T
            polyfile.write(str(points[0,ipt])+', '+str(points[1,ipt])+', '+str(points[2,ipt])+'\n')
            ipt += 1 # update point counter
            ii=ii+1. # update x direction counter
            
    
        #print '\n\n'
        yfs = yfs + dy


    ## Compute the free surface panels
    print ipt



    for j in range(1,ntransomt+1):
        for i in range(1,ntransoml+1):
            panels[0,ipan] = tnpt + i + (j-1)*(ntransoml+1)
            panels[1,ipan] = tnpt + i + j*(ntransoml+1)
            panels[2,ipan] = tnpt + i + 1 + j*(ntransoml+1)
            panels[3,ipan] = tnpt + i + 1 + (j-1)*(ntransoml+1)
            #WRITE(6,'(I4'': '',4(I5,2X))') ipan, panels(:, ipan)
            #polyfile.write(ipan, panels[:, ipan])
            #polyfile.write(str(ipan)+', '+str(panels[0,ipan])+', '+str(panels[1,ipan])+', '+str(panels[2,ipan])++', '+str(panels[3,ipan])'\n')
            polyfile.write(str(panels[0,ipan])+', '+str(panels[1,ipan])+', '+str(panels[2,ipan])+', '+str(panels[3,ipan])+'\n')
            ipan = ipan+1



    #






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
    realFile.write('    '+str(npanels)+'    '+str(nfspanels)+'    '+str(ntransompanels)+'    '+str(npoints)+'\n')
    realFile.write('    '+str(deltax)+'\n')
    realFile.write('    '+str(nfsl)+'\n') 
    realFile.write('    '+str(nfst)+'\n') 
    realFile.write('    '+str(ntransoml)+'\n')
    realFile.write('    '+str(ntransomt)+'\n')
    realFile.write('    '+str(nhl)+'\n')  # num hull panels longitudinaly
    for i in range(npanels+nfspanels+ntransompanels):
        realFile.write(str(panels[0,i])+'   '+str(panels[1,i])+'    '+str(panels[2,i])+'    '+str(panels[3,i])+'\n')
    for i in range(npoints):
        realFile.write(str(points[0,i])+'   '+str(points[1,i])+'    '+str(points[2,i])+'\n')
    for i in range(len(ldx)):
        realFile.write(str(ldx[i])+'\n')
    #if Compute_Curvatures==True:
    #for i in range(len(curvature)):
    #    realFile.write(str(curvature[i])+'\n')
    #for i in range(len(mean_curvature)):
    #    realFile.write(str(mean_curvature[i])+'\n')
    realFile.close()
    #-----------------------------------------------------------------------------------------------------



    return points, panels, Bfs,yfs,npoints,npanels,nfspanels,ntransompanels,ntransomt, ntransoml,tnpt
