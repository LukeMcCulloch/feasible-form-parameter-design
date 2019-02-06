## 20130903
## TLM Python Module to
## compute u,v coordinates
## of the panel centroids

import numpy as np

## - not finished or used at the moment - lets just use u,v here instead

def panelCentroids(corners):


    #1st compute the midpoints of all four panel edges:
    m1=0.5*(corners[:,0]+corners[:,1])
    m2=0.5*(corners[:,1]+corners[:,2])
    m3=0.5*(corners[:,2]+corners[:,3])
    m4=0.5*(corners[:,3]+corners[:,1])

    #2nd compute the coordsys local to the panel
    #iv = m2-m4
    #iv = iv/np.sqrt(sum(iv**2))

    #jvbar=m3-m1

    #Normal Vector:
    #nv = crossproduct(iv,jvbar)
    #nv=nv/np.sqrt(sum(nv**2))

    #Find the second tangent vector ortho to nv & iv
    #jv = crossproduct(nv,iv)

    center = 0.25*(m1+m2+m3+m4) #initial center guess

##    area = 0.
##    momentxi = 0.
##    momenteta = 0.
##
##    for i in range(0,4):
##        dxi  = corners[0,ip1[i]]-corners[0,i]
##        deta = corners[1,ip1[i]]-corners[1,i]
##        area += deta*(corners[0,ip1[i]]+corners[0,i])
##        momentxi += dxi*(corners[1,ip1[i]]**2 + (corners[1,i]*corners[1,ip1[i]]) + corners[0,i]**2)
##        momenteta += deta*( corners[0,ip1[i]]**2 + (corners[0,i]*corners[0,ip1[i]] ) + corners[0,i]**2 )
##
##    area        =  0.5*area
##    momentxi    = -1./6.*momentxi
##    momenteta   =  1./6.*momenteta
##
##    if (area > 0.0000000001):
##        xic    =  momenteta / area
##        etac   =  momentxi  / area
##    else:
##        xic     = 0.
##        etac    = 0.
##
##    center = center + xic*iv + etac*jv
    
    return center
            

        
