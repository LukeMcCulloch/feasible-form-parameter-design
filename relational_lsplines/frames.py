# -*- coding: utf-8 -*-
"""
Created on Wed Dec 02 10:30:55 2015

@author: Luke.McCulloch

http://eli.thegreenplace.net/2015/change-of-basis-in-linear-algebra/
"""
import numpy as np
import cmath as cp

"""
MOSES testing:

At the Disposition menu


freq_response
     rao -period 20 -heading 0
     fr_point 0 0 0  
     fr_point &part(cg  jacket -body)  $$ Georgina says body sys
     fr_point &part(cg  jacket -part)  $$ itr.mac says part sys
     set_variable    hhh -col  3 4 5 6 7 8 9  10  11  12  13  14
end freq_response

&type %hhh%

$ using the part system above gives bad answers

$
$  BUT LOOK AT THIS!
$
&describe part jcg
*jcg &part(cg %par% -p)
&set sacspnts     = %sacspnts  *jcg
&set loc_pt = &token(1 %sacspnts)

freq_response
     rao -period 20 -heading 0
     fr_point %loc_pt
     set_variable    hhh -col  3 4 5 6 7 8 9  10  11  12  13  14
end freq_response

"""
#amp_mask = np.asarray([0==i%2 for i in range(self.dim*2)],dtype=bool)
#phase_mask = np.asarray([1==i%2 for i in range(self.dim*2)],dtype=bool)
def cross(s,r):
    v1 = s[1]*r[2] - r[1]*s[2]
    v2 = s[2]*r[0] - r[2]*s[0]
    v3 = s[0]*r[1] - r[0]*s[1]
    return np.asarray([v1,v2,v3])


def cart2pol(x, y):
    """not used 
    """
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

def pol2cart(rho, phi):
    """not used 
    """
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return(x, y)

polar2z = lambda r,theta: r * np.exp( 1j * theta )
z2polar = lambda z: ( np.abs(z), np.angle(z) )


class Motions(object):
    """Assumes standard basis
    initially
    """
    def __init__(self, rao, origin = None):
        self.rao        = rao
        self.dim        = 6
        if origin is None:
            self.origin = np.asarray([0.,0.,0.])
        else:
            self.origin = origin
        self.amp_mask, self.phase_mask = self.masks()
        self.rao_amp    = rao[self.amp_mask]
        self.rao_phase  = rao[self.phase_mask]
        self.i          = cp.sqrt(-1)
        self.sp         = self.rao_to_cartesian(self.rao_amp, self.rao_phase)
        return
    
    def masks(self):
        amp_mask = np.asarray([0==i%2 for i in range(self.dim*2)],dtype=bool)
        phase_mask = np.asarray([1==i%2 for i in range(self.dim*2)],dtype=bool)
        return amp_mask, phase_mask
        
    def rao_to_cartesian(self, rao_amp, rao_phase):
        """returns the cartesian form of the amp phase
        also changes representation from degrees to radians.
        """
        sp = np.zeros((self.dim),complex)
        sp[0:self.dim/2] = rao_amp[0:self.dim/2]*np.exp( self.i*np.radians(rao_phase[0:self.dim/2]) )
        sp[self.dim/2:self.dim] = np.radians(rao_amp[self.dim/2:self.dim]) \
                                    *np.exp(self.i*np.radians(rao_phase[self.dim/2:self.dim]))
        
        """
        sp = []
        for amp,phase in zip(rao_amp[0:3],rao_phase[0:3]):
            sp.append(amp*np.exp(self.i*np.deg2rad(phase) ) ) #phase*np.pi/180.))#

        for amp,phase in zip(rao_amp[3:],rao_phase[3:]):
            sp.append(np.radians(amp)*np.exp(self.i*np.deg2rad(phase) ) )  #phase*np.pi/180.))#
        sp = np.asarray(sp)
        #"""
        return sp
        
    def rao_to_polar(self, complex_in):
        out_rao = np.zeros(6,float)
        amp     = np.sqrt( (np.real(complex_in))**2 + (np.imag(complex_in))**2)
        phase   = np.arctan2(  np.imag(complex_in)  ,  np.real(complex_in) )
        for i in range(self.dim/2):
            out_rao[i*2] = amp[i]
            out_rao[i*2+1] = phase[i]
        """
        out_rao = []
        for el in complex_in:
            #amp     = np.sqrt(np.conjugate(el)*el)
            amp     = np.sqrt( (np.real(el))**2 + (np.imag(el))**2)
            phase   = np.arctan2(  np.imag(el)  ,  np.real(el) )
            out_rao.append(amp)
            out_rao.append(phase)
        #"""
        return np.asarray(out_rao)
            
    
    def fr_point(self, point):
        r = point - self.origin
        translations = self.rao_to_polar(self.small_angle_transform(self.sp,r))
        rao = np.zeros(12)
        rao[0:6] = translations
        rao[6:] = self.rao[6:]
        for i in [1,3,5]:
            rao[i] = np.rad2deg(rao[i])
        return rao
        
    def small_angle_transform(self, sp, r):
        return sp[:3] + cross(sp[3:],r)
        
    def fr_point_general(self, point, jacket_basis ):
        rjcg                = self.fr_point(point)
        a                   = rjcg[self.amp_mask]
        p                   = rjcg[self.phase_mask]  
        cart                = self.rao_to_cartesian(a,p)
        tr                  = DeepVector( v = cart[0:3], B=E3 )
        rr                  = DeepVector( v = cart[3:6], B=E3 )
        translations        = tr.change_basis(B=jacket_basis)
        rotations           = rr.change_basis(B=jacket_basis)
        polar_translations  = self.rao_to_polar(translations.v)
        polar_rotations     = np.degrees(self.rao_to_polar(rotations.v))
        
        rao = np.zeros(12)
        rao[0:6] = polar_translations
        rao[6:] = polar_rotations
        for i in [1,3,5]:
            rao[i] = np.rad2deg(rao[i])
        return rao


class Frame(object):
    def __init__(self, B):
        self.n = len(B)
        self.E = np.identity(self.n)
        self.B = B
        self.get_transfermatrix_BE()
        return
        
    def __call__(self):
        return self.B
        
    def __str__(self):
        return  str(self.B)
    
    def __repr__(self):
        return 'Frame(B = {})'.format(self.B)
    
    def __getitem__(self, key):
        return self.B[key]
    
    def __setitem__(self, key, value):
        self.B[key] = value
        return
    
    def get_transfermatrix_BE(self):
        self.BtoE = self.B.T
        self.EtoB = np.linalg.inv(self.BtoE)
        return
        
    def make_DCM(self, M=None):
        if M is None:
            M = self.E
        self.DCM = np.outer(self.B,M)
        return
        
        
    
    """ Not finished
    def __add__(self, other):
        if isinstance(other, np.ndarray):
            return Frame(B = self.B + other)
        elif isinstance(other, Frame):
            return Frame(B = self.B + other.B)
        else:
            print 'undefined addition {} + {}'.format(self, other)
            return
        
    def __mult__(self, other):
        if isinstance(other, np.ndarray):
            return Frame(B = np.dot(self.B , other) )
        elif isinstance(other, Frame):
            return Frame(B = np.dot(self.B , other.B) )
        else:
            print 'undefined addition {} + {}'.format(self, other)
            return
    #"""
        

class DeepVector(object):
    def __init__(self, v, B=None):
        self.n = len(v)
        self.v = v
        if B is None:
            self.B = np.identity(self.n)
        else:
            self.B = B
        return
        
    def __call__(self):
        return self.v
        
    def __str__(self):
        return  str(self.v)
    
    def __repr__(self):
        return 'DeepVector(v={}, B = {})'.format(self.v, self.B)
    
    def __getitem__(self, key):
        return self.v[key]
    
    def __setitem__(self, key, value):
        self.v[key] = value
        return
    
    def change_basis(self, B):
        return DeepVector( np.dot(B.EtoB, np.dot(self.B.BtoE,self.v.T)), B  )







    


if __name__ == '__main__':
    testing = True
    
    
    E  = Frame( np.asarray([[1.,0.],[0.,1.]]) )
    B1 = Frame( np.asarray([[2.,3.],[4.,5.]]) )
    B2 = Frame(np.asarray([[-1.,1.],[1.,1.]]) )
    
    
    
    v1_B1 = DeepVector( v = np.asarray([3.,-1.]), B = B1)
    v1_B2 = v1_B1.change_basis(B = B2)
    v1_E = v1_B1.change_basis(B = E)
    
    
    
    E3  = Frame( np.asarray([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]) )
    J3  = Frame( np.asarray([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]) )
    jacket_basis_ = Frame( np.asarray([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]) )
    orient_basis = Frame(np.asarray([[0.9922781, 0.1240329, 0.],[-0.1240329, 0.9922781, 0.],[0., 0., 1.]]))
    #B31 = Frame( np.asarray([[2.,3.,0.],[4.,5.,0.]]) )
    #B32 = Frame(np.asarray([[-1.,1.,0.],[1.,1.,0.]]) )
    
    q = np.asarray([3.,-1.,0.]) #standard basis
    qE = DeepVector( v = q, B=E3 )
    qJ1 = qE.change_basis(B=J3)
    #qB1 = qE.change_basis(B = B31)
    #qB2 = qE.change_basis(B = B32)
    #qBE = qE.change_basis(B = E3)
    
    
    ##
    ##  MOSES Example
    ##
    bodyframe = Frame(np.identity(3))
    
    
    #faltinsen box
    #print 'faltinsen box--------------------------------------------------'
    #Ov = DeepVector(v = np.asarray([0.,0.,0.]) , B = bodyframe)
    #frao = np.asarray([0., 0., 3319., -72., 375., 77., 60134., 106., 7., -162., 8., 52.])
    #self = Motions(frao)
    #print self.fr_point_general(Ov.v, jacket_basis = jacket_basis_)
    #print '---------------------------------------------------------------'
    ##
    ## tow_auto
    ##
    jacketloc = DeepVector(v = np.asarray([225.,0.,50.]) , B = bodyframe)
    rao_origin = [ 0.983108, 116.195, 1.7119E-6, 0., 1.07627, -1.06366, 3.32471E-6, 0., 0.186382, -65.1966,  5.06892E-8, 0.  ]
    
    rao_origin = np.asarray(rao_origin)
    print '\nSingle frequency RAO at the origin'
    print '{} \n'.format(rao_origin)
    
    self = Motions(rao_origin)
    point  = jacketloc.v
    print '\n computing motions from tow_atuto.cif : RAO at point {}\n'.format(point)
    rjcg = self.fr_point(jacketloc.v)
    print 'in the body frame'
    print rjcg
    
    print '\n Give the motion RAO at point {} in the jacket basis:'.format(jacketloc.v)
    rao = self.fr_point_general(jacketloc.v, jacket_basis = jacket_basis_)
    print rao
    
    
    
    ##
    ## Orient
    ##
    quadloc = DeepVector(v = [50.65512, 6.536064, 18.92633], B=jacket_basis_)# B=orient_basis)
    quadloc1 = DeepVector(v = [ 62.25392, -0.2102046, 3.610219], B=jacket_basis_)# B=orient_basis)
    orient_rao_origin = [0.932793, 126.687, 1.12693E-3, 131.412, 1.09329, 4.02739, 5.28684E-5, 0., 0.529304, -54.626, 1.0471E-3, -48.6016  ]
    
    
    orient_rao_origin = np.asarray(orient_rao_origin)
    print '\nSingle frequency RAO at the origin'
    print '{} \n'.format(rao_origin)
    
    print 'ORIENT starts here-----------------------'
    self = Motions(orient_rao_origin)
    print '\n computing motions from orient.cif : RAO at point {}\n'.format(quadloc.v)
    rjcg = self.fr_point(quadloc.v)
    print '\n in the body frame'
    print rjcg
    print '\n at 0 0 0 (also in the body frame)'
    print self.fr_point([0.,0.,0.])
    #print '\nat the cg of the jacket, in the body frame at point {}\n'.format(quadloc1.v)
    #print self.fr_point(quadloc1.v)
    #print '\nat cg of jacket, in jacket frame'
    #print self.fr_point_general(quadloc1.v, jacket_basis = orient_basis)
    print '\n--------------------MOST important Case from Orient-------------------------'
    print '\n Give the motion RAO at point quad cg {} in the orient basis:'.format(quadloc.v)
    rao = self.fr_point_general(quadloc.v, jacket_basis = orient_basis)
    print rao
    
    
    print '\n--------------------faltinsen box-------------------------'
    Ov = DeepVector(v = np.asarray([0.,0.,0.]) , B = bodyframe)
    frao = np.asarray([0., 0., 3319., -72., 375., 77., 60134., 106., 7., -162., 8., 52.])
    self = Motions(frao)
    print self.fr_point_general(Ov.v, jacket_basis = jacket_basis_)
    
    #"""
    fkreal = np.asarray([ -0.3418, 968.16, -821.76, -1.6496E4, -7.0475, 0.87099 ])
    fkimag = np.asarray([ 0.16239, 294.25, 365.44, -2.1862E4, -1.1849, -6.4949 ])
    dreal = np.asarray([ -2.9804E-5, -2.5547E-5, 905.06, 2.5471E-3, -6.8011E-3, 3.9378])
    dimag = np.asarray([-6.8235E-2, -3469.2, 1.6605E-4, 7.969E4,  -0.94477, 1.1074E-2])
    
    freal = fkreal + dreal
    fimag = fkimag+dimag
    complex_in = np.zeros((6),complex)
    complex_in.real[:] = freal[:]
    complex_in.imag[:] = fimag[:]
    out_rao = []
    for el in complex_in:
        #amp     = np.sqrt(np.conjugate(el)*el)
        amp     = np.sqrt( (np.real(el))**2 + (np.imag(el))**2)
        phase   = (np.arctan2(  np.imag(el)  ,  np.real(el) ) )
        out_rao.append(amp)
        out_rao.append(phase)
    for i in [1,3,5]:
        out_rao[i] = np.rad2deg(out_rao[i])
    for i in [7,9,11]:
        out_rao[i] = np.rad2deg(out_rao[i])
    for i in [6,8,10]:
        out_rao[i] = np.deg2rad(np.rad2deg(out_rao[i]))
    #"""
    
    print '\n--------------------Doug-------------------------'
    print '\nDougs RAO'
    og = DeepVector(v = np.asarray([92.8,0.0,12.4]))
    #og = DeepVector(v = np.asarray([0.,0.,0.]))
    dougrao = [0.348, -139.,  0.000, 0.,  0.366,  133.,  0.003,   -26., 0.810,  39., 0.000, 0.]
    dougrao = np.asarray(dougrao)    
    self = Motions(dougrao, origin = og.v)
    #eg = DeepVector(v = [-92.8,0.0,-12.4])
    eg = DeepVector(v = [57.315, 1.546,  28.157 ])
    print self.fr_point(eg.v)
    print '\nCheck it in a new system'
    sos = np.cos(np.radians(45.))
    modbasis = Frame( np.asarray([[sos,-sos,0.],[sos,sos,0.],[0.,0.,1.]]) )
    print self.fr_point_general(eg.v, jacket_basis = modbasis)