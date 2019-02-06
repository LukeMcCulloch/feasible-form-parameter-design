# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 18:54:51 2015

@author: Vanderplas
https://jakevdp.github.io/blog/2012/11/24/simple-3d-visualization-in-matplotlib/

About Quaternions:

    Transformation to and from Axis-Angle Representation:
    
    
    Axis-Angle representation to and from Quaternion:
        GTCG p 897
        
        a = (xo,yo,zo)  : axis of rotation
        theta           : angle of rotation
        
        (1) A-A to Quaternion = w + xi + yj + zk 
            w = cos(theta/2)
            x = a[0]*sin(theta/2)
            y = a[1]*sin(theta/2)
            z = a[2]*sing(theta/2)
            
        (2) Q to Axis-Angle:
            theta   = 2*arccos(w)
            a[0]    = x/sin(theta/2)
            a[1]    = y/sin(theta/2)
            a[2]    = z/sin(theta/2)
"""
import matplotlib.pyplot as plt

import numpy as np



class Quaternion:
    """Quaternions for 3D rotations
    
        self.x = [w, x, y, z] : the Quaternion Quadruple
    """
    def __init__(self, x):
        self.x = np.asarray(x, dtype=float)
        
    @classmethod
    def from_v_theta(cls, v, theta):
        """ Implement (1) above
        Construct quaternion (and normalize it) 
            v       : unit vector
            theta   : rotation angle 
        """
        theta = np.asarray(theta)
        v = np.asarray(v)
        
        s = np.sin(0.5 * theta)
        c = np.cos(0.5 * theta)
        #vnrm = np.sqrt(np.sum(v * v))
        vnrm = np.sqrt(np.dot(v,v))

        q = np.concatenate([[c], s * v / vnrm])
        return cls(q)

    def __repr__(self):
        return "Quaternion:\n" + self.x.__repr__()

    def __mul__(self, other):
        """ multiplication of two quaternions.
        see: Schneider, Eberly
        Geometric Tools for Computer Graphics p. 896"""
        prod = self.x[:, None] * other.x

        return self.__class__([(prod[0, 0] - prod[1, 1]
                                 - prod[2, 2] - prod[3, 3]),
                                (prod[0, 1] + prod[1, 0]
                                 + prod[2, 3] - prod[3, 2]),
                                (prod[0, 2] - prod[1, 3]
                                 + prod[2, 0] + prod[3, 1]),
                                (prod[0, 3] + prod[1, 2]
                                 - prod[2, 1] + prod[3, 0])])

    def as_v_theta(self):
        """Implement (2) above
        Starting wih a normalized Quaternion, Q.x = [w,x,y,z]
        Return the v, theta equivalent representation
        """
        # compute theta
        #norm = np.sqrt((self.x ** 2).sum(0))
        norm = np.sqrt(np.dot(self.x,self.x))
        theta = 2 * np.arccos(self.x[0] / norm)

        # compute the unit vector
        v = np.array(self.x[1:], order='F', copy=True)
        #v /= np.sqrt(np.sum(v ** 2, 0))
        v /= np.sqrt(np.dot(v,v))

        return v, theta

    def as_rotation_matrix(self):
        """Return the rotation matrix of the (normalized) quaternion"""
        v, theta = self.as_v_theta()
        c = np.cos(theta)
        s = np.sin(theta)

        return np.array([[v[0] * v[0] * (1. - c) + c,
                          v[0] * v[1] * (1. - c) - v[2] * s,
                          v[0] * v[2] * (1. - c) + v[1] * s],
                         [v[1] * v[0] * (1. - c) + v[2] * s,
                          v[1] * v[1] * (1. - c) + c,
                          v[1] * v[2] * (1. - c) - v[0] * s],
                         [v[2] * v[0] * (1. - c) - v[1] * s,
                          v[2] * v[1] * (1. - c) + v[0] * s,
                          v[2] * v[2] * (1. - c) + c]])
                          
                          


class CubeAxes(plt.Axes):
    """An Axes for displaying a 3D cube"""
    # fiducial face is perpendicular to z at z=+1
    one_face = np.array([[1, 1, 1], [1, -1, 1], [-1, -1, 1], [-1, 1, 1], [1, 1, 1]])

    # construct six rotators for the face
    x, y, z = np.eye(3)
    rots = [Quaternion.from_v_theta(x, theta) for theta in (np.pi / 2, -np.pi / 2)]
    rots += [Quaternion.from_v_theta(y, theta) for theta in (np.pi / 2, -np.pi / 2)]
    rots += [Quaternion.from_v_theta(y, theta) for theta in (np.pi, 0)]
    
    # colors of the faces
    colors = ['blue', 'green', 'white', 'yellow', 'orange', 'red']
    
    def __init__(self, fig, rect=[0, 0, 1, 1], *args, **kwargs):
        # We want to set a few of the arguments
        kwargs.update(dict(xlim=(-2.5, 2.5), ylim=(-2.5, 2.5), frameon=False,
                           xticks=[], yticks=[], aspect='equal'))
        super(CubeAxes, self).__init__(fig, rect, *args, **kwargs)
        self.xaxis.set_major_formatter(plt.NullFormatter())
        self.yaxis.set_major_formatter(plt.NullFormatter())
        
        # define the current rotation
        self.current_rot = Quaternion.from_v_theta((1, 1, 0), np.pi / 6)
        
    
    def draw_cube(self):
        """draw a cube rotated by theta around the given vector"""
        # rotate the six faces
        Rs = [(self.current_rot * rot).as_rotation_matrix() for rot in self.rots]
        faces = [np.dot(self.one_face, R.T) for R in Rs]
        
        # project the faces: we'll use the z coordinate
        # for the z-order
        faces_proj = [face[:, :2] for face in faces]
        zorder = [face[:4, 2].sum() for face in faces]
        
        # create the polygons if needed.
        # if they're already drawn, then update them
        if not hasattr(self, '_polys'):
            self._polys = [plt.Polygon(faces_proj[i], 
                                       fc=self.colors[i],
                                       alpha=0.9, 
                                       zorder=zorder[i]) 
                                       for i in range(6)]
            for i in range(6):
                self.add_patch(self._polys[i])
        else:
            for i in range(6):
                self._polys[i].set_xy(faces_proj[i])
                self._polys[i].set_zorder(zorder[i])
                
        self.figure.canvas.draw()
        
    def draw_surface(self, surface):
        """TLM methd to draw a Bspline surface
        surface : curve.py Bspline surface
        """
        return
        
        

        
        
class CubeAxesInteractive(CubeAxes):
    """An Axes for displaying an Interactive 3D cube"""
    def __init__(self, *args, **kwargs):
        super(CubeAxesInteractive, self).__init__(*args, **kwargs)
        
        # define axes for Up/Down motion and Left/Right motion
        self._v_LR = (0, 1, 0)  #Left Right
        self._v_UD = (-1, 0, 0) #Up Down
        self._active = False
        self._xy = None

        # connect some GUI events
        self.figure.canvas.mpl_connect('button_press_event',
                                       self._mouse_press)
        self.figure.canvas.mpl_connect('button_release_event',
                                       self._mouse_release)
        self.figure.canvas.mpl_connect('motion_notify_event',
                                       self._mouse_motion)
        self.figure.canvas.mpl_connect('key_press_event',
                                       self._key_press)
        self.figure.canvas.mpl_connect('key_release_event',
                                       self._key_release)
        
        self.figure.text(0.05, 0.05, ("Click and Drag to Move\n"
                                    "Hold shift key to adjust rotation"))

    #----------------------------------------------------------
    # when the shift button is down, change the rotation axis
    # of left-right movement
    def _key_press(self, event):
        """Handler for key press events"""
        if event.key == 'shift':
            self._v_LR = (0, 0, -1)
#        if event.key == 'left':
#            self.rotate(self.y, np.pi / 12)
#        elif event.key == 'right':
#            self.rotate(self.y, -np.pi / 12)
#        elif event.key == 'up':
#            self.rotate(self.x, np.pi / 12)
#        elif event.key == 'down':
#            self.rotate(self.x, -np.pi / 12)
#        fig.canvas.draw()

    def _key_release(self, event):
        """Handler for key release event"""
        if event.key == 'shift':
            self._v_LR = (0, 1, 0)

    #----------------------------------------------------------
    # while the mouse button is pressed, set state to "active"
    # so that motion event will rotate the plot
    def _mouse_press(self, event):
        """Handler for mouse button press"""
        if event.button == 1:
            self._active = True
            self._xy = (event.x, event.y)

    def _mouse_release(self, event):
        """Handler for mouse button release"""
        if event.button == 1:
            self._active = False
            self._xy = None

    def _mouse_motion(self, event):
        """Handler for mouse motion"""
        if self._active:
            dx = event.x - self._xy[0]
            dy = event.y - self._xy[1]
            self._xy = (event.x, event.y)
            
            rot1 = Quaternion.from_v_theta(self._v_UD, 0.01 * dy)
            rot2 = Quaternion.from_v_theta(self._v_LR, 0.01 * dx)

            self.current_rot = (rot2 * rot1 * self.current_rot)
            self.draw_cube()
            
if __name__ =='__main__':
    if True:
        fig = plt.figure(figsize=(4, 4))
        ax = CubeAxesInteractive(fig)
        fig.add_axes(ax)
        ax.draw_cube()
        #display(fig)
        
#    if True:
#        c = Cube()
#        fig, ax = plt.subplots(figsize=(8, 8),
#                               subplot_kw=dict(xticks=[], yticks=[]))
#        c.add_to_ax(ax)
#        c.rotate(c.x - c.y, -np.pi / 6)
#        ax.set_xlim(-10, 10)
#        ax.set_ylim(-10, 10)
#        
#        
#        def zoom(val):
#            c.set_view((0, 0, val))
#            ax.set_xlim(-val, val)
#            ax.set_ylim(-val, val)
#            fig.canvas.draw()
#            
#        slider_ax = fig.add_axes((0.2, 0.05, 0.6, 0.02))
#        slider = Slider(slider_ax, "perspective", 1, 20, valinit=10, color='#AAAAAA')
#        slider.on_changed(zoom)

        