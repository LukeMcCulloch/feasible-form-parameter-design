# -*- coding: utf-8 -*-
"""
Created on Wed Apr 29 18:54:51 2015




@author: 

    JVP's wonder cube is the place from which
    much of the machinery for drawing 3D->2D was taken,
    of all that is in this file:
        https://jakevdp.github.io/blog/2012/11/24/simple-3d-visualization-in-matplotlib/
    
    
Luke McCulloch:
    Curves, and Optimization
    Oh, and I made up the translation operation below so 
    if it's sort of bad, it's my fault.  It worked, so I stopped there.



Matplotlib docs and examples:
    for point picking and zooming

About Quaternions:

    moved to Quaternion.py
"""
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
#from matplotlib.widgets import Button
import numpy as np
import copy

import curve             as spline
import GUIspline         as gs
from   ADILS             import IntervalLagrangeSpline, Lagrangian
from   FormParameter     import FormParameterDict
from   initialValues     import interval_bounds, lagrangian_bounds, InitializeControlVertices

from Quaternion import Quaternion
#from   initialValues     import InitializeControlPoints
#from interval_arithmetic import ia
#from automatic_differentiation import ad



class DrawCurve(plt.Axes):
    """Basics for displaying pickable, movable 3D curves
    that don't make you go blind (zorder!).
    """
    
    def __init__(self, 
                 fig, 
                 rect=[0, 0, 1, 1], 
                 sections_=None, 
                 xkcd=False, 
                 is_3d=False,
                 *args, **kwargs):
        if xkcd:
            plt.xkcd()
        real_faces = None
        sections = sections_.sections
        bow_id         = sections_.bow_id
        for key in sections:
            if key == bow_id:
                if real_faces is None:
                    cv, vt = self.get_bow(sections[key])
                    real_faces  = [cv]
                    vertices    = [vt]
                else:
                    cv, vt = self.get_bow(sections[key])
                    real_faces.append(cv)
                    vertices.append(vt)
            else:
                if real_faces is None:
                    cv, vt = self.get_curve(sections[key])
                    if len(np.unique(vt[:,0]))==1 or \
                        len(np.unique(vt[:,1]))==1 or \
                        len(np.unique(vt[:,2]))==1 :
                        real_faces  = [cv]
                    vertices    = [vt]
                else:
                    cv, vt = self.get_curve(sections[key])
                    if len(np.unique(vt[:,0]))==1 or \
                        len(np.unique(vt[:,1]))==1 or \
                        len(np.unique(vt[:,2]))==1 :
                        real_faces.append( cv)
                    vertices.append(vt)
        self.hull_frame   = sections_
        self.real_faces = real_faces
        self.vertices   = vertices
        self.colors = ['green','blue', 'green', 
                           'white', 'yellow', 'orange', 'red','grey']
        self._show_curvature = False
        
        x, y, z = np.identity(3)
        kwargs.update(dict(xlim=(-25., 25.), 
                           ylim=(-25., 25.), 
                           frameon=True,
                           xticks=[], 
                           yticks=[], 
                           aspect='equal'))
        super(DrawCurve, self).__init__(fig, rect, 
                                        *args, **kwargs)
        self.rots = Quaternion.from_v_theta(x, 0.)
        self.xaxis.set_major_formatter(plt.NullFormatter())
        self.yaxis.set_major_formatter(plt.NullFormatter())
        self.nfaces = len(real_faces)
        # define the current rotation axis and angle
        self._current_rot = Quaternion.from_v_theta((1., 0., 0.), 0.)
    
    def update_face_and_vertices(self, index):
        #print 'update curve',index
        bow_id = self.hull_frame.bow_id
        sections = self.hull_frame.sections
        if bow_id == index:
            cv,vt = self.get_bow(sections[index])
        else:
            cv,vt = self.get_curve(sections[index])
        if len(np.unique(vt[:,0]))==1 or \
            len(np.unique(vt[:,1]))==1 or \
            len(np.unique(vt[:,2]))==1 :
            self.real_faces[index] = cv
        self.vertices[index] = vt
        return
    
    def draw_curves(self):
        """draw a cube rotated by theta around the given vector"""
        #rotation matrix:
        Rs = (self._current_rot * self.rots).as_rotation_matrix()
        
        #translation vector:
        Ts = np.asarray([self._current_trans_x,self._current_trans_y,0.])
        
        # rotate and translate the faces
        faces = [np.dot(face,Rs) + Ts for face in self.real_faces]
        
        #rotate and translate vertices 
        verts = [np.dot(vt,Rs) + Ts for vt in self.vertices]
            
        
        #curvature:
        null_cv = np.zeros([30,2],float)
        if not hasattr(self, '_polys'):
            self.current_curve = 0
            self.cvp = self.hull_frame.sections[0][0].curve.compute_curve_of_curvature()
            cvp = self.get_curvature(self.cvp,0)
            # rather opeque computation 
            # to compound the curve of curvature with the  reversed Bspline curve 
            # to make a contiguous 'polygon surface'.
            cvp = np.asarray( [ np.concatenate( [ cvp[:,el],
                                self.real_faces[self.current_curve][::-1][1:-1,el]  
                                                ] ) for el in [0,1,2] ] ).T
            cvp = np.dot(cvp,Rs)
            cvp += Ts
            cvproj = cvp[:, :2]
            cvorder = cvp[:30,2].sum()
        if self._show_curvature:
            which = self.current_curve
            cvp = self.cvp
            # trick to compound the 0 to 1 curvature with the 1 to 0 curve 
            # for a contiguous 'polygon'.  [::-1] reverses the face points
            cvp = np.asarray( [ 
                                np.concatenate( [cvp[:,el],
                                self.real_faces[which][::-1][1:-1,el] 
                                                ]) for el in [0,1,2] ] ).T
            cvp = np.dot(cvp,Rs) + Ts
            cvproj = cvp[:, :2]
            cvorder = cvp[:30,2].sum()
            
            
        # project the faces: z-order by z coord
        faces_proj = [face[:, :2] for face in faces] #get face[x,y]
        zorder = [face[:32, 2].sum() for face in faces] #sum the z cooridnate
        
        # now the vertices
        pts_proj = [vt[:, :2] for vt in verts]
        zptorder = [vt[:32, 2].sum() for vt in verts] 
        
        
        
        
        # create the polygons if needed.
        # if they're already drawn, then update them
        if not hasattr(self, '_polys'):
            self._polys = [plt.Polygon(faces_proj[i], 
                                       fc=self.colors[0],
                                       alpha=0.2, 
                                       zorder=zorder[i],
                                       picker=True) 
                                       for i in range(self.nfaces)]
            #"""
            # Polygon is, as well as I can tell,
            # required to do picking on the control vertices
            self._controlpoly = [plt.Polygon(pts_proj[i], 
                                       fc=self.colors[0],
                                       alpha=.0, 
                                       zorder=zptorder[i],
                                       fill = False, 
                                       closed = False) 
                                       for i in range(self.nfaces)]
            #"""
            #"""
            # Line2D is, as well as I can tell,
            # required to visualize the control vertices
            self._markers = [Line2D(pts_proj[i][:,0],
                                    pts_proj[i][:,1],
                                       color='blue',
                                       marker='o', 
                                       markerfacecolor='r',
                                       alpha = .4,
                                       zorder=zptorder[i],
                                       picker=True) 
                                       for i in range(self.nfaces)]
            #"""
            #print cvproj
            self._cv_poly = plt.Polygon(cvproj, 
                                       fc=self.colors[7],
                                       alpha=.2, 
                                       zorder=cvorder,
                                       fill = True, 
                                       closed = False) 

            
            for i in range(self.nfaces):
                self.add_patch(self._polys[i])          #show the curve as area
                self.add_patch(self._controlpoly[i])    #picking control vertices
                self.add_line(self._markers[i])         #seeing control vertices
            self.add_patch(self._cv_poly) #curvature
            self._cv_poly.set_xy(null_cv) #curvature stuff
            self.tag_polys()
        else: #just update
            for i in range(self.nfaces):
                self._polys[i].set_xy(faces_proj[i]) #show the curve as area
                self._polys[i].set_zorder(zorder[i]) #ordering
                
                self._controlpoly[i].set_xy(pts_proj[i])        #vpicking control vertices
                self._controlpoly[i].set_zorder(zptorder[i])    #ordering
                
                self._markers[i].set_data(zip(*pts_proj[i])) #seeing control vertices #--inconsistent api?--#
                self._markers[i].set_zorder(zptorder[i]) 
            
            if self._show_curvature:
                #print cvproj
                self._cv_poly.set_xy(cvproj)
                #self._cv_poly.set_zorder(cvorder)
            else:
                self._cv_poly.set_xy(null_cv)
                
        self.figure.canvas.draw()
        
        
    def tag_polys(self):
        for s, poly in zip(self.hull_frame.sections, self._polys ):
            poly.name = s
        return
        
    #def create_graphical_face(self, curve):
    #    return

#    def on_pick_which(self, event):
#        """a failed attempt to find the picked surface!
#        """
#        self.event = event
#        x=self.event.mouseevent.x
#        y=self.event.mouseevent.y
#        for poly in self._polys:
#            if poly.contains(event.mouseevent):
#                #print 'found poly [{}]'.format( poly.name)
#                self.active_poly = poly
#                break
#        return

    def get_surface(self, surf):
        self.surface = surf
        self.plot_wireframe(self.surface[:,:,0], 
                          self.surface[:,:,1], 
                          self.surface[:,:,2], 
                          rstride=3, 
                          cstride=3, 
                          color = 'r',
                          alpha=.3)
        return
        
    def get_curve(self, Lspline_wz, 
                      origin = [0.,0.,0.],
                      is_3dcurve=False):
        """draw a 'surface' 
        defined by a Bspline curve
        
        DONE: accomodate true 3D curves (hacked in)
        """
        Lspline = Lspline_wz[0]
        if Lspline.curve.dim == 2:
            z       = Lspline_wz[1]
            vb = frame.add1D(Lspline.curve.r,
                             which_index=2, 
                             where=z)
            verts3D = frame.add1D(Lspline.curve.vertices, 
                                  which_index=2, 
                                  where=z)
            vbr = list(vb)
            
            #vbr.insert(0,np.asarray([0.,0.,z]))
            #vbr.append(np.asarray([0.,0.,z]))
            
            vbr.insert(0,np.asarray([Lspline.curve.vertices[0,0],
                                     Lspline.curve.vertices[-1,1],
                                     z]))
                                    
            vbr.append(np.asarray([Lspline.curve.vertices[0,0],
                                   Lspline.curve.vertices[-1,1],
                                   z]))
            
            
            vbr         = np.asarray(vbr)
            one_face    = vbr
        elif Lspline.curve.dim==3:
            vb          = Lspline.curve.r
            verts3D     = Lspline.curve.vertices
            vbr  = list(vb)
            vbr.insert(0,np.asarray([Lspline.curve.vertices[0,0],
                                    Lspline.curve.vertices[-1][1],
                                    Lspline.curve.vertices[0,2]]))
                                    
            vbr.append(np.asarray([Lspline.curve.vertices[0,0],
                                    Lspline.curve.vertices[-1,1],
                                    Lspline.curve.vertices[0,2]]))
            
            vbr         = np.asarray(vbr)
            one_face    = vbr
        return one_face, verts3D
    
    def get_bow(self, Lspline_wz, 
                    origin = [0.,0.,0.],is_3dcurve=False):
        """draw a 'surface' 
        defined by a Bspline curve
        
        TODO: accomodate true 3D curves (hacked in)
        """
        Lspline = Lspline_wz[0]
        if Lspline.curve.dim == 2:
            z           = Lspline_wz[1]
            vb = frame.rotate_to_longitudinal(
                            Lspline.curve.r, 
                            where=z)
            verts3D = frame.rotate_to_longitudinal(
                            Lspline.curve.vertices, 
                            where=z)
            vbr = list(vb)
            vbr.insert(0,np.asarray([0.,0.,z]))
            vbr.append(np.asarray([0.,0.,z]))
            vbr = np.asarray(vbr)
            one_face = vbr
        elif Lspline.curve.dim==3:
            vb      = Lspline.curve.r 
            verts3D = Lspline.curve.vertices 
            vbr = list(vb)
            vbr.insert(0,np.asarray([Lspline.curve.vertices[0,0],
                                    Lspline.curve.vertices[-1][1],
                                    Lspline.curve.vertices[0,2]]))
                                    
            vbr.append(np.asarray([Lspline.curve.vertices[0,0],
                                    Lspline.curve.vertices[-1,1],
                                    Lspline.curve.vertices[0,2]]))
            vbr         = np.asarray(vbr)
            one_face    = vbr
        return one_face, verts3D
        
    def get_curvature(self, cvp,index):
        Lspline = self.hull_frame.sections[index][0]
        if Lspline.curve.dim ==2:
            z = self.hull_frame.sections[index][1]
        else:
            z = self.hull_frame.sections[index][1]
            #z = Lspline.curve.vertices[0,2]
        cvp = frame.add1D(cvp, which_index=2, where =z)
        return cvp
    
    
        

        
        
class DrawCurveInteractive(DrawCurve):
    """An Axes for displaying an Interactive 3D cube"""
    def __init__(self, 
                 fig, sections, 
                 xkcd, is_3d=False,
                 *args, **kwargs):
        super(DrawCurveInteractive, self).__init__(fig,
                                                    sections_=sections,
                                                    xkcd=False,
                                                    is_3d=False,
                                                    *args,
                                                    **kwargs)
        self.current_curve  = 0
        self.optimization_curve = self.current_curve 
        self._do_optimization = False
        self.pick_tol       = 5
        self.slow           = 15.
        self.slower         = 20.
        self._ind           = None
        self.active_poly    = None
        self._start_trans_x = 0.
        self._start_trans_y = 0.
        self._start_rot = Quaternion.from_v_theta((1, -1, 0),
                                                  -np.pi / 6)
                                                  
        # define axes for Up/Down motion and Left/Right motion
        self._ax_LR = (0, -1, 0) #Left Right
        self._ax_LR_alt = (0, 0, 1)
        self._step_LR = 0.01
        
        self._ax_UD = (1, 0, 0) #Up Down
        self._step_UD = 0.01

        self._active = False  # true when mouse is over axes
        self._button1 = False  # true when button 1 is pressed
        self._button2 = False  # true when button 2 is pressed
        self._button3 = False  # true when button 2 is pressed
        self._event_xy = None  # store xy position of mouse event
        self._shift = False  # shift key pressed
        
        self._current_trans_x = self._start_trans_x
        self._current_trans_y = self._start_trans_y
        self._current_rot = self._start_rot  

        # connect some GUI events
        self.figure.canvas.mpl_connect('key_press_event',
                                       self._key_press)
        #self.figure.canvas.mpl_connect('pick_event', 
        #                               self.on_pick_which)
        self.figure.canvas.mpl_connect('button_press_event',
                                       self._mouse_press)
        self.figure.canvas.mpl_connect('button_release_event',
                                       self._mouse_release)
        self.figure.canvas.mpl_connect('motion_notify_event',
                                       self._mouse_motion)
        
        self.figure.canvas.mpl_connect('key_release_event',
                                       self._key_release)
        
        self.zoom = self.zoom_factory(base_scale=3.0)
        #self.figure.canvas.mpl_connect('scroll_event',self.zoom_factory())
        
        #self.figure.text(0.05, 0.15, ("Click and Drag to Move\n"
        #                            "Hold shift key to adjust rotation, etc."))
        
        #self.figure.text(0.05, 0.05, ("Right click to zoom\n"
        #                            "shift right click to pan"))
        self.figure.text(0.05, 0.05, ("UNO Naval Architecture"))
    def rotate(self, rotation):
        self._current_rot = self._current_rot * rotation
        
    def translate(self, dx, dy):
        self._current_trans_x = self._current_trans_x + dx
        self._current_trans_y = self._current_trans_y + dy
        pass
        

    #----------------------------------------------------------
    # when the shift button is down, change the rotation axis
    # of left-right movement
    def _key_press(self, event):
        """Handler for key press events"""
        if event.key == 'shift':
            self._shift = True
        elif event.key == 'c':
            if self._show_curvature:
                print 'hide curvature'
                self._show_curvature = False
                self.draw_curves()
            else:
                print 'show curvature'
                self._show_curvature = True
                index = self.current_curve
                cvp = self.hull_frame.sections[index][0].curve.compute_curve_of_curvature()
                z = self.hull_frame.sections[index][1]
                self.cvp = frame.add1D(cvp, which_index=2, where=z)
                self.draw_curves()
        elif event.key == 'o':
            self._do_optimization = True
            #if self._optimize:
            #    self._optimize = False
            #else:
            #    self._optimize = True
            i = self.current_curve
            self.optimization_curve = i
            self.ocurve = gs.setup_optimization(self.hull_frame.sections[i][0], 
                                   optimization_listener = self.update_from_optimization,
                                   mouse_update_listener = self.update_from_manual_vertex_movement)
#            setconstraints = plt.axes([0.50, 0.05, 0.1, 0.075])
#            sconstraints = Button(setconstraints, 'SET')
#            sconstraints.on_clicked(self.ocurve.set_contraints)
#            
#            axconstraints = plt.axes([0.61, 0.05, 0.1, 0.075])
#            bconstraints = Button(axconstraints, 'Form Parm')
#            bconstraints.on_clicked(self.ocurve.add_constraints)
#            
#            axprev = plt.axes([0.72, 0.05, 0.1, 0.075])
#            bprev = Button(axprev, 'Optimize')
#            bprev.on_clicked(self.ocurve.do_optimization)
#            
#            axnext = plt.axes([0.83, 0.05, 0.1, 0.075])
#            bnext = Button(axnext, 'plot')
#            bnext.on_clicked(self.ocurve.plot_contraints)
            
    
            
    def _key_release(self, event):
        """Handler for key release event"""
        if event.key == 'shift':
            self._shift = False
            
    def get_ind_under_point(self, event):
        """get the index of the vertex under point 
        if within tolerance pick_tol
        """
        controlpolys = self._controlpoly
        ind_list = []
        dlist = []
        which_poly = []
        for i,cp in enumerate(controlpolys):
            xyt = cp.get_transform().transform(cp.xy)
            xt, yt = xyt[:, 0], xyt[:, 1]
            d = np.sqrt((xt-event.x)**2 + (yt-event.y)**2)
            indseq = np.nonzero(np.equal(d, np.amin(d)))[0]
            ind_list.append(indseq[0])
            dlist.append(d)
            which_poly.append(i)
        
        ind         = None
        the_poly    = None
        for test, d, ply in zip(ind_list,dlist,which_poly):
            if d[test]>=self.pick_tol:
                pass
            else:
                ind = test
                the_poly = ply
        return ind, the_poly
        
    
    def update_from_optimization(self, verts):
        i = self.optimization_curve
        self.hull_frame.sections[i][0].curve.vertices = verts
        self.hull_frame.sections[i][0].curve.allCurveArray()
        self.update_face_and_vertices(index=i)
        self.draw_curves()
        return
        
    def update_from_manual_vertex_movement(self, x,y,ind):
        i = self.optimization_curve
        self.hull_frame.sections[i][0].curve.vertices[ind,0] = x
        self.hull_frame.sections[i][0].curve.vertices[ind,1] = y
        self.hull_frame.sections[i][0].curve.allCurveArray()
        self.update_face_and_vertices(index=i)
        self.draw_curves()
        return
    
    def update_curve(self, x,y):
        """Use the rotation matrix to ensure
        that things behave intuitively
        """
        i = self.active_poly
        ind = self._ind
        Rs = (self._current_rot * self.rots).as_rotation_matrix()
        xy = np.dot(Rs,[x,y,0.])
        bow_id = self.hull_frame.bow_id
        #if i == bow_id:
        #    self.hull_frame.sections[i][0].curve.vertices[ind,0] += xy[0]/self.slow
        #    self.hull_frame.sections[i][0].curve.vertices[ind,2] += xy[1]/self.slow
        #else:
        self.hull_frame.sections[i][0].curve.vertices[ind,0] += xy[0]/self.slow
        self.hull_frame.sections[i][0].curve.vertices[ind,1] += xy[1]/self.slow
        #end if here
        self.hull_frame.sections[i][0].curve.allCurveArray()
        if self._do_optimization:
            self.current_curve = i
        if self._show_curvature:
            self.current_curve = i
            if i == bow_id:
                cvp = self.hull_frame.sections[i][0].curve.compute_curve_of_curvature()
                z = self.hull_frame.sections[i][1]
                self.cvp = frame.rotate_to_longitudinal(cvp, where=z)
            else:
                cvp = self.hull_frame.sections[i][0].curve.compute_curve_of_curvature()
                z = self.hull_frame.sections[i][1]
                self.cvp = frame.add1D(cvp, which_index=2, where=z)
        self.update_face_and_vertices(index=i)
        return

    def _mouse_press(self, event):
        """Handler for mouse button press"""
        self._event_xy = (event.x, event.y)
        self._ind, self.active_poly = self.get_ind_under_point(event)
        
        if self._ind is not None:
            print '\ncurve = {}'.format(self.active_poly)
            print 'vertex = {}'.format(self._ind)
            if self._show_curvature:
                self.current_curve = self.active_poly
        if event.button == 1: #rotating plot, xor moving vertices
            self._button1 = True
        elif event.button == 2: #mouse wheel.  not doing anything 
            self._button3 = True
        elif event.button == 3: #right click: zoom or translation  
            self._button2 = True

    def _mouse_release(self, event):
        """Handler for mouse button release"""
        self._ind = None
        self.active_poly = None
        self._event_xy = None
        if event.button == 1:
            self._button1 = False
        elif event.button == 2:
            self._button3 = False
        elif event.button == 3:
            self._button2 = False

    def _mouse_motion(self, event):
        """Handler for mouse motion"""
        if self._button1 or self._button2:
            dx = event.x - self._event_xy[0]
            dy = event.y - self._event_xy[1]
            self._event_xy = (event.x, event.y)
            if self._button1: #rotate or move vertex
                if self._ind is not None:
                    poly = self.active_poly
                    self.update_curve(dx,dy)
                    self.draw_curves()
                    return
                if self._shift:
                    ax_LR = self._ax_LR_alt
                else:
                    ax_LR = self._ax_LR
                rot1 = Quaternion.from_v_theta(self._ax_UD,
                                               self._step_UD * dy)
                rot2 = Quaternion.from_v_theta(ax_LR,
                                               self._step_LR * dx)
                self.rotate(rot1 * rot2)
                self.draw_curves()

            if self._button2:
                if self._shift: #translate
                    self.translate(dx/self.slower,dy/self.slower)
                    self.draw_curves()
                else: #zoom
                    factor = 1 - 0.003 * (dx + dy)
                    xlim = self.get_xlim()
                    ylim = self.get_ylim()
                    self.set_xlim(factor * xlim[0], factor * xlim[1])
                    self.set_ylim(factor * ylim[0], factor * ylim[1])
    
                    self.figure.canvas.draw()
#            if self._button3:
#                factor = 1 - 0.003 * (dx + dy)
#                xlim = self.get_xlim()
#                ylim = self.get_ylim()
#                self.set_xlim(factor * xlim[0], factor * xlim[1])
#                self.set_ylim(factor * ylim[0], factor * ylim[1])
#
#                self.figure.canvas.draw()
    ##
    ##quick, zoom!
    ##http://stackoverflow.com/questions/11551049/matplotlib-plot-zooming-with-scroll-wheel
    ##
    def zoom_factory_internet(ax,base_scale = 2.):
        def zoom_fun(event):
            # get the current x and y limits
            cur_xlim = ax.get_xlim()
            cur_ylim = ax.get_ylim()
            cur_xrange = (cur_xlim[1] - cur_xlim[0])*.5
            cur_yrange = (cur_ylim[1] - cur_ylim[0])*.5
            xdata = event.xdata # get event x location
            ydata = event.ydata # get event y location
            if event.button == 'up':
                # deal with zoom in
                scale_factor = 1/base_scale
            elif event.button == 'down':
                # deal with zoom out
                scale_factor = base_scale
            else:
                # deal with something that should never happen
                scale_factor = 1
                print event.button
            # set new limits
            ax.set_xlim([xdata - cur_xrange*scale_factor,
                         xdata + cur_xrange*scale_factor])
            ax.set_ylim([ydata - cur_yrange*scale_factor,
                         ydata + cur_yrange*scale_factor])
            #plt.draw() # force re-draw
            ax.figure.canvas.draw() #TLM draw call
    
        fig = ax.get_figure() # get the figure of interest
        # attach the call back
        fig.canvas.mpl_connect('scroll_event',zoom_fun)
    
        #return the function
        return zoom_fun
        
    def zoom_factory(self,base_scale = 2.):
        """TLM version
        """
        def zoom_fun(event):
            dx = event.x #- self._event_xy[0]
            dy = event.y #- self._event_xy[1]
            self._event_xy = (event.x, event.y)
            factor = 10.*(1 - 0.003 * (dx + dy))
            if event.button == 'up':
                scale_factor = 1./factor
            elif event.button == 'down':
                scale_factor = factor
            else:
                scale_factor = 1.
                print event.button
            
            xlim = self.get_xlim()
            ylim = self.get_ylim()
            self.set_xlim(scale_factor * xlim[0], scale_factor * xlim[1])
            self.set_ylim(scale_factor * ylim[0], scale_factor * ylim[1])

            self.figure.canvas.draw()

    
        #fig = ax.get_figure() # get the figure of interest
        # attach the call back
        self.figure.canvas.mpl_connect('scroll_event',zoom_fun)
    
        #return the function
        return zoom_fun
##
##
##
##---Back to my stuff
##
##  2D to 3D quick and inflexibly
##           
class frame(object):
    def __init__(self):
        return
    @classmethod
    def add1D(cls, array, which_index=2, where = 0.):
        size,dim = np.shape(array) 
        new_array = np.zeros((size,dim+1),float)
        if which_index == 0:
            new_array[:,0]  = where
            new_array[:,1:] = array[:]
        if which_index == 1:
            new_array[:,0]  = array[:,0]
            new_array[:,1]  = where
            new_array[:,2]  = array[:,1]
        if which_index == 2:
            new_array[:,0]  = array[:,0]
            new_array[:,1]  = array[:,1]
            new_array[:,2]  = where
        return new_array
    @classmethod
    def rotate_to_longitudinal(cls, array, where = 0.):
        size,dim = np.shape(array) 
        new_array = np.zeros((size,dim+1),float)
        new_array[:,0]  = 0.
        new_array[:,1]  = array[:,1]
        new_array[:,2]  = where + array[:,0]
        return new_array
    
    @classmethod #Cp
    def xy_to_zyx_elevate(cls, 
                         array, 
                         elevation = 0.):
        size,dim = np.shape(array) 
        new_array = np.zeros((size,dim+1),float)
        new_array[:,0]  = 0.
        new_array[:,1]  = elevation - array[:,1][::-1]
        new_array[:,2]  = array[:,0]
        return new_array
    
    @classmethod #Cp
    def xy_to_zyx_elevate_noflp(cls, 
                         array, 
                         elevation = 0.):
        size,dim = np.shape(array) 
        new_array = np.zeros((size,dim+1),float)
        new_array[:,0]  = 0.
        new_array[:,1]  = elevation - array[:,1]
        new_array[:,2]  = array[:,0]
        return new_array
        
    @classmethod #DWL
    def xy_to_zxy_elevate(cls, 
                         array, 
                         elevation = 0.):
        size,dim = np.shape(array) 
        new_array = np.zeros((size,dim+1),float)
        new_array[:,0]  = array[:,1]
        new_array[:,1]  = elevation
        new_array[:,2]  = array[:,0]
        return new_array
        
    @classmethod
    def translate(cls, *args):
        """redo this to simplify add1D
        """
        return curve
        
        
##
## Quick stations
##
class hull_frame(object):
    """TLM Nov 1 2017
    Longitudinal axis is set to z?
    """
    def __init__(self):
        self.station_id = 0
        self.bow_id = None
        self.sections = {}
        return
        
    def __call__(self,index):
        return self.sections[index]
        
    def add_section(self, Lspline, z):
        self.sections[self.station_id] = [Lspline,z]
        self.station_id +=1
        return
    
    def add_bow(self, Lspline, z):
        self.sections[self.station_id] = [Lspline,z]
        self.bow_id = self.station_id
        self.station_id +=1
        return
    
    def add_3Dsection(self, Lspline):
        z = Lspline.curve.vertices[0,2]
        self.sections[self.station_id] = [Lspline,z]
        self.station_id +=1
        return
    
    def add_3Dbow(self, Lspline):
        z = Lspline.curve.vertices[0,2]
        self.sections[self.station_id] = [Lspline,z]
        self.bow_id = self.station_id
        self.station_id +=1
        return
        

def rotate_curve():
    rt=np.sqrt(2.)/2.
    _current_rot = Quaternion.from_v_theta((rt, 0., rt), 0.)
    
    return

def initial_curve(start=(0.,12.), end=(12.,0.), num=7, k=4,nump=50):
    vertices = linear_vertices(start, end, num)
    k=4
    nump=50
    curve = spline.Bspline(vertices, k, nump)
    return curve    

def linear_vertices(start, end, num):
    start = np.asarray(start)
    end = np.asarray(end)
    D = end - start
    dim = len(D)
    vertices = []
    for i in range(dim):
        vi = np.linspace(start[i],end[i],num)
        vertices.append(vi)
    vertices = np.asarray(vertices)
    vertices = np.transpose(vertices)
    return vertices  
    
if __name__ =='__main__':
    
    
    if True: #old 2D-fake-3D test
        test    = 'SAC_3' 
        HF = hull_frame()
        
        curve = initial_curve((0.,0.,0.),(36.,0.,0.), num=12, k=4,nump=30)
        v = copy.deepcopy(curve.vertices)
        v[4:7,1] = 15.
        
        curve.vertices = v #updates the curve automagically
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = 288.)   
        FPD.add_XcConstraint(kind='equality', value=18.)
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = -0.)
        FPD.add_yPointConstraint(kind='equality', location = 0.333, value = 12.)
        FPD.add_yPointConstraint(kind='equality', location = 0.666, value = 12.)
        FPD.add_E1(kind='LS', weight = .1)
        FPD.add_E2(kind='LS', weight = .5)
        #FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize()
        save1 = copy.deepcopy(Lspline.curve)
        
        curve = Lspline.curve
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = 288.)   
        FPD.add_XcConstraint(kind='equality', value=18.)
        FPD.add_AngleConstraint(kind='equality', location = 0., value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = 1., value = -0.)
        FPD.add_yPointConstraint(kind='equality', location = 0.333, value = 12.)
        FPD.add_yPointConstraint(kind='equality', location = 0.666, value = 12.)
        FPD.add_AngleConstraint(kind='equality', location = .333, value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = .666, value = 0.)
        #FPD.add_CurvatureConstraint(kind='equality', location = 0., value = 0.)
        #FPD.add_CurvatureConstraint(kind='equality', location = 1., value = 0.)
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize()
        
        curve = Lspline.curve
        interval_data, small = interval_bounds(curve)
        FPD = FormParameterDict(curve) 
        FPD.add_AreaConstraint(kind='equality', value = 288.)   
        FPD.add_XcConstraint(kind='equality', value=18.)
        #FPD.add_AngleConstraint(kind='equality', location = 0., value = 0.)
        #FPD.add_AngleConstraint(kind='equality', location = 1., value = -0.)
        FPD.add_yPointConstraint(kind='equality', location = 0.333, value = 12.)
        FPD.add_yPointConstraint(kind='equality', location = 0.666, value = 12.)
        FPD.add_AngleConstraint(kind='equality', location = .333, value = 0.)
        FPD.add_AngleConstraint(kind='equality', location = .666, value = 0.)
        FPD.add_yPointConstraint(kind='equality', location = .5, value = 12.)
        #FPD.add_AngleConstraint(kind='equality', location = .5, value = 0.)
        #FPD.add_CurvatureConstraint(kind='equality', location = .5, value = 0.)
        #FPD.add_CurvatureConstraint(kind='equality', location = .333, value = 0.)
        #FPD.add_CurvatureConstraint(kind='equality', location = .666, value = 0.)
        FPD.add_E1(kind='LS', weight = .5)
        FPD.add_E2(kind='LS', weight = .5)
        FPD.add_E3(kind='LS', weight = .5)
        FPD.add_ArcLengthApprox(kind='LS', weight = .5)
        L = Lagrangian(FPD)
        interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
        Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
        Lspline.optimize()
    
        SAC_curve = copy.deepcopy(Lspline)
        #----------------------------------------------------------------------
        #----------------------------------------------------------------------
        #  Sections...
        #
        #
        if True:
            
            
            #"""
            k=4
            nump=30
            xb = 0.
            yb = 5.#12.
            xe = 5.#12.
            ye = 0.
            
            alphae = 0.
            alphab = 0.
            Cab_given = 0.
            Cae_given = 0
            xc = 2.#4.
            yc = 2.#4.
            curve_area = 12.#72.
            slope = 'down'
            
            ini_v = InitializeControlVertices(alphae = 0., alphab = 0.,Cab_given = 0.,Cae_given =0.)
                
            curve = spline.Bspline(ini_v.vertices, k, nump)  
            ae = alphae
            ab = alphab
            
            ini_v = InitializeControlVertices(xb=0.,yb=5.,xe=5.,ye=0.,
                                              alphae=ae, alphab=ab,
                                              Cab_given=0.,Cae_given=0.,
                                                                  nCV = 7)
            curve = spline.Bspline(ini_v.vertices, k, nump)  
            wi = 3.5 #2.5#
            w  = .5 #.65
            ep = 1.e-10
            sp = 1.e-4
            x = ini_v.vertices[:,0]
            y = ini_v.vertices[:,1]
            interval_data, small = interval_bounds(curve)
            
            #if test == 'basic': 
            FPD = FormParameterDict(curve) 
            FPD.add_AreaConstraint(kind='equality', value = curve_area)#value = AF.fuzzyNumber(60.,70.,80.) )#
            FPD.add_AngleConstraint(kind='equality', location = 0., value = alphab)#AF.fuzzyNumber(-5.,-2.,0.))#
            FPD.add_AngleConstraint(kind='equality', location = 1., value = alphae)
            FPD.add_CurvatureConstraint(kind='equality', location = 0., value = Cab_given)
            FPD.add_CurvatureConstraint(kind='equality', location = 1., value = Cae_given)
            #"""
            FPD.add_E1(kind='LS', weight = 1.)
            FPD.add_E2(kind='LS', weight = 1.)
            FPD.add_E3(kind='LS', weight = 1.)
            #FPD.add_ArcLength(kind='LS', weight = 1.)
            FPD.add_ArcLengthApprox(kind='LS', weight = 1.)
            L = Lagrangian(FPD)
            
            interval_data, small = lagrangian_bounds(L, interval_data, small, 1.e4)
            
            Lspline = IntervalLagrangeSpline(curve, L, data = interval_data)
            Lspline.curve.verbose = True
            Lspline.optimize()
            self = Lspline
            
            #HF.add_section(Lspline, z=0.)
            #HF.add_section(copy.deepcopy(Lspline), z=10.)
            #HF.add_section(copy.deepcopy(Lspline), z=20.)
            #HF.add_bow(copy.deepcopy(Lspline), z=30.)
            
            HF.add_section(copy.deepcopy(Lspline), z=18.)
            HF.add_bow(SAC_curve, z=0.)
            sections = HF.sections
            key=0
            
        if True:
            fig = plt.figure(figsize=(4, 4))
            ax = DrawCurveInteractive(fig, sections=HF, xkcd=False)
            fig.add_axes(ax)
            ax.draw_curves()
            self = ax
            #display(fig)
            
            if False: #test quaternion rotors
                #rotation matrix:
                Rs = (self._current_rot * self.rots).as_rotation_matrix()
                #translation vector:
                Ts = np.asarray([self._current_trans_x,self._current_trans_y,0.])
                # rotate and translate the faces
                faces = [np.dot(face,Rs) + Ts for face in self.real_faces]
                #rotate and translate vertices 
                verts = [np.dot(vt,Rs) + Ts for vt in self.vertices]
                
                pt = faces[0][10]
                q1 = self.rots
                q2 = self._current_rot
                Qrotor = self._current_rot * self.rots
                #Qrotor =  self.rots *self._current_rot 
                print 'rotation method 1: Rotation Matrix'
                print np.dot(pt,Rs)
                print 'rotation method 2: Quaternions'
                a=Qrotor*pt*Qrotor.conjugate()
                b=Qrotor.conjugate()*pt*Qrotor #matches Rs - in w,(x,y,z) form.  Note as_v_theta does not quite!
                c = self._current_rot * self.rots * pt  * self.rots.conjugate() *self._current_rot.conjugate()
                print a
                print a.as_v_theta()[0]
                print b.as_v_theta()[0]
                #order matters!
                b1 = Qrotor.conjugate()
                b2 = Qrotor*pt
                b3 = b2*b1
                #http://stackoverflow.com/questions/4870393/rotating-coordinate-system-via-a-quaternion
                x_axis_unit = (1, 0, 0)
                y_axis_unit = (0, 1, 0)
                z_axis_unit = (0, 0, 1)
                theta = np.pi/2.
                r1 = Quaternion.from_v_theta(x_axis_unit, theta)
                r2 = Quaternion.from_v_theta(y_axis_unit, theta)
                r3 = Quaternion.from_v_theta(z_axis_unit, theta)
                
                v = r1*np.asarray(y_axis_unit)*r2.conjugate() * r3.conjugate()
                v.as_v_theta() 
                
                v = r3*r2*r1*np.asarray(y_axis_unit)*r1.conjugate()*r2.conjugate() * r3.conjugate()
                print v.as_v_theta() 
                # output: (4.930380657631324e-32, 2.220446049250313e-16, -1.0), theta
                
                
                v = r3*r2*r1*np.asarray(x_axis_unit)*r1.conjugate()*r2.conjugate() * r3.conjugate()
                print v.as_v_theta()
                # output: (4.930380657631324e-32, 2.220446049250313e-16, -1.0)
    
        
        print 'Area of section: {}'.format(Lspline.curve.area)
        print 'Value of Objective Function: {}'.format(Lspline.f)
        
        """
        Objective function
        
        Lspline.f
        
        
        Area:
        
        
        Lspline.Lagrangian.equality[0].type   #area
        Lspline.Lagrangian.equality[0].kind #equality
        
        Lspline.Lagrangian.equality[0].computed_value  
        Lspline.Lagrangian.equality[0].computed_value.grad.T
        Lspline.Lagrangian.equality[0].computed_value.hess
        
        """
        """
        How to rotate a curve back in 'real life':
            
        vertices=np.asarray([   [-4.,0.,0.], [-2.,0.,0.],[0.,0.,0.],
                                    [2.,0.,0.],
                                    [3.0,0.,0.],
                                    [6.0,5.,0.],
                                    [9.0,12.,0.],
                                    [11.,12.,0.],
                                    [12.,12.,0.],[14.,12.,0.],[16.,12.,0.]])
        v0 = spline.Bspline(vertices, k, nump)
        
        Rs = (self._current_rot * self.rots).as_rotation_matrix()
     
        # (method 1)
        tv1 = np.zeros((v0.n,3),float)
        for i, vtx in enumerate(v0.vertices):
            tv1[i] = np.dot(vtx,Rs)
            
        
        # (method 2)
        tv2 = np.dot(v0.vertices[:],Rs)
        
        """