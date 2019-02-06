# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 19:05:51 2016

@author: lukemcculloch
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons


class slider_factory(object):
    """
        Accepts: my overloaded subplot from adials_gui
        returns: a slider for a relational interval
    """
    def __init__(self, ax, a=.1,b=10.):
        self.ax = ax
        self.fig = self.ax.figure
        self.plt = plt
        
        t = np.arange(0.0, 1.0, 0.001)
        self.t = t
        a0 = 5
        f0 = 3
        self.a0 = a0
        self.f0 = f0
        s = a0*np.sin(2*np.pi*f0*t)
        self.l, = plt.plot(t, s, lw=2, color='red')
        plt.axis([0, 1, -10, 10])

        self.axcolor = 'lightgoldenrodyellow'
        axcolor = self.axcolor
        self.axfreq = plt.axes([0.25, 0.1, 0.65, 0.03], axisbg=axcolor)
        self.axamp = plt.axes([0.25, 0.15, 0.65, 0.03], axisbg=axcolor)
        
        self.sfreq = Slider(self.axfreq, 'Freq', 0.1, 30.0, valinit=f0)
        self.samp = Slider(self.axamp, 'Amp', a,b, valinit=a0)
        
        
        self.resetax = self.plt.axes([0.8, 0.025, 0.1, 0.04])
        self.button = Button(self.resetax, 'Reset', color=self.axcolor, hovercolor='0.975')
                
        self.rax = plt.axes([0.025, 0.5, 0.15, 0.15], axisbg=axcolor)
        self.radio = RadioButtons(self.rax, ('red', 'blue', 'green'), active=0)
        self.radio.on_clicked(self.colorfunc)
        
        self.update()
        return


    def update_(self, val):
        t=self.t
        amp = self.samp.val
        freq = self.sfreq.val
        self.l.set_ydata(amp*np.sin(2*np.pi*freq*t))
        self.fig.canvas.draw_idle()
        return
        
    def update_freq(self, values):
        """
            TODO: implement
        """
        a = values[0]
        b = values[1]
        
        self.sfreq = Slider(self.axfreq, 'Freq', a,b, 
                            valinit=self.f0)
        self.update()
        self.fig.canvas.draw_idle()
        return
        
    def update_amp(self, values):
        a = values[0]
        b = values[1]
        self.fig.clear()
        self = slider_factory(ax,a,b)
        return
        
    def update(self):
        self.sfreq.on_changed(self.update_)
        self.samp.on_changed(self.update_)
        return


    def reset(self, event):
        self.sfreq.reset()
        self.samp.reset()
    
    def colorfunc(self, label):
        self.l.set_color(label)
        self.fig.canvas.draw_idle()
        return





fig, ax = plt.subplots()
self = slider_factory(ax)

#self.button.on_clicked(self.reset())
plt.show()

#self.update_amp([3.,6.])