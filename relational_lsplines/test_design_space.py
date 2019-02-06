# -*- coding: utf-8 -*-
"""
Created on Sun Jan  1 21:47:38 2017

@author: lukemcculloch
"""
import numpy as np
from design_space import HullSpace
import sqKanren as lp
from extended_interval_arithmetic import ia

from hull_inference_ob_graph import Hull as hullclp
from interval_arithmetic import ia as oia


if __name__ == '__main__':
    self = HullSpace()
    self.ini_coeff()
    self.dsmax = ia(22.,22.)
    self.Cb = ia(.95,.95)
    self.vol = ia(2000.,300000.)
    self.lwl = ia(120.,120.)
    self.bwl = ia(20.,20.) 
    self.draft = ia(20.,20.) 
    #self.AC_revise()
    self.Cwp = ia(.3,1.)
    self.Ccp = ia(0.,1.)
    #self.Cp = ia(.96,.99)
    #self.Cp = ia(.96,.96)
    #self.Cwp = ia(.7,.7)
    #"""
    self.states = self.shape_bow_section_inizers()
    self.states = self.shape_stern_section_inizers()
    self.states = self.SAC_run_area_ini()
    self.states = self.SAC_Xc_init()
    self.AC_revise()
    #self.AC_revise(print_=True)
    #for key in self.
    #for key in self.keys:
    for key in self.Coefficients:
        val = self(key)[0].getpoint(.5)
        #val = ia(val-self.SMALL,val+self.SMALL)
        val = ia(val,val)
        self.__setattr__(key.name,val)
        print key, self(key)
        #self.AC_revise()
    for key in self.Areas:
        val = self(key)[0].getpoint(.5)
        val = ia(val-self.SMALL,val+self.SMALL)
        #val = ia(val,val)
        self.__setattr__(key.name,val)
        print key, self(key)
    #"""
    
    self.LCG = ia(60.,60.)
    self.SAC_fwd_Xc = ia(15.,15.)
    self.SAC_mid_Xc = ia(60.,60.)
                        
    
    """
    h1 = hullclp()
    h1.dsmax = oia(22.,22.)
    h1.Cb = oia(.95,.95)
    #h1.Cp = oia(0.,1.)
    h1.vol = oia(2000.,300000.)
    h1.lwl = oia(120.,120.)
    h1.bwl = oia(20.,20.) 
    h1.draft = oia(19.999,20.001) 
    #"""