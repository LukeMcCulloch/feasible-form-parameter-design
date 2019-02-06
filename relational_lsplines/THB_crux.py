#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 18 19:42:16 2017

@author: luke

Not used!!
"""

def THBasis(self,u,v,ucurve,vcurve):
        
        dummy_bounds = ia(0.,0.)
        ubasis_levels, ulcurve = ucurve.recursive_to_bottom(u)
        vbasis_levels, vlcurve = ucurve.recursive_to_bottom(u)
        ulevel_mats = []
        vlevel_mats = []
        #TODO in THBsurface: make ulcurve.parent make sense.
        #note that it is only needed for basis, -whew-
        while ulcurve.parent is not None: #TODO: make this work for surface curves!
            #child = ulcurve
            ulcurve = ulcurve.parent
            vlcurve = vlcurve.parent
            
            uiBasis =  ubasis_levels[-1]
            uactive_parents, uapb = ulcurve.get_one_level(u)
            uoBasis = ulcurve.TLMBasisFuns(u)
            uchild = ulcurve.children
            
            
            
            viBasis =  vbasis_levels[-1]
            vactive_parents, vapb = vlcurve.get_one_level(u)
            voBasis = vlcurve.TLMBasisFuns(u)
            vchild = vlcurve.children
            
            
            
            uactive_child_bounds_index = uchild.bounds.contains(u)
            ulen_acbi = len(uactive_child_bounds_index)
            assert(len(uactive_child_bounds_index)<2),'error: child bounds active in two spans at once'
            if ulen_acbi==1:
                uchild_bounds = uchild.bounds[uactive_child_bounds_index[0]]
                uinner = uchild.enclosed_basis_in_range(uchild_bounds)
                uouter = uchild.parent.enclosed_basis_in_range(uchild_bounds) #to be eliminated
                
            else:
                uchild_bounds = dummy_bounds
                uinner = []
                uouter = [] #none to be eliminated
                # if none are encl
                #apb=[]
                
            
            
            
            vactive_child_bounds_index = vchild.bounds.contains(u)
            vlen_acbi = len(vactive_child_bounds_index)
            assert(len(vactive_child_bounds_index)<2),'error: child bounds active in two spans at once'
            if vlen_acbi==1:
                vchild_bounds = vchild.bounds[vactive_child_bounds_index[0]]
                vinner = vchild.enclosed_basis_in_range(vchild_bounds)
                vouter = vchild.parent.enclosed_basis_in_range(vchild_bounds) #to be eliminated
                
            else:
                vchild_bounds = dummy_bounds
                vinner = []
                vouter = []
            
        
            if ulcurve.rm is None:
                ulcurve.rm = ulcurve.wgen(ulcurve.p,ulcurve.level+1)#
        
            if vlcurve.rm is None:
                vlcurve.rm = vlcurve.wgen(vlcurve.p,vlcurve.level+1)#
            
            
            
            ulevel_mats.append(ulcurve.rm)
            cl = len(ubasis_levels)
            matrix_reducer = np.matmul(ulevel_mats[-1].T,uiBasis)
            for i in range(cl-1):
                cBasis = ubasis_levels[i]
                local_reducer = np.matmul(ulevel_mats[i].T,cBasis)
                for j in range(i+1,cl):
                    local_reducer = np.matmul(ulevel_mats[j].T,local_reducer)
                matrix_reducer = matrix_reducer + local_reducer
            uoBasis[uapb] = uoBasis[uapb] - matrix_reducer[uapb]
            
            
            uoBasis[uouter] = 0.
            
            uoBasis = ulcurve.clip_basis_outside_refinement(uoBasis, 
                                                    uactive_parents)
            ubasis_levels.append(uoBasis)
            
            
            
            vlevel_mats.append(vlcurve.rm)
            cl = len(vbasis_levels)
            matrix_reducer = np.matmul(vlevel_mats[-1].T,viBasis)
            for i in range(cl-1):
                cBasis = vbasis_levels[i]
                local_reducer = np.matmul(vlevel_mats[i].T,cBasis)
                for j in range(i+1,cl):
                    #transBasis = basis_levels[j]
                    local_reducer = np.matmul(vlevel_mats[j].T,local_reducer)
                matrix_reducer = matrix_reducer + local_reducer
            voBasis[vapb] = voBasis[vapb] - matrix_reducer[vapb]
            
            
            voBasis[vouter] = 0.
            
            voBasis = vlcurve.clip_basis_outside_refinement(voBasis, 
                                                    vactive_parents)
            vbasis_levels.append(voBasis)
            
            
            
        return ubasis_levels, vbasis_levels