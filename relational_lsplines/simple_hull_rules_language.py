#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 6 10:00:00 am 2018

@author: luke

    contains:
        Rules for the Bare Hull.
        
    purpose:
        demonstrate the langauge functionality
        by setting the rule system for describing 
        bare ship hulls in Form Parameter Design (FPD)
        
    TODO: re-write ad hoc meta programming
    so that the hull class will accept state changes
    the old fashioned way!
    
    
    DEV
    ----------
    
        #
        #    SD.hull.SAC.plot3DmultiList(
        #            SD.hull.lcurvenet[:1],
        #            [SD.hull.CProfile])
        #            
        #    
        
        if SD.hull.hullformtype == 'osv':
            thiscurve = SD.hull.DWL_aft
        else:
            thiscurve = SD.hull.CPK_aft
        
        SD.hull.CProfile.plot3DmultiList([SD.hull.CProfile,
                                 SD.hull.CPK_cdr_b],
                                 [SD.hull.bowfairning,
                                  SD.hull.bowtransition,
                                 SD.hull.CPK_fwd,
                                 SD.hull.sterntransition,
                                 SD.hull.sternfairing,
                                 thiscurve,])
        

Reasonable Prior?:
    
A policy that randomly stumbles onto good training examples 
will bootstrap itself much faster than a policy that doesn’t.
 
https://www.alexirpan.com/2018/02/14/rl-hard.html



An aspect of this difficulty involves building an intuition 
for what tool should be leveraged to solve a problem. 

http://ai.stanford.edu/~zayd/why-is-machine-learning-hard.html

"""
##
import opt_simple_hull as opt
#
ShipDesigner = opt.ShipDesigner#Design Helper class
# which abstracts out the various processes of designing a hull
# 'the random way'
#
DesignSpecification = opt.DesignSpecification #tiny class
# which maps from tuples to intervals
# the idea is to instantiate the design space with something familiar (tuple)
#
# personally, I prefer to use the tiny language itself as that's what
# captures the design space intent.
# And that's where the 'real machinary' is.
#
#
#
lp = opt.lp #everything in sqKanren.py
np = lp.np #numpy
#
ia = lp.ia #use ia import from lp instead to match isinstance
#
RulesGraphProcessor = lp.RulesGraphProcessor
#
import random
import sobol #actually using python's random.random at the moment.
# -this was only because I percieved sobol not to be 'random enough'
# -but after some experience with random.random, things seem similar.
import copy 
##
##
##

  
def ini_coeffs(obj):#clist, rgp):
    """initialize coefficients to ia(0.,1.)
    """
    print 'initializing coefficients'
    for co in obj.Coefficients:
        co   = co == ia(0.,1.)
    
        obj.rgp.add_one_rule(co,co.name)
    
        obj.rgp.compute_fresh_rules_graph()
    return obj


class ComplexConsistency(object):
    
    def __init__(self, x0=0.,y0=0.):
        self.x0 = x0
        self.y0 = y0
        return
    
    def est_xc_fwd_sac(self, design_space):
        """
            design as if this curve 
            were flipped -> max y at origin
            
            This gives tightest results.
        """
        x3      = design_space.get_value(
                'SAC_entrance_len')[0] #starting allowable range
        x3.inf = self.x0 
        
        area    = design_space.get_value(
                        'SAC_entrance_area')[0]#SAC_entrance_area
        xb = ia(self.x0,self.x0)
        xe = design_space.get_value(
                'SAC_entrance_len')
        x1 = self.vertices[1,0]
        x2 = self.vertices[2,0]
        x4 = self.vertices[-3,0]
        x5 = self.vertices[-2,0]
        xe = self.vertices[-1,1]
        #
        yb = self.vertices[0,1]
        y1 = self.vertices[1,1]
        y2 = self.vertices[2,1]
        y4 = self.vertices[-3,1]
        y5 = self.vertices[-2,1]
        ye = self.vertices[-1,1]

        #term1 = (1./(area*6.))
        term1 = (area**-1) * (1./6.)
        if isinstance(term1,list):
            term1 = term1[0]
        term2 = 2.*area*(x2+x4) + x1**2*yb - x1*x2*yb - x1*x4*yb
        term3 = x1*x2*y1 - x2*x4*y1 - x1**2*y2 + x1*x4*y2 - x2*x5*y4 + x5**2*y4 + xe**2*y5
        term4 = -xe*x2*y5 - xe*x4*y5 + x2*x4*y5 + xe*x5*y5 - x4*x5*y5
        termA = term1*(term2+term3+term4)
        
        #termB = (1./(6.*area))*(2.*area - x1*yb -x2*y1 + x1*y2 - x4*y2 + x2*y4 - x5*y4 - xe*y5 +x4*y5 )
        
        termB = (2.*area - x1*yb -x2*y1 + \
                     x1*y2 - x4*y2 + x2*y4 - \
                     x5*y4 - xe*y5 +x4*y5 )*term1
        
        xc = termA + termB*x3
        return xc


class HullGeometryGenerator(object):
    """Rules for the Feasible Design
    of a Bare Ship Hull
    
    
    bare hull coordinate system:
    
        y
        |
        |
        o -----x
       /  
      /
     z  
     
     x : transverse
     y : vertical 
     z : longitudinal 
     
     (3D graphics coordiantes... TODO: change this ;)
     
    Notes
    ----------
        New style language based rules and AC_revise
    """
    def __init__(self,rtype='gauss',
                 verbose=False):
        #
        self.rgp = lp.RulesGraphProcessor(
                                verbose=verbose)
        self.tol = 1.e-2#1.e-4
        self.sobol_seq = sobol.sobolSeq([1,1],[1,1])
        self.random_type = {'sobol':self.sobol_seq.next,
                            'gauss':random.random}
        self.verbose = verbose
        self.rtype = rtype
        #toying around with putting a _paver_ back in:
        # but the reason it's not in
        # is that we are only really testing
        # a space by seeing if it becomes empty
        # this is not paving
        # but it could still be used to assemble
        # a tree of design choices along the way
        self.feasible_designs = []
        self.infeasible_designs = []
        
    
    def __str__(self):
        print 'HullGeometryGenerator, self.rgp.env'
        print 'rules graph processor'
        print ' design space (self.rgp.env) is as follows:\n'
        return ''.format(self.rgp.print_state())
    
    def __call__(self,this):
        if isinstance(this, str):
            return self.get_value(this)
        else:
            print 'expects the name of one of the PStates attributes'
            return None
    
    def __setattr__(self, item, value):
        """http://code.activestate.com/recipes/389916-example-setattr-getattr-overloading/
        
        Maps attributes to values.
        Only if they are initialised
        
        This 'interesting' business is so that I can overload the setter
        of attributes that I do not know the names of ahead of time.
        (it's must nicer if I enforce the names ahead of time)
        (also the stupid hull_from_simple_designspace class
        required certain names anyway  -again, not composable, not linguistic-
        so that spoils the fun and we might as well be enforcing names everywhere)
        (but hey, we are on the road to real automation so let's keep our
        half measures for now)
        
        Note that the tame way to go about this is to require
        the use of the standard langaugae construct == to set 
        calss.lp.PStates attributes 
        to map to _whatever_ you set it to 
        in the rules graph processing states list
        
        To think about:
            -use eval and parse = ?  to replace __setattr__ ?
                *overkill
                *complicated
                *gets in the way of = at init => does it even work since this is the case??
                *it should, see sq kanren classes
                
            -drop setattr overloading?
                *loose ability to use self.attr = ia(low,high)
                
        -----------------------------------------------------------------------    
        annotated init hasattr 'helper'long verbose version:
        -----------------------------------------------------------------------
        #print 'type of item is = ',type(item), item
        if hasattr(self,item):
            node = self.__getattribute__(item)
            self.__eq__(node,value)
            #            else:
            #                #print 'PS 2'
            #                node = self.__getattribute__(item)
            #                self.make_attributes_from_rule(node)
            #                node = self.__getattribute__(item)
            #                self.__eq__(node,value)
        else:
            if self.__dict__.has_key(item):
                #print 'using 00'
                #print 'found item',item,value
                self.__eq__(item,value)
            elif not self.__dict__.has_key(item):  # this test allows attributes to be set in the __init__ method
                #print 'using 1, if cannot pass here, do class.eval'
                #try:
                #    self.make_attributes_from_rule(item)
                #    node = self.__getattribute__(item)
                #    self.__eq__(node,value)
                #except:
                #    return dict.__setattr__(self, item, value)
                return dict.__setattr__(self, item, value)
            elif self.__dict__.has_key(item):       # any normal attributes are handled normally
                print 'using 2'
                print 'found item',item,value
                self.__eq__(item,value)
            else:
                print '3'
                #self.make_attributes_from_rule(item)
                #self.__setattr__(item, value)
                self.__eq__(item,value)
        -----------------------------------------------------------------------
        -The idea is to add PStates variables to this class
        automagically.
        -This kind of thing is not strictly necessary but 
        keep cropping up in my designs for the use of sqKanren...
                
        """
    
        if hasattr(self,item):
            node = self.__getattribute__(item)
            self.__eq__(node,value)
            #            if isinstance(node, lp.PStates):
            #                self.__eq__(node,value)
            #            else:
            #                print 'use vanilla setattr'
            #                print 'setting ',node,' => ', value
            #                dict.__setattr__(self, item, value)
        else:
            if self.__dict__.has_key(item):
                self.__eq__(item,value)
            elif not self.__dict__.has_key(item):  # this test allows attributes to be set in the __init__ method
                #print 'using 1, if cannot pass here, do class.eval'
                return dict.__setattr__(self, item, value)
            elif self.__dict__.has_key(item):       # any normal attributes are handled normally
                print 'using 2'
                print 'found item',item,value
                self.__eq__(item,value)
            else:
                print '3'
                self.__eq__(item,value)
    
    def get_value(self, name):
        """
            put in this level of indirection
            because the hull is technically a list of sates
            it's not so wise to [0] everywhere!
            
            TODO: overload getattr to get the vars properly out 
            of the env?
        """
        if isinstance(name, lp.Variable):
            val = self.rgp.env(name)[0]
        elif isinstance(name, str):
            nm,val=self.rgp.get_name_and_value(name)
            val = self.rgp.env(nm)
        return val
    
    def build_list(self, varname, dlist=None):
        """Helper to build a searchable list 
        for constraint filtering via tree_search
        and wrapped search
        
        
        
        hullrulesnet.build_list('vol',dlist)
        hullrulesnet.build_list('lwl',dlist)
        
        
        hullrulesnet.build_list('bwl',dlist)
        
        hullrulesnet.build_list('draft',dlist)
        
        hullrulesnet.build_list('Ccp',dlist)
        
        hullrulesnet.build_list('Cb',dlist)
        
        hullrulesnet.build_list('Cp',dlist)
        
        hullrulesnet.build_list('Awp',dlist)
        
        print dlist
        Out[120]: [PS(lwl), PS(vol), PS(bwl), PS(draft), PS(Ccp), PS(Cb), PS(Cp), PS(Awp)]
        
        hullrulesnet.wrapped_search(dlist, maxinner=3)
        
        dlist = hullrulesnet.get_thick_list()
        hullrulesnet.wrapped_search(dlist, maxinner=3)
        """
        if dlist is None: dlist = []
        if isinstance(varname, lp.PStates):
            var = getattr(self, varname.name)
        else:
            assert(isinstance(varname, str) ),'Error variable is neither PStates nor string'
            var = getattr(self, varname)
        dlist.append(var)
        return dlist
    
    def __eq__(self, var,value,update=True):
        """Use = as a logic variable setter 
        for already ('registered'=existing in the States env) logic variables
        """
        print 'setting item',var
        print 'to value',value
        var = var == value
        self.rgp.add_one_rule(var,var.name)
        if update:
            self.rgp.compute_fresh_rules_graph()
        return

    def get_parameter(self, name):
        """Not used so far
        """
        param = self.rgp(name)
        return param
    
    
    def get_from_rgp(self, name):
        """
        
        dev
        ----------
        -its odd that this works,
        yet dropping the treednode straight into the env
         does nothing
        yet construct returns the right answers
         supposedly without adding a new varialb to env
         
         
        cur_el = self.rgp.env.bind(self.BulbBeam.name)
        >>> self.rgp.env(cur_el)
        [ia(5.0, 8.0)]
        """
        mk_name, mk_value = self.rgp.get_name_and_value(name)
        return mk_name, mk_value
    
    
    def get_for_design(self, var):
        """return a thin value for each FPD parameter
        assumes parameters have already been thinned
        -for use only in taking the midpoint of an interval
        for use in a single hull design using form parameter
        design.
        -the idea here is simply to put a functional layer in between
        the FPD tool and this constraint language hull space class just in  case
        one of them changes.
        """
        mk_name,mk_value = self.get_from_rgp(var.name)
        assert(len(mk_value)==1),"Design space is multivalued"
        x = .5
        return  mk_value[0].getpoint(x)
    
    
    
    #    def smartgetter(self, atstr):
    #        """No Longer Needed
    #        """
    #        val = getattr(self, atstr)
    #        if isinstance(val, lp.PStates):
    #            return val
    #        else:
    #            return None
    
    def tree_search(self,
                    sobol_seq = None,
                    cmax=10,
                    maxACiter=None):
        """Used by ShipDesigner to narrow the design space
        down to 1 thin interval valued set of design parameters
        for use in Form Parameter Design
        
        States tracks it's self history
        This may cause objects never to be 
        garbage collected?
        Anyway thats how design_sq was working in opt_simple_hull
        
        
        dlist = hullrulesnet.get_thick_list()
        hullrulesnet.wrapped_search(dlist, maxinner=1) 
        """
        if sobol_seq is None:sobol_seq=self.sobol_seq
        #
        dlist = self.get_thick_list()
        #
        count = 0
        while dlist and count<cmax:
            print 'beginning tree search\n:'
            print 'count : ',count
            self.wrapped_search(dlist, maxinner=1,
                                maxACrevise_iter=maxACiter) 
            dlist = self.get_thick_list()
            count+=1
        
        #toying around with putting a paver back in:
        self.feasible_designs.append(self.rgp.env.get_latest_child())
                
        return
    
    def wrapped_search(self,
                       klist,
                       maxinner = 1,
                       maxACrevise_iter=None):
        loop = True
        inner = 0
        redolist = klist
        while loop and inner<maxinner:
            redolist = self.search_loop(redolist,
                                        maxiter=maxACrevise_iter)
            if len(redolist) == 0: #now this is redundant
                loop=False
                print '*'
            else:
                print 'redo list contents: ',redolist
            inner +=1
        return
    
    
    def set_random_generator(self, which='gauss'):
        self.rtype = which
        return
    
    def get_next(self, which = None):
        if which is None: which = self.rtype
        return self.random_type[which]()

    def search_loop(self, 
                    inlist,
                    this_iter = 0,
                    maxinner = 1,
                    ith = 0,
                    operator='Default',
                    maxiter=None):
        """
        inlist = search_list
        """
        mc = 0
        print 'using ',self.rtype,' random sequence '
        redo = []
        for i in range(len(inlist)):
            var = inlist[i]
            x = self.get_next(self.rtype)
            if self.verbose: print 'x = ',x
            #
            if operator == 'Default':
                checker,mc = self.set_this(var,x,
                                        maxit=maxiter,
                                        mc=mc)
            elif operator == 'Narrow':
                checker,mc = self.narrow_this(var)
            #
            if not checker:
                redo.append(var)
        if len(redo)==0:
            return redo
        
        else:
            if this_iter>maxinner:
                return redo
            else:
                this_iter += 1
                print 'redo',redo
                if this_iter %2 == 1:
                    return self.search_loop(redo,
                                            this_iter,
                                            maxinner = 1,
                                            ith = 0,
                                            operator='Narrow',
                                            maxiter=maxiter)
                else:
                    return self.search_loop(redo,
                                            this_iter,
                                            maxinner = 1,
                                            ith = 0,
                                            operator='Default',
                                            maxiter=maxiter)
    
    

    
    def set_this(self, var, x, ith=-1,maxit=None,mc=0):
        """
            var = key
            ith=-1
            
            ?Do this the old fashioned way 
            to go back one state in the event
            that an equality rule returns the null state
            
            or, -since that would require adding 
                  def equal()
                  to the code,
                -Which would then mean future users
                  building similar classes would have to
                  be smart about internals...
            
            ?how about copy.copy(self.rgp.env)
            try the equality
            accept new env only if it works.
            
            (only) issue: sometimes backing up
            further is nice
            
            
            NOTE:
                
                
                #-----------------------------------------
                ret_state = copy.copy(self.rgp.env)
                self.rgp.add_one_rule(var,var.name)
                self.rgp.compute_fresh_rules_graph(maxiter=3)#maxiter=8)
                #-----------------------------------------
        """
        mk_name, mk_val = self.get_from_rgp(var) #
        if mk_val[ith] is None:
            return False
        val = mk_val[ith].getpoint(x)
        #------ performance boost ;)
        space = self.tol*.3
        value = ia(val-space,val+space)
        #------ 
        # set new value as this rule:
        var = var == value
        #-----------------------------------------
        ret_state = copy.copy(self.rgp.env)
        #self.lastnode = self.rgp.env.get_latest_child() #toying around with putting a paver back in
        self.rgp.add_one_rule(var,var.name)
        self.rgp.compute_fresh_rules_graph(maxiter=maxit)#maxiter=8)
        #-----------------------------------------
        #
        #now check if ok or not
        #(checking here instead of 
          #RulesGraphProcessor's AC_revise function )
        if len(self.rgp.env.states)==0:
            mc += 1 
            #----------------------------
            if mc>5:
                mc = 0
            #                if ret_state.parent is not None:
            #                    print 'up one level'
            #                    ret_state = ret_state.parent
            #                else:
            #                    return False, mc
            #----------------------------
            self.rgp.reset(ret_state,var)
            #----------------------------
            v1i = mk_val[ith].inf
            v2s = mk_val[ith].sup
            v1 = ia(v1i,val)
            v2 = ia(val,v2s)
            if self.verbose:
                #print 'not splitting ',mk_name
                print ' sending ',mk_name, 'to the redo list'
                #print 'v1 = ',v1
                #print 'v2 = ',v2
            #self.rgp.env = self.rgp.env.split( (mk_val[ith],v1,v2) )
            return False, mc #, ret_state
        else:
            self.rgp.env.parent = ret_state
            ret_state.children.append(self.rgp.env)
            if self.rgp._verbose: print 'done ', mk_name,'=>', value
            self.rgp.varsmap[var.name].flist.pop()
            return True, mc #,  self.rgp.env
    
    

    
    def narrow_this(self, var, ith=-1,mc=0):
        """
        idea
        ----------
            -when a bad value is found, 
            split the design space around it
            -narrow instead of thin
            for actual gains
            
        
        findings
        ----------
            -IN PRACTICE THIN SPLITTING GAINS NOTHING
            -it does not have the geometric 
            power of a split around infinity,
            for instance.  The design space stays 
            compact, or amost compact (math needed
            -> close inspection to check this behavior), 
            even if now technically
            split into 2 pieces
            -this should not be a problem for a narrowing operator
            
        dev miscellaneous
        ----------
            var = key
            ith=-1
        """
        x1=.25
        x2=.75
        mk_name, mk_val = self.get_from_rgp(var) #
        if mk_val[ith] is None:
            return True
        vali = mk_val[ith].getpoint(x1)
        vals = mk_val[ith].getpoint(x2)
        value = ia(vali,vals)
        var = var == value
        #-----------------------------------------
        ret_state = copy.copy(self.rgp.env)
        self.rgp.add_one_rule(var,var.name)
        self.rgp.compute_fresh_rules_graph()#maxiter=8)
        #-----------------------------------------
        if len(self.rgp.env.states)==0:
            mc += 1
            #----------------------------
            if mc>5:
                mc = 0
            #                if ret_state.parent is not None:
            #                    print 'up one level'
            #                    ret_state = ret_state.parent
            #                    mc = 0
            #                    return False, mc #because it makes less rigourous the split logic
            #                else:
            #                    return False, mc
            #----------------------------
            self.rgp.reset(ret_state,var)
            #----------------------------
            v1i = mk_val[ith].inf
            v2s = mk_val[ith].sup
            v1 = ia(v1i,vali)
            v2 = ia(vals,v2s)
            if self.verbose:
                #print 'not splitting ',mk_name
                print ' splitting wide! ',mk_name
                print 'v1 = ',v1
                print 'v2 = ',v2
            self.rgp.env = self.rgp.env.split( (mk_val[ith],v1,v2) )
            return False,mc #, ret_state
        else:
            if self.rgp._verbose: print 'narrowed ', mk_name,'=>', value
            self.rgp.varsmap[var.name].flist.pop()
            return True,mc #,  self.rgp.env
    
    
    def get_thick_list(self):
        todo_list = []
        for var in self.rgp.vars:
            try:
                attr = self.__getattribute__(var.name)
                mk_name, mk_val = self.get_from_rgp(attr)
                # if multivalued, simply pick the first.
                iaval = mk_val[0]
                diff = iaval.sup - iaval.inf
                if np.linalg.norm(diff)<self.tol:
                    pass
                else:
                    todo_list.append(attr)
            except:
                pass
        return todo_list
    
    def get_node_list(self, option=''):
        if option.lower() == 'thick':
            return self.get_thick_list()
        else:
            todo_list = []
            for var in self.rgp.vars:
                try:
                    attr = self.__getattribute__(var.name)
                    mk_name, mk_val = self.get_from_rgp(attr)
                    # if multivalued, simply pick the first.
                    iaval = mk_val[0]
                    diff = iaval.sup - iaval.inf
                    todo_list.append(attr)
                except:
                    pass
            return todo_list
    
    
    
    def get_allnode_list(self, option=''):
        if option.lower() == 'thick':
            return self.get_thick_list()
        else:
            todo_list = []
            for var in self.rgp.vars:
                try:
                    attr = self.__getattribute__(var.name)
                    mk_name, mk_val = self.get_from_rgp(attr)
                    # if multivalued, simply pick the first.
                    todo_list.append(attr)
                except:
                    pass
            return todo_list
    
    
    #def initialize_lists(self):
    #    """Not done as 
    #    the PStates variables and
    #    the rules are made up on the fly
    #    """
    #   return 
    
    def print_hull_state(self):
        ll = self.__dict__.keys()
        self.print_list(ll)
        return
        
    def print_list(self, list_):
        for key in list_:
            print self.get_from_rgp(key)
        return
    
    def print_thick_list(self):
        hl = self.get_thick_list()
        self.print_list(hl)
        return
    
    def make_attributes_from_rule(self, node, name = None,
                                  skipthis=False):
        """local metaprogramming
        to add variables to the HullGeometryGenerator class
        (convieniece in access)
        (does nothing if the rule 'node'is already an attribute)
        (a 'node' is another way to look at a variable)
        
        Notes
        ----------
            automated internal use only (recommended)
        """
        if not skipthis:
            if hasattr(self, node.name):
                pass
            else:
                if name is None:
                    setattr(self, node.name, node)
                else:
                    setattr(self, name, node)
                    
        if node.args:
            for child in node.args:
                if isinstance(child, lp.PStates):
                    if not child.discount:
                        self.make_attributes_from_rule(child)
                    else:
                        self.make_attributes_from_rule(child, skipthis=True)
        return
    
    
    def set_rules(self, *rules):
        """
        Parameters
        ----------
            rule or rules 
                
        Returns
        ----------
            Nothing
            
        Actions
        ----------
        self.make_attributes_from_rule:     local metaprogramming
                                            to add variables to the 
                                            HullGeometryGenerator class
                                            (convieniece in access)
                                            (does nothing if the rule 'node'
                                            is already an attribute)
                                            (a 'node' is another way to look
                                            at a variable)
                                        
        self.rgp.add_one_rule:              really add a rule 
                                            to the database:
                                            
                            1.  construct  -walks the rule graph to find the
                                            logic variables it is using for
                                            each node in the rule.
                                            (each node is a relational composition
                                            of two factors (nodes or constants)
                                            to get a third)
                            1.  compiles  -compiles each rule to a list of 
                                            sqKanren functions.
                                            -adds the functions to each 
                                            involved variable's list of functions
                                            it is involved in.  
                                            (not a tautology, but almost ;)
                                            
                                        
        self.rgp.compute_fresh_rules_graph: update the environment list
                                            (aka the database States)
                                            by filtering on the new rules
                                            furthermore,
                                            call arc consistency (AC revise)
                                            on any variables which get 
                                            updated information.
                                            This ensures propogation 
                                            of information througout
                                            and consistency 
                                            (infeasible intervals are excised)
                                            
            
        Notes
        ----------
            testing for existance not needed
            No matter what the var is going to get updated.
            This does allow the user to make a mistake in naming 
            and have it only show up as a strange result after AC_revision
        """
        for rule in rules:
            self.make_attributes_from_rule(rule)    # local stuff to add an attribute to this outside (non sqKanren) class
            # the purpose of which is to drive the search.
            self.rgp.add_one_rule(rule, rule.name)  # the real thing
        self.rgp.compute_fresh_rules_graph()        # AC revise over what can be different
        return 
    
    



def make_hull_rules_net(hullrulesnet):
    """Primary Bare Hull Rules
    -----------------------------------
    
    -minimal requirements for usage
    -this kind of thing works 'on the fly'
        -what do you mean?
        -how do you demonstrate this?
    
    hullrulesnet = HullGeometryGenerator()
    
    >>>  Awp  =  [ia(2125.0, 4290.0)]
    >>>  bwl  =  [ia(25.0, 33.3333333333)]
    >>>  lfwl  =  [ia(0.0, 171.6)]
    
    
    hullrulesnet.bwl = hullrulesnet.bwl == ia(25.,30.)
    hullrulesnet.rgp.add_one_rule(hullrulesnet.bwl)
    hullrulesnet.rgp.AC_revise()
    
    
    >>>  Awp  =  [ia(2125.0, 3861.0)]
    >>>  bwl  =  [ia(25.0, 30.0)]
    >>>  lfwl  =  [ia(0.0, 154.44)]
    
    hullrulesnet.rgp.print_state()
    
    
    ----------------------------------
    
    hullrulesnet.lfwl = hullrulesnet.lfwl == ia(30.,30.)
    hullrulesnet.rgp.add_one_rule(hullrulesnet.bwl)
    hullrulesnet.rgp.AC_revise()
    
    
    hullrulesnet.bwl = hullrulesnet.bwl == ia(25.,30.)
    hullrulesnet.rgp.add_one_rule(hullrulesnet.bwl)
    hullrulesnet.rgp.AC_revise()
    
    
    >>>  Awp  =  [ia(2125.0, 3861.0)]
    >>>  bwl  =  [ia(25.0, 30.0)]
    >>>  lfwl  =  [ia(30.0, 30.0)]
    ----------------------------------
    
    
    TODO:
    ----------------------------------
        -rule to constrain the total length of the
        flat of side
        
        -how about regulating where the bow and stern
        fairness curves can be?
        
        -how about the fore and aft extents of the flat spot?
        
        DevNotes: somehow DWL and CPKeel (but especially DWL)
        work 'better' than SAC rules, despite having fewer points
    """
    
    #
    #***********************************
    # Instantiate the design space variables
    lwl = lp.PStates(name='lwl')
    bwl = lp.PStates(name='bwl')
    draft = lp.PStates(name='draft')
    vol = lp.PStates(name='vol')
    disp = lp.PStates(name='disp')
    
    Cb = lp.PStates(name='Cb')
    Cp = lp.PStates(name='Cp')
    
    Awp = lp.PStates(name='Awp')
    Cwp = lp.PStates(name='Cwp')
    
    Acp = lp.PStates(name='Acp')
    Ccp = lp.PStates(name='Ccp')
    
    Amsh = lp.PStates(name='Amsh')
    Cmidshp = lp.PStates(name='Cmidshp')
    
    LCG = lp.PStates(name='LCG')
    Clcg = lp.PStates(name='Clcg')
    
    Clb = lp.PStates(name='Clb') #length to beam ratio
    Cdl = lp.PStates(name='Cdl') #displacement to length ratio
    
    
    Cld = lp.PStates(name='Cld') #length to depth ratio (not really used)
    
    
    
    #
    #***********************************
    # Flats
    lfwl = lp.PStates(name='lfwl')      #flat water line
    lfcp = lp.PStates(name='lfcp')      #flat center plane
    lfsac = lp.PStates(name='lfsac')    #flat of SAC
    
    
    #***********************************
    #***********************************
    #***********************************
    #
    # 3 Section SAC Curve Variables
    #
    # 3 lengths, areas, centers
    #
    #***********************************
    # rules:  3-part variables
    SAC_entrance_area = lp.PStates(name='SAC_entrance_area')
    SAC_mid_area = lp.PStates(name='SAC_mid_area')
    SAC_run_area = lp.PStates(name='SAC_run_area')
    
    SAC_entrance_len = lp.PStates(name='SAC_entrance_len')
    #lfsac = lp.PStates(name='lfsac')
    SAC_run_len = lp.PStates(name='SAC_run_len')
    
    SAC_fwd_Xc = lp.PStates(name='SAC_fwd_Xc')
    SAC_mid_Xc = lp.PStates(name='SAC_mid_Xc')
    SAC_run_Xc = lp.PStates(name='SAC_run_Xc')
    
    
    
    
    #
    #***********************************
    # Basic Dimension Relations
    """
    #rule: rename vol to disp
    """
    disp = disp == vol
    
    """-----------------------------------------------
    #rule: block coefficient
    """
    Cb = Cb == vol/(lwl*bwl*draft)
    
    
    """-----------------------------------------------
    #rule: prismatic coefficient
    """
    Cp = Cp == Cb/Cmidshp
    
    
    """-----------------------------------------------
    #rule: waterplane coefficient
    #"""
    Cwp = Cwp == Awp/(lwl*bwl)
    
    
    """-----------------------------------------------
    #rule: Centerplane Coefficient
    #"""
    Ccp = Ccp == Acp/(lwl*draft)
    
    
    """-----------------------------------------------
    #rule: midship coefficient
    #"""
    Cmidshp = Cmidshp == Amsh/(bwl*draft)
    
    
    """-----------------------------------------------
    #rule: LCG Coefficient
    #"""
    LCG = LCG == Clcg*lwl
    
    
    
    
    
    """TODO:  make Coefficients of the flat of curve
    for DWL cLProfile and SAC
    """
    """-----------------------------------------------
    #rule:  flat of WL FOWL
    #"""
    #lfwl = lfwl <= Cp*Awp/bwl
    
    
    """-----------------------------------------------
    #rule: flat of CL Keel FOCP
    #"""
    #lfcp = lfcp <= Cp*Acp/draft
    
    """-----------------------------------------------
    #rule: Len of Flat of SAC (1)
        states = (states * (max_height,mid_len,mid_area))  
        mid_area <= Amsh*lfsac
        lfsac = mid_len
        max_height = Amsh
    #"""
    #lfsac = lfsac <= Cp*vol*ia(.2,.4)/Amsh
    #lfsac = lfsac <= Cp*vol*ia(.3,.5)/Amsh
    lfsac = lfsac <= Cp*vol*ia(0.1,0.8)/Amsh
        
    #
    #***********************************
    # Add flat rules all to the Database:
    
    #hullrulesnet.set_rules(lfwl,
    #                       lfcp,
    #                       lfsac)
    
    hullrulesnet.set_rules(lfsac)
    hullrulesnet.rgp.compute_fresh_rules_graph()
    
    
    
    
    
    
    #
    #***********************************
    # Add the above rules to the Database
    # like so:
    hullrulesnet.set_rules(disp,
                           Cb,
                           Cp,
                           Cwp,
                           Ccp,
                           Cmidshp,
                           LCG)
    #***********************************
    #
    hullrulesnet.rgp.compute_fresh_rules_graph()
    
    
    
    
    
    """-----------------------------------------------
    #rule: Length/beam ratio
    #"""
    Clb = Clb == lwl/bwl
    
    """-----------------------------------------------
    #rule: Length/depth ratio 
    - this would not be correct because draft /= depth
    and more importantly 
    this hull is now 'designed at some other draft'
    i.e. the waterline is no longer the true design waterline
    (not even close)
    #"""
    #Cld = Cld == lwl/draft
    
    """-----------------------------------------------
    #rule: displacement/Length ratio
    #"""
    Cdl = Cdl == vol/(lwl*lwl*lwl)
    #
    #***********************************
    # compile to the relational rules base again:
    hullrulesnet.set_rules(Clb,Cdl)
    #
    """
    #rule: Cbl, the length to beam interval
            
    """
    #Cbl = Cbl == ia(.2,.4)
    Clb = Clb == ia(4.2,5.2)
    #
    #***********************************
    # compile to the relational rules base again:
    hullrulesnet.set_rules(Clb)
    hullrulesnet.rgp.compute_fresh_rules_graph()
    
    
    
    
    
    """-----------------------------------------------
    #rule: RULE: flat of SAC <= lfwl
    """
    lfsac = lfsac <= lfwl#*ia(.5,.7) - bad to modify this rule.  why? 
    # because we already have a rule for the fraction relating these two
    # it's [> lfwl = lfwl <= lfsac*ia(1.0,1.2) <] this guy
    #
    #***********************************
    # Add rule to the Database:
    hullrulesnet.set_rules(lfsac)
    
    hullrulesnet.rgp.compute_fresh_rules_graph()
    #lfsac = lfsac == lwl*ia(.03,.22)
    #
    #lfsac = lfsac == lwl*ia(.05,.1)
    #longer hull seems to cry for longer flat
    #lfsac = lfsac == lwl*ia(.08,.22)
    #
    #or maybe not - push more volume to the ends
    #lfsac = lfsac == lwl*ia(.05,.1)
    #not quite that much
    lfsac = lfsac == lwl*ia(.06,.12)
    
    #to much volume asked for in the outer portions:
    #lfsac = lfsac == lwl*ia(.09,.17)
    #lfsac = lfsac == lwl*ia(.2,.35)
    
    #
    #***********************************
    # Add rule to the Database:
    hullrulesnet.set_rules(lfsac)
    
    hullrulesnet.rgp.compute_fresh_rules_graph()
    
    """-----------------------------------------------
    #rule: flat of SAC <= lfcp
    #"""
    lfsac = lfsac <= lfcp
    
    
    #
    #***********************************
    # Add rule to the Database:
    hullrulesnet.set_rules(lfsac)
    hullrulesnet.rgp.compute_fresh_rules_graph()
    
    
    
    #
    #***********************************
    # more Flat rules
    
    
    
    """
    Do we need these any more?
    lfwl = lfwl <= Cp*Awp/bwl
    
    
    lfcp = lfcp <= Cp*Acp/draft
    """
    
    #rule:  flat of WL FOWL
    #lfwl = lfwl <= lfsac*ia(3.0,10.)
    #lfwl = lfwl == lfsac*ia(3.5,4.)
    
    #try to eliminate the drag on the nose of the hull
    #lfwl = lfwl == lfsac*ia(1.5,2.)
    #stock:
    #lfwl = lfwl == lfsac*ia(1.5,1.8)  
    #new try:
    lfwl = lfwl == lfsac*ia(1.2,1.5)  
    
    
    #rule: flat of CL Keel FOCP
    #lfcp = lfcp <= lfsac*ia(2.,2.5)
    #
    #stock:
    #lfcp = lfcp == lfsac*ia(1.1,1.4)
    #new try:
    #lfcp = lfcp == lfsac*ia(1.0,1.1)
    """
    There is a trick here:
        the cLProfile should be flat to the bulb.
    It was not dsigned with this in mind initially.
    -giving me fits now...
    """
    #lfcp = lfcp == lfsac*ia(1.2,1.8)
    #lfcp = lfcp == lfsac*ia(1.,1.2)
    #lfcp = lfcp == lfsac*ia(1.3,1.8)
    lfcp = lfcp == lfsac*ia(1.,1.1)
        
        
        
    #
    #***********************************
    # Add flat rules all to the Database:
    hullrulesnet.set_rules(lfwl,
                           lfcp)
    hullrulesnet.rgp.compute_fresh_rules_graph()
    
    
    
    
    lfcp = lfcp <= lfwl
    hullrulesnet.set_rules(lfcp)
    hullrulesnet.rgp.compute_fresh_rules_graph()
    
    
    #
    #***********************************
    # Set the design space intervals themselves
    # not done here.  Using something else to do the same
    
    
    #
    #***********************************
    # rule: instantiate Cp with tight limits
    #Cp = Cp == ia(.65,.75) #fail safe  ?
    #Cp = Cp == ia(.4,.65)
    #Cp = Cp == ia(.7,.85)
    #Cp = Cp == ia(.4,.55)
    #Cp = Cp == ia(.6,.7)
    #
    #Cp = Cp == ia(.6,.68)
    Cp = Cp == ia(.68,.8) #june 20
    #Cp = Cp == ia(.8,.88) #July 06
    
    #
    #Cp = Cp == ia(.6,.75)
    #Cp = Cp == ia(.7,.80)
    #Cp = Cp == ia(.75,.88) #postulate:  if Awp goes down, then cp needs to go up to avoid collapsing area
    hullrulesnet.set_rules(Cp)
    hullrulesnet.rgp.compute_fresh_rules_graph()
    
    
    #
    #***********************************
    # NEW rule for OSV DWL area, Awl
    """
    Awp >= bwl*( lfwl + ia(0.5,0.5)*(lwl-lfwl) ) 
    
    because at less than this, there is a good
    chance that the DWL will curve across the Centerplane
    near the bow. - Because its less than the flat plus the triangle area
    with minimum aft (box) area.
    
    -This is an attempt to get a DWL area constraint that is wiser
    about OSV shape - the block area behind the flat of DWL
    will actually require knowing where
    the lfwl actually ends up sitting, 
    unless we just roll both lfwl and rectangular aft shape into
    one extended lfwl.
    -maybe this last bit is the real way to go.
    -the thing is that the 'old' flat of side DWL
    lfwl I mean, is tied to the location of all those curves...
    """
    
    Awp = Awp >= bwl*( lfwl + (lwl-lfwl)*ia(0.75,.75) )#*ia(.5,.5) #do not use 1/2 here.
                                                # design with full beam!
    hullrulesnet.set_rules(Awp)
    hullrulesnet.rgp.compute_fresh_rules_graph()
    
    
    
    
    #""" #NOT USESD
    Awp = Awp >= bwl*( lfwl )#*ia(.5,.5) #do not use 1/2 here.
                        # design with full beam!
    hullrulesnet.set_rules(Awp)
    hullrulesnet.rgp.compute_fresh_rules_graph()
    #"""
    
    
    """ redundant
    Awp = Awp <= bwl*( lwl )#*ia(.5,.5) #do not use 1/2 here.
                        # design with full beam!
    hullrulesnet.set_rules(Awp)
    hullrulesnet.rgp.compute_fresh_rules_graph()
    #"""
    
    
    #Acp = Acp >= draft*( lfcp + (lwl-lfcp)*ia(0.75,0.75) ) 
    Acp = Acp >= draft*( lfcp + (lwl-lfcp)*ia(0.6,0.6) ) 
    hullrulesnet.set_rules(Acp)
    hullrulesnet.rgp.compute_fresh_rules_graph()
    
    
    """ redundant
    Acp = Acp <= draft*( lwl ) 
    hullrulesnet.set_rules(Acp)
    hullrulesnet.rgp.compute_fresh_rules_graph()
    #"""
    
    
    """
#    ###
#    ###
#    ###
#    #***********************************
#    #***********************************
#    #***********************************
#    #
#    # 3 Section SAC Curve Rules
#    #
#    #***********************************
#    # rules:  Xc relationships
#    SAC_fwd_Xc = SAC_fwd_Xc == SAC_entrance_len*ia(.64,.68) # ia(.66,.72) #ia(.66,.72) #ia(.64,.68) #
#    # Add Xc rules all to the Designbase:
#    hullrulesnet.set_rules(SAC_fwd_Xc)
#    hullrulesnet.rgp.compute_fresh_rules_graph()
#    
#    
#    SAC_mid_Xc = SAC_mid_Xc == SAC_entrance_len + lfsac*ia(.5,.5)
#    # Add Xc rules all to the Designbase:
#    hullrulesnet.set_rules(SAC_mid_Xc)
#    hullrulesnet.rgp.compute_fresh_rules_graph()
#    
#    
#    SAC_run_Xc = SAC_run_Xc == SAC_entrance_len + lfsac + SAC_run_len*ia(.32,.36) #ia(.28,.34) #ia(.28,.33) #ia(.32,.36) #
#    # Add Xc rules all to the Designbase:
#    hullrulesnet.set_rules(SAC_run_Xc)
#    hullrulesnet.rgp.compute_fresh_rules_graph()
#    
#    
#    #
#    vol = vol == (SAC_entrance_area*SAC_fwd_Xc +\
#                    SAC_mid_area*SAC_mid_Xc +\
#                        SAC_run_area*SAC_run_Xc)/LCG
#    #
#    #
#    #***********************************
#    # Add area rules all to the Designbase:
#    hullrulesnet.set_rules(vol)
##    hullrulesnet.set_rules(SAC_fwd_Xc,
##                           SAC_mid_Xc,
##                           SAC_run_Xc,
##                           vol)
#    hullrulesnet.rgp.compute_fresh_rules_graph()
#    
#    
#    
#    #
#    #***********************************
#    # rules:  area relationships
#    SAC_entrance_area = SAC_entrance_area == vol*ia(.29,.7) #ia(.29,.6)  #  ia(.25,.4)  # ia(.29,.7) #
#    #SAC_mid_area = SAC_mid_area == vol*ia(.25,.5)
#    #***********************************
#    # Another Area Requirement (mid section box)
#    SAC_mid_area = SAC_mid_area == bwl*draft*lfcp
#    #***********************************
#    SAC_run_area = SAC_run_area == vol*ia(.29,.7) #ia(.29,.6)   # ia(.25,.4)   # ia(.29,.7) #
#    vol = vol == SAC_entrance_area + \
#                    SAC_mid_area + \
#                        SAC_run_area
#    #
#    #***********************************
#    # Add area rules all to the Designbase:
#    hullrulesnet.set_rules(SAC_entrance_area,
#                           SAC_run_area,
#                           SAC_mid_area,
#                           vol)
#    hullrulesnet.rgp.compute_fresh_rules_graph()
#    
#    
#    
#    
#    
#    #
#    #***********************************
#    # rules:  length relationships
#    SAC_entrance_len = SAC_entrance_len == lwl*ia(.29,.7) #ia(.29,.6)  # ia(.36,.48)  # ia(.29,.7) #
#    SAC_run_len = SAC_run_len == lwl*ia(.29,.7) #ia(.29,.6)     # ia(.36,.48)  # ia(.29,.7) #
#    lfsac = lfsac == lwl - (SAC_entrance_len + \
#                                SAC_run_len)
#    #
#    #***********************************
#    # Add length rules all to the Designbase:
#    hullrulesnet.set_rules(SAC_entrance_len,
#                           SAC_run_len,
#                           #lwl,
#                           lfsac)
#    hullrulesnet.rgp.compute_fresh_rules_graph()
    
    ###
    ###
    ###
    #"""
    
    
    
    return hullrulesnet
    



    

    
"""
def get_curve(xb,yb,
              xe,ye,
              alphab=None,alphae=None,
              cb=None,ce=None,
              area=None,
              xc=None,yc=None,
              nCV=None):
##"""
#    """
#        given : the basic form parameters
#        output : a curve satisfying them
#    #"""
"""
    data_store = FPDIdeal(xb,
                             yb,
                             xe,
                             ye,
                             alphab,
                             alphae,
                             cb,
                             ce,
                             area,
                             xc,
                             yc,
                             nCV)
    
    return data_store.make_curve()
"""



def set_design_space():
    """
    USER INPUT HERE:
    
        *Choose your design space input parameters here
        
    
    TODO: add design space rules input hooks for the bbow?
    
    
    use 
    SD.export('designspace')
    SD.export('rhino')
    to save the design
    and the simple curves for ComplexHull maker
    to turn into a 5 part surface for Rhino.
    
    
    1.) narrow the mid-flat zone
        to make more volume go to towards the ends
    """
    return DesignSpecification(
                            #lwl = (80.,130.), #making tubs
                            lwl = (110.,160.),
                            #draft = (12.,23.),
                            draft = (15.,26.),
                            #bwl = (22.,34.), #full width of vessel
                            #bwl = (22.,35.), #full width of vessel
                            bwl = (25.,34),#35.), #full width of vessel
                            vol=(1000.,35000.),
                            #
                            # initial:
                            #LCG = (30.,75.),
                            #Clcg = (.45,.49), #location Coeff for the LCG
                            #
                            #new try (coming up...):
                            #LCG = (35.,75.), #square DWL aft may make 
                            # no more tubs, so we have to adapt LCG (no bearing on tubiness though)
                            LCG = (55.,85.), #square DWL aft may make 
                            #
                            # for more vol than you realize
                            Clcg = (.48,.52), #location Coeff for the LCG
                            # so center up the LCG.
                            #
                            #Cb = (0.5,0.75), #fail safe ?
                            #Cb = (.5,.7),
                            Cb = (.5,.75), #high prismatic, low block.
                            #
                            #Cb = (0.5,0.7),
                            #Cb = (0.65,0.75),
                            #Cb = (0.55,0.8),
                            #Cmidshp = (0.94,.99), #midship coefficient
                            #Cmidshp = (0.84,.89), #midship coefficient
                            ##Fails: Cmidshp = (0.77,.85), #midship coefficient
                            #Cmidshp = (0.84,.94), #midship coefficient
                            #Cmidshp = (0.88,.95), #midship coefficient
                            Cmidshp = (0.9,.95), #midship coefficient
                            #
                            # Why does this tradoff seem
                            # to be in effect?:
                            #
                            # -make everything else pretty?
                            # Cwp = (.6,.95),#water plane coefficient
                            # Ccp = (.6,.95) #centerplane coefficient
                            #
                            # -enforce actual area on the hull curves:
                            
                            # pretty good Cwp:
                            #Cwp = (.5,.88),#water plane coefficient
                            #new try:
                            # (after all, the aft DWL is 'square')
                            
                            #Cwp = (.75,.9),#water plane coefficient
                            # need more here to 
                            # help the fwd fairness cuves stay flat to cL:
                            #Cwp = (.85,.9),#water plane coefficient
                            #Cwp = (.8,.85),#water plane coefficient
                            #Cwp = (.82,.88),#water plane coefficient
                            
                            #Cwp = (.81,.95),#water plane coefficient
                            #Cwp = (.7,.88), #water plane coefficient
                            #Cwp = (.75,.92), #water plane coefficient
                            Cwp = (.89,.94),
                            #Cwp = (.87,.92),
                            #
                            
                            # with increased SAC vol aft
                            # we need to actully cut Ccp down
                            #Ccp = (.76,.85) #centerplane coefficient
                            #Ccp = (.6,.78)
                            #Ccp = (.59,.74)
                            #Ccp = (.65,.85)
                            Ccp = (.78,.87)
                            )
"""Cwp should be > Ccp
since the only difference (besides the 
antisymmetric appearance)
is that Awp is -really square- aft
while Acp is not square fwd.

(At other opposite ends they are similar.)


Cwp should be > Ccp
because Ccp has curveature fore and aft.
Cwp is a box aft and curvy fwd.  ergo Cwp has more area.
"""


def combine_space_with_rules(DS,hullrulesnet):
    """TODO: make ia accept tuples of len 2
    """
    for key in DS.__dict__:
        new_value = DS.__dict__[key]
        var,val = hullrulesnet.rgp.get_name_and_value(key)
        #
        # We could set this as if using '=' on the calss attribute:
        #setattr(hullrulesnet,var.name,ia(new_value[0],new_value[1])) 
        #
        # Or set this using the logic language directly
        #
        # still have to get this class attribute on the fly though, so:
        this_var = getattr(hullrulesnet,key) 
        #
        # here is the rule
        this_var = this_var == ia(new_value[0],new_value[1]) 
        #
        # add it to the database
        hullrulesnet.rgp.add_one_rule(this_var, 
                                      this_var.name) 
    hullrulesnet.rgp.compute_fresh_rules_graph()
    
    
    return hullrulesnet


"""

tcurves = SD.THBhull.child.child.tcurvelist

lcurves = SD.THBhull.child.child.lcurvelist

c = tcurves[0]
c.plot3DmultiList(tcurves,lcurves)


#"""
   
        
if __name__ == "__main__":
    """
    
    u = vertical direction
    v = transverse direction
    u,v = (0.,0.) = keel,fwd
    u,v = (1.,1.) = wL,aft
    
    
    x = transverse direction
    y = vertical direction
    z = longitudinal direction
    """
    
    import_design = False # True # False # 
    
    import relational_lsplines as rlspline 
    
    print '\n--------------------------------------------------------------'
    print 'Start Generation of Hull with Bulbous Bow '
    print 'representation type: THB-spline'
    print '--------------------------------------------------------------\n'
    spp = rlspline.ADILS.SolverPostProcessor
    
    #key demo of the 'domain specific language based' PhD code:
    hullrulesnet = make_hull_rules_net(HullGeometryGenerator(rtype='gauss',
                                                             verbose=True) ) 
    DS = set_design_space()
    
    hullrulesnet = combine_space_with_rules(DS,
                                            hullrulesnet)
    #checking things early on:
    #hullrulesnet.compute_fresh_rules_graph()
    #hullrulesnet.compute_rules_graph()
    #hullrulesnet.tree_search()
    #thicklist = hullrulesnet.get_thick_list()
    #hullrulesnet.tree_search()
    #checkit()
    
    """
    SD.make_bare_hull()
    SD.make_bare_hull()
    SD.bulbous_bow_random_design_selection()
    SD.make_bulbous_bow()
    #"""
    
    def doall():
        SD = ShipDesigner(design_space=hullrulesnet)
        SD.bare_hull_random_design_selection_new()
        SD.bulbous_bow_random_design_selection()
        SD.make_bare_hull()
        SD.make_bulbous_bow()
        #SD.make_THBspline_complex_hull()
        return SD
    #SD = doall()
    
    def dothis():
        SD.bulbous_bow_random_design_selection()
        SD.make_bulbous_bow()
        SD.make_THBspline_complex_hull()
        SD.plot_THB_ship_hull()
        return
    #self = SD.hull
    
    if import_design:
    
        SD = ShipDesigner(design_space=hullrulesnet)
        HyperParameters = SD.hp
        
        importthis = 'designspace_wavy_March1'
        importthis = 'designspace_March3'
        #importthis = 'designspace_March3_2'
        #importthis = 'designspace_toughMarch8'
        #SD.Import('designspace',importthis) #use a prior saved hull design 

        SD.Import('designspace') #use default 
        if True:
            #(or space of designs, 
            #since they have the same representation)
            #
            #*******************************************************
            #make changes to hyperparameters here
            HyperParameters = SD.hyperparameter_basic()
            #
            #*******************************************************
            # new, automatable parameters
            # 
            #**************************************************************
            # Equivalent to the starting 
            # dropout = [3,7] (NONLINEAR)
            # s,s,s good hulls of March 1st 2018:
            #
            HyperParameters.ignore = []# key here
            #
            HyperParameters.dropout = [4,9]
            
            # include transverse vertices in initial
            # linear longitudinal interpolation
            HyperParameters.use_in_linear_interp = [0,1,2,
                                                    3,4,
                                                    6,
                                                    9,10,
                                                    11,12,13] 
            # include transverse vertices in final
            # longitudinal interpolation
            HyperParameters.use_in_longi_solver = [0,1,2,
                                                   3,4,
                                                   6,
                                                   9,10,
                                                   11,12,13]  
            # 
            #************************************************************** 
            #**************************************************************
            # 
            # Pretty good for March 3 tough hull
            #
            HyperParameters.ignore = []
            #
            #HyperParameters.dropout = [4,6,9]
            HyperParameters.dropout = [4,9]
            #HyperParameters.dropout = []  #make lots of transverses to make it hard on the longitudinals
           
            #HyperParameters.dropout = [6]  #hard on the longitudinals, 
            #but maybe avoid the dreaded center unconformity!
            #HyperParameters.dropout = [3] 
            """
            #deliberatly make thing tough!
            #-to provid a nice set of hard curves for the multilevel solver
            #"""
            # include transverse vertices in initial
            # linear longitudinal interpolation
            HyperParameters.use_in_linear_interp = [0,1,2,
                                                    3,4,
                                                    6,
                                                    9,10,
                                                    11,12,13] 
            # include transverse vertices in final
            # longitudinal interpolation
            HyperParameters.use_in_longi_solver = [0,1,2,
                                                   3,4,
                                                   6,
                                                   9,10,
                                                   11,12,13]  
            HyperParameters.fwd_fairness_location_ratio = .3
            HyperParameters.fwd_transition_location_ratio = .6
            HyperParameters.aft_transition_location_ratio = .4
            HyperParameters.aft_fairness_location_ratio = .65
            # 
            #************************************************************** 
            HyperParameters.fwd_fairness_location_ratio     = .3
            HyperParameters.fwd_transition_location_ratio   = .6
            HyperParameters.aft_transition_location_ratio   = .3
            HyperParameters.aft_fairness_location_ratio     = .6
            #**************************************************************
            # 
            #**************************************************************
            #**************************************************************
            # Equivalent to the starting 
            # dropout = [3,7] (NONLINEAR)
            # s,s,s good hulls of yesterday:
            #
            HyperParameters.ignore = []#
            HyperParameters.dropout = []
            HyperParameters.use_in_linear_interp = [0,1,2,
                                                    3,5,
                                                    6,
                                                    8,10,
                                                    11,12,13] 
            HyperParameters.use_in_longi_solver = [0,1,2,
                                                   5,
                                                   8,
                                                   11,12,13]  
            HyperParameters.fwd_fairness_location_ratio     = .45
            HyperParameters.fwd_transition_location_ratio   = .7
            HyperParameters.aft_transition_location_ratio   = .35
            HyperParameters.aft_fairness_location_ratio     = .6
    
            #
            
            # 
            #************************************************************** 
            #set hyperparameters 
            # this way:
            SD.bare_hull_random_design_selection_new(
                                hyperparameters=HyperParameters)
            SD.bulbous_bow_random_design_selection()
            #more directly:
            #SD.thin_hull_design.change_hyperparameters(hyperparameters=None)
            #
            print '----------'
            SD.thin_hull_design.print_hyperparameters()
            print '----------'
            #
            #or over here:
            SD.make_bare_hull(hyperparameters = None)
            
            #passing None will not effect hyperparameters 
            # that have already been set
            
            # using the same space over and over again to 
            # get a baseline for hyperparameters
            # - only for use in tuning the tuning!
            # - or in passing around bare hull designs in really succinct form
            SD.plot_bare_hull_results()
            SD.make_THB_bare_hull()
            SD.THBhull.re_sample(50,50)
            SD.THBhull.compute_surface()
            SD.THBhull.plotSurface(Shade=True)
            #SD.THBhull.plotSurface(Rstride=10,
            #                       Cstride=10,
            #                       Shade=False)
            
            def check_match():
                import utility_optimization as uopt
                uopt.THBspline_to_Bspline_ckdiff_evaluator(SD.THBhull)
                return
            """
            
            print SD.hull.fwd_fairness_location_ratio
            print SD.hull.fwd_transition_location_ratio
            print SD.hull.aft_transition_location_ratio
            print SD.hull.aft_fairness_location_ratio
            print SD.hull.dropout
            print SD.hull.midship_design_vector
            #"""
            ##
            ###
            ##
            #"""
            #
            SD.make_bulbous_bow()
            SD.make_THBspline_complex_hull()
            #does this:
            #  SD.make_THB_bare_hull()
            SD.THBhull.compute_surface()
            SD.plot_THB_ship_hull()
            #
            #"""
            #
            #
            #
            print SD.hull.lcurvenet[4].vertices
            SD.THBhull.plotSurface(colorize_curves=True)
    else:
        SD = ShipDesigner(design_space=hullrulesnet)
        #
        #*******************************************************
        #make changes to hyperparameters here
        HyperParameters = SD.hyperparameter_basic()
        #**************************************************************
        #**************************************************************
        # Equivalent to the starting 
        # dropout = [3,7] (NONLINEAR)
        # s,s,s good hulls of yesterday:
        #
        HyperParameters.ignore = []#
        #
        HyperParameters.dropout = [4,6,9]
        #HyperParameters.dropout = [4,6]
        #HyperParameters.dropout = []  #make lots of transverses to make it hard on the longitudinals
       
        #HyperParameters.dropout = [6]  #hard on the longitudinals, 
        #but maybe avoid the dreaded center unconformity!
        HyperParameters.dropout = []# [3] 
#        HyperParameters.use_in_linear_interp = [0,1,2,
#                                                3,4,
#                                                6,
#                                                9,10,
#                                                11,12,13] 
#        HyperParameters.use_in_longi_solver = [0,1,2,
#                                               3,4,
#                                               6,
#                                               9,10,
#                                               11,12,13]  
        HyperParameters.use_in_linear_interp = [0,1,2,
                                                3,5,
                                                6,
                                                8,10,
                                                11,12,13] 
        HyperParameters.use_in_longi_solver = [0,1,2,
                                               5,
                                               8,
                                               11,12,13]  
        
        
        #HyperParameters.fwd_fairness_location_ratio     = .55 #try to avoid bow slope.
        HyperParameters.fwd_fairness_location_ratio     = .36
        # getting farther and farther from the original intention 
        # with all of this infrastructure. :(
        HyperParameters.fwd_transition_location_ratio   = .75
        #HyperParameters.aft_transition_location_ratio   = .35
        #HyperParameters.aft_fairness_location_ratio     = .6
        HyperParameters.aft_transition_location_ratio   = .3
        HyperParameters.aft_fairness_location_ratio     = .65
        #
        #
        # New Hull Design Parameters (Summer 2018):
        HyperParameters.aft_drop = .25
        #
        #
        #HyperParameters.bulb_volume_fraction = .2##.25##.3#.35#.275#.25#.127  #actually section area at bow fairness
        #HyperParameters.bulb_length_fraction = .5##.25##.3#.25 #portion of fwd fairness to 0.0 station
        #
        # good blend, to fat
        #HyperParameters.bulb_volume_fraction = .21
        #HyperParameters.bulb_length_fraction = .43
        #
        # 
        HyperParameters.bulb_volume_fraction = .2
        HyperParameters.bulb_length_fraction = .43
        #
        # .25 for bulb vol fraction & | bulb length fraction where near the limits of the rules as written.
        #  there were long struggles to find a consistent candidate at those settings.
        #
        
        #
        #
        #**************************************************
        # END of Architectural Choices
        #**************************************************    
        #
        
        #
        #1.) make bare hull and derive bulbous bow from that:
        #
        # print '----------'
        #SD.bare_hull_random_design_selection_new(
        #                hyperparameters=HyperParameters)
        #
        
        #2.) OR choose bare hull 
        #       and bulbous bow form parameters
        #       first and 
        #       incorporate that into the bare hull:
        #       integrate bulb FPs with BareHull FPs:
        print '----------'
        SD.complex_hull_design_and_generation(
                        hyperparameters=HyperParameters)
        
        print '----------'
        SD.thin_hull_design.print_hyperparameters()
        print '----------'
        SD.make_bare_hull()
        print '----------'
        
        #"""
        #
        SD.hull.plot_primary_curves() #SD.plot_bare_hull_results()
        #
        #print '----------'
        #SD.bulbous_bow_random_design_selection() #use with option 1.) above
        print '----------'
        SD.make_bulbous_bow(dmy=True) #skip tree_search as we already did it
        # to make an integrated "bare hull" (no longer so bare)...
        
        
        
        
        
        #
        #"""
        SD.make_THBspline_complex_hull() #removed 90 nose
        SD.plot_THB_ship_hull()
        #"""
        #SD.hull.plot_hull_transverse_curves_and_longitudinal_curves()
        
        spp = rlspline.ADILS.SolverPostProcessor
        try:
            SACpp = SD.hull.Lsplines.SAC.postprocess
        except:
            print 'SAC failed to solve'
            Lspline = SD.hull.Lsplines.SAC
            SACpp = spp(Lspline)
        try:
            sfcpp = SD.hull.Lsplines.sternfairing.postprocess
        except:
            print 'stern fairness curve failed to solve'
            sfcpp = spp(SD.hull.Lsplines.sternfairing)
        # removed long_keel: doing knot insert delete to match old CProfile split portions
        #    try:
        #        longkeel = spp(SD.hull.Lsplines.long_keel)
        #        print '\n Checking Longitudinal Keel Split Curve:'
        #        print 'split from the CProfile curve and reparameterized'
        #        print longkeel.badparameters
        #    except:
        #        print 'longitudinal keel split from CProfile failed to solve'
        #        sfcpp = spp(SD.hull.Lsplines.sternfairing)
        
        
        self = SD.hull
        """
        import thbsurface as thbspline
        from myspline import rbspline
        self = SD.THBhull
        tlist = SD.THBhull.tcurvelist
        cv = tlist[0]
        cv.plot3DmultiList([],tlist)
        cv.plot3DmultiList(SD.hull.tcurvenet,tlist)
        SD.THBhull.plotSurface()
        #"""
        SD.THBhull.plotSurface(colorize_curves=True)
        #
        SD.THBhull.plotSurface(colorize_curves=False,
                               simplefancy=False,
                               fancy=False)
        #
        SD.THBhull.plotSurface(colorize_curves=True,
                               simplefancy=False,
                               fancy=False)
        #
        SD.THBhull.plotSurface(fancy=True,
                               simplefancy=True)
        
    
    ax=SD.hull.SAC.plotcurve_detailed(normalize=True)
    ax=SD.hull.DWL2D.plotcurve_detailed(normalize=True,
                                        scale_curve=SD.hull.SAC,
                                        canvas=ax)
    ax=SD.hull.CProfile2D.plotcurve_detailed(normalize=True,
                                             scale_curve=SD.hull.SAC,
                                             canvas=ax)
    
    
    
    #ax = SD.hull.SAC.plotcurve_detailed(normalize=True)
    
    ax=SD.hull.Bulb_SAC.plotcurve_detailed(canvas=ax,
                                           normalize=True,
                                           scale_curve=SD.hull.SAC)