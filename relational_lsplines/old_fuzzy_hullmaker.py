# -*- coding: utf-8 -*-
"""
Created on Sat Sep 20 13:33:12 2014

@author: lukemcculloch
"""
#python standards:
import numpy as np
np.set_printoptions(precision=3)
np.set_printoptions(linewidth=300)
import matplotlib.pyplot as plt
import copy
#from mayavi.mlab import * 



# animated fuzzy logic:
import sys
import math
from operator import itemgetter
from functools import partial


#import TLM B-Spline Classes
from curve import Bspline               as Bspline     
from curve import interpolatedBspline   as interpolatedBspline
from curve import curveSet
from curve import BsplineSurface
from automatic_differentiation import ad as scalarAD
from initialValues          import InitializeControlPoints, InitializeControlVertices
from Equalityconstraints    import FormParameter
from InequalityConstraints  import InequalityConstraint
#from interfaceGUItoOpt      import interface_ADLsplineOpt
from ADL_LS_OptBspline      import LagrangeSpline#, ADLsplineOpt
from CustomRange            import customrange
# import TLM 2D optimization with Least Squares Min included as Optimum Criteria:


#point interpolation routines:
from minDistance import naive_minimizer

#import TLM fuzzy logic
#from fich import get_TSK_SAC



def range(start, stop, step=1.):
    """Replacement for built-in range function.

    :param start: Starting value.
    :type start: number
    :param stop: End value.
    :type stop: number
    :param step: Step size.
    :type step: number
    :returns: List of values from `start` to `stop` incremented by `size`.
    :rtype: [float]
    """

    start = float(start)
    stop = float(stop)
    step = float(step)

    result = [start]
    current = start
    while current < stop:
        current += step
        result.append(current)
    return result


# membership functions
def square(a,b,x):
    a = float(a)
    b = float(b)
    x = float(x)
    if x < a:
        return 0.0
    elif x <= b:
        return x
    return 0.
    
def up(a, b, x):
    a = float(a)
    b = float(b)
    x = float(x)
    if x < a:
        return 0.0
    if x < b:
        return (x - a) / (b - a)
    return 1.0


def down(a, b, x):
    return 1. - up(a, b, x)


def tri(a, b, x):
    a = float(a)
    b = float(b)
    m = (a + b) / 2.
    first = (x - a) / (m - a)
    second = (b - x) / (b - m)
    return max(min(first, second), 0.)


def trap(a, b, c, d, x):
    offset = 0.
    first = offset+(x - a) / (b - a)
    second = offset+(d - x) / (d - c)
    return max(min(first, offset+1., second), offset)


# hedges
def hedge(p, mvalue):
    """Generic definition of a function that alters a given membership function
    by intensifying it in the case of *very* of diluting it in the case of
    *somewhat*.  """

    mvalue = float(mvalue)
    if not p:
        return 0.0
    return math.pow(mvalue, p)

very = partial(hedge, 2.)
extermely = partial(hedge, 3.)
extremely = partial(hedge, 3.)
somewhat = partial(hedge, 0.5)
slightly = partial(hedge, 1. / 3.)


def fuzziness(domain, func):
    """The fuzziness of a fuzzy subset is the degree to which the values
    of its membership function cluster around 0.5

    >>> fuzziness(range(-10, 30, 1), profitable)
    0.182114

    :param domain: the domain of the function
    :type domain: list
    :param func: membership function
    :type func: function
    :returns: fuzziness value
    :rtype: float
    """
    domain_size = float(len(domain))
    delta = lambda x: x if (x < 0.5) else (1.0 - x)
    result = (2. / domain_size) * sum([delta(func(val)) for val in domain])
    return result


def approximate(fuzz, n, domain):
    hw = fuzz * (max(domain) - min(domain))
    return partial(tri, n - hw, n + hw)

# fuzzy database queries example
companies = [
    ('a', 500, 7), ('b', 600, -9), ('c', 800, 17),
    ('d', 850, 12), ('e', 900, -11), ('f', 1000, 15),
    ('g', 1100, 14), ('h', 1200, 1), ('i', 1300, -2),
    ('j', 1400, -6), ('k', 1500, 12)
]

profit = itemgetter(2)
sales = itemgetter(1)

percentages = map(float, range(-10, 30, 1))
profitable = partial(up, 0., 15.)
high = partial(up, 600., 1150.)

fand = min


def ffilter(predicate, items):
    snd = itemgetter(1)
    return filter(
        lambda x: snd(x) != 0.0,
        map(predicate, items)
    )


def p1(company):
    value = profitable(profit(company))
    return (company, fand(value, 1))


def p2(company):
    a = profitable(profit(company))
    b = high(sales(company))
    return (company, fand(a, b))


def p3(company):
    a = somewhat(profitable(profit(company)))
    b = very(high(sales(company)))
    return (company, fand(a, b))


# shoe example

sizes = range(4, 13, 0.5)

short = partial(down, 1.5, 1.625)
medium = partial(tri, 1.525, 1.775)
tall = partial(tri, 1.675, 1.925)
very_tall = partial(up, 1.825, 1.95)


#small = partial(down, 4., 6.)
def small(size):
    return down(4., 6., size)


#average = partial(tri, 5., 9.)
def average(size):
    return tri(5., 9., size)


#big = partial(tri, 8., 12.)
def big(size):
    return tri(8., 12., size)


#very_big = partial(up, 11., 13.)
def very_big(size):
    return up(11., 13., size)


#fl.near(20, fl.range(0, 40, 1))(17.5)
near = partial(approximate, 0.125)
around = partial(approximate, 0.25)
roughly = partial(approximate, 0.375)


rules = [
    (short, small),
    (medium, average),
    (tall, big),
    (very_tall, very_big)
]


def updated_func(val, func, size):
    first = func(size)
    return (val * first)


def rulebase(height):
    updated = []
    for input_func, output_func in rules:
        val = input_func(height)
        updated.append(
            partial(updated_func, val, output_func)
        )

    rulebase_function = lambda s: sum([r(s) for r in updated])
    return rulebase_function


def centroid(domain, membership_function):
    fdom = map(membership_function, domain)
    first = sum([a * b for (a, b) in zip(domain, fdom)])
    second = sum(fdom)
    return first / second


def shoe_example(h):
    result = centroid(sizes, rulebase(h))
    return result


def centroid_example():
    domain = map(float, range(0, 10))
    membership_function = partial(trap, 2, 3, 6, 9)
    return centroid(domain, membership_function)


def mand(funcs, val):
    return min([func(val) for func in funcs])


def price_example(man_costs=13.25, comp_price=29.99):
    """
    Pricing goods (Cox, 1994).
    The price should be as high as possible to maximize takings but as low as
    possible to maximize sales. We also want to make a healthy profit (100%
    mark-up on the cost price). We also want to consider what the competition
    is charging.

    rule1: our price must be high
    rule2: our price must be low
    rule3: our price must be around twice the manufacturing costs.
    rule4: if the competition price is not very high then our price must be
           around the competition price.
    """

    prices = range(15., 35., 0.5)
    high = partial(up, 15., 35.)
    low = lambda p: 1 - high(p)
    not_very = lambda v: 1 - very(high(v))

    our_price1 = centroid(prices, partial(mand, [high, low]))
    our_price2 = centroid(
        prices,
        partial(mand, [high, low, around(2.0 * man_costs, prices)]),
    )
    our_price3 = centroid(
        prices,
        partial(
            mand, [high, low, around(2.0 * man_costs, prices), lambda p: not_very(comp_price) * around(comp_price, prices)(p)
            ]
        )
    )

    print our_price1, our_price2, our_price3




class ship_example(object):
    """
        Question:
            What about fuzzy resistance curves instead of just fuzzy membership?
             I don't know what the function looks like, I just know, e.g.
                 increasing, decreasing, etc....
    
            Ship design will start with simple guesses for
            the physical relationships amongst parameters.
            
            As knowledge of the design space grows,
            the functional relationships will change and grow as well.
        """
    """
        Ship Design Rules of Thumb:
        
            1. The ship should be as long as possible to be fast
            
            2.  The ship should be as big as possible to be profitable
            
            3.  The ship should be as small as possible to be strong
            
            4.  The ship should be as thin as possible to be fast
            
            5.  The ship should be as wide as possible to be stable
            
    """
    def __init__(self):
        self.length = np.asarray(range(50.,200.,10.0))
        self.width = np.asarray(range(10.,50.,5.0))
        self.depth = np.asarray(range(10.,30.,5.0))
        
        # simple function modifiers
        self.quickly            = partial(hedge, 3.)
        self.somewhat_quickly   = partial(hedge, 2.)
        self.moderately         = partial(hedge, 1.)
        self.somewhat_slowly    = partial(hedge, 0.5)
        self.slowly             = partial(hedge, 1. / 3.)
        
        return
        
    def ship_wave_resistance(self, length):   
        """
            Wave resistance is some decreasing function of length
                could use partial:
                    #decreasing_w_length = partial(down, self.length[0], self.length[-1])
        """
        minimum_wave_resistance = 100.
        decreasing_w_length = lambda length : minimum_wave_resistance + minimum_wave_resistance*down(self.length[0], self.length[-1], length)
        return self.quickly(decreasing_w_length(length)) #a number
        
    def ship_frictional_resistance(self, length):
        """
            Frictional resistance is some increasing function of length
                could use partial:
                #increasing_w_length = partial(up, self.length[0], self.length[-1])
        """
        minimum_frictional_resistance = 88.
        increasing_w_length = lambda length : minimum_frictional_resistance+ 100.*up(self.length[0], self.length[-1], length)
        return self.quickly(increasing_w_length(length)) #a number
       
    def ship_cost_vs_length(self, length):
        """
            Cost is an increasing function of length
        """
        fixed_cost = 100.
        cost_per_unit_length = 10.*(self.length[-1]-self.length[0])
        increasing_w_length = lambda length : fixed_cost+ cost_per_unit_length*up(self.length[0], self.length[-1], length)
        return self.somewhat_quickly(increasing_w_length(length))
        
    
        
    ##
    ## Plotting
    ##
    
    ## Cost Plots
    def graph_ship_cost_vs_length(self, plot=True, title='Trend, Ship Cost vs Length'):
        self.cstvslgth = []
        for el in self.length:
            self.cstvslgth.append(self.ship_cost_vs_length(el))
        self.cstvslgth = np.asarray(self.cstvslgth)
        if plot==True:
            plt.title(title)
            plt.xlabel('Length')
            plt.ylabel('Cost')
            plt.plot(self.length, self.cstvslgth)
            plt.show()
        return
       
       
    ## Resistance Plots
    def graph_ship_wave_resistance(self, plot=True, title='Trend, Wave Resistance vs Length'):
        self.wr = []
        for el in self.length:
            self.wr.append(self.ship_wave_resistance(el))
        self.wr = np.asarray(self.wr)
        if plot==True:
            plt.xlabel('Length')
            plt.ylabel('Resistance')
            plt.plot(self.length, self.wr)
            plt.title(title)
            plt.show()
        return        
    def graph_ship_frictional_resistance(self, plot=True, title='Trend, Frictional Resistance vs Length'):
        self.fr = []
        for el in self.length:
            self.fr.append(self.ship_frictional_resistance(el))
        self.fr = np.asarray(self.fr)
        if plot==True:
            plt.xlabel('Length')
            plt.ylabel('Resistance')
            plt.plot(self.length, self.fr)
            plt.title(title)
            plt.show()
        return
    def graph_ship_resistance(self, title = 'Predicted Wave and Frictional Resistance'):
        try:
            plt.xlabel('Length')
            plt.ylabel('Resistance')
            plt.plot(self.length, self.wr, label='Wave Resistance')
            plt.plot(self.length, self.fr, label='Frictional Resistance')
            plt.legend()
            plt.title(title)
            plt.show()
        except:
            self.graph_ship_frictional_resistance(plot=False)
            self.graph_ship_wave_resistance(plot=False)
            plt.xlabel('Length')
            plt.ylabel('Resistance')
            plt.plot(self.length, self.wr, label='Wave Resistance')
            plt.plot(self.length, self.fr, label='Frictional Resistance')
            plt.legend()
            plt.title(title)
            plt.show()
        return
        
    
    
        
    
    def ship_length(self):
        """
            1. the ship must be as long as possible to minimize wave drag
            1. the ship must be as short as possible to minimize frictinal resistance
            3. 
        """
        longship = partial(up, 50.,200.)
        shortship = lambda p: 1 - high(p)
        not_very = lambda v: 1 - very(high(v))
        
        self.L1 = centroid(self.length, partial(mand, [longship, shortship]))
        self.L2 = centroid(self.length, partial(mand, [longship, 
                                                       shortship, 
                                                       around(10.0 * man_costs, prices)]),)
        return
    
    def print_generic(self, func, interval):
        f=[]
        for i in interval:
            f.append(func(i))
        plt.plot(interval, f)
        plt.show()
        return

self = ship_example()
self.graph_ship_resistance()
self.graph_ship_cost_vs_length()
#self.print_generic( partial(down,0.,10.), [1.,2.,3.,4.,5.])

class Hull(object):
    
    def __init__(self, displacement, speed, guessLength = 100., Cb = .75, design_units='metric', length_type='considerProductionCost'):
    #def __init__(self, displacement, speed, guessLength = 100., Cb = .75, design_units='metric', length_type='Posdunine'):
    #def __init__(self, displacement, speed, guessLength = 100., Cb = .75, design_units='metric', length_type='classic_empirical_formula'):
        
        self.design_units               = design_units
        self.units_dict                 = {}
        self.units_dict['g']            = {'english':32.17,'metric':9.81}
        self.g                          = self.units_dict['g'][self.design_units]
        self.displacement               = displacement
        self.speed                      = speed
        self.Cb                         = Cb
        self.tol                        = 1.e-7
        Guess_Length                    = self.displacement/100.
        Posdunine_length_coefficient    = 10.
        self.setup_design_database(Guess_Length, Posdunine_length_coefficient)
        self.Lpp, self.Fn = self.find_length_between_perpendiculars(length_type)
        
        #Cb = property(get_Cb, set_Cb)  #update on call...
        return
        
    def setup_design_database(self, Guess_Length, Posdunine_length_coefficient):
        self.DesignDatabase         = {}
        self.DesignDatabase['Lpp']  =   {'guess_length':Guess_Length,
                                        'Posdunine_length_coefficient':Posdunine_length_coefficient}
        return
    
    def find_length_between_perpendiculars(self,kind):
        if kind == 'considerProductionCost':
            Lpp, Fn = self.length_for_minCost(self.DesignDatabase['Lpp']['guess_length'])
        elif kind == 'Posdunine':
            Lpp, Fn = self.length_by_simple_coefficient(self.DesignDatabase['Lpp']['Posdunine_length_coefficient'])
        elif kind == 'classic_empirical_formula':
            Lpp, Fn = self.length_by_Ayre(self.DesignDatabase['Lpp']['guess_length'])
        return Lpp, Fn
        
    def length_for_minCost(self, initial_guess):
        """
           Schneekluth’s Formula
           Schneekluth and Bertram, page 2 in the pdf
           
               Depends on:
                   Cb
                   Fn(L) - trencendental
           
               applicability:
                   for bulbous bow ships
                   disp >= 1,000 tons
                   0.16 <= Fn <= 0.32
                   
                Notes:
                    Ships optimized for yeild are about 10%
                    longer than ships optimized using this equation.
                    
                    More recent ships are often shorter than this formula for Lpp
                    
        """
        #initial_guess = self.displacement/100.
        Disp    = self.displacement
        V       = self.speed
        g       = self.g
        Cb      = self.Cb
        tol     = self.tol
        #define design variable:
        L   = scalarAD(initial_guess,1.,0.)
        Fn  = lambda L : V/np.sqrt(L*g)
        obj = lambda L : L - 3.2 * Disp**0.3 * V**0.3 * (Cb+0.5)/( (0.145/Fn(L))+0.5 )        
        niter = 0
        while abs(obj(L).value) > tol and niter < 50:
            F = obj(L)
            dL = F.value/F.der
            L = L - dL
            niter +=1
        return L.value, Fn(L).value
    
    def length_by_simple_coefficient(self, C):
        """
            Posdunine's Formula
            see Martec Proceedings 2010
            Optimization of Ship Hull Prameter of Inland Vessel
            with Respect to Regression Based Resistance Analysis
            
            Saha and Sarker
            
            Depends on C -  an empirical fit coefficient
        """
        Disp    = self.displacement
        V       = self.speed  #in knots
        g       = self.g
        #define design variable:        
        L = C*((V/(V+2.))**2)*((Disp)**(1./3.))
        fn = V/np.sqrt(L*g)
        return L, fn
        
    def length_by_Ayre(self, initial_guess):
        """
            Depends on Fn(L) - trancendental!
        """
        Disp    = self.displacement
        V       = self.speed
        g       = self.g
        tol     = self.tol
        #define design variable:
        L   = scalarAD(initial_guess,1.,0.)
        Fn  = lambda L : V/np.sqrt(L*g)
        obj = lambda L : L - (Disp**(1./3.))*(3.33 + 10.2 * Fn(L)) 
        niter = 0
        while abs(obj(L).value) > tol and niter < 50:
            F = obj(L)
            dL = F.value/F.der
            L = L - dL
            niter +=1
        return L.value, Fn(L).value
        
if __name__ =='__main__':
    displacement    = 10000.
    speed           = 10.
    trad_hull = Hull(displacement, speed)
    print 'Lpp = {},  Fn = {}'.format(trad_hull.Lpp, trad_hull.Fn)