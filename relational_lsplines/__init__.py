"""
pip usage from shell: cd to the directory where this file is found

and type:
    if brand new:
        pip install . 
    if updating:
        pip install . --upgrage
"""

import opt_simple_hull
#
# this is not the way:
#import opt_simple_hull.DesignSpecification
#import opt_simple_hull.ShipDesigner

from opt_simple_hull import DesignSpecification
from opt_simple_hull import ShipDesigner
#from opt_simple_hull import DesignSpace
#from opt_simple_hull import DesignParams
#from opt_simple_hull import TestParameters

#
from  automatic_differentiation import ad
lp = opt_simple_hull.lp #my microKanren like classes and extensions thereof
ia = lp.ia #extended interval arithmetic 
# interval isinstance needs to use the same ia(todo: fix?)
#
import  curve as spline
from   ADILS         import IntervalLagrangeSpline, Lagrangian
from   FormParameter import FormParameterDict

from simple_hull_rules_language import HullGeometryGenerator
