import copy
from  interval_arithmetic import ia

def vector_AND_(xi, delta):
    return_values = copy.deepcopy(xi)
    for i in range(len(xi)):
        
        #element of xi
        a = xi[i].inf
        c = xi[i].sup
        x_simple = ia(a,c)
        
        #element of delta
        d = delta[i].inf
        f = delta[i].sup
        dz_simple = ia(d,f)
        
        new_values = x_simple & dz_simple
        
        
        return_values[i].inf = new_values.inf
        return_values[i].sup = new_values.sup
        return_values[i].isempty = new_values.isempty
        
    return return_values


def vector_nonAD_AND_(xi, delta):
    return_values = copy.deepcopy(xi)
    for i in range(len(xi)):
        #element of xi
        a = xi[i].inf
        c = xi[i].sup
        x_simple = ia(a,c)
        
        #element of delta
        d = delta[i].inf
        f = delta[i].sup
        dz_simple = ia(d,f)
        
        new_values = x_simple & dz_simple
        
        
        return_values[i].inf = new_values.inf
        return_values[i].sup = new_values.sup
        return_values[i].isempty = new_values.isempty
        
    return return_values


