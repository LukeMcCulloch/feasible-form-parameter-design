# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 14:36:48 2016

@author: lukemcculloch
"""

import numpy as np
import curve as spline


class Signal(object):
    def __init__(self):
        self._handlers = []

    def connect(self, handler):
        self._handlers.append(handler)

    def fire(self, *args):
        for handler in self._handlers:
            handler(*args)
            
class Subject(object):
    def __init__(self):
        self._observers = []

    def attach(self, observer):
        if not observer in self._observers:
            self._observers.append(observer)

    def detach(self, observer):
        try:
            self._observers.remove(observer)
        except ValueError:
            pass
    
    def notify(self, modifier=None, *args):
        print 'notify'
        for observer in self._observers:
            if modifier != observer:
                observer.update(self, *args)

# Example usage
class Data(Subject):
    def __init__(self, name=''):
        super(Data, self).__init__()
        self.name = name
        self._data = 0

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, value):
        print 'setter', value
        self._data = value
        self.notify()
        
class v1(Data):
    default_v1 = 0
    def __init__(self, name=''):
        super(v1, self).__init__(name)
        #self.relations = {}
    
    def update(self, subject):
        print 'v1 update callback', subject
        #if subject.name in relations:
        #    relations[subject.name]._data
        print self._observers[0].name, '<=>',self._observers[0].data
        return
    
class v2(Data):
    default_v2 = 0
    def __init__(self, name=''):
        super(v2, self).__init__(name)
        
    def update(self, subject):
        print 'v2 update callback', subject
        print self._observers[0].name, '<=>',self._observers[0].data
        return
        
class rule1(Data):
    default_r1 = True#False
    def __init__(self, name=''):
        super(rule1, self).__init__(name)
        
    def update(self, subject):
        print 'rule1 update callback', subject
        print self._observers[0].name, '<=>',self._observers[0].data
        return

class HexViewer(object):
    def update(self, subject):
        print 'HexViewer: Subject %s has data 0x%x' % (subject.name, subject.data)


class DecimalViewer(object):  
    def update(self, subject):
        print 'DecimalViewer: Subject %s has data %d' % (subject.name, subject.data)


class hull_data(object):
    def __init__(self):
        return


# Example usage...
def main():
    data1 = Data('Data 1')
    print "Setting Data 1 = 10"
    data1.data = 10

if __name__ == "__main__": 
    main()
    a = v1('var 1')
    b = v2('var 2')
    a.attach(b)
    b.attach(a)
    
    a.data = 10
    b.data = 33