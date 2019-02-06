##
##  Unit testing AD IA
##
## Luke McCulloch
## September 25 2015
##
#import numpy as np
import unittest
#import pytest


from time_analysis import Timer

from interval_analysis import ia 
from automatic_differentiation import ad


class AdditionTests(unittest.TestCase):

    def testOne(self):
        self.failUnless(IsOdd(1))

    def testTwo(self):
        self.failIf(IsOdd(2))

def main():
    unittest.main()

if __name__ == '__main__':
    main()