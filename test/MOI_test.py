
import unittest
from main.MOI import *
from math import pi

class TestMoI(unittest.TestCase):
        

    def test_calculate_inertia_rotated_rectangle(self):
        aileron_obj = aileron()
        self.assertEqual(calculate_inertia_rotated_rectangle(1.0,1.0,pi/2), 1/12.)

    #def test_calculate_inertia_circular_skin(self):
    #    self.assertEqual(calculate_inertia_circular_skin(),0)

    def test_calculate_rotated_inertia(self):
        self.assertEqual(calculate_rotated_inertia(2, 3, 1, 45), (1.5, 3.5, -0.5))

if __name__ == '__main__':
    unittest.main()