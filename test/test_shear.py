import unittest
from forces_moments_func import int_shear, int_moment


class TestShearAndMoment(unittest.TestCase):
    def test_zero(self):
        s_f = int_shear([0, 10, 20], [0, 0 ,0], 0)
        m_f = int_moment(s_f)
        for x in range(20):
            self.assertEqual(s_f(x), 0)
            self.assertEqual(m_f(x), 0)

    def test_fy(self):
        s_f = int_shear([0, 10, 20], [10, 20, -30], 0)
        m_f = int_moment(s_f)
        ans_shear = [10,10,30,30,0]
        ans_m = [0,50,100,250,400]
        for i, x in enumerate([0, 5, 10, 15, 20]):
            self.assertAlmostEqual(s_f(x), ans_shear[i])
            self.assertAlmostEqual(m_f(x), ans_m[i])

    def test_q(self):
        s_f = int_shear([0, 10, 20], [-10, -20, -70], 5)
        m_f = int_moment(s_f)
        ans_shear = [-10, 15, 20, 45, 0]
        ans_m = [0, -10+15*3/2, -10+40*8/2, -10 + 40*8/2 + (20+45)/2*5, -10 + 40*8/2 + (20+70)/2*10]
        for i, x in enumerate([0, 5, 10, 15, 20]):
            self.assertAlmostEqual(s_f(x), ans_shear[i])
            self.assertAlmostEqual(m_f(x), ans_m[i])