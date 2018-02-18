import unittest
from main.MOI import example


class TestMoI(unittest.TestCase):
    def test_example(self):
        self.assertEqual(example(2), 4)

    def test_example2(self):
        self.assertEqual(example(0), 0)