#!/usr/bin/env python

from greengenes.metrics import calc_invariant, calc_nonACGT
from cogent.util.unit_test import TestCase, main

class MetricsTests(TestCase):
    def setUp(self):
        pass

    def test_calc_invariant(self):
        """calc the invariants!!"""
        invariants = [('a','AANNNTNNTNTCCGGCCCNNACN', 14),
                      ('b','AANNTTNNTNCTCGGCCCNNACN', 15)]
        input_sq = "AATACATTCTTCCGGCCCA---T"
        exp = 10.0 / 14.0
        invar_len = 14.0

        obs = calc_invariant(input_sq, invariants)
        self.assertFloatEqual(obs,exp)

    def test_calc_nonACGT(self):
        """Calculate the non-acgt percentage off of nast alignment"""
        a = "--A--T--NN-CC---G"
        b = "--A--T--GG-CC---G"
        c = "--A--T--CC-CC---G"
        d = "--A--T--XYZCC---G"
        e = "--A--N--AAACC---G"
        
        exp_a = 2.0 / 7.0
        exp_b = 0.0 / 7.0
        exp_c = 0.0 / 7.0
        exp_d = 3.0 / 8.0
        exp_e = 1.0 / 8.0

        obs_a = calc_nonACGT(a)
        obs_b = calc_nonACGT(b)
        obs_c = calc_nonACGT(c)
        obs_d = calc_nonACGT(d)
        obs_e = calc_nonACGT(e)

        self.assertEqual(obs_a, exp_a)
        self.assertEqual(obs_b, exp_b)
        self.assertEqual(obs_c, exp_c)
        self.assertEqual(obs_d, exp_d)
        self.assertEqual(obs_e, exp_e)

if __name__ == '__main__':
    main()
