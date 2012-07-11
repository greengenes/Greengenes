#!/usr/bin/env python

from numpy import array
from greengenes.masking import inflate_by_mask, img_best_16s_per_genome
from cogent.util.unit_test import TestCase, main

class MaskTests(TestCase):
    def setUp(self):
        pass

    def test_inflate_by_mask(self):
        """Inflate a sequence by a mask"""
        mask = array([False, False, True, True, False, True])
        seq = 'ATG'
        exp = '--AT-G'
        obs = inflate_by_mask(seq, mask)
        self.assertEqual(obs,exp)

        self.assertRaises(IndexError, inflate_by_mask, 'AT', mask)
        self.assertRaises(ValueError, inflate_by_mask, 'ATTTT', mask)

    def test_img_best_16s_per_genome(self):
        """the best, ONLY THE BEST"""
        masked = {'1|a':'aatt--ggcc', '2|a':'atgc---.-a',
                  '1|b':'aattcc'}
        unmasked = {'1|a':'aaattt--gggccc','2|a':'atgcccc-.-a',
                    '1|b':'aattcc'}
        exp = {'1|a\tmasked_length=8\tunmasked_length=12':'aatt--ggcc',
               '1|b\tmasked_length=6\tunmasked_length=6':'aattcc'}
        obs = img_best_16s_per_genome(masked,unmasked)
        self.assertEqual(obs, exp)

if __name__ == '__main__':
    main()
