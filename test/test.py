from psl_gwas import utility
import unittest


class Test(unittest.TestCase):
    def test_complement(self):
        self.assertEqual(utility.complement('ATCG'), 'TAGC')

if __name__ == '__main__':
    unittest.main()
