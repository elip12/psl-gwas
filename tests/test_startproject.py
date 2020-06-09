import sys
import subprocess
import unittest
import shutil
from os.path import join

class E2E(unittest.TestCase):
    def test_start_project(self):
        rawsrc = join('tests', 'data', 'raw')
        rawdst = join('testproject', 'data')

        #completed_process = subprocess.run(['./bin/startproject.sh', 'testproject'], input=b'y', stdout=subprocess.DEVNULL)
        #self.assertEqual(completed_process.returncode, 0)
        #subprocess.run(['cp', '-r', rawsrc, rawdst]) 
        


    @classmethod
    def tearDownClass(cls):
        pass
        #subprocess.run(['rm', '-rf', 'testproject'])

if __name__ == '__main__':
    unittest.main()
