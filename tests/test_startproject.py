import sys
import subprocess
import unittest

class E2E(unittest.TestCase):
    def test_start_project(self):
        completed_process = subprocess.run(['./bin/startproject.sh', 'testproject'], input=b'y', stdout=subprocess.DEVNULL)

    @classmethod
    def tearDownClass(cls):
        subprocess.run(['rm', '-rf', 'testproject'])

if __name__ == '__main__':
    unittest.main()
