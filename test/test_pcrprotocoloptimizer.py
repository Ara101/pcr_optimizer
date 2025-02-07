"""
Created Fall 2022
Authors: Lily Torp, K. Lionel Tukei
"""

import unittest 
import os
import sys

# Get the current file's directory
current_dir = os.path.dirname(os.path.abspath(__file__))
print(current_dir)
# Get the root directory (assuming your script is two levels deep in the directory structure)
root_dir = os.path.abspath(os.path.join(current_dir, '..'))
print(root_dir)
# Add the root directory to the Python path
sys.path.append(root_dir)

# Now you can import your module
from pcr_optimizer.pcrprotocoloptimizer import pcr


class Testpcr(unittest.TestCase):

    def setUp(self):
        # Annealing positions are standard
        gene = "atggagacagacacactcctgctatgggtactgctgctctgggttccaggttccactggtgacacaagtttgtacaaaaaagttggcaccaagtcgatcctagatggccttgcagataccaccttccgcaccatcaccactgacctcctgtacgtgggctcaaatgacattcagtacgaagacatcaaaggtgacatggcatccaaattagggtacttcccacagaaattccctttaacttcctttaggggaagtcccttccaagagaagatgactgcgggagacaacccccagctagtcccagcagaccaggtgaacattacagaattttacaacaagtctctctcgtccttcaaggagaatgaggagaacatccagtgtggggagaacttcatggacatagagtgtttcatggtcctgaaccccagccagcagctggccattgcagtcctgtccctcacgctgggcaccttcacggtcctggagaacctcctggtgctgtgcgtcatcctccactcccgcagcctccgctgcaggccttcctaccacttcatcggcagcctggcggtggcagacctcctggggagtgtcatttttgtctacagcttcattgacttccacgtgttccaccgcaaagatagccgcaacgtgtttctgttcaaactgggtggggtcacggcctccttcactgcctccgtgggcagcctgttcctcacagccatcgacaggtacatatccattcacaggcccctggcctataagaggattgtcaccaggcccaaggccgtggtggcgttttgcctgatgtggaccatagccattgtgatcgccgtgctgcctctcctgggctggaactgcgagaaactgcaatctgtttgctcagacattttcccacacattgatgaaacctacctgatgttctggatcggggtcaccagcgtactgcttctgttcatcgtgtatgcgtacatgtatattctctggaaggctcacagccacgccgtccgcatgattcagcgtaccgacgcgctggacctggaggagggaggaaacgtctatatcaaggccgacaagcagaagaacggcatcaaggcgaacttctgcatccgccacaacatcgaggacggcggcgtgcagctcgcctaccactaccagcagaacacccccatcggcgacggccccgtgctgctgcccgacaaccactacctgagcgtgcagtccaaactttcgaaagaccccaacgagaagcgcgatcacatggtcctgctggagttcgtgaccgccgccgggatcactttcggcatggacgagctgtacaagggcggtaccggagggagcatggtgagaaagggcgaggagctgttcaccggggtggtgcccatcctggtcgagctggacggcgacgtaaacggccacaagttcagcgtgggcggcgagggtgagggcgatgccaccgttggcaagctgaccctgaagttcatctgcaccaccggcaagctgcccgtgccctggcccaccctcgtgaccaccctgacctacggcgtgcagtgcttcagccgctaccccgaccacatgaagcagcacgacttcttcaagtccgccatgcccgaaggctacatccaggagcgcaccatcttcttcaaggacgacggcaactacaagacccgcgccgaggtgaagttcgagggcgacaccctggtgaaccgcatcgagctgaagggcatcgacttcaaggaggacggcaacatcctggggcacaagctggagtacaacaccggagcagcagcacgctggcgcgggcggcgcatggacattaggttagccaagaccctggtcctgatcctggtggtgttgatcatctgctggggccctctgcttgcaatcatggtgtatgatgtctttgggaagatgaacaagctcattaagacggtgtttgcattctgcaccatgctctgcctgctgaactccaccgtgaaccccatcatctatgctctgaggagtaaggacctgcgacacgctttccggagcatgtttccctcttgtgaaggcactgcgcagcctctggataacagcatgggggactcggactgcctgcacaaacacgcaaacaatgcagccagtgttcacagggccgcagaaagctgcatcaagagcacggtcaagattgccaaggtaaccatgtctgtgtccacagacacgtctgccgaggctctg"
        forward_primer = "atggagacagacacactcctgctatgg"
        reverse_primer = "cagagcctcggcagacgtgt"

        # Annealing positions are different | startr = 10, stopr = 31 , startf = 79, stopf = 60
        gene0 = "caatatggtgagcaatcttttgctcttggggattgaccacagtcgtcgacgcgatgcttgccctagtaagcttggcgtacaattacaagtgataagaggt"
        forward_primer0 = "GAGCAATCTTTTGCTCTTGGGG" 
        reverse_primer0 = "TACGCCAAGCTTACTAGGGC"

        # Forward Primer incompatible, Reverse Primer incompatible
        gene1 = "acggagatttcccccttcctagttccctaatggataagatgtttaagatgttcaaacaaggatattcactgtgccagccctgtggtcgggaagtcgaata"
        forward_primer1 = "ACGGAGATTTCCCCCTTCCTAK" 
        reverse_primer1 = "TTATTCGACTTCCCGACCACAG" 

        self.entry = pcr(gene,forward_primer,reverse_primer)
        self.entry0 = pcr(gene0,forward_primer0,reverse_primer0)
        self.entry1 = pcr(gene1,forward_primer1,reverse_primer1)
        
    def test_countGCcontent(self):
        result = self.entry.countGCcontent()
        self.assertIsNotNone(result)
        self.assertIsInstance(result,tuple)
        self.assertIsInstance(result[0],str)
        self.assertIsInstance(result[1],float)
        self.assertEqual(result[1],56.72990063233966)
        
    def test_check(self):
        result = self.entry.check()
        result0 = self.entry0.check(startr = 10, stopr = 31 , startf = 79, stopf = 60)
        result1 = self.entry1.check()
        positive_result = ["Gene looks good!","Forward primer looks good!","Reverse primer looks good!","Primers and Gene are compatible!"]
        negative_result = ["Gene looks good!","Unacceptable character in forward primer. Check sequence","Reverse primer looks good!","Both Primers and gene are incompatible. Check annealing location for both primers"]
        self.assertIsNotNone(result)
        self.assertEqual(result.iloc[:,0].tolist(),positive_result) 
        self.assertEqual(result0.iloc[:,0].tolist(),positive_result)
        self.assertEqual(result1.iloc[:,0].tolist(),negative_result)

    def test_recommend(self):
        result = self.entry.recommend(factor="cost")
        result0 = self.entry0.recommend(factor="time")
        result1 = self.entry1.recommend()
        self.assertIsNotNone(result)
        self.assertIsInstance(result1,str)
        self.assertEqual(result1,"Unacceptable character in forward primer. Check sequence. Consider use check function to troubleshoot entire selection")
        self.assertEqual(result.iloc[1,:].tolist(),[0.1, 0.18, 0.25, 0.45, 0.5, 0.89, 1.0, 1.78])
        self.assertEqual(result0.iloc[:,0].tolist(),['56.73 degrees Celcius','30 seconds or 0.5 minutes','1.5 seconds or 0.025 minutes','34.71 minutes or 0.58 hours'])
        
if __name__ == '__main__':
    unittest.main()