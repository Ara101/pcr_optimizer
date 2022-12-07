"""
Created Fall 2022
Authors: Lily Torp, K. Lionel Tukei
"""

import unittest 
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
        self.assertIsInstance(result,tuple)
        self.assertEqual(result,None) 
        self.assertEqual(result0,None)
        self.assertEqual(result1,None,)

    def test_recommend(self):
        result0 = self.entry.recommend(factor="cost")
        result = self.entry.recommend(factor="time")
        self.assertEqual(len(self.entry.recommend()),17)
        self.assertIsNotNone(result)
        self.assertIsInstance(result,str)
        self.assertIsNotNone(result0)
        self.assertIsInstance(result0,str)
        
if __name__ == '__main__':
    unittest.main()