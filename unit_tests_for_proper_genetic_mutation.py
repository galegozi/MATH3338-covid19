import unittest
import proper_genetic_mutation as P
import random as R
import Levenshtein as L
class TestProperGeneticMutation(unittest.TestCase):
    def test_mutate(self):
        alphabet = ['A', 'C', 'G', 'T']
        for _ in range(100):
            length = R.randint(1, 100000)
            test_str = ''.join(R.choices(population=alphabet, k=length))
            mutated = P.mutate(test_str)
            self.assertEqual(L.distance(mutated, test_str), 1)
            self.assertEqual(len(mutated), len(test_str))
            for char in mutated:
                self.assertTrue(char in alphabet)
    
    def test_padding_worker_fxn(self):
        for _ in range(100):
            length = R.randint(1, 100000)
            padding = P.padding_worker_fxn()


if __name__ == "__main__":
    unittest.main()