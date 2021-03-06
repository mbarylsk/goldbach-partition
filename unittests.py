#
# Copyright (c) 2016 - 2018, Marcin Barylski
# All rights reserved.

# Redistribution and use in source and binary forms, with or without modification, 
# are permitted provided that the following conditions are met:
# 
# 1. Redistributions of source code must retain the above copyright notice, 
#    this list of conditions and the following disclaimer.
# 
# 2. Redistributions in binary form must reproduce the above copyright notice, 
#    this list of conditions and the following disclaimer in the documentation 
#    and/or other materials provided with the distribution.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
# IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
# INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, 
# OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY 
# OF SUCH DAMAGE.
# 

import sys
import unittest
import goldbach
sys.path.insert(0, '..\\primes\\')
import primes
import dataprocessing

#############################################################
# Unit tests
#############################################################

class TestMethods(unittest.TestCase):
    def test_isprime(self):
        p = primes.Primes(False)
        self.assertTrue(p.is_prime(2))
        self.assertTrue(p.is_prime(3))
        self.assertTrue(p.is_prime(5))
        self.assertTrue(p.is_prime(7))
        self.assertTrue(p.is_prime(11))

    def test_isprime_cuda(self):
        p = primes.Primes(False)
        self.assertTrue(p.is_prime_cuda(2))
        self.assertTrue(p.is_prime_cuda(3))
        self.assertTrue(p.is_prime_cuda(5))
        self.assertTrue(p.is_prime_cuda(7))
        self.assertTrue(p.is_prime_cuda(11))

    def test_isnotprime(self):
        p = primes.Primes(False)
        self.assertFalse(p.is_prime(1))
        self.assertFalse(p.is_prime(4))
        self.assertFalse(p.is_prime(6))
        self.assertFalse(p.is_prime(8))
        self.assertFalse(p.is_prime(10))
        self.assertFalse(p.is_prime(3379995))

    def test_isnotprime_cuda(self):
        p = primes.Primes(False)
        self.assertFalse(p.is_prime_cuda(1))
        self.assertFalse(p.is_prime_cuda(4))
        self.assertFalse(p.is_prime_cuda(6))
        self.assertFalse(p.is_prime_cuda(8))
        self.assertFalse(p.is_prime_cuda(10))
        self.assertFalse(p.is_prime_cuda(3379995))

    def test_get_ith_prime(self):
        p = primes.Primes(False)
        self.assertEqual(p.get_ith_prime(0), 2)
        self.assertEqual(p.get_ith_prime(1), 3)
        self.assertEqual(p.get_ith_prime(2), 5)

    def test_get_ith_twinprime(self):
        p = primes.Primes(False)
        self.assertEqual(p.get_ith_twinprime(0), 3)
        self.assertEqual(p.get_ith_twinprime(1), 5)
        self.assertEqual(p.get_ith_twinprime(2), 7)
        self.assertEqual(p.get_ith_twinprime(3), 11)
        self.assertEqual(p.get_ith_twinprime(4), 13)

    def test_search_for_partition_a1(self):
        p = primes.Primes(False)
        gp = goldbach.GoldbachPartition (p)
        (p1, p2, duration, iterations) = gp.search_for_partition (5, 7, lambda iteration: gp.delta_constant_minus_2(iteration))
        self.assertEqual (p1, 5)
        self.assertEqual (p2, 7)
        self.assertEqual (iterations, 1)
        (p1, p2, duration, iterations) = gp.search_for_partition (7, 9, lambda iteration: gp.delta_constant_minus_2(iteration))
        self.assertEqual (p1, 5)
        self.assertEqual (p2, 11)
        self.assertEqual (iterations, 2)

    def test_search_for_partition_a2(self):
        p = primes.Primes(False)
        gp = goldbach.GoldbachPartition (p)
        (p1, p2, duration, iterations) = gp.search_for_partition (3, 7, lambda iteration: gp.delta_constant_plus_2(iteration))
        self.assertEqual (p1, 3)
        self.assertEqual (p2, 7)
        self.assertEqual (iterations, 1)
        (p1, p2, duration, iterations) = gp.search_for_partition (3, 9, lambda iteration: gp.delta_constant_plus_2(iteration))
        self.assertEqual (p1, 5)
        self.assertEqual (p2, 7)
        self.assertEqual (iterations, 2)

    def test_search_for_partition_a6(self):
        p = primes.Primes(False)
        gp = goldbach.GoldbachPartition (p)
        p.add_to_primes_set(2)
        p.add_to_primes_set(3)
        p.add_to_primes_set(5)
        p.add_to_primes_set(7)
        p.sort_primes_set()
        (p1, p2, duration, iterations) = gp.search_for_partition (3, 5, lambda iteration: gp.delta_prime(iteration))
        self.assertEqual (p1, 3)
        self.assertEqual (p2, 5)
        self.assertEqual (iterations, 1)
        (p1, p2, duration, iterations) = gp.search_for_partition (3, 9, lambda iteration: gp.delta_prime(iteration))
        self.assertEqual (p1, 5)
        self.assertEqual (p2, 7)
        self.assertEqual (iterations, 2)

    def test_search_for_sym_primes_a1(self):
        p = primes.Primes(False)
        gp = goldbach.GoldbachPartition (p)
        (p1, p2, duration, iterations) = gp.search_for_sym_primes (2, lambda iteration: gp.delta_constant_minus_1(iteration))
        self.assertEqual (p1, 2)
        self.assertEqual (p2, 2)
        self.assertEqual (iterations, 1)
        (p1, p2, duration, iterations) = gp.search_for_sym_primes (4, lambda iteration: gp.delta_constant_minus_1(iteration))
        self.assertEqual (p1, 3)
        self.assertEqual (p2, 5)
        self.assertEqual (iterations, 2)

    def test_sum_of_prime_numbers(self):
        p = primes.Primes(False)
        gp = goldbach.GoldbachPartition (p)
        self.assertEqual(gp.find_sum_of_prime_numbers(4), [(2,2)])
        self.assertEqual(gp.find_sum_of_prime_numbers(6), [(3,3)])
        self.assertEqual(gp.find_sum_of_prime_numbers(8), [(3,5)])
        self.assertEqual(gp.find_sum_of_prime_numbers(10), [(3,7),(5,5)])
        self.assertEqual(gp.find_sum_of_prime_numbers(22), [(3,19),(5,17),(11,11)])

    def test_search_for_difference(self):
        p = primes.Primes(False)
        gp = goldbach.GoldbachPartition (p)
        p.sort_primes_set()
        (p1, p2, duration, iterations) = gp.search_for_difference (2)
        self.assertEqual (p1, 3)
        self.assertEqual (p2, 5)
        self.assertEqual (iterations, 1)
        (p1, p2, duration, iterations) = gp.search_for_difference (4)
        self.assertEqual (p1, 3)
        self.assertEqual (p2, 7)
        self.assertEqual (iterations, 1)
        (p1, p2, duration, iterations) = gp.search_for_difference (6)
        self.assertEqual (p1, 5)
        self.assertEqual (p2, 11)
        self.assertEqual (iterations, 2)

    def search_for_6kpm1_in_goldbach_neg(self):
        p = primes.Primes(False)
        gp = goldbach.GoldbachPartition (p)
        p.sort_primes_set()
        self.assertFalse (gp.check_for_6kpm1_in_partition (p, 4, 2, 2))
        self.assertFalse (gp.check_for_6kpm1_in_partition (p, 6, 3, 3))
        self.assertFalse (gp.check_for_6kpm1_in_partition (p, 8, 3, 5))
        self.assertFalse (gp.check_for_6kpm1_in_partition (p, 14, 3, 11))

    def search_for_6kpm1_in_goldbach(self):
        p = primes.Primes(False)
        gp = goldbach.GoldbachPartition (p)
        p.sort_primes_set()
        self.assertFalse (gp.check_for_6kpm1_in_partition (p, 10, 5, 5))
        self.assertFalse (gp.check_for_6kpm1_in_partition (p, 12, 5, 7))
        self.assertFalse (gp.check_for_6kpm1_in_partition (p, 14, 7, 7))

    def test_reduce_prime_for_goldbach_min(self):
        p = primes.Primes(False)
        gp = goldbach.GoldbachPartition (p)
        p.sort_primes_set()
        (n, q, res) = gp.reduce_prime_for_goldbach (5, False)
        self.assertEqual (n, 1)
        self.assertEqual (q, 3)
        self.assertTrue (res)
        (n, q, res) = gp.reduce_prime_for_goldbach (7, False)
        self.assertEqual (n, 1)
        self.assertEqual (q, 5)
        self.assertTrue (res)
        (n, q, res) = gp.reduce_prime_for_goldbach (11, False)
        self.assertEqual (n, 2)
        self.assertEqual (q, 7)
        self.assertTrue (res)

    def test_reduce_prime_for_goldbach_max(self):
        p = primes.Primes(False)
        gp = goldbach.GoldbachPartition (p)
        p.sort_primes_set()
        (n, q, res) = gp.reduce_prime_for_goldbach (5, True)
        self.assertEqual (n, 1)
        self.assertEqual (q, 3)
        self.assertTrue (res)
        (n, q, res) = gp.reduce_prime_for_goldbach (7, True)
        self.assertEqual (n, 2)
        self.assertEqual (q, 3)
        self.assertTrue (res)
        (n, q, res) = gp.reduce_prime_for_goldbach (11, True)
        self.assertEqual (n, 3)
        self.assertEqual (q, 3)
        self.assertTrue (res)

    def test_dictionary_cleanup(self):
        dp = dataprocessing.DataProcessing ()
        d = {0:2, 6:1, 2:3}
        print (d)
        res = dp.dictionary_cleanup(d)
        self.assertEqual(len(res), 4)
        print (res)

    def test_divide(self):
        dp = dataprocessing.DataProcessing ()
        res = dp.divide_list_into_chunks(range(10), 5)
        self.assertEqual(len(res), 2)
        self.assertEqual(res[0], range(0,5))
        self.assertEqual(res[1], range(5,10))

    def test_get_max_factor(self):
        dp = dataprocessing.DataProcessing ()
        self.assertEqual(dp.get_max_factor([(3,5)]), 5)
        self.assertEqual(dp.get_max_factor([(5,3)]), 5)
        self.assertEqual(dp.get_max_factor([(3,3)]), 3)

    def test_get_min_factor(self):
        dp = dataprocessing.DataProcessing ()
        self.assertEqual(dp.get_min_factor([(3,5)]), 3)
        self.assertEqual(dp.get_min_factor([(5,3)]), 3)
        self.assertEqual(dp.get_min_factor([(3,3)]), 3)
        self.assertEqual(dp.get_min_factor([(3,3),(5,7)]), 3)

    def test_get_avg_factor(self):
        dp = dataprocessing.DataProcessing ()
        self.assertEqual(dp.get_avg_factor([(3,5)]), 4)
        self.assertEqual(dp.get_avg_factor([(5,3)]), 4)
        self.assertEqual(dp.get_avg_factor([(3,3)]), 3)
        self.assertEqual(dp.get_avg_factor([(7,5),(2,11)]), 6.25)

    def test_get_max_ratio_in_factors(self):
        dp = dataprocessing.DataProcessing ()
        self.assertEqual(dp.get_max_ratio_in_factors([(3,5)]), 0.6)
        self.assertEqual(dp.get_max_ratio_in_factors([(3,7), (5,5)]), 1)
        self.assertEqual(dp.get_max_ratio_in_factors([(5,5), (3,7)]), 1)

    def test_get_perc_max_saturation_from_factors(self):
        dp = dataprocessing.DataProcessing ()
        self.assertEqual(dp.get_perc_max_saturation_from_factors(4, [(2,2)]), 50)
        self.assertEqual(dp.get_perc_max_saturation_from_factors(8, [(3,5)]), 62.5)

    def test_get_max_value_from_list(self):
        dp = dataprocessing.DataProcessing ()
        self.assertEqual(dp.get_max_value_from_list([1, 2, 3]), 3)
        self.assertEqual(dp.get_max_value_from_list([3, 2, 1]), 3)
        self.assertEqual(dp.get_max_value_from_list([1, 3, 2]), 3)
        self.assertEqual(dp.get_max_value_from_list([32, 32, 32]), 32)

    def test_get_min_value_from_list(self):
        dp = dataprocessing.DataProcessing ()
        self.assertEqual(dp.get_min_value_from_list([1, 2, 3]), 1)
        self.assertEqual(dp.get_min_value_from_list([2, 1, 3]), 1)
        self.assertEqual(dp.get_min_value_from_list([3, 2, 1]), 1)

    def test_get_avg_value_from_list(self):
        dp = dataprocessing.DataProcessing ()
        self.assertEqual(dp.get_avg_value_from_list([1, 2, 3]), 2)
        self.assertEqual(dp.get_avg_value_from_list([2, 1, 3]), 2)
        self.assertEqual(dp.get_avg_value_from_list([3, 2, 1]), 2)
        self.assertEqual(dp.get_avg_value_from_list([5]), 5)
        self.assertEqual(dp.get_avg_value_from_list([3, 2, 1, 10]), 4)

    def test_get_number_of_pairs(self):
        dp = dataprocessing.DataProcessing ()
        self.assertEqual(dp.get_number_of_pairs([(2,3)]), 1)
        self.assertEqual(dp.get_number_of_pairs([]), 0)
        self.assertEqual(dp.get_number_of_pairs([(2,3), (5,7)]), 2)

    def test_get_diff_in_factors(self):
        dp = dataprocessing.DataProcessing ()
        self.assertEqual(dp.get_diff_in_factors([(2,3)]), [1])
        self.assertEqual(dp.get_diff_in_factors([(2,5),(2,3)]), [3,1])
        self.assertEqual(dp.get_diff_in_factors([(2,5),(3,2)]), [3,-1])

    def test_read_num_from_line (self):
        dp = dataprocessing.DataProcessing ()
        self.assertEqual(dp.read_num_from_line ("Number: 10"), 10)
        self.assertEqual(dp.read_num_from_line ("Number: 34453432"), 34453432)

    def test_read_sums_from_line (self):
        dp = dataprocessing.DataProcessing ()
        self.assertEqual(dp.read_sums_from_line (" Pairs: [(3, 7), (5, 5)]"), [(3, 7), (5, 5)])
        self.assertEqual(dp.read_sums_from_line (" Pairs: [(7, 83), (11, 79), (17, 73), (19, 71), (23, 67), (29, 61), (31, 59), (37, 53), (43, 47)]"), [(7, 83), (11, 79), (17, 73), (19, 71), (23, 67), (29, 61), (31, 59), (37, 53), (43, 47)])

#############################################################
# Main - run unit tests
#############################################################

unittest.main()
