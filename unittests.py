#
# Copyright (c) 2016 - 2017, Marcin Barylski
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
        p.add_to_prime_set(2)
        p.add_to_prime_set(3)
        p.add_to_prime_set(5)
        p.add_to_prime_set(7)
        p.sort_prime_set()
        self.assertEqual(p.get_ith_prime(0), 2)
        self.assertEqual(p.get_ith_prime(1), 3)
        self.assertEqual(p.get_ith_prime(2), 5)

    def test_search_for_partition_a1(self):
        p = primes.Primes(False)
        gp = goldbach.GoldbachPartition (p)
        (p1, p2, duration, iterations) = gp.search_for_partition (5, 7, lambda iteration: gp.delta_constant_minus(iteration))
        self.assertEqual (p1, 5)
        self.assertEqual (p2, 7)
        self.assertEqual (iterations, 1)
        (p1, p2, duration, iterations) = gp.search_for_partition (7, 9, lambda iteration: gp.delta_constant_minus(iteration))
        self.assertEqual (p1, 5)
        self.assertEqual (p2, 11)
        self.assertEqual (iterations, 2)

    def test_search_for_partition_a2(self):
        p = primes.Primes(False)
        gp = goldbach.GoldbachPartition (p)
        (p1, p2, duration, iterations) = gp.search_for_partition (3, 7, lambda iteration: gp.delta_constant_plus(iteration))
        self.assertEqual (p1, 3)
        self.assertEqual (p2, 7)
        self.assertEqual (iterations, 1)
        (p1, p2, duration, iterations) = gp.search_for_partition (3, 9, lambda iteration: gp.delta_constant_plus(iteration))
        self.assertEqual (p1, 5)
        self.assertEqual (p2, 7)
        self.assertEqual (iterations, 2)

    def test_search_for_partition_a6(self):
        p = primes.Primes(False)
        gp = goldbach.GoldbachPartition (p)
        p.add_to_prime_set(2)
        p.add_to_prime_set(3)
        p.add_to_prime_set(5)
        p.add_to_prime_set(7)
        p.sort_prime_set()
        (p1, p2, duration, iterations) = gp.search_for_partition (3, 5, lambda iteration: gp.delta_prime(iteration))
        self.assertEqual (p1, 3)
        self.assertEqual (p2, 5)
        self.assertEqual (iterations, 1)
        (p1, p2, duration, iterations) = gp.search_for_partition (3, 9, lambda iteration: gp.delta_prime(iteration))
        self.assertEqual (p1, 5)
        self.assertEqual (p2, 7)
        self.assertEqual (iterations, 2)

    def test_sumofprimenumbers(self):
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
        p.sort_prime_set()
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

    def test_divide(self):
        p = primes.Primes(False)
        gp = goldbach.GoldbachPartition (p)
        res = gp.divide_into_chunks(range(10), 5)
        self.assertEqual(len(res), 2)
        self.assertEqual(res[0], range(0,5))
        self.assertEqual(res[1], range(5,10))

#############################################################
# Main - run unit tests
#############################################################

unittest.main()
