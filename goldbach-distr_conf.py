#
# Copyright (c) 2016, Marcin Barylski
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

import math
import unittest
import os
import time
import numpy as np

#############################################################
# Settings - configuration
#############################################################

# Run unit tests
#   o True  - run unit tests
#   o False - run main program
run_unit_tests = False

# Minimal even number checked against Goldbach conjecture
minimum_n = 6
# Maximum even number checked against Goldbach conjecture
maximum_n = 10000000
# Chunk size for distributed computation
max_chunk_size = 100000

# Caching previous primality results
#   o True  - auxilary sets of primes and composite numbers will grow
#             it will speed up further primality tests but more RAM will
#             be occupied
#   o False - do not cache new primality test results
caching_primality_results = False

# Helper files
#   o file_input_primes - contains prime numbers
#   o file_input_nonprimes - contains composite numbers
file_input_primes = 't_prime_numbers.txt'
file_input_nonprimes = 't_nonprime_numbers.txt'

#############################################################
# Business logic
#############################################################

set_prime = set ()
set_nonprime = set ()

list_sorted_prime = []

def init_set (filename, is_prime):
    global list_sorted_prime
    if os.path.exists(filename):
        f = open(filename, "r")
        lines = f.readlines()
        for line in lines:
            line = line.replace('[', '')
            line = line.replace(']', '')
            numbers = line.split(',')
            for number in numbers:
                if is_prime:
                    add_to_prime_set(int(number))
                else:
                    add_to_nonprime_set(int(number))

def is_in_prime_set (n):
    global set_prime
    if n in set_prime:
        return True
    else:
        return False

def is_in_nonprime_set (n):
    global set_nonprime
    if n in set_nonprime:
        return True
    else:
        return False

def add_to_prime_set (n):
    global set_prime
    set_prime.add(n)

def add_to_nonprime_set (n):
    global set_nonprime
    set_nonprime.add(n)

def is_prime (n):
    if n == 1:
        return False
    elif n == 2 or n == 3:
        return True
    elif is_in_prime_set (n):
        return True
    elif is_in_nonprime_set (n):
        return False
    elif n % 2 == 0 or n % 3 == 0:
        return False
    result = True
    i = 5
    while i*i <= n:
        if n %  i == 0 or n % (i + 2) == 0:
            result = False
            break
        i += 6

    if caching_primality_results:
        if result:
            add_to_prime_set (n)
        else:
            add_to_nonprime_set (n)
    return result

def get_ith_prime (i):
    global list_sorted_prime
    return list_sorted_prime[i]

def delta_prime (iteration):
    if iteration == 0:
        delta = 0
    else:
        delta = get_ith_prime(iteration + 1) - get_ith_prime(iteration)
    return delta

def search_for_partition (p1, p2, delta):
    found = False
    iteration = 0

    startTime = time.time()
    while (not found):
        iteration += 1
        if (is_prime (p1) and is_prime(p2)):
            found = True
        if (not found):
            p1 = p1 + delta (iteration)
            p2 = p2 - delta (iteration)
        if p2 < 2 or p1 < 2:
            raise ("Could not find GP")
    duration = time.time() - startTime
    return p1, p2, duration, iteration

def verify (chunk):
    for n in chunk:
        p1 = 3
        p2 = n - 3
        (p1, p2, d, i) = search_for_partition (p1, p2, lambda i: delta_prime(i))
      
def divide (r, size):
    out = []
    last = 0
    while last < len(r):
        out.append (r[int(last):int(last + size)])
        last += size
    return out

#############################################################
# Unit tests
#############################################################

class TestPrimeMethods(unittest.TestCase):
    def test_isprime(self):
        self.assertTrue(is_prime(2))
        self.assertTrue(is_prime(3))
        self.assertTrue(is_prime(5))
        self.assertTrue(is_prime(7))
        self.assertTrue(is_prime(11))

    def test_isnotprime(self):
        self.assertFalse(is_prime(1))
        self.assertFalse(is_prime(4))
        self.assertFalse(is_prime(6))
        self.assertFalse(is_prime(8))
        self.assertFalse(is_prime(10))
        self.assertFalse(is_prime(3379995))

    def test_get_ith_prime(self):
        global set_prime
        global list_sorted_prime
        set_prime.add (2)
        set_prime.add (3)
        set_prime.add (5)
        set_prime.add (7)
        list_sorted_prime = sorted (set_prime)
        self.assertEqual(get_ith_prime(0), 2)
        self.assertEqual(get_ith_prime(1), 3)
        self.assertEqual(get_ith_prime(2), 5)

    def test_search_for_partition_a6(self):
        global set_prime
        global list_sorted_prime
        set_prime.add (2)
        set_prime.add (3)
        set_prime.add (5)
        set_prime.add (7)
        list_sorted_prime = sorted (set_prime)
        (p1, p2, duration, iterations) = search_for_partition (3, 5, lambda iteration: delta_prime(iteration))
        self.assertEqual (p1, 3)
        self.assertEqual (p2, 5)
        self.assertEqual (iterations, 1)
        (p1, p2, duration, iterations) = search_for_partition (3, 9, lambda iteration: delta_prime(iteration))
        self.assertEqual (p1, 5)
        self.assertEqual (p2, 7)
        self.assertEqual (iterations, 2)

    def test_divide(self):
        res = divide(range(10), 5)
        self.assertEqual(len(res), 2)
        self.assertEqual(res[0], range(0,5))
        self.assertEqual(res[1], range(5,10))

#############################################################
# Main - unit tests
#############################################################

if run_unit_tests:
    unittest.main()

#############################################################
# Main
#############################################################

print ("Loading helper sets...")
init_set(file_input_primes, True)
init_set(file_input_nonprimes, False)
print ("DONE")

print ("Sorting primes...")
list_sorted_prime = sorted(list(set_prime))
print ("DONE")

print ("Verification for all even numbers from", minimum_n, "to", maximum_n, "started ...")
i = 1
n_of_chunks = (maximum_n - minimum_n) / 2 / max_chunk_size
chunks = divide (range(minimum_n, maximum_n, 2), max_chunk_size)
for chunk in chunks:
    verify (chunk)
    perc = int (100 * i / n_of_chunks)
    print (" Chunk #", i, ":", chunk, "verified (", perc, "% completed )")
    i+= 1
print ("DONE")
