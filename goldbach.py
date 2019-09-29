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

import math
import unittest
import sys
import numpy as np
import os
from datetime import datetime
import time
sys.path.insert(0, '..\\primes\\')
import primes

class GoldbachPartition:

    # object for carrying out operations on prime numbers
    primes = ""

    # currently analyzed number
    number = 0

    def __init__(self, p):
        self.primes = p
        self.number = 0

    def set_number (self, n):
        self.number = n

    def add_to_primes_set_to_be_excluded (self, p1):
        self.primes.add_to_primes_set_to_be_excluded(p1)
    
    def delta_constant_plus_2 (self, iteration):
        return 2

    def delta_constant_minus_2 (self, iteration):
        return -2

    def delta_constant_minus_1 (self, iteration):
        return -1

    def delta_variable (self, iteration):
        if iteration == 0:
            delta = 0
        elif iteration == 1:
            delta = 2
        elif iteration == 2:
            delta = -4
        elif iteration > 2:
            delta = 2
        return delta

    def delta_variable_3x (self, iteration):
        if iteration < 2:
            delta = 0
        else:
            delta = 3*iteration
        return delta

    def delta_variable_3xplus1 (self, iteration):
        if iteration < 2:
            delta = 0
        else:
            delta = 3*iteration + 1
        return delta

    def delta_variable_3xminus1 (self, iteration):
        if iteration < 2:
            delta = 0
        else:
            delta = 3*iteration - 1
        return delta

    # skip first prime: 2
    def delta_prime (self, iteration):
        if iteration == 0:
            delta = 0
        else:
            delta = self.primes.get_ith_prime(iteration + 1) - self.primes.get_ith_prime(iteration)
        return delta

    def delta_twinprime (self, iteration):
        if iteration == 0:
            delta = 0
        else:
            delta = self.primes.get_ith_twinprime(iteration) - self.primes.get_ith_twinprime(iteration - 1)
        return delta

    def delta_3kmin1 (self, iteration):
        delta= self.number - self.number - 1 - 3*iteration
        if iteration == 0:
            delta = 0
        else:
            delta = self.primes.get_ith_prime(iteration + 1) - self.primes.get_ith_prime(iteration)
        return delta

    def search_for_partition (self, p1, p2, delta):
        found = False
        iteration = 0

        startTime = time.time()
        while not found:
            iteration += 1
            if (self.primes.is_prime (p1) and self.primes.is_prime (p2)):
                found = True
            if not found:
                p1 = p1 + delta (iteration)
                p2 = p2 - delta (iteration)
            if p2 < 2 or p1 < 2:
                raise Exception ("Could not find GP for ", p1 + p2)
        duration = time.time() - startTime
        return p1, p2, duration, iteration

    def search_for_sym_primes (self, num, delta):
        found = False
        iteration = 0

        startTime = time.time()
        p1 = num
        p2 = num
        while (not found):
            iteration += 1
            if (self.primes.is_prime (p1) and self.primes.is_prime (p2)):
                found = True
            if (not found):
                p1 = p1 + delta (iteration)
                p2 = p2 - delta (iteration)
            if p2 < 2 or p1 < 2:
                raise Exception ("Could not find symmetrical primes")
        duration = time.time() - startTime
        return p1, p2, duration, iteration

    def find_sum_of_prime_numbers (self, n):
        min_i = 2
        max_i = int(n / 2) + 1
        factors = []
        for i in range (min_i, max_i):
            p1 = i
            p2 = n - p1
            if self.primes.is_prime(p1) and self.primes.is_prime (p2):
                pair = (p1, p2)
                factors.append (pair)
        return factors

    def search_for_difference (self, num):
        found = False
        iteration = 0
    
        startTime = time.time()
        while (not found):
            iteration += 1
            p1 = self.primes.get_ith_prime (iteration) # skip 2; start from 3
            p2 = p1 + num
            if self.primes.is_prime (p2):
                found = True
        duration = time.time() - startTime
        return p1, p2, duration, iteration

    # returns True if 6kpm1 hypothesis is fulfilled for num and its
    # Goldbach Partition p1+p2:
    # - if n mod 3 = 0, then GSC (2n, 6a-1, 6b+1)
    # - if n mod 3 = 1, then GSC (2n, 6a+1, 6b+1)
    # - if n mod 3 = 2, then GSC (2n, 6a-1, 6b-1)
    # Otherwise returns False.
    def check_for_6kpm1_in_partition (self, p, num, p1, p2):
        if num/2 % 3 == 0:
            if (p.is_6km1 (p1) and p.is_6kp1 (p2)) or (p.is_6km1 (p2) and p.is_6kp1 (p1)):
                return True
        elif num/2 % 3 == 1:
            if (p.is_6kp1 (p1) and p.is_6kp1 (p2)):
                return True
        else:
            if (p.is_6km1 (p1) and p.is_6km1 (p2)):
                return True
        return False

    def reduce_prime_for_goldbach (self, p, look_for_max):
        if look_for_max:
            n = int(math.log(p, 2))
            diff = -1
            min_num = 0
        else:
            n = 1
            diff = 1
            min_num = 2
        q = p
        found = False
        while (q > min_num):
            q = p - 2**n
            if self.primes.is_prime(q):
                found = True
                break
            else:
                n += diff
            if n < 1:
                break
        return (n, q, found)
