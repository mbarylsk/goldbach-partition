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
import os
import numpy
import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule

#############################################################
# OpenCL kernels
#############################################################

kernel = SourceModule("""
    __global__ void first_factor(long long n, int *value, int *primes){
        int i = threadIdx.x;
        int p = primes[i];
        if(n % p == 0){
            value[i] = p;
        }else{
            value[i] = 0;
        }
    }
""")

#############################################################
# Class
#############################################################

class Primes:

    set_prime = set ()
    set_nonprime = set ()
    list_sorted_prime = []

    # Caching previous primality results
    #   o True  - auxilary sets of primes and composite numbers will grow
    #             it will speed up further primality tests but more RAM will
    #             be occupied
    #   o False - do not cache new primality test results
    caching_primality_results = False

    def __init__(self, cache_results):
        self.caching_primality_results = cache_results

    def init_set (self, filename, is_prime):
        if os.path.exists(filename):
            f = open(filename, "r")
            lines = f.readlines()
            for line in lines:
                line = line.replace('[', '')
                line = line.replace(']', '')
                numbers = line.split(',')
                for number in numbers:
                    if is_prime:
                        self.add_to_prime_set(int(number))
                    else:
                        self.add_to_nonprime_set(int(number))

    def is_in_prime_set (self, n):
        if n in self.set_prime:
            return True
        else:
            return False

    def is_in_nonprime_set (self, n):
        if n in self.set_nonprime:
            return True
        else:
            return False

    def add_to_prime_set (self, n):
        self.set_prime.add(n)

    def sort_prime_set (self):
        self.list_sorted_prime = sorted (self.set_prime)

    def add_to_nonprime_set (self, n):
        self.set_nonprime.add(n)

    def is_prime (self, n):
        if n == 1:
            return False
        elif n == 2 or n == 3:
            return True
        elif self.is_in_prime_set (n):
            return True
        elif self.is_in_nonprime_set (n):
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

        if self.caching_primality_results:
            if result:
                self.add_to_prime_set (n)
            else:
                self.add_to_nonprime_set (n)
        return result

    def is_prime_cuda (self, n):
        def min2(list, bound=0):
            for item in list:
                if item > bound:
                    return item
            return None

        first_factor = kernel.get_function('first_factor')

        if n == 1:
            return False
        if n == 2 or n == 3:
            return True

        allPrimes = [2,3]
        numPrimes = len(allPrimes)
        numThreads = 384

        result = numpy.zeros(numPrimes, numpy.int32)
        primes = numpy.copy(allPrimes).astype(numpy.int32)

        while True:
            first_factor(numpy.int64(n), cuda.InOut(result), cuda.In(primes), block=(numThreads, 1, 1))

            prime = min2(result, 1)
            if prime == None:
                break
            return False
        return True

    def get_ith_prime (self, i):
        return self.list_sorted_prime[i]

    
