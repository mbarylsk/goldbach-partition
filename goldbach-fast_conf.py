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
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.mlab as mlab
import numpy as np
import os
from datetime import datetime
import time

#############################################################
# Settings - configuration
#############################################################

run_unit_tests = False
be_verbose = False

min_num = 7
max_num = 200000000
step_factor = 2
checkpoint_value = 200000
big_prime_threshold = 300
min_prime_count_threshold = 500
file_input_primes = 't_prime_numbers.txt'
file_input_nonprimes = 't_nonprime_numbers.txt'

#############################################################
# Settings - output directory and files
#############################################################

if not run_unit_tests:
	directory = str(step_factor*max_num)
	if not os.path.exists(directory):
		os.makedirs(directory)

#############################################################
# Business logic
#############################################################

set_prime = set ()
set_nonprime = set ()

list_nums = []
list_sorted_prime = []

list_checkpoints_duration_a2 = []
list_checkpoints_duration_a6 = []
list_checkpoints_iters_a2 = []
list_checkpoints_iters_a6 = []
list_checkpoints = []

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
                    set_prime.add(int(number))
                else:
                    set_nonprime.add(int(number))

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

    if result:
        add_to_prime_set (n)
    else:
        add_to_nonprime_set (n)
    return result

def find_sum_of_prime_numbers (n):
    min_i = 2
    max_i = int(n / 2) + 1
    factors = []
    for i in range (min_i, max_i):
        p1 = i
        p2 = n - p1
        if is_prime(p1) and is_prime (p2):
            pair = (p1, p2)
            factors.append (pair)
    return factors

def get_diff_in_factors (factors):
    diffs = []
    for pair in factors:
        (p1, p2) = pair
        diff = p2 - p1
        diffs.append(diff)
    return diffs

def get_max_factor (factors):
    maxv = 0
    for pair in factors:
        (p1, p2) = pair
        if (p1 > maxv):
            maxv = p1
        if (p2 > maxv):
            maxv = p2
    return maxv

def get_avg_factor (factors):
    sumv = 0
    count = 0
    for pair in factors:
        (p1, p2) = pair
        count += 2
        sumv += p1 + p2
    avgv = sumv / count
    return avgv

def get_min_factor (factors):
    minv = sys.maxsize
    for pair in factors:
        (p1, p2) = pair
        if (p1 < minv):
            minv = p1
        if (p2 < minv):
            minv = p2
    return minv

def get_number_of_pairs (factors):
    return len(factors)

def get_perc_max_saturation_from_factors (n, factors):
    max_p = get_max_factor (factors)
    return max_p / n * 100

def get_perc_min_saturation_from_factors (n, factors):
    min_p = get_min_factor (factors)
    return min_p / n * 100

def get_perc_avg_saturation_from_factors (n, factors):
    avg_p = get_avg_factor (factors)
    return avg_p / n * 100

def get_max_value_from_list (l):
    return max(l)

def get_min_value_from_list (l):
    return min(l)

def get_avg_value_from_list (l):
    return (sum(l)/len(l))

def get_ith_prime (i):
    global list_sorted_prime
    return list_sorted_prime[i]

def search_for_partition (p1, p2):
    found = False
    iterations = 0
    
    startTime = time.time()
    while (not found):
        iterations += 1
        if (is_prime (p1) and is_prime(p2)):
            found = True
        if (not found):
                p2 = p2 - 2
                p1 = p1 + 2
        if p2 < 2 or p1 < 2:
                raise ("Could not find GP")
    duration = time.time() - startTime
    return p1, p2, duration, iterations

def search_for_partition2 (p1, p2):
    found = False
    iterations = 0

    startTime = time.time()
    while (not found):
        if iterations == 0:
                delta = 0
        elif iterations == 1:
                delta = -2
        elif iterations == 2:
                delta = 4
        elif iterations > 2:
                delta = 2
                
        p2 = p2 - delta
        p1 = p1 + delta

        iterations += 1
        if (is_prime (p1) and is_prime(p2)):
            found = True

        if p1 < 2 or p2 < 2:
                raise ("Could not find GP")
            
    duration = time.time() - startTime
    return p1, p2, duration, iterations

def search_for_partition3 (p1, p2):
    found = False
    iterations = 0
    n = p1 + p2
    startTime = time.time()
    while (not found):
        iterations += 1
        if (is_prime (p1) and is_prime(p2)):
            found = True
        if (not found):
            p1 = get_ith_prime (iterations)
            p2 = n - p1
        if p1 < 2 or p2 < 2:
            raise ("Could not find GP")
    duration = time.time() - startTime
    return p1, p2, duration, iterations

#############################################################
# Presentation
#############################################################

def write_results_to_figures (directory, last_loop):
    # results - figures    
    plt.figure(1)
    plt.plot(list_checkpoints, list_checkpoints_duration_a6, 'b.', ms=2)
    plt.plot(list_checkpoints, list_checkpoints_duration_a2, 'r.', ms=2)
    red_patch = mpatches.Patch(color='red', label='A2')
    blue_patch = mpatches.Patch(color='blue', label='A6')
    plt.legend(handles=[red_patch, blue_patch], loc='upper right', bbox_to_anchor=(0.4, 0.8))
    plt.xlabel('Number')
    plt.ylabel('Time [s]')
    plt.title('Duration of total calculations')
    plt.grid(True)
    plt.savefig(directory + "/f_checkpoint_duration_alg.png")

    plt.figure(2)
    plt.plot(list_checkpoints, list_checkpoints_iters_a6, 'b.', ms=2)
    plt.plot(list_checkpoints, list_checkpoints_iters_a2, 'r.', ms=2)
    red_patch = mpatches.Patch(color='red', label='A2')
    blue_patch = mpatches.Patch(color='blue', label='A6')
    plt.legend(handles=[red_patch, blue_patch], loc='upper right', bbox_to_anchor=(0.4, 0.8))
    plt.xlabel('Number')
    plt.ylabel('Iterations')
    plt.title('Total iterations')
    plt.grid(True)
    plt.savefig(directory + "/f_checkpoint_iters_alg.png")

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

    def test_sumofprimenumbers(self):
        self.assertEqual(find_sum_of_prime_numbers(4), [(2,2)])
        self.assertEqual(find_sum_of_prime_numbers(6), [(3,3)])
        self.assertEqual(find_sum_of_prime_numbers(8), [(3,5)])
        self.assertEqual(find_sum_of_prime_numbers(10), [(3,7),(5,5)])
        self.assertEqual(find_sum_of_prime_numbers(22), [(3,19),(5,17),(11,11)])

    def test_get_max_factor(self):
        self.assertEqual(get_max_factor([(3,5)]), 5)
        self.assertEqual(get_max_factor([(5,3)]), 5)
        self.assertEqual(get_max_factor([(3,3)]), 3)

    def test_get_min_factor(self):
        self.assertEqual(get_min_factor([(3,5)]), 3)
        self.assertEqual(get_min_factor([(5,3)]), 3)
        self.assertEqual(get_min_factor([(3,3)]), 3)
        self.assertEqual(get_min_factor([(3,3),(5,7)]), 3)

    def test_get_avg_factor(self):
        self.assertEqual(get_avg_factor([(3,5)]), 4)
        self.assertEqual(get_avg_factor([(5,3)]), 4)
        self.assertEqual(get_avg_factor([(3,3)]), 3)
        self.assertEqual(get_avg_factor([(7,5),(2,11)]), 6.25)

    def test_get_perc_max_saturation_from_factors(self):
        self.assertEqual(get_perc_max_saturation_from_factors(4, [(2,2)]), 50)
        self.assertEqual(get_perc_max_saturation_from_factors(8, [(3,5)]), 62.5)

    def test_get_max_value_from_list(self):
        self.assertEqual(get_max_value_from_list([1, 2, 3]), 3)
        self.assertEqual(get_max_value_from_list([3, 2, 1]), 3)
        self.assertEqual(get_max_value_from_list([1, 3, 2]), 3)
        self.assertEqual(get_max_value_from_list([32, 32, 32]), 32)

    def test_get_min_value_from_list(self):
        self.assertEqual(get_min_value_from_list([1, 2, 3]), 1)
        self.assertEqual(get_min_value_from_list([2, 1, 3]), 1)
        self.assertEqual(get_min_value_from_list([3, 2, 1]), 1)

    def test_get_avg_value_from_list(self):
        self.assertEqual(get_avg_value_from_list([1, 2, 3]), 2)
        self.assertEqual(get_avg_value_from_list([2, 1, 3]), 2)
        self.assertEqual(get_avg_value_from_list([3, 2, 1]), 2)
        self.assertEqual(get_avg_value_from_list([5]), 5)
        self.assertEqual(get_avg_value_from_list([3, 2, 1, 10]), 4)

    def test_get_number_of_pairs(self):
        self.assertEqual(get_number_of_pairs([(2,3)]), 1)
        self.assertEqual(get_number_of_pairs([]), 0)
        self.assertEqual(get_number_of_pairs([(2,3), (5,7)]), 2)

    def test_get_diff_in_factors(self):
        self.assertEqual(get_diff_in_factors([(2,3)]), [1])
        self.assertEqual(get_diff_in_factors([(2,5),(2,3)]), [3,1])
        self.assertEqual(get_diff_in_factors([(2,5),(3,2)]), [3,-1])

    def test_get_ith_prime(self):
        self.assertEqual(get_ith_prime(0), 2)
        self.assertEqual(get_ith_prime(1), 3)
        self.assertEqual(get_ith_prime(2), 5)

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

dt_start = datetime.now()
dt_current_previous = dt_start

dt_diff2 = 0
dt_diff6 = 0

dt_iter2 = 0
dt_iter6 = 0

# new calculations
for k in range (min_num, max_num):
    num = step_factor*k

    list_nums.append (num)

    # algorithm 2
    # start from the biggest possible difference between primes
    p1 = 3
    p2 = num - 3

    (p1, p2, duration, iterations) = search_for_partition (p1, p2)
    dt_diff2 += duration
    dt_iter2 += iterations

    if p1 + p2 != num:
        print ("Alg #2: violation of sum for p1=", p1, "p2=", p2, "n=", num)


    # algorithm 6
    # A2 but next p1 is always prime
    p1 = 3
    p2 = num - 3

    (p1, p2, duration, iterations) = search_for_partition3 (p1, p2)
    dt_diff6 += duration
    dt_iter6 += iterations

    if p1 + p2 != num:
        print ("Alg #6: violation of sum for p1=", p1, "p2=", p2, "n=", num)

    
    # checkpoint - partial results
    if num % checkpoint_value == 0:
        dt_current = datetime.now()
        dt_diff_current = (dt_current - dt_current_previous).total_seconds()
        list_checkpoints.append(num)
        list_checkpoints_duration_a2.append(dt_diff2)
        list_checkpoints_duration_a6.append(dt_diff6)
        list_checkpoints_iters_a2.append(dt_iter2)
        list_checkpoints_iters_a6.append(dt_iter6)
        print ("Checkpoint", k, "of total", max_num, "took", dt_diff_current, "seconds")
        
        # remember results so far
        write_results_to_figures (directory, False)

dt_end = datetime.now()

# final results - time of processing
dt_diff = dt_end - dt_start
print ("Total calculations lasted:", dt_diff)

# final results - figures
write_results_to_figures (directory, True)
