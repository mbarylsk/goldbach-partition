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

# Run unit tests
#   o True  - run unit tests
#   o False - run main program
run_unit_tests = False

# Minimal even number checked against Goldbach conjecture
#   o number = min_num * step_factor
min_num = 7
step_factor = 2
# Maximum even number checked against Goldbach conjecture
#   o number = max_num * step_factor
max_num = 200000

# Checkpoint value when partial results are drawn/displayed
checkpoint_value = 2000

# Caching previous primality results
#   o True  - auxilary sets of primes and composite numbers will grow
#             it will speed up further primality tests but more RAM will
#             be occupied
#   o False - do not cache new primality test results
caching_primality_results = False

# Algorithms to be checked
algo_to_check = {'a1', 'a2', 'a3', 'a4', 'a5', 'a6'}

# Helper files
#   o file_input_primes - contains prime numbers
#   o file_input_nonprimes - contains composite numbers
file_input_primes = 't_prime_numbers.txt'
file_input_nonprimes = 't_nonprime_numbers.txt'

#############################################################
# Settings - output directory and files
#############################################################

if not run_unit_tests:
    directory = str(step_factor*max_num)
    if not os.path.exists(directory):
        os.makedirs(directory)
    algo = ""
    if 'a1' in algo_to_check:
        algo += "a1"
    if 'a2' in algo_to_check:
        algo += "a2"
    if 'a3' in algo_to_check:
        algo += "a3"
    if 'a4' in algo_to_check:
        algo += "a4"
    if 'a5' in algo_to_check:
        algo += "a5"
    if 'a6' in algo_to_check:
        algo += "a6"
    file_output_iters_alg = directory + "/f_checkpoint_iters_alg" + algo + ".png"
    file_output_duration_alg = directory + "/f_checkpoint_duration_alg" + algo + ".png"

#############################################################
# Business logic
#############################################################

set_prime = set ()
set_nonprime = set ()

list_nums = []
list_sorted_prime = []

dt_diff = [0, 0, 0, 0, 0, 0]
dt_iter = [0, 0, 0, 0, 0, 0]

list_checkpoints_duration = [[], [], [], [], [], []]
list_checkpoints_iters = [[], [], [], [], [], []]
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

def delta_constant_plus (iteration):
    return 2

def delta_constant_minus (iteration):
    return -2

def delta_variable (iteration):
    if iteration == 0:
        delta = 0
    elif iteration == 1:
        delta = 2
    elif iteration == 2:
        delta = -4
    elif iteration > 2:
        delta = 2
    return delta

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

#############################################################
# Presentation
#############################################################

def write_results_to_figures (directory):
    # results - figures    
    plt.figure(1)
    g_patch = mpatches.Patch(color='green', label='A1')
    r_patch = mpatches.Patch(color='red', label='A2')
    b_patch = mpatches.Patch(color='blue', label='A3')
    y_patch = mpatches.Patch(color='yellow', label='A4')
    c_patch = mpatches.Patch(color='cyan', label='A5')
    m_patch = mpatches.Patch(color='magenta', label='A6')

    list_of_handles = []
    if 'a1' in algo_to_check:
        plt.plot(list_checkpoints, list_checkpoints_duration[0], 'g.', ms=2)
        list_of_handles.append(g_patch)
    if 'a2' in algo_to_check:
        plt.plot(list_checkpoints, list_checkpoints_duration[1], 'r.', ms=2)
        list_of_handles.append(r_patch)
    if 'a3' in algo_to_check:
        plt.plot(list_checkpoints, list_checkpoints_duration[2], 'b.', ms=2)
        list_of_handles.append(b_patch)
    if 'a4' in algo_to_check:
        plt.plot(list_checkpoints, list_checkpoints_duration[3], 'y.', ms=2)
        list_of_handles.append(y_patch)
    if 'a5' in algo_to_check:
        plt.plot(list_checkpoints, list_checkpoints_duration[4], 'c.', ms=2)
        list_of_handles.append(c_patch)
    if 'a6' in algo_to_check:
        plt.plot(list_checkpoints, list_checkpoints_duration[5], 'm.', ms=2)
        list_of_handles.append(m_patch)

    plt.legend(handles=list_of_handles, loc='upper right', bbox_to_anchor=(0.4, 0.8))
    plt.xlabel('Number')
    plt.ylabel('Time [s]')
    plt.title('Duration of total calculations')
    plt.grid(True)
    plt.savefig(file_output_duration_alg)

    plt.figure(2)
    if 'a1' in algo_to_check:
        plt.plot(list_checkpoints, list_checkpoints_iters[0], 'g.', ms=2)
    if 'a2' in algo_to_check:
        plt.plot(list_checkpoints, list_checkpoints_iters[1], 'r.', ms=2)
    if 'a3' in algo_to_check:
        plt.plot(list_checkpoints, list_checkpoints_iters[2], 'b.', ms=2)
    if 'a4' in algo_to_check:
        plt.plot(list_checkpoints, list_checkpoints_iters[3], 'y.', ms=2)
    if 'a5' in algo_to_check:
        plt.plot(list_checkpoints, list_checkpoints_iters[4], 'c.', ms=2)
    if 'a6' in algo_to_check:
        plt.plot(list_checkpoints, list_checkpoints_iters[5], 'm.', ms=2)

    plt.legend(handles=list_of_handles, loc='upper right', bbox_to_anchor=(0.4, 0.8))
    plt.xlabel('Number')
    plt.ylabel('Iterations')
    plt.title('Total iterations')
    plt.grid(True)
    plt.savefig(file_output_iters_alg)

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

    def test_search_for_partition_a1(self):
        (p1, p2, duration, iterations) = search_for_partition (5, 7, lambda iteration: delta_constant_minus(iteration))
        self.assertEqual (p1, 5)
        self.assertEqual (p2, 7)
        self.assertEqual (iterations, 1)
        (p1, p2, duration, iterations) = search_for_partition (7, 9, lambda iteration: delta_constant_minus(iteration))
        self.assertEqual (p1, 5)
        self.assertEqual (p2, 11)
        self.assertEqual (iterations, 2)

    def test_search_for_partition_a2(self):
        (p1, p2, duration, iterations) = search_for_partition (3, 7, lambda iteration: delta_constant_plus(iteration))
        self.assertEqual (p1, 3)
        self.assertEqual (p2, 7)
        self.assertEqual (iterations, 1)
        (p1, p2, duration, iterations) = search_for_partition (3, 9, lambda iteration: delta_constant_plus(iteration))
        self.assertEqual (p1, 5)
        self.assertEqual (p2, 7)
        self.assertEqual (iterations, 2)

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

# new calculations
for k in range (min_num, max_num):
    num = step_factor*k

    list_nums.append (num)

    # algorithm 1
    # start from half of n - the smallest possible difference between primes
    if 'a1' in algo_to_check:
        p1 = num / 2
        if p1 % 2 == 0:
            p1 = p1 - 1
            p2 = p1 + 2
        else:
            p2 = p1

        (p1, p2, duration, iterations) = search_for_partition (p1, p2, lambda iteration: delta_constant_minus(iteration))
        dt_diff[0] += duration
        dt_iter[0] += iterations

        if p1 + p2 != num:
            print ("Alg #1: violation of sum for p1=", p1, "p2=", p2, "n=", num)

    # algorithm 2
    # start from the biggest possible difference between primes
    if 'a2' in algo_to_check:
        p1 = 3
        p2 = num - 3

        (p1, p2, duration, iterations) = search_for_partition (p1, p2, lambda iteration: delta_constant_plus(iteration))
        dt_diff[1] += duration
        dt_iter[1] += iterations

        if p1 + p2 != num:
            print ("Alg #2: violation of sum for p1=", p1, "p2=", p2, "n=", num)

    # algorithm 3
    # start from 1/3 of n
    if 'a3' in algo_to_check:
        p1 = int(num / 3)
        p2 = num - p1
        if p1 % 2 == 0:
            p1 = p1 - 1
            p2 = p2 + 1

        (p1, p2, duration, iterations) = search_for_partition (p1, p2, lambda iteration: delta_constant_minus(iteration))
        dt_diff[2] += duration
        dt_iter[2] += iterations

        if p1 + p2 != num:
            print ("Alg #3: violation of sum for p1=", p1, "p2=", p2, "n=", num)

    # algorithm 4
    # start from the biggest possible difference between primes
    if 'a4' in algo_to_check:
        p1 = 5
        p2 = num - 5

        (p1, p2, duration, iterations) = search_for_partition (p1, p2, lambda iteration: delta_constant_plus(iteration))
        dt_diff[3] += duration
        dt_iter[3] += iterations

        if p1 + p2 != num:
            print ("Alg #4: violation of sum for p1=", p1, "p2=", p2, "n=", num)

    # algorithm 5
    # start from the biggest possible difference between primes
    if 'a5' in algo_to_check:
        p1 = 5
        p2 = num - 5

        (p1, p2, duration, iterations) = search_for_partition (p1, p2, lambda iteration: delta_variable(iteration))
        dt_diff[4] += duration
        dt_iter[4] += iterations

        if p1 + p2 != num:
            print ("Alg #5: violation of sum for p1=", p1, "p2=", p2, "n=", num)

    # algorithm 6
    # A2 but next p1 is always prime
    if 'a6' in algo_to_check:
        p1 = 3
        p2 = num - 3

        (p1, p2, duration, iterations) = search_for_partition (p1, p2, lambda iteration: delta_prime(iteration))
        dt_diff[5] += duration
        dt_iter[5] += iterations

        if p1 + p2 != num:
            print ("Alg #6: violation of sum for p1=", p1, "p2=", p2, "n=", num)
    
    # checkpoint - partial results
    if num % checkpoint_value == 0:
        dt_current = datetime.now()
        dt_diff_current = (dt_current - dt_current_previous).total_seconds()
        list_checkpoints.append(num)
        if 'a1' in algo_to_check:
            list_checkpoints_duration[0].append(dt_diff[0])
            list_checkpoints_iters[0].append(dt_iter[0])
        if 'a2' in algo_to_check:
            list_checkpoints_duration[1].append(dt_diff[1])
            list_checkpoints_iters[1].append(dt_iter[1])
        if 'a3' in algo_to_check:
            list_checkpoints_duration[2].append(dt_diff[2])
            list_checkpoints_iters[2].append(dt_iter[2])
        if 'a4' in algo_to_check:
            list_checkpoints_duration[3].append(dt_diff[3])
            list_checkpoints_iters[3].append(dt_iter[3])
        if 'a5' in algo_to_check:
            list_checkpoints_duration[4].append(dt_diff[4])
            list_checkpoints_iters[4].append(dt_iter[4])
        if 'a6' in algo_to_check:
            list_checkpoints_duration[5].append(dt_diff[5])
            list_checkpoints_iters[5].append(dt_iter[5])

        perc_completed = str(int(k * 100 / max_num))
        print ("Checkpoint", k, "of total", max_num, "took", dt_diff_current, "seconds. (" + perc_completed + "% completed)")
        
        # remember results so far
        write_results_to_figures (directory)

dt_end = datetime.now()

# final results - time of processing
dt_diff = dt_end - dt_start
print ("Total calculations lasted:", dt_diff)

# final results - figures
write_results_to_figures (directory)
