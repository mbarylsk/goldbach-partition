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
#   o True - run unit tests
#   o False - run main program
run_unit_tests = False

# Minimal even number checked against Goldbach conjecture
#   o number = min_num * step_factor
min_num = 7
step_factor = 2
# Maximum even number checked against Goldbach conjecture
#   o number = max_num * step_factor
max_num = 10000000

# Checkpoint value when partial results are drawn/displayed
checkpoint_value = 20000
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

#############################################################
# Business logic
#############################################################

set_prime = set ()
set_nonprime = set ()

list_nums = []
list_sorted_prime = []

list_checkpoints_duration_a1 = []
list_checkpoints_duration_a2 = []
list_checkpoints_duration_a3 = []
list_checkpoints_duration_a4 = []
list_checkpoints_duration_a5 = []
list_checkpoints_duration_a6 = []
list_checkpoints_iters_a1 = []
list_checkpoints_iters_a2 = []
list_checkpoints_iters_a3 = []
list_checkpoints_iters_a4 = []
list_checkpoints_iters_a5 = []
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

def search_for_partition2 (p1, p2):
    found = False
    iteration = 0
    n = p1 + p2
    startTime = time.time()
    while (not found):
        iteration += 1
        if (is_prime (p1) and is_prime(p2)):
            found = True
        if (not found):
            p1 = get_ith_prime (iteration)
            p2 = n - p1
        if p1 < 2 or p2 < 2:
            raise ("Could not find GP")
    duration = time.time() - startTime
    return p1, p2, duration, iteration

#############################################################
# Presentation
#############################################################

def write_results_to_figures (directory, last_loop):
    # results - figures    
    plt.figure(1)
    plt.plot(list_checkpoints, list_checkpoints_duration_a1, 'g.', ms=2)
    plt.plot(list_checkpoints, list_checkpoints_duration_a2, 'r.', ms=2)
    plt.plot(list_checkpoints, list_checkpoints_duration_a2, 'b.', ms=2)
    plt.plot(list_checkpoints, list_checkpoints_duration_a4, 'y.', ms=2)
    plt.plot(list_checkpoints, list_checkpoints_duration_a5, 'c.', ms=2)
    plt.plot(list_checkpoints, list_checkpoints_duration_a6, 'm.', ms=2)
    g_patch = mpatches.Patch(color='green', label='A1')
    r_patch = mpatches.Patch(color='red', label='A2')
    b_patch = mpatches.Patch(color='blue', label='A3')
    y_patch = mpatches.Patch(color='yellow', label='A4')
    c_patch = mpatches.Patch(color='cyan', label='A5')
    m_patch = mpatches.Patch(color='magenta', label='A6')
    plt.legend(handles=[g_patch, r_patch, b_patch, y_patch, c_patch, m_patch], loc='upper right', bbox_to_anchor=(0.4, 0.8))
    plt.xlabel('Number')
    plt.ylabel('Time [s]')
    plt.title('Duration of total calculations')
    plt.grid(True)
    plt.savefig(directory + "/f_checkpoint_duration_alg.png")

    plt.figure(2)
    plt.plot(list_checkpoints, list_checkpoints_iters_a1, 'g.', ms=2)
    plt.plot(list_checkpoints, list_checkpoints_iters_a2, 'r.', ms=2)
    plt.plot(list_checkpoints, list_checkpoints_iters_a3, 'b.', ms=2)
    plt.plot(list_checkpoints, list_checkpoints_iters_a4, 'y.', ms=2)
    plt.plot(list_checkpoints, list_checkpoints_iters_a5, 'c.', ms=2)
    plt.plot(list_checkpoints, list_checkpoints_iters_a6, 'm.', ms=2)
    g_patch = mpatches.Patch(color='green', label='A1')
    r_patch = mpatches.Patch(color='red', label='A2')
    b_patch = mpatches.Patch(color='blue', label='A3')
    y_patch = mpatches.Patch(color='yellow', label='A4')
    c_patch = mpatches.Patch(color='cyan', label='A5')
    m_patch = mpatches.Patch(color='magenta', label='A6')
    plt.legend(handles=[g_patch, r_patch, b_patch, y_patch, c_patch, m_patch], loc='upper right', bbox_to_anchor=(0.4, 0.8))
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

dt_diff1 = 0
dt_diff2 = 0
dt_diff3 = 0
dt_diff4 = 0
dt_diff5 = 0
dt_diff6 = 0

dt_iter1 = 0
dt_iter2 = 0
dt_iter3 = 0
dt_iter4 = 0
dt_iter5 = 0
dt_iter6 = 0

# new calculations
for k in range (min_num, max_num):
    num = step_factor*k

    list_nums.append (num)

    #print ("A1")
    # algorithm 1
    # start from half of n - the smallest possible difference between primes
    p1 = num / 2
    if p1 % 2 == 0:
        p1 = p1 - 1
        p2 = p1 + 2
    else:
        p2 = p1

    (p1, p2, duration, iterations) = search_for_partition (p1, p2, lambda iteration: delta_constant_plus(iteration))
    dt_diff1 += duration
    dt_iter1 += iterations

    if p1 + p2 != num:
        print ("Alg #1: violation of sum for p1=", p1, "p2=", p2, "n=", num)

    #print ("A2")
    # algorithm 2
    # start from the biggest possible difference between primes
    p1 = 3
    p2 = num - 3

    (p1, p2, duration, iterations) = search_for_partition (p1, p2, lambda iteration: delta_constant_plus(iteration))
    dt_diff2 += duration
    dt_iter2 += iterations

    if p1 + p2 != num:
        print ("Alg #2: violation of sum for p1=", p1, "p2=", p2, "n=", num)

    #print ("A3")
    # algorithm 3
    # start from 1/3 of n
    p1 = int(num / 3)
    p2 = num - p1
    if p1 % 2 == 0:
        p1 = p1 - 1
        p2 = p2 + 1

    (p1, p2, duration, iterations) = search_for_partition (p1, p2, lambda iteration: delta_constant_plus(iteration))
    dt_diff3 += duration
    dt_iter3 += iterations

    if p1 + p2 != num:
        print ("Alg #3: violation of sum for p1=", p1, "p2=", p2, "n=", num)

    #print ("A4")
    # algorithm 4
    # start from the biggest possible difference between primes
    p1 = 5
    p2 = num - 5

    (p1, p2, duration, iterations) = search_for_partition (p1, p2, lambda iteration: delta_constant_plus(iteration))
    dt_diff4 += duration
    dt_iter4 += iterations

    if p1 + p2 != num:
        print ("Alg #4: violation of sum for p1=", p1, "p2=", p2, "n=", num)

    #print ("A5")
    # algorithm 5
    # start from the biggest possible difference between primes
    p1 = 5
    p2 = num - 5

    (p1, p2, duration, iterations) = search_for_partition (p1, p2, lambda iteration: delta_variable(iteration))
    dt_diff5 += duration
    dt_iter5 += iterations

    if p1 + p2 != num:
        print ("Alg #5: violation of sum for p1=", p1, "p2=", p2, "n=", num)

    #print ("A6")
    # algorithm 6
    # A2 but next p1 is always prime
    p1 = 3
    p2 = num - 3

    (p1, p2, duration, iterations) = search_for_partition2 (p1, p2)
    dt_diff6 += duration
    dt_iter6 += iterations

    if p1 + p2 != num:
        print ("Alg #6: violation of sum for p1=", p1, "p2=", p2, "n=", num)
    
    # checkpoint - partial results
    if num % checkpoint_value == 0:
        dt_current = datetime.now()
        dt_diff_current = (dt_current - dt_current_previous).total_seconds()
        list_checkpoints.append(num)
        list_checkpoints_duration_a1.append(dt_diff1)
        list_checkpoints_duration_a2.append(dt_diff2)
        list_checkpoints_duration_a3.append(dt_diff3)
        list_checkpoints_duration_a4.append(dt_diff4)
        list_checkpoints_duration_a5.append(dt_diff5)
        list_checkpoints_duration_a6.append(dt_diff6)
        list_checkpoints_iters_a1.append(dt_iter1)
        list_checkpoints_iters_a2.append(dt_iter2)
        list_checkpoints_iters_a3.append(dt_iter3)
        list_checkpoints_iters_a4.append(dt_iter4)
        list_checkpoints_iters_a5.append(dt_iter5)
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
