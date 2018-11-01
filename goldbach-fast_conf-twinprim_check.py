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
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.mlab as mlab
import os
from datetime import datetime
import time
import numpy
import pickle
import goldbach
sys.path.insert(0, '..\\primes\\')
import primes

#############################################################
# Settings - configuration
#############################################################

# Minimal even number checked against Goldbach conjecture
#   o number = min_num * step_factor + step_sum
min_num = 1
step_factor = 6
step_sum = -2
# Maximum even number checked against Goldbach conjecture
#   o number = max_num * step_factor + step_sum
max_num = 4000000

# Checkpoint value when partial results are drawn/displayed
checkpoint_value = 10000

# Caching previous primality results
#   o True  - auxilary sets of primes and composite numbers will grow
#             it will speed up further primality tests but more RAM will
#             be occupied
#   o False - do not cache new primality test results
caching_primality_results = False

# Helper files
#   o file_input_primes - contains prime numbers
#   o file_input_nonprimes - contains composite numbers
file_input_primes = '..\\primes\\t_prime_numbers.txt'
file_input_twinprimes = '..\\primes\\t_twinprime_numbers.txt'
file_input_nonprimes = '..\\primes\\t_nonprime_numbers.txt'

#############################################################
# Settings - output directory and files
#############################################################

directory = "results/" + str(step_factor*max_num)
if not os.path.exists(directory):
    os.makedirs(directory)
file_output_pickle = directory + "/objs_lessertwinprimes_in_6km2.pickle"

#############################################################
# Results of calculations
#############################################################

list_A321221 = []
k_current = 0

#############################################################
# Presentation
#############################################################

def save_current_results (file_output_pickle):
    global k_current, list_A321221
    with open(file_output_pickle, 'wb') as f:
        pickle.dump([k_current, list_A321221], f)

def restore_previous_results (file_output_pickle):
    global k_current, list_A321221
    if os.path.exists(file_output_pickle):
        with open(file_output_pickle, 'rb') as f:
            k_current, list_A321221 = pickle.load(f)

#############################################################
# Main
#############################################################

print ("Initialize objects...")
p = primes.Primes(caching_primality_results)
gp = goldbach.GoldbachPartition (p)
print ("DONE")
print ("Loading helper sets...")
p.init_set(file_input_primes, 1)
p.init_set(file_input_twinprimes, 2)
p.init_set(file_input_nonprimes, 3)
print ("DONE")
print ("Sorting primes...")
p.sort_primes_set()
p.sort_twinprimes_set()
print ("DONE")
print ("Restoring previous results...")
restore_previous_results (file_output_pickle)
if k_current > 0:
    min_num = k_current
    print ("Resuming calculations at", min_num)
    print (" A321221 = ", list_A321221)
print ("DONE")

dt_start = datetime.now()
dt_current_previous = dt_start

# new calculations
for k in range (min_num, max_num):
    num = k*step_factor + step_sum
    
    factors = gp.find_sum_of_prime_numbers (num)

    found= False
    for (p1, p2) in factors:
        if p.is_lesser_twin_prime(p1) and p.is_lesser_twin_prime(p2):
            found = True
            break
    if not found:
        print ("NOT FOUND for:", num)
        list_A321221.append(num)

    # checkpoint - partial results
    if num % checkpoint_value == 0:
        dt_current = datetime.now()
        dt_diff_current = (dt_current - dt_current_previous).total_seconds()

        perc_completed = str(int(k * 100 / max_num))
        print ("Checkpoint", k, "of total", max_num, "took", dt_diff_current, "seconds. (" + perc_completed + "% completed)")

        # remember results so far
        k_current = k
        save_current_results(file_output_pickle)


dt_end = datetime.now()

# final results - figures
k_current = max_num
save_current_results(file_output_pickle)
# final results - time of processing
dt_diff = dt_end - dt_start

print ("A321221 = ", list_A321221)
print ("Total calculations lasted:", dt_diff)

