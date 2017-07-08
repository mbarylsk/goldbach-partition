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

import math
import sys
import numpy as np
import os
from datetime import datetime
import time
import goldbach
sys.path.insert(0, '..\\primes\\')
import primes

#############################################################
# Settings - configuration
#############################################################

# Minimal even number checked against Goldbach conjecture
#   o number = min_num * step_factor
min_num = 7
step_factor = 2
# Maximum even number checked against Goldbach conjecture
#   o number = max_num * step_factor
max_num = 20000

# Searching for max reduction
#   o True  - maximum 2^n
#   o False - minimal 2^n
search_for_max_reduction = False

# Caching previous primality results
#   o True  - auxilary sets of primes and composite numbers will grow
#             it will speed up further primality tests but more RAM will
#             be occupied
#   o False - do not cache new primality test results
caching_primality_results = False

# Checkpoint value when partial results are drawn/displayed
checkpoint_value = 2000

# Helper files
#   o file_input_primes - contains prime numbers
#   o file_input_nonprimes - contains composite numbers
file_input_primes = '..\\primes\\t_prime_numbers.txt'
file_input_nonprimes = '..\\primes\\t_nonprime_numbers.txt'

#############################################################
# Settings - output directory and files
#############################################################

directory = "results/" + str(step_factor*max_num)
if not os.path.exists(directory):
    os.makedirs(directory)
file_output_data_goldbach_reduction = directory + "/f_goldbach_reduction.csv"

#############################################################
# Presentation
#############################################################

def file_init_header ():
    header = "n,hypothesis_fulfilled?,p1,p2,n1,n2,q1,q2"
    f = open(file_output_data_goldbach_reduction, "w+")
    f.write (header + "\n")
    f.close ()

def file_write_line (results):
    f = open(file_output_data_goldbach_reduction, "a+")
    f.write (results + "\n")
    f.close ()

#############################################################
# Main
#############################################################

print ("Initialize objects...")
p = primes.Primes(caching_primality_results)
gp = goldbach.GoldbachPartition (p)
print ("DONE")
print ("Loading helper sets...")
p.init_set(file_input_primes, True)
p.init_set(file_input_nonprimes, False)
print ("DONE")
print ("Preparing result file...")
file_init_header ()
print ("DONE")

dt_start = datetime.now()
dt_current_previous = dt_start

# new calculations
for k in range (min_num, max_num):
    num = step_factor*k
    
    factors = gp.find_sum_of_prime_numbers (num)

    found = False
    max_q = 0
    max_q1 = 0
    max_q2 = 0
    max_n1 = 0
    max_n2 = 0
    max_p1 = 0
    max_p2 = 0
    for (p1, p2) in factors:
        (n1, q1, res) = gp.reduce_prime_for_goldbach (p1, search_for_max_reduction)
        (n2, q2, res) = gp.reduce_prime_for_goldbach (p2, search_for_max_reduction)

        if search_for_max_reduction:
            if q1 + q2 > max_q:
                max_q1 = q1
                max_q2 = q2
                max_n1 = n1
                max_n2 = n2
                max_p1 = p1
                max_p2 = p2
        else:
            if (q1 >= 2) and (q2 >= 2):
                found = True
                break

    if found:
        if search_for_max_reduction:
            results = str(num) + ",True," + str(max_p1) + "," + str(max_p2) + "," + str(max_n1) + "," + str(max_n2) + "," + str(max_q1) + "," + str(max_q2)
        else:
            results = str(num) + ",True," + str(p1) + "," + str(p2) + "," + str(n1) + "," + str(n2) + "," + str(q1) + "," + str(q2)

    else:
        Exception ("Could not find reduction")

    file_write_line (results)

    # checkpoint - partial results
    if num % checkpoint_value == 0:
        dt_current = datetime.now()
        dt_diff_current = (dt_current - dt_current_previous).total_seconds()

        perc_completed = str(int(k * 100 / max_num))
        print ("Checkpoint", k, "of total", max_num, "took", dt_diff_current, "seconds. (" + perc_completed + "% completed)")
        
dt_end = datetime.now()

# final results - time of processing
dt_diff = dt_end - dt_start
print ("Total calculations lasted:", dt_diff)
