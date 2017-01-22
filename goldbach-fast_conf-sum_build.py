#
# Copyright (c) 2017, Marcin Barylski
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

from datetime import datetime
import goldbach
import primes

#############################################################
# Settings - configuration
#############################################################

min_prime_index = 1
max_prime_index = 1000

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
# Results of calculations
#############################################################

# skipping 2 - because it is a part of 1 GP only: 4 = 2 + 2
# starting from 3
local_primes = []
lp_sum = 0

verified_intervals = set()

#############################################################
# Presentation
#############################################################

def add_to_verified_intervals (num):
    global verified_intervals
    verified_intervals.add(num)

def get_data_from_verified_intervals ():
    global verified_intervals
    missing = set()
    item_previous = 4
    list_verified_intervals = sorted (verified_intervals)
    for item in list_verified_intervals:
        finished = False
        k = 2
        while not finished:
            if item_previous + k < item:
                missing.add(item_previous + k)
                k += 2
            else:
                finished = True
        item_previous = item
    return (missing, min(list_verified_intervals), max(list_verified_intervals))

def print_stats ():
    (miss_nums, num_min, num_max) = get_data_from_verified_intervals()
    print ("Primes used to generate GP:", local_primes)
    print (" Verified min :", num_min)
    print (" Verified max :", num_max)
    if len(miss_nums) == 0:
        print (" Missing nums : -")
    else:
        print (" Missing nums :", miss_nums)

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
print ("Sorting primes...")
p.sort_prime_set()
print ("DONE")

dt_start = datetime.now()
dt_current_previous = dt_start

# new calculations
for k in range (min_prime_index, max_prime_index):

    # add new prime to GP build base
    local_primes.append(p.get_ith_prime(k))

    # build all possible GPs
    for lp in local_primes:
        lp_sum = lp + max(local_primes)
        add_to_verified_intervals (lp_sum)

dt_end = datetime.now()

# final results - time of processing
dt_diff = dt_end - dt_start
print_stats ()
print ("Total calculations lasted:", dt_diff)
