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
max_prime_index = 100000
checkpoint_value = 500

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

local_primes = []
lp_sum = 0

verified_intervals = set()

to_be_verified = set()
already_verified = set()
num_where_all_verified = 4

#############################################################
# Presentation
#############################################################

def add_to_verified_intervals (num):
    global verified_intervals
    verified_intervals.add(num)

def add_nums_to_be_verified (num):
    global to_be_verified, num_where_all_verified
    if num > num_where_all_verified:
        for k in range (num_where_all_verified + 2, num, 2):
            if k not in already_verified:
                to_be_verified.add(k)

def remove_nums_to_be_verified (num, max_num):
    global to_be_verified, num_where_all_verified, already_verified

    already_verified.add(num)
    
    if num in to_be_verified:
        to_be_verified.remove(num)
    if (len(to_be_verified) > 0):
        min_num = min(to_be_verified)
        if min_num > num_where_all_verified + 2:
            num_where_all_verified = min_num
    else:
        num_where_all_verified = max_num

    copy_already_verified = already_verified.copy()
    for n in copy_already_verified:
        if n < num_where_all_verified:
            already_verified.remove(n)
        
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
    print ("================================================================================")
    print ("Primes used to generate GP:", local_primes)
    print (" Verified min :", num_min)
    print (" Verified max :", num_max)
    if len(miss_nums) == 0:
        print (" Missing nums : -")
    else:
        print (" Missing nums :", miss_nums)

def print_stats_2 (i):
    global to_be_verified, num_where_all_verified
    print ("===============================================================================")
    print ("Iteration:", i)
    print (" Even numbers to be verified:", to_be_verified)
    print (" # of numbers to be verified:", len(to_be_verified))
    print (" All even numbers verified up to:", num_where_all_verified)
    print (" Spare verified numbers:", already_verified)
    print (" # of spare verified numbers:", len(already_verified))

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

### new calculations
##for k in range (min_prime_index, max_prime_index):
##
##    # add new prime to GP build base
##    local_primes.append(p.get_ith_prime(k))
##
##    # build all possible GPs
##    for lp in local_primes:
##        lp_sum = lp + max(local_primes)
##        add_to_verified_intervals (lp_sum)

k = 0
max_num = max_prime_index - min_prime_index
for ip1 in range (min_prime_index, max_prime_index):

    p1 = p.get_ith_prime(ip1)
    add_nums_to_be_verified (2*p1)
    #print ("DEBUG1:",ip1,p1)
    
    for ip2 in range (1, ip1+1):
        p2 = p.get_ith_prime(ip2)
        #print ("DEBUG2:",ip2,p2)
        num = p1 + p2
        #print ("DEBUG3:", p1, p2, num)
        remove_nums_to_be_verified (num, 2*p1)

    # checkpoint - partial results
    if ip1 % checkpoint_value == 0:
        dt_current = datetime.now()
        dt_diff_current = (dt_current - dt_current_previous).total_seconds()

        perc_completed = str(int(k * 100 / max_num))
        print ("Checkpoint", k, "of total", max_num, "took", dt_diff_current, "seconds. (" + perc_completed + "% completed)")
        
        print_stats_2 (ip1)
        
    k += 1

dt_end = datetime.now()

# final results - time of processing
dt_diff = dt_end - dt_start
print_stats_2 (ip1)
#print_stats ()
print ("Total calculations lasted:", dt_diff)
