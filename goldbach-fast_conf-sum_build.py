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
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
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

# Methods:
#   o method 1 - check all possible sums of primes, build all possible GPs
#   o method 2 - check all possible sums of primes, optimized, show min number where all GPs verified
method = 2

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

list_num_where_all_verified = []
list_to_be_verified = []
list_already_verified = []
list_num_max = []
list_max_actual_diff = []
list_max_actual_diff_perc = []
list_iters = []
list_checkpoints = []

#############################################################
# Business logic
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

#############################################################
# Presentation
#############################################################

def print_stats_1 (i):
    (miss_nums, num_min, num_max) = get_data_from_verified_intervals()
    print ("================================================================================")
    print ("Iteration:", i)
    print (" Primes used to generate GP:", local_primes)
    print ("  Verified min :", num_min)
    print ("  Verified max :", num_max)
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

    fig1 = plt.figure(1)
    r_patch = mpatches.Patch(color='red', label='all verified')
    g_patch = mpatches.Patch(color='green', label='theoretical max')
    list_of_handles = []
    list_of_handles.append(r_patch)
    list_of_handles.append(g_patch)
    plt.legend(handles=list_of_handles, loc='upper right', bbox_to_anchor=(0.4, 0.8))
    plt.plot(list_checkpoints, list_num_where_all_verified, 'r-', ms=2)
    plt.plot(list_checkpoints, list_num_max, 'g-', ms=2)
    plt.xlabel('round of algorithm')
    plt.ylabel('number')
    plt.title('Min number where all numbers below are verified from GP standpoint')
    plt.grid(True)
    file_output_fig1 = "results/fig1_effectiveness.png"
    plt.savefig(file_output_fig1)
    plt.close(fig1)

    fig2 = plt.figure(1)
    r_patch = mpatches.Patch(color='red', label='to be verified')
    g_patch = mpatches.Patch(color='green', label='already verified')
    list_of_handles = []
    list_of_handles.append(r_patch)
    list_of_handles.append(g_patch)
    plt.legend(handles=list_of_handles, loc='upper right', bbox_to_anchor=(0.4, 0.8))
    plt.plot(list_checkpoints, list_to_be_verified, 'r-', ms=1)
    plt.plot(list_checkpoints, list_already_verified, 'g-', ms=1)
    plt.xlabel('iteration')
    plt.ylabel('number of elements')
    plt.title('Spare verified numbers and numbers to be verified')
    plt.grid(True)
    file_output_fig2 = "results/fig2_spare_and_to_be_verified.png"
    plt.savefig(file_output_fig2)
    plt.close(fig2)

    fig3 = plt.figure(1)
    r_patch = mpatches.Patch(color='red', label='diff')
    g_patch = mpatches.Patch(color='green', label='delta in iters')
    list_of_handles = []
    list_of_handles.append(r_patch)
    list_of_handles.append(g_patch)
    plt.legend(handles=list_of_handles, loc='upper right', bbox_to_anchor=(0.4, 0.8))
    plt.plot(list_checkpoints, list_max_actual_diff, 'r-', ms=1)
    plt.plot(list_checkpoints, list_max_actual_diff_perc, 'g-', ms=1)
    plt.xlabel('iteration')
    plt.ylabel('difference')
    plt.title('Difference between actual and theoretical minimum')
    plt.grid(True)
    file_output_fig3 = "results/fig3_diff.png"
    plt.savefig(file_output_fig3)
    plt.close(fig3)

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

k = 0
max_num = max_prime_index - min_prime_index
diff_previous = 0
iterations = 0

# new calculations
if method == 1: 
    for ip in range (min_prime_index, max_prime_index):

        # add new prime to GP build base
        local_primes.append(p.get_ith_prime(ip))
        
        # build all possible GPs
        for lp in local_primes:
            lp_sum = lp + max(local_primes)
            add_to_verified_intervals (lp_sum)

        # checkpoint - partial results
        if k % checkpoint_value == 0:
            dt_current = datetime.now()
            dt_diff_current = (dt_current - dt_current_previous).total_seconds()

            perc_completed = str(int(k * 100 / max_num))
            print ("Checkpoint", k, "of total", max_num, "took", dt_diff_current, "seconds. (" + perc_completed + "% completed)")
        
            print_stats_1 (k)


        k += 1
    
elif method == 2:
    
    for ip1 in range (min_prime_index, max_prime_index):

        p1 = p.get_ith_prime(ip1)
        add_nums_to_be_verified (2*p1)
    
        for ip2 in range (1, ip1+1):
            p2 = p.get_ith_prime(ip2)
            num = p1 + p2
            remove_nums_to_be_verified (num, 2*p1)

        # checkpoint - partial results
        if k % checkpoint_value == 0:
            dt_current = datetime.now()
            dt_diff_current = (dt_current - dt_current_previous).total_seconds()

            perc_completed = str(int(k * 100 / max_num))
            print ("Checkpoint", k, "of total", max_num, "took", dt_diff_current, "seconds. (" + perc_completed + "% completed)")
        
            print_stats_2 (k)

        list_checkpoints.append (k)
        list_num_where_all_verified.append (num_where_all_verified)
        list_num_max.append (num)
        diff_now = num - num_where_all_verified
        list_max_actual_diff.append (diff_now)
        list_max_actual_diff_perc.append (abs(diff_now - diff_previous))
        if num == num_where_all_verified:
            print ("Verified all till max for iteration", k)
        list_to_be_verified.append(len(to_be_verified))
        list_already_verified.append(len(already_verified))
        list_iters.append(iterations)
        
        k += 1
        diff_previous = diff_now

dt_end = datetime.now()

# final results - time of processing
dt_diff = dt_end - dt_start
if method == 1:
    print_stats_1 (k)
elif method == 2:
    print_stats_2 (k)

print ("Total calculations lasted:", dt_diff)
