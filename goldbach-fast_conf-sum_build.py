#
# Copyright (c) 2017 - 2018, Marcin Barylski
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

import sys
import os
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import goldbach
import pickle
sys.path.insert(0, '..\\primes\\')
import primes
import dataprocessing

#############################################################
# Settings - configuration
#############################################################

min_prime_index = 1
max_prime_index = 100000
max_num = max_prime_index - min_prime_index
checkpoint_value = 5000

# Caching previous primality results
#   o True  - auxilary sets of primes and composite numbers will grow
#             it will speed up further primality tests but more RAM will
#             be occupied
#   o False - do not cache new primality test results
caching_primality_results = False

# Detecting terms of http://oeis.org/A301776
#   o True  - detect terms
#   o False - do not detect terms
detect_no_backlog = True

# Methods:
#   o method 1 - check all possible sums of primes, build all possible GPs
#   o method 2 - check all possible sums of primes, optimized, show min number where all GPs verified
method = 2

# Helper files
#   o file_input_primes - contains prime numbers
#   o file_input_nonprimes - contains composite numbers
file_input_primes = '..\\primes\\t_prime_numbers.txt'
file_input_nonprimes = '..\\primes\\t_nonprime_numbers.txt'

#############################################################
# Settings - output directory and files
#############################################################

directory = "results/" + str(max_num)
if not os.path.exists(directory):
    os.makedirs(directory)
file_output_fig1 = directory + "/f_sumbuild_effectiveness.png"
file_output_fig2 = directory + "/f_sumbuild_spare_and_to_be_verified.png"
file_output_fig3 = directory + "/f_sumbuild_diff_and_delta.png"
file_output_fig4 = directory + "/f_sumbuild_diff.png"
file_output_fig5 = directory + "/f_sumbuild_ratio_to_be_ver.png"
file_output_fig6 = directory + "/f_sumbuild_ratio_alr_ver.png"
file_output_fig7 = directory + "/f_sumbuild_ratio_to_be_ver_diff.png"
file_output_fig8 = directory + "/f_sumbuild_ratio_alr_ver_diff.png"
file_output_pickle = directory + "/objs_sym_primes.pickle"

#############################################################
# Results of calculations
#############################################################

local_primes = []
lp_sum = 0

verified_intervals = set()

to_be_verified = set()
already_verified = set()
num_where_all_verified = 4

list_num_where_all_verified =   []
list_to_be_verified =           []
list_already_verified =         []
list_to_be_verified_avg =       []
list_already_verified_avg =     []
list_num_max =                  []
# [0] - actual diff, [1] - delta between two consecutive diffs
list_diff =                     [[],[]]
list_checkpoints =              []
# [0] - diff/# of to be verified, [1] - diff/# of already verified
list_ratio =                    [[],[]]
# [0] - average of list_ratio[0], [1] - average of list_ratio[1]
list_ratio_avg =                [[],[]]
list_ratio_avg_diff =           [[],[]]
# http://oeis.org/A301776
list_A301776 =                  [2]

k_current = 0
k = 0
diff_previous = 0

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

def save_current_results (file_output_pickle):
    global k_current, to_be_verified, already_verified, num_where_all_verified, list_num_where_all_verified, list_to_be_verified, list_already_verified, list_to_be_verified_avg, list_already_verified_avg, list_num_max, list_diff, list_ratio, list_ratio_avg, list_checkpoints, diff_previous
    global local_primes, lp_sum, verified_intervals, list_A301776
    with open(file_output_pickle, 'wb') as f:
        pickle.dump([k_current, list_A301776, to_be_verified, already_verified, num_where_all_verified, list_num_where_all_verified, list_to_be_verified, list_already_verified, list_to_be_verified_avg, list_already_verified_avg, list_num_max, list_diff, list_ratio, list_ratio_avg, list_checkpoints, diff_previous, local_primes, lp_sum, verified_intervals], f)

def restore_previous_results (file_output_pickle):
    global k_current, to_be_verified, already_verified, num_where_all_verified, list_num_where_all_verified, list_to_be_verified, list_already_verified, list_to_be_verified_avg, list_already_verified_avg, list_num_max, list_diff, list_ratio, list_ratio_avg, list_checkpoints, diff_previous
    global local_primes, lp_sum, verified_intervals, list_A301776
    if os.path.exists(file_output_pickle):
        with open(file_output_pickle, 'rb') as f:
            k_current, list_A301776, to_be_verified, already_verified, num_where_all_verified, list_num_where_all_verified, list_to_be_verified, list_already_verified, list_to_be_verified_avg, list_already_verified_avg, list_num_max, list_diff, list_ratio, list_ratio_avg, list_checkpoints, diff_previous, local_primes, lp_sum, verified_intervals = pickle.load(f)

#############################################################
# Presentation
#############################################################

def print_results_1 (i):
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

def print_results_2 (i):
    global to_be_verified, num_where_all_verified, already_verified, list_diff
    print ("===============================================================================")
    print ("Iteration:", i)
    print (" Even numbers to be verified:", to_be_verified)
    print (" # of numbers to be verified:", len(to_be_verified))
    print (" All even numbers verified up to:", num_where_all_verified)
    print (" Spare verified numbers:", already_verified)
    print (" # of spare verified numbers:", len(already_verified))
    print (" Current vs. theoretical max (diff):", list_diff[0][i])
    print (" Delta in diff between two last iterations:", list_diff[1][i])
    print (" A301776:", list_A301776)
    print (" diff / # of numbers to be verified:", list_ratio_avg[0][i])
    print (" diff / # of numbers to already verified:", list_ratio_avg[1][i])


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
    plt.savefig(file_output_fig1)
    plt.close(fig1)

    fig2 = plt.figure(1)
    r_patch = mpatches.Patch(color='red', label='# to be verified')
    g_patch = mpatches.Patch(color='green', label='# already verified')
    k_patch = mpatches.Patch(color='black', label='avg # to be verified')
    b_patch = mpatches.Patch(color='blue', label='avg # already verified')
    list_of_handles = []
    list_of_handles.append(r_patch)
    list_of_handles.append(g_patch)
    list_of_handles.append(k_patch)
    list_of_handles.append(b_patch)
    plt.legend(handles=list_of_handles, loc='upper right', bbox_to_anchor=(0.4, 0.8))
    plt.plot(list_checkpoints, list_to_be_verified, 'r-', ms=1)
    plt.plot(list_checkpoints, list_already_verified, 'g-', ms=1)
    plt.plot(list_checkpoints, list_to_be_verified_avg, 'k-', ms=1)
    plt.plot(list_checkpoints, list_already_verified_avg, 'b-', ms=1)
    plt.xlabel('iteration')
    plt.ylabel('number of elements')
    plt.title('Spare verified numbers and numbers to be verified')
    plt.grid(True)
    plt.savefig(file_output_fig2)
    plt.close(fig2)

    fig3 = plt.figure(1)
    r_patch = mpatches.Patch(color='red', label='diff')
    g_patch = mpatches.Patch(color='green', label='delta in iters')
    list_of_handles = []
    list_of_handles.append(r_patch)
    list_of_handles.append(g_patch)
    plt.legend(handles=list_of_handles, loc='upper right', bbox_to_anchor=(0.4, 0.8))
    plt.plot(list_checkpoints, list_diff[0], 'r-', ms=1)
    plt.plot(list_checkpoints, list_diff[1], 'g-', ms=1)
    plt.xlabel('iteration')
    plt.ylabel('difference')
    plt.title('Difference between actual and theoretical max (and delta)')
    plt.grid(True)
    plt.savefig(file_output_fig3)
    plt.close(fig3)

    fig4 = plt.figure(1)
    b_patch = mpatches.Patch(color='blue', label='diff')
    list_of_handles = []
    list_of_handles.append(b_patch)
    plt.legend(handles=list_of_handles, loc='upper right', bbox_to_anchor=(0.4, 0.8))
    plt.plot(list_checkpoints, list_diff[0], 'b.', ms=1)
    plt.xlabel('iteration')
    plt.ylabel('difference')
    plt.title('Difference between actual and theoretical max')
    plt.grid(True)
    plt.savefig(file_output_fig4)
    plt.close(fig4)

    fig5 = plt.figure(1)
    b_patch = mpatches.Patch(color='blue', label='diff/# to be verified')
    r_patch = mpatches.Patch(color='red', label='avg')
    list_of_handles = []
    list_of_handles.append(b_patch)
    list_of_handles.append(r_patch)
    plt.legend(handles=list_of_handles, loc='upper right', bbox_to_anchor=(0.4, 0.8))
    plt.plot(list_checkpoints, list_ratio[0], 'b-', ms=1)
    plt.plot(list_checkpoints, list_ratio_avg[0], 'r-', ms=1)
    plt.xlabel('iteration')
    plt.ylabel('ratio')
    plt.title('Ratio')
    plt.grid(True)
    plt.savefig(file_output_fig5)
    plt.close(fig5)

    fig6 = plt.figure(1)
    g_patch = mpatches.Patch(color='green', label='diff/# already verified')
    r_patch = mpatches.Patch(color='red', label='avg')
    list_of_handles = []
    list_of_handles.append(g_patch)
    list_of_handles.append(r_patch)
    plt.legend(handles=list_of_handles, loc='upper right', bbox_to_anchor=(0.4, 0.8))
    plt.plot(list_checkpoints, list_ratio[1], 'g-', ms=1)
    plt.plot(list_checkpoints, list_ratio_avg[1], 'r-', ms=1)
    plt.xlabel('iteration')
    plt.ylabel('ratio')
    plt.title('Ratio')
    plt.grid(True)
    plt.savefig(file_output_fig6)
    plt.close(fig6)

#############################################################
# Main
#############################################################

print ("Initialize objects...")
p = primes.Primes(caching_primality_results)
gp = goldbach.GoldbachPartition (p)
dp = dataprocessing.DataProcessing()
print ("DONE")
print ("Loading helper sets...")
p.init_set(file_input_primes, 1)
p.init_set(file_input_nonprimes, 3)
print ("DONE")
print ("Sorting primes...")
p.sort_primes_set()
print ("DONE")
print ("Restoring previous results...")
restore_previous_results (file_output_pickle)
if k_current > 0:
    min_prime_index = k_current
    k = k_current
    print ("Resuming calculations at", min_prime_index)
print ("DONE")

dt_start = datetime.now()
dt_current_previous = dt_start

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
        
            print_results_1 (k)
            k_current = k
            save_current_results(file_output_pickle)

        k += 1
    
elif method == 2:
    
    for ip1 in range (min_prime_index, max_prime_index):

        # check all possible sums
        p1 = p.get_ith_prime(ip1)
        add_nums_to_be_verified (2*p1)
        for ip2 in range (1, ip1+1):
            p2 = p.get_ith_prime(ip2)
            num = p1 + p2
            remove_nums_to_be_verified (num, 2*p1)

        # special case - max == current (no backlog)
        if detect_no_backlog:
            if len(to_be_verified) == 0:
                print ("+-----------------------------------------------------------------------------------------------------")
                print ("| Interesting - For iteration", k, "(this is prime", p1, ") I have verified all possible numbers till", num)
                print ("+-----------------------------------------------------------------------------------------------------")
                list_A301776.append (p1)

        # checkpoint - partial results
        if k % checkpoint_value == 0:
            dt_current = datetime.now()
            dt_diff_current = (dt_current - dt_current_previous).total_seconds()

            perc_completed = str(int(k * 100 / max_num))
            print ("Checkpoint", k, "of total", max_num, "took", dt_diff_current, "seconds. (" + perc_completed + "% completed)")
        
            print_results_2 (k)
            k_current = k
            save_current_results(file_output_pickle)

        list_checkpoints.append (k)
        list_num_where_all_verified.append (num_where_all_verified)
        list_num_max.append (num)
        diff_now = num - num_where_all_verified
        list_diff[0].append (diff_now)
        list_diff[1].append (abs(diff_now - diff_previous))
        if num == num_where_all_verified:
            print ("Verified all till max for iteration", k)
        list_to_be_verified.append(len(to_be_verified))
        list_already_verified.append(len(already_verified))
        if len(to_be_verified) > 0:
            list_ratio[0].append (diff_now/len(to_be_verified))
        else:
            list_ratio[0].append (0)
        if len(already_verified) > 0:
            list_ratio[1].append (diff_now/len(already_verified))
        else:
            list_ratio[1].append (0)
        list_ratio_avg[0].append(dp.get_avg_value_from_list (list_ratio[0]))
        list_ratio_avg[1].append(dp.get_avg_value_from_list (list_ratio[1]))
        list_to_be_verified_avg.append(dp.get_avg_value_from_list(list_to_be_verified))
        list_already_verified_avg.append(dp.get_avg_value_from_list(list_already_verified))
        
        k += 1
        diff_previous = diff_now

dt_end = datetime.now()

# final results - time of processing
dt_diff = dt_end - dt_start
if method == 1:
    print_results_1 (k)
elif method == 2:
    print_results_2 (k)
k_current = max_num
save_current_results(file_output_pickle)

print ("Total calculations lasted:", dt_diff)

threshold_max = 700
threshold_min = 50
neighbour = 20
for x in range (1, k):
    if (list_diff[0][x] > list_diff[0][x-neighbour] + threshold_max) and (list_diff[0][x] > list_diff[0][x+neighbour] + threshold_max) and list_diff[0][x] > 2000:
        print ("MAX: k=", x, "previous=", list_diff[0][x-neighbour], "current=", list_diff[0][x], "next=", list_diff[0][x+neighbour])

for x in range (0, k):
    if list_diff[0][x] < threshold_min:
        print ("MIN: k=", x, "current=", list_diff[0][x])

list_ratio_avg_diff[0].append(0)
list_ratio_avg_diff[1].append(0)
for x in range (1, len(list_ratio_avg[0])):
    list_ratio_avg_diff[0].append(list_ratio_avg[0][x-1] - list_ratio_avg[0][x])
    list_ratio_avg_diff[1].append(list_ratio_avg[1][x-1] - list_ratio_avg[1][x])

fig7 = plt.figure(1)
b_patch = mpatches.Patch(color='blue', label='diff')
list_of_handles = []
list_of_handles.append(b_patch)
plt.legend(handles=list_of_handles, loc='upper right', bbox_to_anchor=(0.4, 0.8))
plt.plot(list_checkpoints, list_ratio_avg_diff[0], 'b-', ms=1)
plt.xlabel('iteration')
plt.ylabel('diff')
plt.ylim(ymax = 0.02, ymin = -0.02)
plt.title('Average diff')
plt.grid(True)
plt.savefig(file_output_fig7)
plt.close(fig7)

fig8 = plt.figure(1)
g_patch = mpatches.Patch(color='green', label='diff')
list_of_handles = []
list_of_handles.append(g_patch)
plt.legend(handles=list_of_handles, loc='upper right', bbox_to_anchor=(0.4, 0.8))
plt.plot(list_checkpoints, list_ratio_avg_diff[1], 'g-', ms=1)
plt.ylim(ymax = 0.001, ymin = -0.001)
plt.xlabel('iteration')
plt.ylabel('diff')
plt.title('Average diff')
plt.grid(True)
plt.savefig(file_output_fig8)
plt.close(fig8)
