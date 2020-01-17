#
# Copyright (c) 2019 - 2020, Marcin Barylski
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
import numpy as np
import os
from datetime import datetime
import goldbach
sys.path.insert(0, '..\\primes\\')
import primes
import dataprocessing

#############################################################
# Settings - configuration
#############################################################

# set to True if you want to see more output during calculations
be_verbose = False

# Caching previous primality results
#   o True  - auxilary sets of primes and composite numbers will grow
#             it will speed up further primality tests but more RAM will
#             be occupied
#   o False - do not cache new p    rimality test results
caching_primality_results = False

min_num = 2
max_num = 12500
step_factor = 2
checkpoint_value = 500
file_input_primes = '..\\primes\\t_prime_numbers.txt'
file_input_nonprimes = '..\\primes\\t_nonprime_numbers.txt'

#############################################################
# Settings - output directory and files
#############################################################

directory = "results/" + str(step_factor*max_num)
if not os.path.exists(directory):
    os.makedirs(directory)

file_output_list_of_eliminated_primes = directory + "/t_list_of_eliminated_primes_" + str(step_factor*max_num) + ".txt"
file_output_list_of_primes_in_partitions = directory + "/t_list_of_primes_in_partitions_" + str(step_factor*max_num) + ".txt"
file_output_list_of_sets_required_factors = directory + "/t_list_of_sets_required_factors_" + str(step_factor*max_num) + ".txt"

#############################################################
# Business logic
#############################################################

list_of_sets_required_factors = [{2}]
temp_list_of_sets_required_factors = []
list_of_primes_in_partitions = []
list_of_eliminated_primes = [] 
list_no_of_eliminated_primes = []
list_no_of_eliminated_primes_to_no_of_partitions = []
list_no_of_eliminated_primes_to_no_of_partitions_avg = []

list_nums = []
list_no_of_required_factors = []
list_no_of_primes = []
list_no_of_primes_half = []
list_no_of_primes_half_n = []
list_no_of_sets_required_factors = []

def update_metrics (p, dp, num, factors):

    global list_nums, list_of_eliminated_primes, list_no_of_sets_required_factors, list_of_primes_in_partitions
    global list_no_of_primes_half, list_no_of_primes, list_no_of_primes_half_n, list_no_of_required_factors, list_no_of_eliminated_primes

    list_nums.append (num)

    min_lenght = sys.maxsize
    for set_required_factors in list_of_sets_required_factors:
        if len(set_required_factors) < min_lenght:
            min_lenght = len(set_required_factors)

    no_primes = p.get_all_primes_leq(num)
    no_primes_half_n = p.get_all_primes_leq(math.floor(num / 4))
    list_no_of_primes_half.append (math.ceil( no_primes / 4))
    list_no_of_primes.append(no_primes)
    list_no_of_primes_half_n.append(no_primes_half_n)
    list_no_of_required_factors.append (min_lenght)

    list_no_of_sets_required_factors.append (len(list_of_sets_required_factors))

    list_of_primes_in_partitions.sort()

    list_of_eliminated_primes = []
    for p in list_of_primes_in_partitions:
        list_of_eliminated_primes.append (p)

    for set_required_factors in list_of_sets_required_factors:
        for q in set_required_factors:
            if q in list_of_eliminated_primes:
                list_of_eliminated_primes.remove (q)

    list_no_of_eliminated_primes.append (len(list_of_eliminated_primes))

    list_no_of_eliminated_primes_to_no_of_partitions.append(len(list_of_eliminated_primes)/dp.get_number_of_pairs (factors))
    list_no_of_eliminated_primes_to_no_of_partitions_avg.append(dp.get_avg_value_from_list (list_no_of_eliminated_primes_to_no_of_partitions))

#############################################################
# Presentation
#############################################################

def write_results_to_figures (directory):

    global list_nums, list_no_of_primes, list_no_of_primes_half, list_no_of_primes_half_n, list_no_of_required_factors, list_no_of_sets_required_factors, list_no_of_eliminated_primes

    plt.figure(1)
    r_patch = mpatches.Patch(color='red', label='ceil (pi(n)/4)')
    g_patch = mpatches.Patch(color='green', label='pi(n)')
    m_patch = mpatches.Patch(color='magenta', label='pi(n/4)')
    b_patch = mpatches.Patch(color='blue', label='# of req primes for GSC')
    list_of_handles = []
    list_of_handles.append(g_patch)
    list_of_handles.append(r_patch)
    list_of_handles.append(m_patch)
    list_of_handles.append(b_patch)
    plt.legend(handles=list_of_handles, loc='upper right', bbox_to_anchor=(0.5, 0.8))
    plt.plot(list_nums, list_no_of_primes, 'g.', ms=2)
    plt.plot(list_nums, list_no_of_primes_half, 'r.', ms=2)
    plt.plot(list_nums, list_no_of_primes_half_n, 'm.', ms=2)
    plt.plot(list_nums, list_no_of_required_factors, 'b.', ms=2)
    plt.xlabel('n')
    plt.ylabel('Count')
    plt.title('How big are the sets?')
    plt.grid(True)
    plt.savefig(directory + "/f_required_primes.png")

    plt.figure(2)
    plt.plot(list_nums, list_no_of_sets_required_factors, 'b-', ms=1)
    plt.xlabel('n')
    plt.ylabel('Count')
    plt.title('How many sets?')
    plt.grid(True)
    plt.savefig(directory + "/f_number_of_possible_sets.png")

    plt.figure(3)
    plt.plot(list_nums, list_no_of_eliminated_primes, 'b-', ms=1)
    plt.xlabel('n')
    plt.ylabel('Count')
    plt.title('How many eliminated primes?')
    plt.grid(True)
    plt.savefig(directory + "/f_number_of_eliminated_primes.png")

    plt.figure(4)
    blue_patch = mpatches.Patch(color='blue', label='ratio')
    red_patch = mpatches.Patch(color='red', label='avg')
    plt.legend(handles=[red_patch, blue_patch], loc='upper right', bbox_to_anchor=(0.9, 0.9))
    plt.plot(list_nums, list_no_of_eliminated_primes_to_no_of_partitions, 'b-', ms=1)
    plt.plot(list_nums, list_no_of_eliminated_primes_to_no_of_partitions_avg, 'r-', ms=1)
    plt.xlabel('n')
    plt.ylabel('Ratio')
    plt.title('How many eliminated primes to partitions?')
    plt.grid(True)
    plt.savefig(directory + "/f_eliminated_primes_to_partitions.png")

def write_results_to_files (directory):
    global file_output_list_of_eliminated_primes, file_output_list_of_primes_in_partitions, file_output_list_of_sets_required_factors
    
    f = open(file_output_list_of_eliminated_primes, "w")
    f.write (str(list_of_eliminated_primes))
    f.close ()

    f = open(file_output_list_of_primes_in_partitions, "w")
    f.write (str(list_of_primes_in_partitions))
    f.close ()

    f = open(file_output_list_of_sets_required_factors, "w")
    f.write (str(list_of_sets_required_factors))
    f.close ()

#############################################################
# Main - Phase 1
# Preload files & restore previous calculations
#############################################################

print ("---------------------------------------------------")
print ("Initialize objects...")
p = primes.Primes(caching_primality_results)
gp = goldbach.GoldbachPartition (p)
dp = dataprocessing.DataProcessing()
print ("DONE")
print ("Loading helper sets...")
p.init_set(file_input_primes, True)
p.init_set(file_input_nonprimes, False)
print ("DONE")
print ("Sorting primes...")
p.sort_primes_set()
print ("DONE")
print ("Output result folder: ", directory)
print ("---------------------------------------------------")

#############################################################
# Main - Phase 2
# New calculations
#############################################################

dt_start = datetime.now()
dt_current_previous = dt_start

for k in range (min_num, max_num):
    num = step_factor*k

    if be_verbose:
        print ("=============")
        print ("num=", num)
    
    factors = gp.find_sum_of_prime_numbers (num)
    if be_verbose:
        print ("current factors:", factors)

    # step 1:
    # check #1: maybe set_required_factors already contains required factors?
    fullfiled = False
    for pair in factors:
        (p1, p2) = pair
        
        # remember all primes present in partions
        if p1 not in list_of_primes_in_partitions:
            list_of_primes_in_partitions.append (p1)
        if p2 not in list_of_primes_in_partitions:
            list_of_primes_in_partitions.append (p2)
        
        for set_required_factors in list_of_sets_required_factors:
            if not fullfiled and p1 in set_required_factors and p2 in set_required_factors:
                fullfiled = True
                if be_verbose:
                    print ("num=", num, "is fullfiled")
                    print (list_of_sets_required_factors)

    # check #2: p1 or p2 is not on a set_required_factors
    if not fullfiled:

        if be_verbose:
            print ("not fullfiled")
        
        for pair in factors:
            (p1, p2) = pair
            for my_set_required_factors in list_of_sets_required_factors:
                my_temp_set = my_set_required_factors.copy()
                if p1 in my_temp_set and p2 not in my_temp_set:
                    my_temp_set.add (p2)
                    if be_verbose:
                        print ("case 1", my_temp_set)
                    temp_list_of_sets_required_factors.append (my_temp_set)
                    my_temp_set = {}

        for pair in factors:
            (p1, p2) = pair
            for my_set_required_factors in list_of_sets_required_factors:
                my_temp_set = my_set_required_factors.copy()
                if p1 not in my_temp_set and p2 in my_temp_set:
                    my_temp_set.add (p1)
                    if be_verbose:
                        print ("case 2", my_temp_set)
                    temp_list_of_sets_required_factors.append (my_temp_set)
                    my_temp_set = {}

        for pair in factors:
            (p1, p2) = pair
            for my_set_required_factors in list_of_sets_required_factors:
                my_temp_set = my_set_required_factors.copy()
                if p1 not in my_temp_set and p2 not in my_temp_set:
                    my_temp_set.add (p1)
                    my_temp_set.add (p2)
                    if be_verbose:
                        print ("case 3", my_temp_set)
                    temp_list_of_sets_required_factors.append (my_temp_set)
                    my_temp_set = {}
    else:
        for set_required_factors in list_of_sets_required_factors:
            condition_met = False
            for pair in factors:
                (p1, p2) = pair
                if p1 in set_required_factors and p2 in set_required_factors:
                    condition_met = True

            if not condition_met:
                list_of_sets_required_factors.remove (set_required_factors)
                if be_verbose:
                    print ("Removed", set_required_factors, "from list_of_sets_required_factors")
                    print ("Now list_of_sets_required_factors", list_of_sets_required_factors)

    # step 2: cleanup of set_temp_list_required_factors
    if len(temp_list_of_sets_required_factors) > 0:
        min_lenght = sys.maxsize
        for temp_set_required_factors in temp_list_of_sets_required_factors:
            if be_verbose:
                print (temp_set_required_factors)
            if len(temp_set_required_factors) < min_lenght:
                min_lenght = len(temp_set_required_factors)

        list_of_sets_required_factors = []
        for temp_set_required_factors in temp_list_of_sets_required_factors:
            if len(temp_set_required_factors) == min_lenght:
                list_of_sets_required_factors.append (temp_set_required_factors)

        temp_list_of_sets_required_factors = []

    update_metrics (p, dp, num, factors)

    # checkpoint - partial results
    if num % checkpoint_value == 0:
        dt_current = datetime.now()
        dt_diff_current = (dt_current - dt_current_previous).total_seconds()
        print ("Iteration", k, "of total", max_num, "took", dt_diff_current, "seconds")
        # remember results so far
        write_results_to_figures (directory)
        write_results_to_files (directory)

    if be_verbose:
        print ("final result - currently required factors", list_of_sets_required_factors)

dt_end = datetime.now()
write_results_to_figures (directory)

# final results
dt_diff = dt_end - dt_start
print ("Total calculations lasted:", dt_diff)

print (" + Max examined number:", num)
print (" + currently required factors:", list_of_sets_required_factors)
print (" + primes present in partitions:", list_of_primes_in_partitions)
print (" + number of required factors:", list_no_of_required_factors)
print (" + eliminated primes:", list_of_eliminated_primes)

write_results_to_files (directory)
