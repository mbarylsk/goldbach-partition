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
#   o False - do not cache new primality test results
caching_primality_results = False

# Check basic stats
check_regular = False
# Check stats related to twin primes
check_twins = False
check_twins_extended = False

# Check hypothesis related to n = p1 + p2
# For every even n>=10
#  if n= 6k then p1= 6x-1 and p2= 6y+1
#  if n= 6k+2 then p1= 6x+1 and p2= 6y+1
#  if n= 6k+4 then p1= 6x-1 and p2= 6y-1
check_6kpm1_hypothesis = True

min_num = 2
max_num = 1000000
step_factor = 2
checkpoint_value = 10000
big_prime_threshold = 100
min_prime_count_threshold = 500
file_input_primes = '..\\primes\\t_prime_numbers.txt'
file_input_nonprimes = '..\\primes\\t_nonprime_numbers.txt'
file_input_oldpairs = 'input\\t_prime_sum_pairs.txt'

#############################################################
# Settings - output directory and files
#############################################################

directory = "results/" + str(step_factor*max_num)
if not os.path.exists(directory):
    os.makedirs(directory)

file_output_primes_lower_than = directory + "/t_prime_numbers_lower_than_" + str(step_factor*max_num) + ".txt"
file_output_nonprimes_lower_than = directory + "/t_nonprime_numbers_lower_than_" + str(step_factor*max_num) + ".txt"
file_output_primes_frequency_lower_than = directory + "/t_min_prime_numbers_freq_lower_than_" + str(step_factor*max_num) + ".txt"
file_output_prime_sum_pairs_below = directory + "/t_prime_sum_pairs_below_" + str(step_factor*max_num) + ".txt"
file_output_big_smallest_primes_below = directory + "/t_big_smallest_primes_below_" + str(step_factor*max_num) + ".txt"

#############################################################
# Business logic
#############################################################

list_min_sats = []
list_max_sats = []
list_avg_sats = []
list_num_distinct_pairs = []
list_max_diff_in_pairs = []
list_min_diff_in_pairs = []
list_avg_diff_in_pairs = []
list_max_diff_in_pairs_trend = []
list_min_diff_in_pairs_trend = []
list_avg_diff_in_pairs_trend = []
list_min_prime = []
list_nums = []
list_checkpoints_duration = []
list_checkpoints = []
prev_max_diff = 0
prev_min_diff = 0
prev_avg_diff = 0
max_diff_trend_factor = 0
min_diff_trend_factor = 0
avg_diff_trend_factor = 0
dict_min_primes_count = dict()
min_prime = 0
list_num_twin_greater = []
list_num_twin_lesser = []
list_num_twin_all = []
list_num_twin_lesser_distinct = []
list_num_twin_greater_distinct = []
list_num_diff_twin_greater_lesser_distinct = []
list_num_diff_twin_greater_lesser_distinct_avg = []
list_num_ratio_twin_gp = []
list_num_ratio_twin_gp_avg = []
list_zero_twin_lesser = []
list_zero_twin_greater = []
counter_zero_twin_lesser = 0
counter_zero_twin_greater = 0
list_zero_twin_lesser_count = []
list_zero_twin_greater_count = []
list_count_6kpm1 = []
list_count_6kpm1_ratio = []

def calculate_metrics (num, factors, dp, p):
    global prev_max_diff, max_diff_trend_factor, prev_min_diff, min_diff_trend_factor, prev_avg_diff, avg_diff_trend_factor, min_prime
    global counter_zero_twin_lesser, counter_zero_twin_greater

    list_nums.append(num)
    num_of_pairs = dp.get_number_of_pairs (factors)

    if check_regular:
        min_prime = dp.get_min_factor (factors)
        if min_prime not in dict_min_primes_count:
            dict_min_primes_count[min_prime] = 1
        else:
            counter = dict_min_primes_count[min_prime]
            counter = counter + 1
            dict_min_primes_count[min_prime] = counter
        
        diffs_in_factors = dp.get_diff_in_factors (factors)
        
        min_diff = dp.get_min_value_from_list (diffs_in_factors)
        max_diff = dp.get_max_value_from_list (diffs_in_factors)
        avg_diff = dp.get_avg_value_from_list (diffs_in_factors)

        max_sat = dp.get_perc_max_saturation_from_factors (num, factors)
        min_sat = dp.get_perc_min_saturation_from_factors (num, factors)
        avg_sat = dp.get_perc_avg_saturation_from_factors (num, factors)
        
        list_max_sats.append(max_sat)
        list_avg_sats.append(avg_sat)
        list_min_sats.append(min_sat)
        list_min_prime.append(min_prime)
        list_num_distinct_pairs.append(num_of_pairs)
        list_max_diff_in_pairs.append(max_diff)
        list_min_diff_in_pairs.append(min_diff)
        list_avg_diff_in_pairs.append(avg_diff)

        if prev_max_diff > max_diff:
            max_diff_trend_factor = max_diff_trend_factor - 1
        elif prev_max_diff < max_diff:
            max_diff_trend_factor = max_diff_trend_factor + 1
        list_max_diff_in_pairs_trend.append(max_diff_trend_factor)

        if prev_min_diff > min_diff:
            min_diff_trend_factor = min_diff_trend_factor - 1
        elif prev_max_diff < min_diff:
            min_diff_trend_factor = min_diff_trend_factor + 1
        list_min_diff_in_pairs_trend.append(min_diff_trend_factor)

        if prev_avg_diff > avg_diff:
            avg_diff_trend_factor = avg_diff_trend_factor - 1
        elif prev_max_diff < avg_diff:
            avg_diff_trend_factor = avg_diff_trend_factor + 1
        list_avg_diff_in_pairs_trend.append(avg_diff_trend_factor)

        prev_max_diff = max_diff
        prev_min_diff = min_diff
        prev_avg_diff = avg_diff

    if check_twins:
        number_of_twins_lesser = 0
        number_of_twins_greater = 0
        number_of_twins = 0
        min_twin_lesser = 0
        min_twin_greater = 0
        set_twins_lesser = set()
        set_twins_greater = set()

        for (p1, p2) in factors:
        
            if p.is_lesser_twin_prime (p1):
                if p1 not in set_twins_lesser:
                    set_twins_lesser.add (p1)
                number_of_twins += 1
            if p.is_lesser_twin_prime (p2):
                if p2 not in set_twins_lesser:
                    set_twins_lesser.add (p2)
                number_of_twins += 1
            if p.is_greater_twin_prime (p1):
                if p1 not in set_twins_greater:
                    set_twins_greater.add (p1)
                number_of_twins += 1
            if p.is_greater_twin_prime (p2):
                if p2 not in set_twins_greater:
                    set_twins_greater.add (p2)
                number_of_twins += 1
            
            if (p1 != p2):
                # lesser twins
                if p.is_lesser_twin_prime (p1):
                    number_of_twins_lesser += 1
                    if number_of_twins_lesser == 1:
                        min_twin_lesser = p1
                    elif p1 < min_twin_lesser:
                        min_twin_lesser = p1
                if p.is_lesser_twin_prime (p2):
                    number_of_twins_lesser += 1
                    if number_of_twins_lesser == 1:
                        min_twin_lesser = p2
                    elif p2 < min_twin_lesser:
                        min_twin_lesser = p2
                # greater twins
                if p.is_greater_twin_prime (p1):
                    number_of_twins_greater += 1
                    if number_of_twins_greater == 1:
                        min_twin_greater = p1
                    elif p1 < min_twin_greater:
                        min_twin_greater = p1
                if p.is_greater_twin_prime (p2):
                    number_of_twins_greater += 1
                    if number_of_twins_greater == 1:
                        min_twin_greater = p2
                    elif p2 < min_twin_greater:
                        min_twin_greater = p2
            else:
                # lesser twins
                if p.is_lesser_twin_prime (p1):
                    number_of_twins_lesser += 1
                    if number_of_twins_lesser == 1:
                        min_twin_lesser = p1
                    elif p1 < min_twin_lesser:
                        min_twin_lesser = p1
                # greater twins
                if p.is_greater_twin_prime (p1):
                    number_of_twins_greater += 1
                    if number_of_twins_greater == 1:
                        min_twin_greater = p1
                    elif p1 < min_twin_greater:
                        min_twin_greater = p1

        list_num_twin_all.append (number_of_twins)
        list_num_ratio_twin_gp.append(number_of_twins/(num_of_pairs*2))
        list_num_ratio_twin_gp_avg.append(dp.get_avg_value_from_list (list_num_ratio_twin_gp))
        
    if check_twins_extended:
        number_of_twins_lesser_distinct = len(set_twins_lesser)
        number_of_twins_greater_distinct = len(set_twins_greater)
        list_num_twin_lesser.append (number_of_twins_lesser)
        list_num_twin_greater.append (number_of_twins_greater)
        list_num_twin_lesser_distinct.append (number_of_twins_lesser_distinct)
        list_num_twin_greater_distinct.append (number_of_twins_greater_distinct)
        list_num_diff_twin_greater_lesser_distinct.append (number_of_twins_greater_distinct-number_of_twins_lesser_distinct)
        list_num_diff_twin_greater_lesser_distinct_avg.append (dp.get_avg_value_from_list (list_num_diff_twin_greater_lesser_distinct))
        
        if number_of_twins_lesser_distinct == 0:
            list_zero_twin_lesser.append(num)
            counter_zero_twin_lesser += 1
            list_zero_twin_lesser_count.append(counter_zero_twin_lesser)
        if number_of_twins_greater_distinct == 0:
            list_zero_twin_greater.append(num)
            counter_zero_twin_greater += 1
            list_zero_twin_greater_count.append(counter_zero_twin_greater)

    if check_6kpm1_hypothesis:
        count = 0
        if num > 8:
            for (p1, p2) in factors:
                if num/2 % 3 == 0:
                    if (p.is_6km1 (p1) and p.is_6kp1 (p2)) or (p.is_6km1 (p2) and p.is_6kp1 (p1)):
                        count += 1
                if num/2 % 3 == 1:
                    if (p.is_6kp1 (p1) and p.is_6kp1 (p2)):
                        count += 1
                if num/2 % 3 == 2:
                    if (p.is_6km1 (p1) and p.is_6km1 (p2)):
                        count += 1

        list_count_6kpm1.append(count)
        if (count > 0):
            list_count_6kpm1_ratio.append(num_of_pairs/count)
        else:
            list_count_6kpm1_ratio.append(0)
        if count == 0:
            print ("WARNING: 6kpm1 Condition not met for", num)

#############################################################
# Presentation
#############################################################

def print_iteration_output ():
    global num, factors, max_sat, min_sat, diffs_in_factors, num_of_pairs, prev_max_diff, prev_min_diff, prev_avg_diff, max_diff_trend_factor, min_diff_trend_factor, avg_diff_trend_factor
    print ("=============================================================")
    print ("n to be analyzed: %d" % num)
    print ("n = p1 + p2: ", factors)
    print ("max saturation: ", max_sat)
    print ("min saturation: ", min_sat)
    print ("diff in factors: ", diffs_in_factors)
    print ("# of pairs:", num_of_pairs)
    print (prev_max_diff, prev_min_diff, prev_avg_diff)
    print (max_diff_trend_factor, min_diff_trend_factor, avg_diff_trend_factor)

def write_results_to_files (directory, p):
    global file_output_primes_lower_than, file_output_nonprimes_lower_than, file_output_primes_frequency_lower_than
    
    f = open(file_output_primes_lower_than, "w")
    f.write (str(p.get_list_sorted_prime()))
    f.close ()

    f = open(file_output_nonprimes_lower_than, "w")
    f.write (str(get_list_sorted_nonprime()))
    f.close ()

    f = open(file_output_primes_frequency_lower_than, "w")
    for key in sorted(dict_min_primes_count.keys()):
        ikey = int(key)
        ivalue = int(dict_min_primes_count[key])
        f.write (str(ikey) + "->" + str(ivalue) + "\n")
    f.close ()

def write_auxiliary_results_to_file (directory, num, factors):
    global big_prime_threshold, min_prime, file_output_prime_sum_pairs_below, file_output_big_smallest_primes_below

    # prime numbers
    f = open(file_output_prime_sum_pairs_below, "a+")
    f.write ("Number: " + str(num) + "\n")
    f.write (" Pairs: " + str(factors) + "\n")
    f.close ()

    # search for big minimal primes
    if min_prime > big_prime_threshold:
        f = open(file_output_big_smallest_primes_below, "a+")
        f.write ("Number: " + str(num) + "\n")
        f.write (" Minimum prime: " + str(min_prime) + "\n")
        f.close ()

def read_results_from_file (file_input_oldpairs, dp, p):
    
    f = open(file_input_oldpairs, "r")
    lines = f.readlines()
    read_pairs = False
    for line in lines:
        if read_pairs:
            read_pairs = False

            # read factors from file and calculate metrics
            factors_from_file = dp.read_sums_from_line (line)
            calculate_metrics (num_from_file, factors_from_file, dp, p)

            write_auxiliary_results_to_file (directory, num_from_file, factors_from_file)         
        else:
            read_pairs = True

            # read number from file
            num_from_file = dp.read_num_from_line (line)

            if num_from_file % checkpoint_value == 0:
                print (".")

    return num_from_file

def write_results_to_figures (directory, last_loop):
    # results - figures

    if check_regular:
        plt.figure(1)
        plt.plot(list_nums, list_max_sats, 'b.', ms=2)
        plt.plot(list_nums, list_min_sats, 'r.', ms=2)
        plt.plot(list_nums, list_avg_sats, 'g.', ms=2)
        blue_patch = mpatches.Patch(color='blue', label='max(p1, p2)/number*100%')
        red_patch = mpatches.Patch(color='red', label='min(p1, p2)/number*100%')
        green_patch = mpatches.Patch(color='green', label='avg(p1, p2)/number*100%')
        plt.legend(handles=[red_patch, blue_patch, green_patch], loc='upper right', bbox_to_anchor=(0.8, 0.8))
        plt.xlabel('Number')
        plt.ylabel('Percentage of number')
        plt.title('Participation of max/min/avg prime from Goldbach partition in number')
        plt.grid(True)
        plt.savefig(directory + "/f_saturation.png")

        plt.figure(2)
        plt.plot(list_nums, list_num_distinct_pairs, 'bo', ms=2)
        plt.xlabel('Number')
        plt.ylabel('Number of Goldbach partitions for a given number')
        plt.title('Number of Goldbach partitions')
        plt.grid(True)
        plt.savefig(directory + "/f_pairs.png")

        plt.figure(3)
        plt.plot(list_nums, list_max_diff_in_pairs, 'bo', ms=2)
        plt.xlabel('Number')
        plt.ylabel('Maximum difference in Goldbach partition')
        plt.title('Maximum difference in Goldbach partition for a given number')
        plt.grid(True)
        plt.savefig(directory + "/f_max_diff_pairs.png")

        plt.figure(4)
        plt.plot(list_nums, list_min_diff_in_pairs, 'bo', ms=2)
        plt.xlabel('Number')
        plt.ylabel('Minimum difference in Goldbach partition')
        plt.title('Minimum difference in Goldbach partition for a given number')
        plt.grid(True)
        plt.savefig(directory + "/f_min_diff_pairs.png")

        plt.figure(5)
        plt.plot(list_nums, list_avg_diff_in_pairs, 'bo', ms=2)
        plt.xlabel('Number')
        plt.ylabel('Average difference in Goldbach partition')
        plt.title('Average difference in Goldbach partition for a given number')
        plt.grid(True)
        plt.savefig(directory + "/f_avg_diff_pairs.png")

        plt.figure(6)
        plt.plot(list_nums, list_max_diff_in_pairs_trend, 'b.', ms=1, label='Trend in max diff')
        plt.plot(list_nums, list_min_diff_in_pairs_trend, 'r.', ms=1, label='Trend in min diff')
        plt.plot(list_nums, list_avg_diff_in_pairs_trend, 'g.', ms=1, label='Trend in avg diff')
        blue_patch = mpatches.Patch(color='blue', label='Trend in max diff')
        red_patch = mpatches.Patch(color='red', label='Trend in min diff')
        green_patch = mpatches.Patch(color='green', label='Trend in avg diff')
        plt.legend(handles=[red_patch, blue_patch, green_patch])
        plt.xlabel('Number')
        plt.ylabel('Change in trend')
        plt.title('Trends in Goldbach partition')
        plt.grid(True)
        plt.savefig(directory + "/f_max_diff_pairs_trend.png")

        plt.figure(7)
        plt.plot(list_checkpoints, list_checkpoints_duration)
        plt.xlabel('Checkpoint [iteration]')
        plt.ylabel('Time [s]')
        plt.title('Duration of total calculations between checkpoints in seconds')
        plt.grid(True)
        plt.savefig(directory + "/f_checkpoint_duration.png")

        plt.figure(8)
        plt.plot(list_nums, list_min_prime, 'bo', ms=2)
        plt.xlabel('Number')
        plt.ylabel('Minimum prime number in Goldbach partition')
        plt.title('Minimum prime number in Goldbach partition for a given number')
        plt.grid(True)
        plt.savefig(directory + "/f_min_prime_in_sum.png")

        plt.figure(9)
        list_primes = []
        list_primes_count = []
        for key in sorted(dict_min_primes_count.keys()):
            ikey = int(key)
            ivalue = int(dict_min_primes_count[key])
            list_primes.append(ikey)
            list_primes_count.append(ivalue)
            if last_loop:
                if ivalue > min_prime_count_threshold:
                    plt.annotate(str(ikey)+":"+str(ivalue), xy=(ikey,ivalue), xytext=(5,5), textcoords='offset points')
        plt.plot(list_primes, list_primes_count, 'bo-', ms=3)
        plt.xlabel('Prime number')
        plt.ylabel('Apperances')
        plt.title('Frequency of minimum primes in Goldbach partition')
        plt.grid(True)
        plt.savefig(directory + "/f_min_prime_in_sum_histogram.png")

    if check_twins_extended:
        plt.figure(10)
        plt.plot(list_nums, list_num_twin_lesser, 'b.', ms=1)
        plt.xlabel('n')
        plt.ylabel('Number of lesser twin primes')
        plt.title('Number of lesser twin primes in distinct Goldbach partition of n')
        plt.grid(True)
        plt.savefig(directory + "/f_twin_primes_lesser.png")

        plt.figure(11)
        plt.plot(list_nums, list_num_twin_greater, 'b.', ms=1)
        plt.xlabel('n')
        plt.ylabel('Number of greater twin primes')
        plt.title('Number of greater twin primes in Goldbach partition of n')
        plt.grid(True)
        plt.savefig(directory + "/f_twin_primes_greater.png")

        plt.figure(13)
        plt.plot(list_nums, list_num_twin_lesser_distinct, 'b.', ms=1)
        plt.xlabel('n')
        plt.ylabel('Number of lesser twin primes')
        plt.title('Number of distinct lesser twin primes in Goldbach partition of n')
        plt.grid(True)
        plt.savefig(directory + "/f_twin_primes_lesser_distinct.png")

        plt.figure(14)
        plt.plot(list_nums, list_num_twin_greater_distinct, 'b.', ms=1)
        plt.xlabel('n')
        plt.ylabel('Number of greater twin primes')
        plt.title('Number of distinct greater twin primes in Goldbach partition of n')
        plt.grid(True)
        plt.savefig(directory + "/f_twin_primes_greater_distinct.png")

        fig = plt.figure(15)
        ax1 = fig.add_subplot(111)
        ax1.bar(list_nums, list_num_twin_greater_distinct, color='blue', edgecolor='blue', bottom=list_num_twin_greater_distinct)
        ax1.bar(list_nums, list_num_twin_lesser_distinct, color='red', edgecolor='red')
        plt.xlabel('n')
        plt.ylabel('Number of primes')
        plt.title('Number of distinct twin primes in Goldbach partition of n')
        plt.grid(True)
        plt.savefig(directory + "/f_twin_primes_lesser_greater_distinct.png")

        plt.figure(17)
        plt.plot(list_zero_twin_lesser, list_zero_twin_lesser_count, 'b.', ms=1)
        plt.xlabel('n')
        plt.ylabel('Count')
        plt.title('Count of even numbers without lesser twin prime in GP')
        plt.grid(True)
        plt.savefig(directory + "/f_counter_no_lesser_twin_prime.png")

        plt.figure(18)
        plt.plot(list_zero_twin_greater, list_zero_twin_greater_count, 'r.', ms=1)
        plt.xlabel('n')
        plt.ylabel('Count')
        plt.title('Count of even numbers without greater twin prime in GP')
        plt.grid(True)
        plt.savefig(directory + "/f_counter_no_greater_twin_prime.png")

        plt.figure(19)
        plt.bar(list_zero_twin_lesser, [1]*len(list_zero_twin_lesser), color='blue')
        plt.bar(list_zero_twin_greater, [1]*len(list_zero_twin_greater), color='red')
        blue_patch = mpatches.Patch(color='blue', label='no of lesser twin primes')
        red_patch = mpatches.Patch(color='red', label='no of greater twin primes')
        plt.legend(handles=[red_patch, blue_patch], loc='upper right', bbox_to_anchor=(0.8, 0.8))
        plt.xlabel('n')
        plt.ylabel('Count')
        plt.title('Numbers without lesser and greater twin prime in GP')
        plt.grid(True)
        plt.savefig(directory + "/f_no_twin_prime.png")

        plt.figure(20)
        plt.plot(list_nums, list_num_diff_twin_greater_lesser_distinct, 'b.', ms=1)
        plt.plot(list_nums, list_num_diff_twin_greater_lesser_distinct_avg, 'r.', ms=1)
        blue_patch = mpatches.Patch(color='blue', label='diff')
        red_patch = mpatches.Patch(color='red', label='avg')
        plt.legend(handles=[red_patch, blue_patch], loc='upper right', bbox_to_anchor=(0.9, 0.9))
        plt.xlabel('n')
        plt.ylabel('Diff')
        plt.title('Diff between distinct greater and lesser twin primes in GP')
        plt.grid(True)
        plt.savefig(directory + "/f_diff_greater_lesser_twin_primes.png")

    if check_twins:
        plt.figure(12)
        plt.plot(list_nums, list_num_twin_all, 'b.', ms=1)
        plt.xlabel('n')
        plt.ylabel('Number of twin primes')
        plt.title('Number of twin primes in Goldbach partition of n')
        plt.grid(True)
        plt.savefig(directory + "/f_twin_primes_all.png")

        plt.figure(16)
        plt.plot(list_nums, list_num_ratio_twin_gp, 'b.', ms=1)
        plt.plot(list_nums, list_num_ratio_twin_gp_avg, 'r.', ms=1)
        blue_patch = mpatches.Patch(color='blue', label='ratio')
        red_patch = mpatches.Patch(color='red', label='avg')
        plt.legend(handles=[red_patch, blue_patch], loc='upper right', bbox_to_anchor=(0.9, 0.9))
        plt.xlabel('n')
        plt.ylabel('Ratio')
        plt.title('Ratio: # of twin primes / # of primes in Goldbach partions')
        plt.grid(True)
        plt.savefig(directory + "/f_ratio_twin_primes_to_all_primes.png")

    if check_6kpm1_hypothesis:
        plt.figure(21)
        plt.plot(list_nums, list_count_6kpm1, 'b.', ms=1)
        plt.xlabel('n')
        plt.ylabel('Count of matching pairs')
        plt.title('Hypothesis for 6k+-1')
        plt.grid(True)
        plt.savefig(directory + "/f_hypo_6kpm1.png")

        plt.figure(22)
        plt.plot(list_nums, list_count_6kpm1_ratio, 'b-', ms=1)
        plt.xlabel('n')
        plt.ylabel('Ratio')
        plt.title('Ratio of GB pairs to pairs matching 6k+-1 hypothesis')
        plt.grid(True)
        plt.savefig(directory + "/f_hypo_6kpm1_ratio.png")

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

if os.path.exists(file_input_oldpairs):
    print ("Restoration of previous calculations started ...")
    
    num_from_file = read_results_from_file (file_input_oldpairs, dp, p)
    min_num = int(num_from_file / step_factor)
    print ("Restoration of previous", min_num, "calculations has been completed.")

#############################################################
# Main - Phase 2
# New calculations
#############################################################

dt_start = datetime.now()
dt_current_previous = dt_start

for k in range (min_num, max_num):
    num = step_factor*k
    
    factors = gp.find_sum_of_prime_numbers (num)

    calculate_metrics (num, factors, dp, p)
    
    # checkpoint - partial results
    if num % checkpoint_value == 0:
        dt_current = datetime.now()
        dt_diff_current = (dt_current - dt_current_previous).total_seconds()
        list_checkpoints.append(num)
        list_checkpoints_duration.append(dt_diff_current)
        print ("Iteration", k, "of total", max_num, "took", dt_diff_current, "seconds")
        # remember results so far
        write_results_to_figures (directory, False)

    write_auxiliary_results_to_file (directory, num, factors)

    if be_verbose:
        print_iteration_output ()

dt_end = datetime.now()

# final results
dt_diff = dt_end - dt_start
print ("Total calculations lasted:", dt_diff)

#write_results_to_files (directory)
write_results_to_figures (directory, False)
