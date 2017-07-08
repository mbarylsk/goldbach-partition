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
#   o number = min_num * step_factor
min_num = 1
step_factor = 2
# Maximum even number checked against Goldbach conjecture
#   o number = max_num * step_factor
max_num = 200000000

# Checkpoint value when partial results are drawn/displayed
# should be greater than zero
checkpoint_value = 200000

# Type of output figures
#   o True  - both basic and advanced figures
#   o False - basic figures only
create_detailed_figures = False

# Caching previous primality results
#   o True  - auxilary sets of primes and composite numbers will grow
#             it will speed up further primality tests but more RAM will
#             be occupied
#   o False - do not cache new primality test results
caching_primality_results = False

# Algorithms to be checked
algo_to_check = {'a1'}

# Helper files
#   o file_input_primes - contains prime numbers
#   o file_input_nonprimes - contains composite numbers
file_input_primes = '..\\primes\\t_prime_numbers.txt'
file_input_nonprimes = '..\\primes\\t_nonprime_numbers.txt'

#############################################################
# Settings - output directory and files
#############################################################

directory = str(step_factor*max_num)
if not os.path.exists(directory):
    os.makedirs(directory)
algo = ""
if 'a1' in algo_to_check:
    algo += "a1"
file_output_iters_alg = directory + "/f_checkpoint_diff_iters_alg" + algo + ".png"
file_output_iters_alg_stats = directory + "/f_checkpoint_diff_iters_stats_alg" + algo + ".png"
file_output_ind_iters_alg = directory + "/f_checkpoint_diff_ind_iters_alg" + algo + ".png"
file_output_max_iters_alg = directory + "/f_checkpoint_diff_max_iters_alg" + algo + ".png"
file_output_duration_alg = directory + "/f_checkpoint_diff_duration_alg" + algo + ".png"
file_output_pickle = directory + "/objs.pickle"

#############################################################
# Results of calculations
#############################################################

list_nums = []
list_iterations = []
list_iterations_max = []
list_iterations_avg = []
list_iterations_diff = []

dt_diff = [0, 0, 0, 0, 0, 0]
dt_iter = [0, 0, 0, 0, 0, 0]

list_checkpoints_duration = [[], [], [], [], [], []]
list_checkpoints_iters = [[], [], [], [], [], []]
list_checkpoints_iters_max = [[], [], [], [], [], []]
list_checkpoints_iters_avg = [[], [], [], [], [], []]
list_checkpoints = []

max_iterations = 0
max_iteration_details = ""
avg_iterations = 0
diff_iterations = 0
previous_iterations = 0
loops = 0
k_current = 0

#############################################################
# Presentation
#############################################################

def write_results_to_figures (directory):
    # results - figures    
    plt.figure(1)
    g_patch = mpatches.Patch(color='blue', label='A1')

    list_of_handles = []
    if 'a1' in algo_to_check:
        plt.plot(list_checkpoints, list_checkpoints_duration[0], 'b.', ms=2)
        list_of_handles.append(g_patch)

    plt.legend(handles=list_of_handles, loc='upper right', bbox_to_anchor=(0.4, 0.8))
    plt.xlabel('n')
    plt.ylabel('Time [s]')
    plt.title('Duration of total calculations')
    plt.grid(True)
    plt.savefig(file_output_duration_alg)

    plt.figure(2)
    if 'a1' in algo_to_check:
        plt.plot(list_checkpoints, list_checkpoints_iters[0], 'b.', ms=2)

    plt.legend(handles=list_of_handles, loc='upper right', bbox_to_anchor=(0.4, 0.8))
    plt.xlabel('n')
    plt.ylabel('Sum(I(n))')
    plt.title('Total iterations')
    plt.grid(True)
    plt.savefig(file_output_iters_alg)

    plt.figure(3)
    if 'a1' in algo_to_check:
        plt.plot(list_checkpoints, list_checkpoints_iters_max[0], 'r.', ms=2)
        plt.plot(list_checkpoints, list_checkpoints_iters_avg[0], 'g.', ms=2)

    ri_patch = mpatches.Patch(color='red', label='max')
    gi_patch = mpatches.Patch(color='green', label='avg')
    list_of_handles_iter = []
    list_of_handles_iter.append(ri_patch)
    list_of_handles_iter.append(gi_patch)

    plt.legend(handles=list_of_handles_iter, loc='upper right', bbox_to_anchor=(0.4, 0.8))
    plt.xlabel('n')
    plt.ylabel('f(I(n))')
    plt.title('Max and average iterations')
    plt.grid(True)
    plt.savefig(file_output_iters_alg_stats)

    if create_detailed_figures:
        plt.figure(4)
        if 'a1' in algo_to_check:
            plt.plot(list_nums, list_iterations, 'b.', ms=2)

        plt.legend(handles=list_of_handles, loc='upper right', bbox_to_anchor=(0.4, 0.8))
        plt.xlabel('n')
        plt.ylabel('I(n)')
        plt.title('Iterations per number')
        plt.grid(True)
        plt.savefig(file_output_ind_iters_alg)

        plt.figure(5)
        if 'a1' in algo_to_check:
            plt.plot(list_nums, list_iterations_diff, 'g-', ms=1)
            plt.plot(list_nums, list_iterations_max, 'r.', ms=2)
            plt.plot(list_nums, list_iterations_avg, 'b.', ms=2)

        ri_patch = mpatches.Patch(color='red', label='max')
        bi_patch = mpatches.Patch(color='blue', label='avg')
        gi_patch = mpatches.Patch(color='green', label='delta')

        list_of_handles_iter = []
        list_of_handles_iter.append(ri_patch)
        list_of_handles_iter.append(bi_patch)
        list_of_handles_iter.append(gi_patch)

        plt.legend(handles=list_of_handles_iter, loc='upper right', bbox_to_anchor=(0.4, 0.8))
        plt.xlabel('n')
        plt.ylabel('f(I(n))')
        plt.title('Iteration statistics')
        plt.grid(True)
        plt.savefig(file_output_max_iters_alg)

def save_current_results (file_output_pickle):
    global k_current, max_iterations, max_iteration_details, avg_iterations, diff_iterations, previous_iterations, loops, dt_diff, dt_iter, list_checkpoints, list_checkpoints_iters, list_checkpoints_iters_max, list_checkpoints_iters_avg, list_checkpoints_duration
    with open(file_output_pickle, 'wb') as f:
        pickle.dump([k_current, max_iterations, max_iteration_details, avg_iterations, diff_iterations, previous_iterations, loops, dt_diff, dt_iter, list_checkpoints, list_checkpoints_iters, list_checkpoints_iters_max, list_checkpoints_iters_avg, list_checkpoints_duration], f)

def restore_previous_results (file_output_pickle):
    global k_current, max_iterations, max_iteration_details, avg_iterations, diff_iterations, previous_iterations, loops, dt_diff, dt_iter, list_checkpoints, list_checkpoints_iters, list_checkpoints_iters_max, list_checkpoints_iters_avg, list_checkpoints_duration
    if os.path.exists(file_output_pickle):
        with open(file_output_pickle, 'rb') as f:
            k_current, max_iterations, max_iteration_details, avg_iterations, diff_iterations, previous_iterations, loops, dt_diff, dt_iter, list_checkpoints, list_checkpoints_iters, list_checkpoints_iters_max, list_checkpoints_iters_avg, list_checkpoints_duration = pickle.load(f)

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
print ("Restoring previous results...")
restore_previous_results (file_output_pickle)
if k_current > 0:
    min_num = k_current
    print ("Resuming calculations at", min_num)
print ("DONE")

dt_start = datetime.now()
dt_current_previous = dt_start

# new calculations
for k in range (min_num, max_num):
    num = step_factor*k
    loops += 1

    if create_detailed_figures:
        list_nums.append (num)

    # algorithm 1
    if 'a1' in algo_to_check:
        (p1, p2, duration, iterations) = gp.search_for_difference (num)
        dt_diff[0] += duration
        dt_iter[0] += iterations

        avg_iterations = int (dt_iter[0]/loops)
        diff_iterations = abs (previous_iterations - iterations)
        if iterations > max_iterations:
            max_iteration_details = str (num) + "=" +str(p2) + "-" + str(p1) + "(max:" + str(iterations) + ", avg:" + str(avg_iterations) + ")"
            max_iterations = iterations

        if create_detailed_figures:
            list_iterations.append(iterations)
            list_iterations_max.append(max_iterations)
            list_iterations_avg.append(avg_iterations)
            list_iterations_diff.append(diff_iterations)

        if p2 - p1 != num:
            print ("Alg #1: violation of diff for n=", num, "p2=", p2, "p1=", p1)

        previous_iterations = iterations

    # checkpoint - partial results
    if num % checkpoint_value == 0:
        dt_current = datetime.now()
        dt_diff_current = (dt_current - dt_current_previous).total_seconds()

        list_checkpoints.append(num)
        if 'a1' in algo_to_check:
            list_checkpoints_duration[0].append(dt_diff[0])
            list_checkpoints_iters[0].append(dt_iter[0])
            list_checkpoints_iters_avg[0].append(avg_iterations)
            list_checkpoints_iters_max[0].append(max_iterations)

        perc_completed = str(int(k * 100 / max_num))
        print ("Checkpoint", k, "of total", max_num, "took", dt_diff_current, "seconds. (" + perc_completed + "% completed)", max_iteration_details)
        
        # save results collected so far
        write_results_to_figures (directory)
        k_current = k
        save_current_results(file_output_pickle)

dt_end = datetime.now()

# final results - time of processing
dt_diff = dt_end - dt_start
print ("Total calculations lasted:", dt_diff)

# final results
write_results_to_figures (directory)
k_current = max_num
save_current_results(file_output_pickle)
