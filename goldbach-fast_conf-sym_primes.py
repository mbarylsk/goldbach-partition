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
min_num = 3
step_factor = 1
# Maximum even number checked against Goldbach conjecture
#   o number = max_num * step_factor
max_num = 40000

# Checkpoint value when partial results are drawn/displayed
checkpoint_value = 500

# Caching previous primality results
#   o True  - auxilary sets of primes and composite numbers will grow
#             it will speed up further primality tests but more RAM will
#             be occupied
#   o False - do not cache new primality test results
caching_primality_results = False

# Algorithms to be checked
algo_to_check = {'a1', 'a2', 'a3'}

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
algo = ""
if 'a1' in algo_to_check:
    algo += "a1"
if 'a2' in algo_to_check:
    algo += "a2"
if 'a3' in algo_to_check:
    algo += "a3"
file_output_iters_alg = directory + "/f_checkpoint_sym_primes_iters_alg" + algo + ".png"
file_output_iters_prob_alg = directory + "/f_checkpoint_sym_primes_iters_prob_alg" + algo + ".png"
file_output_duration_alg = directory + "/f_checkpoint_sym_primes_duration_alg" + algo + ".png"
file_output_pickle = directory + "/objs_sym_primes.pickle"

#############################################################
# Results of calculations
#############################################################

dt_diff = [0, 0, 0]
dt_iter = [0, 0, 0]

list_checkpoints_duration = [[], [], []]
list_checkpoints_iters = [[], [], []]
list_checkpoints = []

k_current = 0

#############################################################
# Presentation
#############################################################
            
def update_algo_results (idc, duration, iterations):
    global dt_diff, dt_iter, dt_iter_which
    dt_diff[idc] += duration
    dt_iter[idc] += iterations

def write_results_to_figures (directory):
    # results - figures    
    fig1 = plt.figure(1)
    g_patch = mpatches.Patch(color='green', label='A1')
    r_patch = mpatches.Patch(color='red', label='A2')
    b_patch = mpatches.Patch(color='blue', label='A3')

    list_of_handles = []
    if 'a1' in algo_to_check:
        plt.plot(list_checkpoints, list_checkpoints_duration[0], 'g.', ms=2)
        list_of_handles.append(g_patch)
    if 'a2' in algo_to_check:
        plt.plot(list_checkpoints, list_checkpoints_duration[1], 'r.', ms=2)
        list_of_handles.append(r_patch)
    if 'a3' in algo_to_check:
        plt.plot(list_checkpoints, list_checkpoints_duration[2], 'b.', ms=2)
        list_of_handles.append(b_patch)

    plt.legend(handles=list_of_handles, loc='upper right', bbox_to_anchor=(0.4, 0.8))
    plt.xlabel('n')
    plt.ylabel('time [s]')
    plt.title('Duration of total calculations')
    plt.grid(True)
    plt.savefig(file_output_duration_alg)
    plt.close(fig1)

    fig2 = plt.figure(2)
    if 'a1' in algo_to_check:
        plt.plot(list_checkpoints, list_checkpoints_iters[0], 'g.', ms=2)
    if 'a2' in algo_to_check:
        plt.plot(list_checkpoints, list_checkpoints_iters[1], 'r.', ms=2)
    if 'a3' in algo_to_check:
        plt.plot(list_checkpoints, list_checkpoints_iters[2], 'b.', ms=2)

    plt.legend(handles=list_of_handles, loc='upper right', bbox_to_anchor=(0.4, 0.8))
    plt.xlabel('n')
    plt.ylabel('sum(I(n))')
    plt.title('Total iterations')
    plt.grid(True)
    plt.savefig(file_output_iters_alg)
    plt.close(fig2)

def save_current_results (file_output_pickle):
    global k_current, dt_diff, dt_iter, list_checkpoints_duration, list_checkpoints_iters, list_checkpoints
    with open(file_output_pickle, 'wb') as f:
        pickle.dump([k_current, dt_diff, dt_iter, list_checkpoints_duration, list_checkpoints_iters, list_checkpoints], f)

def restore_previous_results (file_output_pickle):
    global k_current, dt_diff, dt_iter, list_checkpoints_duration, list_checkpoints_iters, list_checkpoints
    if os.path.exists(file_output_pickle):
        with open(file_output_pickle, 'rb') as f:
            k_current, dt_diff, dt_iter, list_checkpoints_duration, list_checkpoints_iters, list_checkpoints = pickle.load(f)

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
    
    # algorithm 1
    # delta is constant: +/- 1
    if 'a1' in algo_to_check:
        (p1, p2, duration, iterations) = gp.search_for_sym_primes (num, lambda iteration: gp.delta_constant_minus_1(iteration))
        update_algo_results (0, duration, iterations)

        if (p1 + p2)/2 != num:
            print ("Alg #1: violation of sum for p1=", p1, "p2=", p2, "n=", num)

    # algorithm 2
    # if n is even, delta is: +/- 1, otherwise: +/- 2
    if 'a2' in algo_to_check:
        if num % 2 == 0:
            (p1, p2, duration, iterations) = gp.search_for_sym_primes (num, lambda iteration: gp.delta_constant_minus_1(iteration))
        else:
            (p1, p2, duration, iterations) = gp.search_for_sym_primes (num, lambda iteration: gp.delta_constant_minus_2(iteration))
        update_algo_results (1, duration, iterations)

        if (p1 + p2)/2 != num:
            print ("Alg #2: violation of sum for p1=", p1, "p2=", p2, "n=", num)

##    # algorithm 3
##    # next candidate based on form of num 6a+i
##    if 'a3' in algo_to_check:
##        if num == 2:
##            p1 = 2
##            p2 = 2
##            iterations = 0
##            duration = 0
##        elif num == 3:
##            p1 = 3
##            p2 = 3
##            iterations = 0
##            duration = 0
##        elif num == 4:
##            p1 = 3
##            p2 = 5
##            iterations = 0
##            duration = 0
##        elif num == 5:
##            p1 = 5
##            p2 = 5
##            iterations = 0
##            duration = 0
##        elif num % 6 == 0 or num % 6 == 3:
##            (p1, p2, duration, iterations) = gp.search_for_sym_primes (num, lambda iteration: gp.delta_constant(iteration))
##        elif num % 6 == 1 or num % 6 == 4:
##            (p1, p2, duration, iterations) = gp.search_for_sym_primes (num, lambda iteration: gp.delta_constant(iteration))
##        else:
##            (p1, p2, duration, iterations) = gp.search_for_sym_primes (num, lambda iteration: gp.delta_constant(iteration))
##        update_algo_results (2, duration, iterations)
##        
##        if (p1 + p2)/2 != num:
##            print ("Alg #3: violation of sum for p1=", p1, "p2=", p2, "n=", num)
    
    # checkpoint - partial results
    if num % checkpoint_value == 0:
        dt_current = datetime.now()
        dt_diff_current = (dt_current - dt_current_previous).total_seconds()
        list_checkpoints.append(num)
        if 'a1' in algo_to_check:
            list_checkpoints_duration[0].append(dt_diff[0])
            list_checkpoints_iters[0].append(dt_iter[0])
        if 'a2' in algo_to_check:
            list_checkpoints_duration[1].append(dt_diff[1])
            list_checkpoints_iters[1].append(dt_iter[1])
        if 'a3' in algo_to_check:
            list_checkpoints_duration[2].append(dt_diff[1])
            list_checkpoints_iters[2].append(dt_iter[1])

        perc_completed = str(int(k * 100 / max_num))
        print ("Checkpoint", k, "of total", max_num, "took", dt_diff_current, "seconds. (" + perc_completed + "% completed)")
        
        # remember results so far
        write_results_to_figures (directory)
        k_current = k
        save_current_results(file_output_pickle)

dt_end = datetime.now()

# final results - time of processing
dt_diff = dt_end - dt_start
print ("Total calculations lasted:", dt_diff)

# final results - figures
write_results_to_figures (directory)
k_current = max_num
save_current_results(file_output_pickle)
