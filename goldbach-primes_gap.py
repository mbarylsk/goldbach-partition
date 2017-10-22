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

import math
import sys
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.mlab as mlab
import os
from datetime import datetime
import time
import numpy
sys.path.insert(0, '..\\primes\\')
import primes
import goldbach
import dataprocessing

#############################################################
# Settings - configuration
#############################################################

min_num = 2
max_num = 1000000

# Checkpoint value when partial results are drawn/displayed
checkpoint_value = 5000

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
# Settings - output directory and files
#############################################################

directory = "results/" + str(max_num)
if not os.path.exists(directory):
    os.makedirs(directory)

file_output_gap_num_of_sym_primes = directory + "/f_checkpoint_gap_num_of_sym_primes.png"
file_output_gap_ratio_to_goldbach = directory + "/f_checkpoint_gap_ratio_to_goldbach.png"
file_output_gap_ratio_avg_to_goldbach = directory + "/f_checkpoint_gap_ratio_avg_to_goldbach.png"
file_output_max_diff_in_sym_primes = directory + "/f_checkpoint_gap_max_diff_in_sym_primes.png"
file_output_min_index_sym_primes = directory + "/f_checkpoint_gap_min_index_sym_primes.png"
file_output_max_index_sym_primes = directory + "/f_checkpoint_gap_max_index_sym_primes.png"

############################################################
# Results of calculations
#############################################################

dt_diff = [0]
dt_iter = [0]

list_checkpoints = []
list_checkpoints_primes = []
list_checkpoints_ratio = []
list_checkpoints_ratio_avg = []
list_checkpoints_max_diff = []
list_checkpoints_min_index = []
list_checkpoints_max_index = []

#############################################################
# Presentation
#############################################################

def write_results_to_figures (directory):
    # results - figures
    fig1 = plt.figure(1)
    plt.plot(list_checkpoints, list_checkpoints_primes, 'g.', ms=2)
    plt.xlabel('n')
    plt.ylabel('Count')
    plt.title('Number of symmetric primes')
    plt.grid(True)
    plt.savefig(file_output_gap_num_of_sym_primes)
    plt.close(fig1)

    fig2 = plt.figure(2)
    plt.plot(list_checkpoints, list_checkpoints_ratio, 'g.', ms=2)
    plt.xlabel('n')
    plt.ylabel('Ratio')
    plt.title('Ratio: # of GPs / # of symmetric primes')
    plt.grid(True)
    plt.savefig(file_output_gap_ratio_to_goldbach)
    plt.close(fig2)

    fig3 = plt.figure(3)
    plt.plot(list_checkpoints, list_checkpoints_max_diff, 'g.', ms=2)
    plt.xlabel('n')
    plt.ylabel('Delta')
    plt.title('Maximum difference in index in symmetric primes')
    plt.grid(True)
    plt.savefig(file_output_max_diff_in_sym_primes)
    plt.close(fig3)

    fig4 = plt.figure(4)
    plt.plot(list_checkpoints, list_checkpoints_min_index, 'g.', ms=2)
    plt.xlabel('n')
    plt.ylabel('Index')
    plt.title('Minimal index in symmetric primes')
    plt.grid(True)
    plt.savefig(file_output_min_index_sym_primes)
    plt.close(fig4)

    fig5 = plt.figure(5)
    plt.plot(list_checkpoints, list_checkpoints_max_index, 'g.', ms=2)
    plt.xlabel('n')
    plt.ylabel('Index')
    plt.title('Maximum index in symmetric primes')
    plt.grid(True)
    plt.savefig(file_output_max_index_sym_primes)
    plt.close(fig5)

    fig6 = plt.figure(6)
    plt.plot(list_checkpoints, list_checkpoints_ratio_avg, 'g.', ms=2)
    plt.xlabel('n')
    plt.ylabel('Avg(Ratio)')
    plt.title('Average ratio: # of GPs / # of symmetric primes')
    plt.grid(True)
    plt.savefig(file_output_gap_ratio_avg_to_goldbach)
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
p.init_set(file_input_primes, True)
p.init_set(file_input_nonprimes, False)
print ("DONE")
print ("Sorting primes...")
p.sort_prime_set()
print ("DONE")

dt_start = datetime.now()
dt_current_previous = dt_start

# new calculations
for num in range (min_num, max_num):

    primes_found = []
    index_found = []
    for i in range (0, num - 1):
        (both_primes, k1, k2) = p.is_symmetric_prime (num, i)
        if both_primes:
            primes_found.append(k1)
            index_found.append(i)
            if k1 != k2:
                primes_found.append(k2)
            
    primes_found.sort()
    count = len(primes_found)
    max_index = dp.get_max_value_from_list (index_found)
    min_index = dp.get_min_value_from_list (index_found)
    
    list_checkpoints.append (num)
    list_checkpoints_primes.append (count)
    list_checkpoints_max_diff.append (max_index - min_index)
    list_checkpoints_min_index.append (min_index)
    list_checkpoints_max_index.append (max_index)

    factors = gp.find_sum_of_prime_numbers (2 * num)
    list_checkpoints_ratio.append (dp.get_number_of_pairs (factors) / count)
    avg_ratio = dp.get_avg_value_from_list (list_checkpoints_ratio)
    list_checkpoints_ratio_avg.append (avg_ratio)

    if count == 0:
        print ("WARNING: For n=", n, "no symmetric primes found")

    if num % checkpoint_value == 0:
        write_results_to_figures (directory)
        dt_current = datetime.now()
        dt_diff_current = (dt_current - dt_current_previous).total_seconds()
        perc_completed = str(int(num * 100 / max_num))
        print ("Checkpoint", num, "of total", max_num, "took", dt_diff_current, "seconds. (" + perc_completed + "% completed)")
        print (" Average ratio:", avg_ratio)

dt_end = datetime.now()

# final results - time of processing
dt_diff = dt_end - dt_start
print ("Total calculations lasted:", dt_diff)

write_results_to_figures (directory)
