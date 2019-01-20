#
# Copyright (c) 2019, Marcin Barylski
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

min_num_def = 2
min_num = min_num_def
max_num_def = 10000
max_num = max_num_def
step_factor = 2
checkpoint_value = 100
file_input_primes = '..\\primes\\t_prime_numbers.txt'
file_input_nonprimes = '..\\primes\\t_nonprime_numbers.txt'

#############################################################
# Settings - output directory and files
#############################################################

directory = "results/" + str(step_factor*max_num)
if not os.path.exists(directory):
    os.makedirs(directory)

#############################################################
# Business logic
#############################################################

list_primes = []
list_primes_nom = []
list_primes_den = []
list_primes_factor_perc = []

def update_data (p, nominator, denominator):
    global list_primes, list_primes_nom, list_primes_den
    index = 0
    found = False
    for q in list_primes:
        if q == p:
            list_primes_nom[index] = nominator
            list_primes_den[index] = denominator
            found = True
            break
        index += 1
    if not found:
        list_primes.append(p)
        list_primes_nom.append(nominator)
        list_primes_den.append(denominator)

def set_prime_as_confimed (p):
    update_data (p, 1, 1)

def get_data (p):
    index = 0
    for q in list_primes:
        if q == p:
            return (list_primes_nom[index], list_primes_den[index])
            break
        index += 1
    return (0, 0)

def is_prime_confirmed (p):
    index = 0
    for q in list_primes:
        if q == p:
            if list_primes_nom[index] == 1 and list_primes_den[index] == 1:
                return True
            else:
                return False
        index += 1
    # not found - also return false
    return False

def clear_not_confirmed_primes ():
    index = 0
    for p in list_primes:
        if not is_prime_confirmed (p):
            list_primes_nom[index] = 0
            list_primes_den[index] = 0
        index += 1

def update_metrics ():
    index = 0
    list_primes_factor_perc.clear()
    for p in list_primes:
        if list_primes_den[index] > 0:
            list_primes_factor_perc.append (list_primes_nom[index]/list_primes_den[index])
        else:
            list_primes_factor_perc.append (0)
        index += 1

def print_candidates (threshold):
    index = 0
    for p in list_primes:
        if list_primes_factor_perc[index] > threshold:
            print ("Prime", p, "has value:", list_primes_factor_perc[index])
        index += 1

# returns True if new msut have primes found - this means that calculations must be restarted
def calculate_metrics (num, factors, dp, p):

    found = False
    no_of_partitions = dp.get_number_of_pairs (factors)
    
    if no_of_partitions == 1:
        for (p1, p2) in factors:
            if is_prime_confirmed (p1) and is_prime_confirmed (p2):
                found = True
                if be_verbose:
                    print ("Number", num, " Primes:", p1, "and", p2, "were already present as required.")
                return False
            else:
                set_prime_as_confimed (p1)
                set_prime_as_confimed (p2)
                found = True
                if be_verbose:
                    print ("Number", num, " Primes:", p1, "and", p2, "were added as required.")
                return True
    else:
        for (p1, p2) in factors:
            if is_prime_confirmed (p1) and is_prime_confirmed (p2):
                found = True
                if be_verbose:
                    print ("Number", num, " Primes:", p1, "and", p2, "are already confirmed as required.")
                break
    if not found:
        for (p1, p2) in factors:
            (n1, d1) = get_data (p1)
            (n2, d2) = get_data (p2)
            if not is_prime_confirmed (p1):
                update_data (p1, n1 + 1, d1 + no_of_partitions)
                if be_verbose:
                    print ("Number", num, " Prime:", p1, "will have nom=", n1 + 1, "den=", d1 + no_of_partitions)
            if not is_prime_confirmed (p2):
                update_data (p2, n2 + 1, d2 + no_of_partitions)
                if be_verbose:
                    print ("Number", num, " Prime:", p2, "will have nom=", n2 + 1, "den=", d2 + no_of_partitions)

    if be_verbose:   
        print (list_primes)
        print (list_primes_nom)
        print (list_primes_den)
        print ("END")

    return False

#############################################################
# Presentation
#############################################################

def write_results_to_figures (directory):

    plt.figure(1)
    plt.plot(list_primes, list_primes_factor_perc, 'b.', ms=2)
    plt.xlabel('Prime')
    plt.ylabel('Factor [%]')
    plt.title('Factor for primes - is it required for GSC?')
    plt.grid(True)
    plt.savefig(directory + "/f_required_primes.png")

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

completed_all = False
need_to_restart = True

while not completed_all:

    print ("---------------------------")
    print ("New calcuations started ...")
    print ("Current must-have primes:", list_primes)
    if need_to_restart:
        min_num = min_num_def
        clear_not_confirmed_primes ()
        need_to_restart = False

    for k in range (min_num, max_num_def):
        num = step_factor*k
    
        factors = gp.find_sum_of_prime_numbers (num)

        if calculate_metrics (num, factors, dp, p):
            print ("Need to restart because new must-have prime found ...")
            need_to_restart = True
            break

    # checkpoint - partial results
    if num % checkpoint_value == 0:
        dt_current = datetime.now()
        dt_diff_current = (dt_current - dt_current_previous).total_seconds()
        list_checkpoints.append(num)
        list_checkpoints_duration.append(dt_diff_current)

        update_metrics ()
        write_results_to_figures (directory)
        
        print ("Iteration", k, "of total", max_num, "took", dt_diff_current, "seconds")

    if num > max_num_def - 1:
        completed_all = True

dt_end = datetime.now()

update_metrics ()
write_results_to_figures (directory)
print (list_primes)
print (list_primes_nom)
print (list_primes_den)

print_candidates (0.3)

# final results
dt_diff = dt_end - dt_start
print ("Total calculations lasted:", dt_diff)
