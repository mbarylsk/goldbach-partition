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

import matplotlib.pyplot as plt
import primes
import os
import numpy as np
import pickle

#############################################################
# Settings - configuration
#############################################################

# Minimal even number checked against Goldbach conjecture
#   o number = min_num * step_factor
min_num = 1
step_factor = 2
# Maximum even number checked against Goldbach conjecture
#   o number = max_num * step_factor
max_num = 2000000

# Checkpoint value when partial results are drawn/displayed
# should be greater than zero
checkpoint_value = 100000

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

color_no_turn = 0
color_turn = 50

#############################################################
# Settings - output directory and files
#############################################################

directory = str(step_factor*max_num)
if not os.path.exists(directory):
    os.makedirs(directory)
file_output_shape_1 = directory + "/f_shape_1.png"
file_output_shape_2 = directory + "/f_shape_2.png"
file_output_shape_3 = directory + "/f_shape_3.png"
file_output_shape_4 = directory + "/f_shape_4.png"
file_output_pickle = directory + "/objs_shape.pickle"

#############################################################
# Results of calculations
#############################################################

new_x = [0, 0, 0, 0]
new_y = [0, 0, 0, 0]
delta_x = [1, 1, 1, 1]
delta_y = [0, 0, 0, 0]
num_current = [0, 0, 0, 0]

datax = [[0], [0], [0], [0]]
datay = [[0], [0], [0], [0]]
colors = [[0], [0], [0], [0]]

is_previous_prime = [False, False, False, False]
sign = 1
k_current = 0

#############################################################
# Presentation
#############################################################

def write_results_to_figures():
    area = np.pi * 2
    fig = plt.figure(1)
    plt.scatter(datax[0], datay[0], s=area, c=colors[0], alpha=0.5)
    title = "n=2k+1 (k=1,2,3...); n from 3 to " + str(num_current[0])
    fig.suptitle(title, fontsize=10)
    plt.savefig(file_output_shape_1)

    fig = plt.figure(2)
    plt.scatter(datax[1], datay[1], s=area, c=colors[1], alpha=0.5)
    title = "n=6k+1 (k=1,2,3...); n from 7 to " + str(num_current[1])
    fig.suptitle(title, fontsize=10)
    plt.savefig(file_output_shape_2)

    fig = plt.figure(3)
    plt.scatter(datax[2], datay[2], s=area, c=colors[2], alpha=0.5)
    title = "n=6k-1 (k=1,2,3...); n from 5 to " + str(num_current[2])
    fig.suptitle(title, fontsize=10)
    plt.savefig(file_output_shape_3)

    fig = plt.figure(4)
    plt.scatter(datax[3], datay[3], s=area, c=colors[3], alpha=0.5)
    title = "n = 6k-+1 (k=1,2,3...); n from 5 to " + str(num_current[3])
    fig.suptitle(title, fontsize=10)
    plt.savefig(file_output_shape_4)

def next_delta (delta_x, delta_y, turn):
    if turn:
        if delta_x == 1 and delta_y == 0:
            delta_x = 0
            delta_y = 1
            return (delta_x, delta_y)
        if delta_x == 0 and delta_y == 1:
            delta_x = -1
            delta_y = 0
            return (delta_x, delta_y)
        if delta_x == -1 and delta_y == 0:
            delta_x = 0
            delta_y = -1
            return (delta_x, delta_y)
        if delta_x == 0 and delta_y == -1:
            delta_x = 1
            delta_y = 0
            return (delta_x, delta_y)
    return (delta_x, delta_y)

def next_turn (p, num, is_previous_prime):
    turn = False
    if p.is_prime(num) and not is_previous_prime:
        turn = True
        is_previous_prime = True
    if not p.is_prime(num) and is_previous_prime:
        turn = True
        is_previous_prime = False
    return (turn, is_previous_prime)

def next_sign (sign):
    if sign == 1:
        return -1
    else:
        return 1

def save_current_results (file_output_pickle):
    global k_current, new_x, new_y, delta_x, delta_y, num_current, datax, datay, colors, is_previous_prime, sign
    with open(file_output_pickle, 'wb') as f:
        pickle.dump([k_current, new_x, new_y, delta_x, delta_y, num_current, datax, datay, colors, is_previous_prime, sign], f)

def restore_previous_results (file_output_pickle):
    global k_current, new_x, new_y, delta_x, delta_y, num_current, datax, datay, colors, is_previous_prime, sign
    if os.path.exists(file_output_pickle):
        with open(file_output_pickle, 'rb') as f:
           k_current, new_x, new_y, delta_x, delta_y, num_current, datax, datay, colors, is_previous_prime, sign = pickle.load(f)

#############################################################
# Main
#############################################################

print ("Initialize objects...")
p = primes.Primes(caching_primality_results)
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

# new calculations
for k in range (min_num, max_num):

    # case 1: subsequent odd numbers
    num = k*step_factor + 1
    num_current[0] = num

    (turn, is_previous_prime[0]) = next_turn (p, num, is_previous_prime[0])
    (delta_x[0], delta_y[0]) = next_delta (delta_x[0], delta_y[0], turn)

    new_x[0]+= delta_x[0]
    new_y[0]+= delta_y[0]
    datax[0].append(new_x[0])
    datay[0].append(new_y[0])
    if turn:
        colors[0].append(color_turn)
    else:
        colors[0].append(color_no_turn)

    # case 2:
    num = k*3*step_factor + 1
    num_current[1] = num

    (turn, is_previous_prime[1]) = next_turn (p, num, is_previous_prime[1])
    (delta_x[1], delta_y[1]) = next_delta (delta_x[1], delta_y[1], turn)

    new_x[1]+= delta_x[1]
    new_y[1]+= delta_y[1]
    datax[1].append(new_x[1])
    datay[1].append(new_y[1])
    if turn:
        colors[1].append(color_turn)
    else:
        colors[1].append(color_no_turn)

    # case 3:
    num = k*3*step_factor - 1
    num_current[2] = num

    (turn, is_previous_prime[2]) = next_turn (p, num, is_previous_prime[2])
    (delta_x[2], delta_y[2]) = next_delta (delta_x[2], delta_y[2], turn)

    new_x[2]+= delta_x[2]
    new_y[2]+= delta_y[2]
    datax[2].append(new_x[2])
    datay[2].append(new_y[2])
    if turn:
        colors[2].append(color_turn)
    else:
        colors[2].append(color_no_turn)

    # case 4:
    num = k*3*step_factor - 1*sign
    num_current[3] = num

    (turn, is_previous_prime[3]) = next_turn (p, num, is_previous_prime[3])
    sign = next_sign (sign)
    (delta_x[3], delta_y[3]) = next_delta (delta_x[3], delta_y[3], turn)

    new_x[3]+= delta_x[3]
    new_y[3]+= delta_y[3]
    datax[3].append(new_x[3])
    datay[3].append(new_y[3])
    if turn:
        colors[3].append(color_turn)
    else:
        colors[3].append(color_no_turn)

    if num % checkpoint_value == 1:

        perc_completed = str(int(k * 100 / max_num))
        print ("Checkpoint", k, "of total", max_num, "(" + perc_completed + "% completed)")
        
        # save results collected so far
        write_results_to_figures ()
        k_current = k
        save_current_results(file_output_pickle)


# final results
write_results_to_figures ()
k_current = max_num
save_current_results(file_output_pickle)

