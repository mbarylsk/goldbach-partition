#
# Copyright (c) 2016 - 2019, Marcin Barylski
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

class DataProcessing:

    def get_diff_in_factors (self, factors):
        diffs = []
        for pair in factors:
            (p1, p2) = pair
            diff = p2 - p1
            diffs.append(diff)
        return diffs

    def get_max_factor (self, factors):
        maxv = 0
        for pair in factors:
            (p1, p2) = pair
            if (p1 > maxv):
                maxv = p1
            if (p2 > maxv):
                maxv = p2
        return maxv

    def get_avg_factor (self, factors):
        sumv = 0
        count = 0
        for pair in factors:
            (p1, p2) = pair
            count += 2
            sumv += p1 + p2
        avgv = sumv / count
        return avgv

    def get_min_factor (self, factors):
        minv = sys.maxsize
        for pair in factors:
            (p1, p2) = pair
            if (p1 < minv):
                minv = p1
            if (p2 < minv):
                minv = p2
        return minv

    def get_unique_factors (self, factors):
        uf = []
        for pair in factors:
            (p1, p2) = pair
            if p1 not in uf:
                uf.append (p1)
            if p1 not in uf:
                uf.append (p2)
        return uf

    def get_max_ratio_in_factors (self, factors):
        maxr = 0
        for pair in factors:
            (p1, p2) = pair
            ratio = p1 / p2
            if ratio <=1 and ratio > maxr:
                maxr = ratio
        return maxr

    def get_number_of_pairs (self, factors):
        return len(factors)

    def get_perc_max_saturation_from_factors (self, n, factors):
        max_p = self.get_max_factor (factors)
        return max_p / n * 100

    def get_perc_min_saturation_from_factors (self, n, factors):
        if n == 0:
            return 0
        else:
            min_p = self.get_min_factor (factors)
            return min_p / n * 100

    def get_perc_avg_saturation_from_factors (self, n, factors):
        if n == 0:
            return 0
        else:
            avg_p = self.get_avg_factor (factors)
            return avg_p / n * 100

    def get_max_value_from_list (self, l):
        return max(l)

    def get_min_value_from_list (self, l):
        return min(l)

    def get_avg_value_from_list (self, l):
        if (len(l) > 0):
            return (sum(l)/len(l))
        else:
            return 0

    def dictionary_cleanup (self, d):
        max_x = 0
        for x, y in d.items():
            if x > max_x:
                max_x = x
        new_l = []
        for z in range (0, max_x+1, 2):
            if z in d:
                new_l.append (d[z])
            else:
                new_l.append (0)
        return new_l

    def divide_list_into_chunks (self, r, size):
        out = []
        last = 0
        while last < len(r):
            out.append (r[int(last):int(last + size)])
            last += size
        return out

    def read_sums_from_line (self, line): 
        line = line.replace('Pairs: ', '')
        line = line.replace('[', '')
        line = line.replace(']', '')
        line = line.replace('(', '')
        line = line.replace(')', '')
        line = line.replace('\n', '')
        line = line.replace(' ', '')
        numbers = line.split(',')
        read_second_number = False
        factors_from_file = []
        for number in numbers:
            if read_second_number:
                p2 = int(number)
                pair = (p1, p2)
                factors_from_file.append (pair)
                read_second_number = False
            else:
                p1 = int(number)
                read_second_number = True
        return factors_from_file

    def read_num_from_line (self, line):
        line = line.replace('Number: ', '')
        num_from_file = int(line)
        return num_from_file
