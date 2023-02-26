"""
BSD 3-Clause License

Copyright (c) 2023, Özlem Salman
Copyright (c) 2023, Armağan Salman

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its
   contributors may be used to endorse or promote products derived from
   this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""


#(
# Python system modules:
import itertools
#)

#(
# 3rd party modules:
import Bio
from Bio import Align
#)

#(
# Local modules:
import csv_file_io as Csvio
#)


def get_grouped_data(data_file_path, group_key_func = None):
#(
	if group_key_func == None:
	#(
		# By default, use activation class as group key
		group_key_func = lambda x: x[2]
	#)
	data_with_header = Csvio.read_data_from_csv(data_file_path)
	data = data_with_header[1:] # Skip data header
	
	# activation class string might have spaces. strip them:
	data_iter = map(lambda r: [r[0], r[1], r[2].strip()] , data)
	
	# Sorting is necessary for itertools.groupby
	sorted_data = sorted(data_iter, key = lambda x: x[2]) # x[2] = activation class
	
	return itertools.groupby(sorted_data, key = group_key_func)
#)


def select_rows_by_min_seqlen(data_rows, desired_min_len, desired_max_len = None):
#(
    """ row = id, sequence, activation class string
    """
    
    selected_rows = []
    for rw in data_rows:
    #(
        seq = rw[1]
        if len(seq) >= desired_min_len and (desired_max_len == None or len(seq) <= desired_max_len):
            selected_rows.append(rw)
        
    #)      
    return selected_rows
#)


def cutoff_seq_suffixes(data_rows, cutoff_len):
#(
    """ row = id, sequence, activation class string
    """
    
    selected_rows = []
    for rw in data_rows:
    #(
        seq = rw[1]
        
        if len(seq) >= cutoff_len:
        #(
            new_row = [rw[0], seq[:cutoff_len], rw[2]]
            selected_rows.append(new_row)
        #)
    #)      
    return selected_rows
#)


def PSS_same_len(first_seq, second_seq):
#(
    #pairwise_alignment_score_same_len = PSS
    if len(first_seq) == len(second_seq):
    #(
        aligner = Align.PairwiseAligner()
        aligner.open_gap_score = -0.5
        aligner.extend_gap_score = -0.1
        aligner.target_end_gap_score = 0.0
        aligner.query_end_gap_score = 0.0
        
        # return aligner.align(first_seq, second_seq)
        return aligner.score(first_seq, second_seq)
    #)
    else:
    #(
        raise Exception("Length of the sequences must be the same.")
    #)
#)


def pairwise_align_two_datasets(csv_file):
#(
    #mod. active ile inactive exp gruplarından 80'er tane row alınıp, her row karşısındaki row ile karşılaştırılacak
    get_class_value = lambda x: x[2]
    grouped_class = get_grouped_data(csv_file, group_key_func = get_class_value)

    A_set = []
    B_set = []
    
    for x in grouped_class:
    #(
        if x[0] == "inactive - exp":
            A_set = list(x[1])
            
        if x[0] == "mod. active":
            B_set = list(x[1])
    #)
    
    min_seqlen = 12
    
    A_set_12_len = select_rows_by_min_seqlen(A_set, min_seqlen)
    B_set_12_len = select_rows_by_min_seqlen(B_set, min_seqlen) 
    
    score = pairwise_alignment_score_same_len("abc", "abc")
    print(f"Alignment score: {score}")
    
    for x in A_set_12_len:
        print(x)
    
    print("~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~")
    print("A_set_cut_12:")
    
    A_set_cut_12 = cutoff_seq_suffixes(A_set_12_len, min_seqlen)
    
    for x in A_set_cut_12:
        print(x)
    #
    
#)
def delta_func(i, j):
       
    if i == j:
        
        return 1
    else:
        
        return 0  

        

def total_similarity_scores(A_set, B_set):
#(
    A_set_len = len(A_set)
    B_set_len = len(B_set)
        
            
    #PSS_ij = PSS_same_len(A_set[i], B_set[j])
    double_sum = 0
    for i in range(A_set_len):
    #(
        for j in range(B_set_len):
        #(
            A_seq = A_set[i]
            B_seq = B_set[j]
           
            PSS_ij = PSS_same_len(A_seq, B_seq)
            delta_ij = delta_func(i, j)
            delta_AB = delta_func(A_set, B_set)

            inner_equation = 1 - delta_ij * delta_AB
            product = PSS_ij * inner_equation
            
            double_sum += product  #hesaplananları burda tutuyor
        #)   
    #) 
        
    # product içteki sum kısmını hesaplıyor
    first_equation = 1 / A_set_len * (B_set_len - delta_AB)
    TSS = first_equation * double_sum
    return TSS
        
        
        
    

#)


def main(args):
#(
    #mod. active ile inactive exp gruplarından 80'er tane row alınıp, her row karşısındaki row ile karşılaştırılacak

    pairwise_align_two_datasets("ACPs_Breast_cancer.csv")
    """
    first_file = args["first_sequence_file"]
    second_file = args["second_sequence_file"]

    get_class_value = lambda x: x[2]
    grouped_class = get_grouped_data(first_file, group_key_func = get_class_value)

    for key, group in grouped_class:
    #(
        print("~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~")
        print(key)
        for x in group:
        #(
            print(x)
        #)
        print("~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~")
    #)
    """


    print("Executed main.")
#)


if __name__ == "__main__":
#(
	args = dict()
	
	f = "ACPs_Breast_cancer.csv"
	sf = "ACPs_Lung_cancer.csv"
	assert(f != sf)
	
	args["first_sequence_file"] = f
	args["second_sequence_file"] = sf
	
	main(args)
#)



"""
select_by_seqlen("data.csv", 12)
select_by_seqlen("data.csv", 15)
select_by_seqlen("data.csv", 20)
select_by_seqlen("data.csv", 3)
"""

"""
#peptitleri uzunluğuna göre alınacak, aşağıdaki de dosya okumak için
    data_with_header = Csvio.read_data_from_csv(csv_file_path)
	rows = data_with_header[1:]
"""