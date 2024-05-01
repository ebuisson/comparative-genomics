import csv
import matplotlib.pyplot as plt
from collections import Counter
import math

'''
Function to identify conserved regions in a Multiple Sequence Alignment file with Shannon Entropy method
input = MSA file in .aln format
output = list with Shannon entropy values for each position in the MSA file
'''

def shannon_entropy(aa_counts):
    # gaps included as a 21st AA
    # positions with >50% gaps are ignored (output = -1)
    # pi = fraction of residues of each aa type
    # M = number of amino acid types
    H = 0
    residues_count = int(sum(aa_counts.values()))
    if residues_count <= 13:
        H = 1
    else:
        for value in aa_counts.values():
            pi = int(value) / 27
            H = H + (pi*math.log(pi,2))
        if residues_count < 27:
            gap_pi = (27-residues_count)/ 27
            H = H + (gap_pi*math.log(gap_pi,2))
    return -H

#included to compare to Shannon entropy scores
def variability(N,k,freq):
    # excluding gaps
    # N = the number of sequences in the alignment
    # k = the number of different amino acids at a given position
    # n = the frequency of the most common amino acid at that position
    return((N*k)/freq)

#make a dictionary with key=species and value=string of protein sequence with gaps
def make_dict(file):
    dict = {}
    with open(file) as my_file:
        tsv_reader = csv.reader(my_file, delimiter=" ", skipinitialspace=True)
        next(tsv_reader)
        next(tsv_reader)
        next(tsv_reader)
        for row in tsv_reader: 
            if len(row)==0 or len(row[0])==0 or row[0].startswith(tuple([' ','.',':','*'])):
                continue
            if dict.get(row[0]) is None:
                dict.update({row[0]:''})
            dict[row[0]]=dict[row[0]] + row[1]
    return dict

#to plot resulting entropy values; values between 0 and 1 (inclusive) will be plotted in red
def barplot(y_values):
    print(y_values)
    x_axis = range(len(y_values))
    y_axis = y_values
    col = []
    for num in y_axis:
        if num < 0 or num > 1:
            col.append('blue')
        else:
            col.append('red')
    plt.bar(x_axis, y_axis, color=col)
    plt.xlabel('MSA index')
    plt.ylabel('Shannon entropy')
    plt.show()

def functional_regions(dict):
    len_seq = int(sum(len(x) for x in dict.values())/len(dict))
    var = []
    shannon = []
    for i in range(len_seq):
        str = ''
        for value in dict.values(): 
            x = list(value)
            if x[i] == '-':
                continue
            else:
                str = str + x[i]
        s = set(str)
        c = Counter(str)
        v = variability(len(str), len(s), (c.most_common(1)[0][1])/len(str))
        var.append(v)
        se = shannon_entropy(c)
        shannon.append(se)
        str=''
    return(shannon)

    
    
    


        
