def readFasta( inFileName ):
    "Read fasta files"
    from Bio import SeqIO
    import sys

    headerList = []
    seqList = []
    
    fasta_sequences = SeqIO.parse(open(inFileName),'fasta')
    for each_sequence in fasta_sequences:
        headerList.append(each_sequence.id)
        seqList.append(str(each_sequence.seq))
    return headerList, seqList

def findMutations( seqList ):
    "Find unique mutation codons and counts"
    from operator import itemgetter

    # Turn seqList into numbers to save memory
    numList = []
    for seq in seqList:
        num = [ord(x) for x in seq]
        numList.append(num)
    
    germline = numList[0]
    mutations = []

    for i in range(1,len(numList)):
        for j in range(len(numList[i])):
            if germline[j]!=numList[i][j] and germline[j]!=ord('N') and numList[i][j]!=ord('N'):
                if j%3 == 0:
                    nt_idx = [j,j+1,j+2]
                elif j%3 == 1:
                    nt_idx = [j-1,j,j+1]
                elif j%3 == 2:
                    nt_idx = [j-2,j-1,j]
                mutations.append([''.join(itemgetter(*nt_idx)(seqList[0])), str(j+1), ''.join(itemgetter(*nt_idx)(seqList[i]))]) # Create codons based on position
    
    # Claculate unique mutations and counts
    unique_mutations = []
    unique_count = []
    for mutation in mutations:
        if mutation not in unique_mutations:
            unique_mutations.append(mutation)
            unique_count.append(1)
        elif mutation in unique_mutations:
            unique_count[unique_mutations.index(mutation)] += 1
    return unique_mutations, unique_count
            
def divideMutation( unique_mutations,unique_count ):
    "devide mutations into single mutation (one mutation in a codon) and multiple mutation (more than one mutations in a codon)"
    from itertools import izip

    single_mutations = []
    single_count = []
    multiple_mutations = []
    multiple_count = []
    
    for mutation in unique_mutations:
        if sum(c1 != c2 for c1,c2 in izip(mutation[0],mutation[2])) == 1:
            single_mutations.append(mutation)
            single_count = unique_count[unique_mutations.index(mutation)]
        else:
            multiple_mutations.append(mutation)
            multiple_count = unique_count[unique_mutations.index(mutation)]

    return single_mutations, single_count, multiple_mutations, multiple_count 
       

def groupMultipleMutation( multiple_mutations ):
    "Group mutations from the same codon"
    from math import ceil

    multiple_group = [0]*len(multiple_mutations)

    for i in range(len(multiple_mutations)):
        if multiple_group[i]==0:
            multiple_group[i] = i+1
            for j in range(i+1,len(multiple_mutations)):
                if multiple_group[j]==0:
                    if multiple_mutations[i][0]==multiple_mutations[j][0] and multiple_mutations[i][2]==multiple_mutations[j][2] and ceil(float(multiple_mutations[i][1])/3)==ceil(float(multiple_mutations[j][1])/3):
                        multiple_group[j] = i+1

    return multiple_group


def mutationType(single_mutations):
    "Find mutations type (R/S) for single mutation"
    from Bio import Seq
    
    for i in range(len(single_mutations)):
        germline = single_mutations[i][0]
        mutated = single_mutations[i][2]
        if '-' not in germline and '-' not in mutated:
            if Seq.translate(germline) == Seq.translate(mutated):
                single_mutations[i].append('silent')
            else:
                single_mutations[i].append('replacement')
        else:
            sinlge_mutations[i].append('unknown')

        return sinlge_mutations


def applyBias( multiple_mutations,multiple_group,bias ):
    ""
    from itertools import permutations
    from Bio import Seq
    from collections import Counter
    
    counted = [0]*len(multiple_mutations)
    for mutation in multiple_mutations:
        if counted[multiple_mutations.index(mutation)] != 1:
            mismatch_positions = [i for i in range(len(mutation[0])) if mutation[0][i]!=mutation[2][i]]
            p = list(permutations(mismatch_positions))
            type_list = []
            type_count = []
            for i in range(len(p)):
                types = []
                germline = mutation[0]
                for j in range(len(p[i])):
                    mutated = germline[:p[i][j]] + mutation[2][p[i][j]] + germline[p[i][j]+1:]
                    if 'N' not in germline and 'N' not in mutated:
                        if Seq.translate(germline) == Seq.translate(mutated):
                            types.append('silent')
                        else:
                            types.append('replacement')
                    else:
                        types.append('unkown')   
                    germline = mutated
                type_list.append(types)
                type_frequency = Counter(types)
                type_count.append(type_frequency[bias])
            type_list = type_list[type_count.index(max(type_count))]
                         
            indices = [i for i, x in enumerate(multiple_group) if x == multiple_group[multiple_mutations.index(mutation)]]
            for idx in indices:
                counted[idx] = 1
                multiple_mutations[idx].append(type_list[indices.index(idx)])
    return multiple_mutations          
            
