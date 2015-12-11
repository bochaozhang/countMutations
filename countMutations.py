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
                mutations.append([''.join(itemgetter(*nt_idx)(seqList[0])), j+1, ''.join(itemgetter(*nt_idx)(seqList[i]))]) # Create codons based on position
    
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
    "devide mutations into single mutation (one mutation in a codon), double mutation (two mutations in a codon) and triple mutation (three mutations in a codon)"
    from itertools import izip

    single_mutations = []
    single_count = []
    double_mutations = []
    double_count = []
    triple_mutations = []
    triple_count = []
    
    for mutation in unique_mutations:
        if sum(c1 != c2 for c1,c2 in izip(mutation[0],mutation[2])) == 1:
            single_mutations.append(mutation)
            single_count = unique_count[unique_mutations.index(mutation)]
        elif sum(c1 != c2 for c1,c2 in izip(mutation[0],mutation[2])) == 2:
            double_mutations.append(mutation)
            double_count = unique_count[unique_mutations.index(mutation)]
        if sum(c1 != c2 for c1,c2 in izip(mutation[0],mutation[2])) == 3:
            triple_mutations.append(mutation)
            triple_count = unique_count[unique_mutations.index(mutation)

    return single_mutations, single_count, double_mutations, double_count, triple_mutations, triple_count 


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


def filterMutations( previous_mutations,previous_count,multiple_mutations,multiple_count ):
    "Try to explain multiple mutation by previous mutation"
    from itertools import permutations

    getting_removed = []
    for i in range(len(multiple_mutations)):
        germline = multiple_mutations[i][0]
        mutated = multiple_mutations[i][2]
        positions = [j for j in range(len(germline)) if germline[j]!=mutatied[j]]          
        for position in positions:
            mutation = []
            
        indices = [j for j, x in enumerate(multiple_linkage) if x == 1]
