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
            

        
