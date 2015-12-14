import countMutations
headerList,seqList = countMutations.readFasta(‘test.fasta’)
unique_mutations,unique_count = countMutations.findMutations(seqList)
single_mutations,single_count,multiple_mutations,multiple_count = countMutations.divideMutation(unique_mutations,unique_count)
multiple_group = countMutations.groupMultipleMutation(multiple_mutations)
multiple_mutations = countMutations.applyBias(multiple_mutations,multiple_count,multiple_group,bias)