import sys
from numpy import *

def fixStableGene(list_gene, probab, num_selected):
	sort_idx = sorted(range(len(probab)), key=lambda x:probab[x], reverse=True)
	fixed_gene = []
	for i in range(num_selected):
		fixed_gene.append(list_gene[sort_idx[i]])

	return fixed_gene

def crossover(gene1, gene2, num_cut):
	len_gene1 = len(gene1)
	len_gene2 = len(gene2)

	if len_gene1 != len_gene2:
		print 'length of two genes is different'
	else:
		child_gene = []
		prev_num = 1
		for i in range(num_cut):
			num_ran = random.random_integers(prev_num + 1, len_gene1 - (num_cut - i - 1))
#			print num_ran
			if i%2 == 0:
				child_gene = child_gene + gene1[prev_num-1:num_ran-1]
			else:
				child_gene = child_gene + gene2[prev_num-1:num_ran-1]
			prev_num = num_ran

		if i%2 == 1:
			child_gene = child_gene + gene1[prev_num-1:]
                else:
                        child_gene = child_gene + gene2[prev_num-1:]

		return child_gene

def mutation(gene, wid_range):
	spin = [-1, 1]

	len_gene = len(gene)
	num_ran = random.random_integers(0, len_gene-1); #print num_ran

	child_gene = []
	if num_ran + wid_range > len_gene-1:
		rand_gene = ndarray.tolist(random.choice(spin, wid_range))
		child_gene = rand_gene[:num_ran+wid_range-len_gene] + gene[num_ran+wid_range-len_gene:num_ran] + rand_gene[num_ran+wid_range-len_gene:]
	else:
		rand_gene = ndarray.tolist(random.choice(spin, wid_range))
		child_gene = gene[0:num_ran] + rand_gene + gene[num_ran+wid_range:]

	child_gene = array(child_gene)
	if child_gene[0] == -1:
		child_gene = -1 * child_gene

	return ndarray.tolist(child_gene)

def randomGen(num_at):
	spin = [-1,1]
	rand_gene = random.choice(spin, num_at)
	if rand_gene[0] == -1:
		rand_gene = -1 * rand_gene

	return ndarray.tolist(rand_gene)

#print crossover([1,2,3,4,5,6,7,8,9], [11,12,13,14,15,16,17,18,19], 3)
#print randomGen(9)
#print mutation([1,2,3,4,5,6,7,8,9], 4)
#a = ['a', 'b', 'c', 'd','e','f']
#pro = [0.3, 0.2, 0.5, 0.01, 0.4, 0.5]
#print fixStableGene(a, pro, 4)
