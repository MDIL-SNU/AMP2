# This is a package of modules for genetic algorithm.
import sys,os,yaml
from numpy import *
import genetic_operator as go
from module_vector import *
from copy import deepcopy as dcopy

# This function is for generating next generation in genetic algorithm.
def genNextGeneration(size_gene, list_gene, probab, portion, mut_range):
	newlist_gene = []

	if probab != [] and list_gene != []:
		if portion[0] != 0:
			newlist_gene = newlist_gene + go.fixStableGene(list_gene, probab, portion[0])

		if portion[1] != 0:
			i = 0
			while i < portion[1]:
				parent_gene = selectParents(probab, 2)
#				print parent_gene
				child_candidate = go.crossover(list_gene[parent_gene[0]], list_gene[parent_gene[1]], 1)
				if not child_candidate in newlist_gene:
					newlist_gene.append(child_candidate)
					i = i + 1

		if portion[2] != 0:
			i = 0
			while i < portion[2]:
				parent_gene = selectParents(probab, 1)
				child_candidate = go.mutation(list_gene[parent_gene[0]], mut_range)
				if not child_candidate in newlist_gene:
					newlist_gene.append(child_candidate)
					i = i + 1

	i = 0
	while i < portion[3]:
		child_candidate = go.randomGen(size_gene)
		if not child_candidate in newlist_gene:
			newlist_gene.append(child_candidate)
			i = i + 1

	return newlist_gene

# This function is for select parents making child configuration
def selectParents(probab, s_num):
	list_selected = []
	probab_copy = copy(probab)
	probab_sum = sum(probab)

	for i in range(s_num):
		num_ran = random.random_sample() * probab_sum
		tot_sum = 0
		for j,item in enumerate(probab_copy):
			tot_sum = tot_sum + item
			if num_ran < tot_sum:
				list_selected.append(j)
				probab_copy[j] = 0.0
				probab_sum = sum(probab_copy)
				break

	return list_selected

# This function is for converting from energy to probability
def energyToProbab(energy, gamma):
	proba = []
	E_max = max(energy); E_min = min(energy)
	for E in energy:
		proba.append((E_max-E) + (E_max-E_min)*gamma/(gamma-1))

	return proba

# This function is for calculating energy after switching spin
def getEnergy(list_gene, info_pair, para_J):
	list_E = []
	for gene in list_gene:
		E_sum = 0
		for pair_stat in info_pair:
			E_sum = E_sum + gene[pair_stat[2]] * gene[pair_stat[3]] * para_J[pair_stat[5]]
		list_E.append(E_sum/2.0)
	
	return list_E

# This function is for making pair list
def pairCategorize(atom_pos, axis, mag_true, cutoff):
	list_pair = []
	mag_atom_pos = []
	for i in range(len(atom_pos)):
		if mag_true[i] == 1:
			mag_atom_pos.append(atom_pos[i])
	for i in range(len(mag_atom_pos)):
		for j in range(len(mag_atom_pos)):
			tmp_dist = calcDistwithPBC(axis, mag_atom_pos[i][0:3], mag_atom_pos[j][0:3], cutoff)
			for dist in tmp_dist:
				list_pair.append([mag_atom_pos[i][4],mag_atom_pos[j][4],i,j,dist,mag_atom_pos[i][4]+mag_atom_pos[j][4]+str(dist)])

	return list_pair
			
# This function is for calculating distance considered periodic boundary condition
def calcDistwithPBC(lparam, ptA, ptB, cutoff):
        dist = []
        vecAB = array(ptA) - array(ptB)

        for i in [-2,-1,0,1,2]:
                for j in [-2,-1,0,1,2]:
                        for k in [-2,-1,0,1,2]:
                                vec_cur = copy(vecAB)
                                vec_cur = vec_cur + dot(array([i,j,k]), array(lparam))
                                dist_cur = round(linalg.norm(vec_cur), 2)

                                if  dist_cur < cutoff and dist_cur > 0.1:
                                        dist.append(dist_cur)

        return dist

def listsearch(list, item):
        check = 0
        for i in range(len(list)):
                if round((list[i] - item), 5) == 0.0:
                        check = 1
                        break

        if check == 1:
                return i
        else:
                return -1
