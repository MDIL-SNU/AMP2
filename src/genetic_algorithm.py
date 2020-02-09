####################################
# date : 2019-04-08                #
# Author : yybbyb@snu.ac.kr        #
####################################
# This is for performing genetic algorithm to find the most stable magnetic spin ordering
import sys,os,yaml
from numpy import *
import genetic_operator as go
from module_AF import *
from module_amp2_input import *
from module_GA import *
from module_vector import *
from copy import deepcopy as dcopy

## read input
inp_file = sys.argv[1]
with open(inp_file,'r') as f:
	inp_yaml = yaml.safe_load(f)

mag_atom_list = inp_yaml['mag_atom_list']
[axis,atom_pos] = read_poscar(sys.argv[2])
atom_pos_cart = []
for i in range(len(atom_pos)):
	atom_pos_cart.append(dir_to_cart(atom_pos[i][0:3],axis)+atom_pos[i][3:])
mag_true = []
mag_sum = 0
for i in range(len(atom_pos)):
	if atom_pos[i][4] in mag_atom_list:
		mag_true.append(1)
		mag_sum = mag_sum + 1
	else:
		mag_true.append(0)

mut_length = 6
if mut_length*2 > mag_sum:
	mut_length = mag_sum/2

# get tot_population
population = [inp_yaml['population']['best'],inp_yaml['population']['crossover'],inp_yaml['population']['mutation'],inp_yaml['population']['random']]
tot_population = sum(population)
if tot_population > 2**(mag_sum-1):
	for i in range(4):
		population[i] = int(floor(float(population[i])/float(tot_population)*(2.0**(mag_sum-1))))
	tot_population = sum(population)

num_iter = mag_sum*10
if inp_yaml['min_iteration'] > num_iter:
	num_iter = inp_yaml['min_iteration']
elif inp_yaml['max_iteration'] > inp_yaml['min_iteration'] and inp_yaml['max_iteration'] < num_iter:
	num_iter = inp_yaml['max_iteration']

if mag_sum < 9:
	tot_population = 2**(mag_sum-1)
	num_iter = 1

info_pair = pairCategorize(atom_pos_cart, axis, mag_true, inp_yaml['cutoff'])
gene_list = genNextGeneration(mag_sum, [], [], [0,0,0,tot_population], mut_length)

list_portion = [[], [], []]

param_list = open(inp_yaml['pair_coeff'],'r').readlines()
#param_val = []
#param_len = []
para_J = {}
# set para_J
for ll in param_list:
	line = ll.split()
	param_len = around(float(line[0]),decimals=2)
	param_val = float(line[1])
	if line[2] == line[3]:
		para_J[line[2]+line[3]+str(param_len)] = param_val
		para_J[line[2]+line[3]+str(param_len+0.01)] = param_val
		para_J[line[2]+line[3]+str(param_len-0.01)] = param_val
	else:
		para_J[line[2]+line[3]+str(param_len)] = param_val
		para_J[line[3]+line[2]+str(param_len)] = param_val
		para_J[line[2]+line[3]+str(param_len+0.01)] = param_val
		para_J[line[3]+line[2]+str(param_len+0.01)] = param_val
		para_J[line[2]+line[3]+str(param_len-0.01)] = param_val
		para_J[line[3]+line[2]+str(param_len-0.01)] = param_val

ERES = open('energy.dat', 'w')
#print para_J
energy_gene = getEnergy(gene_list, info_pair, para_J)
#print energy_gene
ERES.write(str(min(energy_gene)) + '\n')
prob_list = energyToProbab(energy_gene, inp_yaml['selection_coeff'])


###### run_GA ##########
for i in range(1, num_iter):
	gene_list = genNextGeneration(mag_sum, gene_list, prob_list, population, mut_length)
	energy_gene = getEnergy(gene_list,info_pair, para_J)

	ERES.write(str(min(energy_gene)) + '\n')

	prob_list = energyToProbab(energy_gene, 2)

ERES.close()
####################

num_min_gene = 10
if len(gene_list) < num_min_gene:
	num_min_gene = len(gene_list)

min_gene = [[x,y] for (y,x) in sorted(zip(energy_gene, gene_list))][0:num_min_gene]
energy_tolerance = float(sys.argv[3])*mag_sum # eV
loop_num = 0
for kk in range(num_min_gene):
	if abs(round(min_gene[kk][1]-min_gene[0][1],5)) < energy_tolerance:
		loop_num = kk+1

# make new poscar from the result of G.A.
for kk in range(loop_num):
	idx = 0
	tmp_up = []
	tmp_down = []
	new_atom_pos = []
	if mag_true[0] == 1:
		if min_gene[kk][0][idx] == 1:
			tmp = dcopy(atom_pos[0])
			tmp[4] = tmp[4]+'_up'
			tmp_up.append(tmp)
		else:
			tmp = dcopy(atom_pos[0])
			tmp[4] = tmp[4]+'_down'
			tmp_down.append(tmp)
		idx = idx+1
	else:
		new_atom_pos.append(atom_pos[0])
	prev_type = atom_pos[0][4]
	for i in range(1,len(atom_pos)):
		if not prev_type == atom_pos[i][4] and mag_true[i-1] == 1:
			for line in tmp_up:
				new_atom_pos.append(line)
			for line in tmp_down:
				new_atom_pos.append(line)
			tmp_up = []
			tmp_down = []

		if mag_true[i] == 1:
			if min_gene[kk][0][idx] == 1:
				tmp = dcopy(atom_pos[i])
				tmp[4] = tmp[4]+'_up'
				tmp_up.append(tmp)
			else:
				tmp = dcopy(atom_pos[i])
				tmp[4] = tmp[4]+'_down'
				tmp_down.append(tmp)

			idx = idx+1
		else:
			new_atom_pos.append(atom_pos[i])
	
		prev_type = atom_pos[i][4]

	if mag_true[-1] == 1:
		for line in tmp_up:
			new_atom_pos.append(line)
		for line in tmp_down:
			new_atom_pos.append(line)
		tmp_up = []
		tmp_down = []

	write_poscar(axis,new_atom_pos,'POSCAR_spin'+str(kk),'New poscar from G.A.')
#	print min_gene[kk]

