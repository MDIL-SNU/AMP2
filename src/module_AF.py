####################################
# date : 2019-04-08                #
# Author : yybbyb@snu.ac.kr        #
####################################
import os,sys,yaml,subprocess
import numpy as np
from module_vector import *
from module_vasprun import *
from module_amp2_input import *

def check_spin(ref_dir,inp_pos,min_mom):
	[axis,atom_pos] = read_poscar(inp_pos)
	magnet_out = pygrep('magnetization (x)',ref_dir+'/OUTCAR',0,len(atom_pos)+3).splitlines()[-len(atom_pos):]
#	magnet_out = subprocess.check_output(['grep','-A'+str(len(atom_pos)+3),'magnetization (x)',ref_dir+'/OUTCAR']).splitlines()[-len(atom_pos):]
	magnet_list = [float(x.split()[-1]) for x in magnet_out]
	type_name = []
#	type_num = []
	mag_val = {}
	for i in range(len(atom_pos)):
		over = 0
		for j in range(len(type_name)):
			if type_name[j] == atom_pos[i][4]:
				over = 1
				break
		if over == 0:
			type_name.append(atom_pos[i][4])
			mag_val[atom_pos[i][4]] = magnet_list[i]
	new_type_name = []
	for typ in type_name:
		if abs(mag_val[typ]) > min_mom:
			new_type_name.append(typ)
		else:
			mag_val[typ] = 0
	return new_type_name,mag_val

def find_pair(poscar,cutoff,mag_atom_list,mag_val,tol):
	[axis,atom_pos] = read_poscar(poscar)
	sole_list = []
	pair_list = []
	mag_list = []
	for i in range(len(atom_pos)):
		if not atom_pos[i][4] in mag_atom_list:
			mag_list.append(0)
		else:
			mag_list.append(mag_val[atom_pos[i][4]])
			if len(sole_list) == 0:
				sole_list.append([atom_pos[i][4],i])
			else:
				if not atom_pos[i][4] in [x[0] for x in sole_list]:
					sole_list.append([atom_pos[i][4],i])
			for j in range(i,len(atom_pos)):
				if atom_pos[j][4] in mag_atom_list:
					length = short_dist(atom_pos[i][0:3],atom_pos[j][0:3],axis)
					if length <= cutoff and length > 0.1:
						if len(pair_list) == 0:
							pair_list.append([length,atom_pos[i][4],atom_pos[j][4],i,j])
						else:
							dup = 0
							for k in range(len(pair_list)):
								if abs(length-pair_list[k][0]) < tol:
									if atom_pos[i][4] == pair_list[k][1] and atom_pos[j][4] == pair_list[k][2] :
										dup = 1
										break
									elif atom_pos[i][4] == pair_list[k][2] and atom_pos[j][4] == pair_list[k][1] :
										dup = 1
										break
 							if dup == 0:
								pair_list.append([length,atom_pos[i][4],atom_pos[j][4],i,j])
	return [pair_list,sole_list,mag_list]

def write_pair_list(pair_list,sole_list,mag_list,target):
	with open(target+'/pair_list_table.dat','w') as plw:
		for line in pair_list:
			plw.write('\t'.join([str(x) for x in line])+'\n')
	tot_mag_list = []
	[spin,num_spin] = reduce_mag_list(mag_list)
	tot_mag_list.append(' '.join([str(num_spin[x])+'*'+str(spin[x]) for x in range(len(spin))]))
	for line in sole_list:
		new_mag_list = mag_list[0:line[1]]+[-mag_list[line[1]]]+mag_list[line[1]+1:]
		[spin,num_spin] = reduce_mag_list(new_mag_list)
		tot_mag_list.append(' '.join([str(num_spin[x])+'*'+str(spin[x]) for x in range(len(spin))]))
	for line in pair_list:
		new_mag_list = mag_list[0:line[3]]+[-mag_list[line[3]]]+mag_list[line[3]+1:line[4]]+[-mag_list[line[4]]]+mag_list[line[4]+1:]
		[spin,num_spin] = reduce_mag_list(new_mag_list)
		tot_mag_list.append(' '.join([str(num_spin[x])+'*'+str(spin[x]) for x in range(len(spin))]))
	with open(target+'/mag_list.dat','w') as mlw:
		for line in tot_mag_list:
			mlw.write(line+'\n')
	return tot_mag_list

def reduce_mag_list(mag_list):
	spin = []
	num_spin = []
	for i in mag_list:
		if len(spin) == 0:
			spin.append(i)
			num_spin.append(1)
		else:
			if i == spin[-1]:
				num_spin[-1] = num_spin[-1]+1
			else:
				spin.append(i)
				num_spin.append(1)
	return [spin,num_spin]

def calc_ising_param(sole_list,pair_list,ene_file):
	with open(ene_file,'r') as enef:
		lines = enef.readlines()
	ENE = {}
	ENE['Ferro'] = float(lines[0].split()[1])
	for i in range(len(sole_list)):
		ENE[sole_list[i][0]] = float(lines[1+i].split()[1])
	ene_pair = []
	JJ = []
	for i in range(len(pair_list)):
		pair_ene = float(lines[len(sole_list)+i+1].split()[1])
		JJ.append(pair_list[i][0:3]+[(ENE['Ferro']+pair_ene-ENE[pair_list[i][1]]-ENE[pair_list[i][2]])/4])
	return JJ

def write_ising_param(JJ,target):
	with open(target+'/Pair_coeff.dat','w') as pcw:
		for line in JJ:
			pcw.write(str(line[0])+'\t'+str(line[3])+'\t'+line[1]+'\t'+line[2]+'\n')

def write_inp_for_GA(mag_atom_list,inp_af,target):
	inp_yaml = {}
	for key in inp_af['genetic_algorithm'].keys():
		inp_yaml[key] = inp_af['genetic_algorithm'][key]
	inp_yaml['mag_atom_list'] = mag_atom_list
	inp_yaml['cutoff'] = inp_af['cutoff_for_parameter']
	inp_yaml['pair_coeff'] = target+'/Pair_coeff.dat'
	with open(target+'/input_GA.yaml','w') as inp_GA:
		yaml.dump(inp_yaml,inp_GA,default_flow_style = False)

def count_mag_atom_num(poscar,mag_atom_list):
	[axis,atom_pos] = read_poscar(poscar)
	mag_atom_num = 0
	for line in atom_pos:
		if line[4] in mag_atom_list:
			mag_atom_num = mag_atom_num + 1
	return mag_atom_num

def set_supercell_list(poscar,mag_atom_list):
	sup_size_list = [1,2]
	supercell_list = []
	mag_atom_num = count_mag_atom_num(poscar,mag_atom_list)
	for i in sup_size_list:
		for j in sup_size_list:
			for k in sup_size_list:
				if i*j*k*mag_atom_num % 2 == 0:
					supercell_list.append([i,j,k,i*j*k])
	sort_supercell_list = sorted(supercell_list,key=lambda x : x[3])
	supercell_list = []
	for line in sort_supercell_list:
		supercell_list.append(line[0:3])
	return supercell_list

def get_bond_matrix(poscar):
	[axis,atom_pos] = read_poscar(poscar)
	bond_matrix = []
	for i in range(len(atom_pos)):
		tmp = []
		for j in range(len(atom_pos)):
			tmp.append([round(short_dist(atom_pos[i][0:3],atom_pos[j][0:3],axis),6),atom_pos[i][4],atom_pos[j][4]])
		bond_matrix.append(sorted(sorted(tmp,key = lambda x:x[0]),key = lambda x:x[2]))

	sorted_bond_matrix = sorted(sorted(bond_matrix,key = lambda x:x[0][0]),key = lambda x:x[0][1])
	
	return sorted_bond_matrix

def resized_kpoints(ref_path,targ_path):
	sym = int(open(targ_path+'/sym','r').readline())
	KPT_ref = [float(x) for x in open(ref_path+'/KPOINTS','r').readlines()[3].split()]
	axis = poscar_to_axis(ref_path+'/POSCAR')
	recipro_latt = reciprocal_lattice(axis)
	l = []
	for i in range(3):
		l.append((recipro_latt[i][0]**2.+recipro_latt[i][1]**2.+recipro_latt[i][2]**2.)**0.5)
	min_dk = max([x/y for x,y in zip(l,KPT_ref)])
	KP = []

	axis = poscar_to_axis(targ_path+'/POSCAR')
	recipro_latt = reciprocal_lattice(axis)
	l = []
	for i in range(3):
		l.append((recipro_latt[i][0]**2.+recipro_latt[i][1]**2.+recipro_latt[i][2]**2.)**0.5)
	
	kpoint = open(targ_path+'/KPOINTS','w')
	# Gamma-centred mesh for hexagoanl symmetry
	if sym==12 or sym==13 or sym==14:
		KPset = 'Gamma-centered'
	else :
		KPset = 'Monk-horst'
	for i in range(3) :
		if sym in [5,6,10]: # symmetry BCT
			KP.append(str(max([int(round(x/min_dk)) for x in l])))
		else:
			KP.append(str(int(round(l[i]/min_dk))))
		if KP[-1] == '0' :
			KP[-1] = '1'
	# Write KPOINTS
	kpoint.write("Auto k-point\n 0\n"+KPset+"\n  "+KP[0]+"  "+KP[1]+"  "+KP[2]+"\n  0  0  0\n")
	kpoint.close()

def set_mag_info(atom_pos,target,mag_val):
	magmom = []
	for line in atom_pos:
		if len(line[4].split('_')) > 1:
			sline = line[4].split('_')
			mag_type = '_'.join(sline[:-1])
			if sline[-1] == 'up':
				magmom.append(mag_val[mag_type])
			elif sline[-1] == 'down': 
				magmom.append(-1.0*mag_val[mag_type])
			else:
				magmom.append(mag_val[line[4]])
		else:
			magmom.append(mag_val[line[4]])
	[spin,num_spin] = reduce_mag_list(magmom)
	with open(target+'/mag_info','w') as inp:
		inp.write(' '.join([str(num_spin[x])+'*'+str(spin[x]) for x in range(len(spin))]))
