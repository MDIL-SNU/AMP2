import os, sys, subprocess, yaml,shutil,math
from module_log import *
from module_vasprun import *
from module_converge import *
from module_relax import *
from module_band import *
from module_effm import *
from module_AF import *
from module_amp2_input import *
import spglib
import numpy as np

def calc_effm_new(target,carrier_type,Temp,weight):
	import scipy.constants as sc
	scp = sc.physical_constants
	with open(target+'/inp_grid') as inp:
		number_of_pocket = int(inp.readline())
		grid_size = float(inp.readline())
#		max_E_diff = 0.1
		max_E_diff = float(inp.readline())
#		max_E_diff = 0.2
		max_E_diff = calc_max_E_diff(0.90,300)
		pocket_kpt = []
		pocket_inform = []
		num_kp = []
		for k in range(number_of_pocket):
			pocket_kpt.append([float(x) for x in inp.readline().split()])
			num_kp.append([])
			for i in range(3):
				num_kp[k].append(sum([int(x) for x in inp.readline().split()],1)) # num_negative_direction + num_positive_direction + 1
			pocket_inform.append([int(x) for x in inp.readline().split()])
#		new_line = inp.readline().split()
#		for k in range(number_of_pocket):
#			pocket_inform.append([int(new_line[0]),0])
		extreme_E = float(inp.readline())

	Band = PROCAR_to_array(target+'/PROCAR') #Band[band_idx][kpt_idx][spin_idx]
	diff_K = grid_size * np.pi * 2.0

	EN_tot = []
	Edk2_tot = []
	wei_tot = []
	shift = 0
	for k in range(number_of_pocket):
		if carrier_type == 'hole':
			min_band_idx = 0
			max_band_idx = pocket_inform[k][0]+1
		else:
			min_band_idx = pocket_inform[k][0]
			max_band_idx = len(Band)

		spin_idx = pocket_inform[k][1]
		for n in range(min_band_idx,max_band_idx):
			for x_id in range(1,num_kp[k][0]-3):
				for y_id in range(1,num_kp[k][1]-3):
					for z_id in range(1,num_kp[k][2]-3):
						energy = Band[n][idx_change(x_id,y_id,z_id,num_kp[k],shift)][spin_idx]
						if abs(energy-extreme_E) < max_E_diff:
							EN_tot.append(energy)
							Edk2_tot.append(cal_dk2(Band,n,x_id,y_id,z_id,num_kp[k],shift,spin_idx,diff_K))
							wei_tot.append(weight[k])
		shift = shift + num_kp[k][0] * num_kp[k][1] * num_kp[k][2]

	sum_fermi = 0
	sum_deriv = [0,0,0,0,0,0]
	for i in range(len(EN_tot)):
		sum_fermi = sum_fermi + Fermi_dist(EN_tot[i],extreme_E,carrier_type,Temp) * wei_tot[i]
		for j in range(6):
			sum_deriv[j] = sum_deriv[j] + Edk2_tot[i][j] * Fermi_dist(EN_tot[i],extreme_E,carrier_type,Temp) * wei_tot[i]

	deriv_tensor = np.array([[sum_deriv[0],sum_deriv[3],sum_deriv[4]],[sum_deriv[3],sum_deriv[1],sum_deriv[5]],[sum_deriv[4],sum_deriv[5],sum_deriv[2]]])/sum_fermi
	deriv_mat = np.linalg.inv(deriv_tensor)
	deriv_dia = np.linalg.eigvals(deriv_tensor)

	effm_dia = []
	effm = []
	for i in range(3):
		effm_dia.append(scp['natural unit of action in eV s'][0]**2.0*scp['atomic unit of charge'][0]*1.0e20/(deriv_dia[i])/sc.m_e)
		effm.append([scp['natural unit of action in eV s'][0]**2.0*scp['atomic unit of charge'][0]*1.0e20*(deriv_mat[i][x])/sc.m_e for x in range(3)])

	return [effm_dia,effm]



def vec_mod(vec):
	new_vec = []
	for xx in vec:
		new_vec.append(round(xx-np.floor(xx+0.5),7))
	return new_vec

def vec_same(vec1,vec2):
	summ = sum([abs(vec1[x]-vec2[x]) for x in range(len(vec1))])
	if summ < 0.00001:
		return 1
	else:
		return 0

def dir_to_cart(vec,axis):
	cart_vec = []
	for i in range(3):
		cart_vec.append(vec[0]*axis[0][i]+vec[1]*axis[1][i]+vec[2]*axis[2][i])
	return cart_vec


def dup_pos(pos):
	irre_pos = []
	made_pos = []
	idx = 0
	dic = {}

	for i in range(len(pos)):
		for j in range(len(oper)):
			made_pos.append(vec_mod(dir_to_cart(pos[i],oper[j])))
	for i in range(len(made_pos)):
		if len(irre_pos) == 0:
			irre_pos.append(made_pos[i])
		else:
			sw = 0
			for j in range(len(irre_pos)):
				if vec_same(made_pos[i],irre_pos[j]):
					sw = 1
					break
			if sw == 0:
				irre_pos.append(made_pos[i])
	return irre_pos


dir_effm = sys.argv[2]
band_path = sys.argv[3]
pocket_path = dir_effm+'/pocket.log'


[axis,atom_pos] = read_poscar(sys.argv[1])

L = np.mat(axis)
pos = []
atom_type = []
atom_dic = {}
type_dic = {}
index = 0
for line in atom_pos:
	pos.append(line[0:3])
	if not line[4] in atom_dic.keys():
		atom_dic[line[4]] = index
		type_dic[index] = [line[3],line[4]]
		index = index+1
	atom_type.append(atom_dic[line[4]])
D = np.mat(pos)
Cell = (L,D,atom_type)

rot_oper = spglib.get_symmetry(Cell)['rotations']
oper = np.ndarray.tolist(rot_oper)

print calc_max_E_diff(0.98,300)
if not os.path.isfile(pocket_path):
	pocket(dir_effm,band_path,'hole',calc_max_E_diff(0.99,300),0.1)

with open(pocket_path,'r') as f:
	lines = f.readlines()

poc_kp = []
for line in lines:
	poc_kp.append([[float(x) for x in line.split()[0:3]]])

weight = []
for poc in poc_kp:
	prev_poc_len = len(poc)
	while 1:
		prev_poc_len = len(poc)
		poc = dup_pos(poc)
		if len(poc) == prev_poc_len:
			break
		else:
			prev_poc_len = len(poc)

	weight.append(1)
#	weight.append(len(poc))

[effm_dia,effm] = calc_effm_new(dir_effm,'hole',300,weight)
for i in range(3):
	print str(effm[i][0])+'\t'+str(effm[i][1])+'\t'+str(effm[i][2])
print str(effm_dia[0])+'\t'+str(effm_dia[1])+'\t'+str(effm_dia[2])
