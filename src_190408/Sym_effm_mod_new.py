import os, sys, subprocess, yaml,shutil,math
from module_log import *
from module_vasprun import *
from module_vector import *
from module_converge import *
from module_relax import *
from module_band import *
from module_effm import *
from module_AF import *
from module_amp2_input import *
import spglib
import numpy as np
def calc_effm_new(target,carrier_type,Temp,oper):
	import scipy.constants as sc
	scp = sc.physical_constants
	pos_file = target+'/POSCAR'
	rec_axis = reciprocal_lattice(poscar_to_axis(pos_file))
	with open(target+'/inp_grid') as inp:
		number_of_pocket = int(inp.readline())
		grid_size = float(inp.readline())
		max_E_diff = float(inp.readline())
		pocket_kpt = []
		frac_pocket_kpt = []
		pocket_inform = []
		num_kp = []
		num_kp_sep = []
		for k in range(number_of_pocket):
			pocket_kpt.append([float(x) for x in inp.readline().split()])
			frac_pocket_kpt.append(cart_to_dir(pocket_kpt[-1],rec_axis))
			num_kp.append([])
			num_kp_sep.append([])
			for i in range(3):
				num_kp_sep[k].append([int(x) for x in inp.readline().split()])
				num_kp[k].append(1+num_kp_sep[k][i][0]+num_kp_sep[k][i][1]) # num_negative_direction + num_positive_direction + 1
			pocket_inform.append([int(x) for x in inp.readline().split()])
		extreme_E = float(inp.readline())

	Band = PROCAR_to_array(target+'/PROCAR') #Band[band_idx][kpt_idx][spin_idx]
	diff_K = grid_size * np.pi * 2.0

	EN_tot = []
	Edk2_tot = []
	wei_tot = []
	shift = 0
	unique_kpt_idx = []
	for k in range(number_of_pocket):
		if carrier_type == 'hole':
			min_band_idx = 0
			max_band_idx = pocket_inform[k][0]+1
		else:
			min_band_idx = pocket_inform[k][0]
			max_band_idx = len(Band)

		# check pocket overlap
		sw_overlap = 0
		if not k == 0:
			sym_kpts = symmetrical_kpts(frac_pocket_kpt[k],oper)
			for new_kpt in sym_kpts:
				for i in unique_kpt_idx:
					if valid_pocket(dir_to_cart(new_kpt,rec_axis),pocket_kpt[i],num_kp_sep[i],grid_size,rec_axis) == 0:
						sw_overlap = 1
						break
				if sw_overlap == 1:
					break
	
		if sw_overlap == 0:
			unique_kpt_idx.append(k)
			sym_kpts = symmetrical_kpts(frac_pocket_kpt[k],oper)
			num_overlap_sym_kpt = 0
			for new_kpt in sym_kpts:
				if valid_pocket(dir_to_cart(new_kpt,rec_axis),pocket_kpt[k],num_kp_sep[k],grid_size,rec_axis) == 0:
					num_overlap_sym_kpt = num_overlap_sym_kpt +1
			weight = float(len(sym_kpts))/float(num_overlap_sym_kpt)
			print weight
			spin_idx = pocket_inform[k][1]
			for n in range(min_band_idx,max_band_idx):
				for x_id in range(1,num_kp[k][0]-1):
					for y_id in range(1,num_kp[k][1]-1):
						for z_id in range(1,num_kp[k][2]-1):
							energy = Band[n][idx_change(x_id,y_id,z_id,num_kp[k],shift)][spin_idx]
							sw_for_kp_in_cell = kpt_frac_in_cell(x_id-num_kp_sep[k][0][0],y_id-num_kp_sep[k][1][0],z_id-num_kp_sep[k][2][0],grid_size,rec_axis)
							if abs(energy-extreme_E) < max_E_diff and sw_for_kp_in_cell == 1:
								EN_tot.append(energy)
								Edk2_tot.append(cal_dk2(Band,n,x_id,y_id,z_id,num_kp[k],shift,spin_idx,diff_K))
								wei_tot.append(weight)

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


poscar_type_path = sys.argv[1]
dir_effm = sys.argv[2]

oper = read_operation(poscar_type_path)

[effm_dia,effm] = calc_effm_new(dir_effm,'hole',300,oper)
for i in range(3):
	print str(effm[i][0])+'\t'+str(effm[i][1])+'\t'+str(effm[i][2])
print str(effm_dia[0])+'\t'+str(effm_dia[1])+'\t'+str(effm_dia[2])
