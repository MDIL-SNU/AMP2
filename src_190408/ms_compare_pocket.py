import os,sys,glob
from module_effm import *
import subprocess
from module_band import EIGEN_to_array
from module_vector import *
from module_vasprun import poscar_to_axis,pygrep,pyhead,pytail
import numpy as np


def pocket_check(target,band_path,carrier_type,E_width,search_space):
	spin = pygrep('ISPIN',band_path+'/OUTCAR',0,0).split()[2]
#	spin = subprocess.check_output(['grep','ISPIN',band_path+'/OUTCAR']).split()[2]
	ncl = pygrep('NONCOL',band_path+'/OUTCAR',0,0).split()[2]
#	ncl = subprocess.check_output(['grep','NONCOL',band_path+'/OUTCAR']).split()[2]
	[KPT,Band,nelect] = EIGEN_to_array(band_path+'/EIGENVAL',spin)
	axis = poscar_to_axis(band_path+'/POSCAR')
	rec_lat = reciprocal_lattice(axis)
	with open(band_path+'/Band_gap.log','r') as inp:
		lines = inp.readlines()
		VBM_k = lines[2].split()[1:4]
		VBM_E = float(lines[2].split()[5])
		VBM_n = int(lines[-2].split()[1])
		VBM_s = int(lines[-2].split()[3])

		CBM_k = lines[3].split()[1:4]
		CBM_E = float(lines[3].split()[5])
		CBM_n = int(lines[-1].split()[1])
		CBM_s = int(lines[-1].split()[3])

	targ_n = [0,0]
	if carrier_type == 'hole':
		sign = 1
		targ_E = VBM_E
		targ_n[VBM_s-1] = VBM_n-1
		pocket  = [[VBM_k,VBM_n,VBM_s,VBM_E]]
		if not VBM_s == CBM_s:
			targ_n[CBM_s-1] = CBM_n-2
		elif spin == '2' or not nelect%2 == 0 or ncl == 'T': # spin exists.
			for i in range(len(Band)):
				if Band[i][0][2-VBM_s] > targ_E:
					targ_n[2-VBM_s] = i-1
					break
	else:
		sign = -1
		targ_E = CBM_E
		targ_n[CBM_s-1] = CBM_n-1
		pocket  = [[CBM_k,CBM_n,CBM_s,CBM_E]]
		if not CBM_s == VBM_s:
			targ_n[VBM_s-1] = VBM_n
		elif spin == '2' or not nelect%2 == 0 or ncl == 'T': # spin exists.
			for i in range(len(Band)):
				if Band[i][0][2-CBM_s] > targ_E:
					targ_n[2-CBM_s] = i
					break
	KPT_num = [[],[]] # [spin_idx][kpt, energy, # of overlap kpt, index of overlap kpt]
	for i in range(int(spin)):
		for k in range(len(KPT)):
			swit = 0
			if abs(Band[targ_n[i]][k][i]-targ_E) < E_width:
				for j in range(len(KPT_num[i])):
					if KPT[k] == KPT_num[i][j][0]:
						KPT_num[i][j][2] = KPT_num[i][j][1] + 1
						KPT_num[i][j][3].append(k)
						swit = 1
						break
				if swit == 0:
					KPT_num[i].append([KPT[k],Band[targ_n[i]][k][i],1,[k]])


	extreme_log = [[1],[1]] # [spin idx][kpt_num idx]
	extreme = []
	for i in range(int(spin)):
		for k in range(1,len(KPT_num[i])):
			extreme_log[i].append(1)
			for j in range(k):
				if dist_vec(KPT_num[i][k][0],KPT_num[i][j][0],rec_lat) < search_space:
					if sign*(KPT_num[i][k][1] - KPT_num[i][j][1]) > 0:
						extreme_log[i][j] = 0
					else:
						extreme_log[i][k] = 0
		for k in range(len(KPT_num[i])):
			if extreme_log[i][k] == 1:
				extreme.append(KPT_num[i][k]+[i+1])

	for i in range(len(extreme)):
		swit = 1
		for j in range(len(pocket)):
			if dist_vec(extreme[i][0],pocket[j][0],rec_lat) < search_space:
				swit = 0
				break
		if swit == 1:
			pocket.append([extreme[i][0],targ_n[extreme[i][4]-1]+1,extreme[i][4],extreme[i][1]])

	with open(target+'/pocket_check.log','w') as out:
		for i in range(len(pocket)) :
			out.write(' '.join([str(x) for x in pocket[i][0]])+' :\t'+' '.join([str(x) for x in pocket[i][1:4]])+'\n')








path_list = glob.glob(sys.argv[1]+'/*')
for mat_path in path_list:
	if os.path.isdir(mat_path+'/effm_GGA/hole'):
#		if not os.path.isfile(mat_path+'/effm_GGA/hole/pocket_check.log'):
		dir_effm = mat_path+'/effm_GGA/hole'
		band_path = mat_path+'/band_GGA'
		pocket_check(dir_effm,band_path,'hole',calc_max_E_diff(0.99,300),0.1)

		pocket_path = mat_path+'/effm_GGA/hole/pocket_check.log'

		icsd = mat_path.split('_')[-1]
		with open(pocket_path,'r') as f:
			lines = f.readlines()
		if not len(lines) == 1:
			del_e = float(lines[0].split()[-1]) - max([float(lines[x+1].split()[-1]) for x in range(len(lines)-1)])
		else:
			del_e = 0.0
		with open(mat_path+'/effm_GGA/hole/inp_grid','r') as f:
			new_len = len(f.readlines())
		print mat_path+'\t'+icsd+'\t'+str(len(lines))+'\t'+str(del_e)+'\t'+str(new_len)
