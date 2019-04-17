###########################################
### Date: 2019-01-22			###
### yybbyb@snu.ac.kr			###
###########################################
import subprocess,os
from module_band import EIGEN_to_array
from module_vector import *
from module_vasprun import poscar_to_axis
import numpy as np

def calc_max_E_diff(Fermi_cut,temp):
	import scipy.constants as sc
	scp = sc.physical_constants
	max_E_diff = -scp['Boltzmann constant in eV/K'][0]*temp*np.log(1./Fermi_cut-1.)
	return max_E_diff

def pocket(target,band_path,carrier_type,E_width,search_space):
	spin = subprocess.check_output(['grep','ISPIN',band_path+'/OUTCAR']).split()[2]
	ncl = subprocess.check_output(['grep','NONCOL',band_path+'/OUTCAR']).split()[2]
	[KPT,Band,nelect] = EIGEN_to_array(band_path+'/EIGENVAL',spin)
	axis = poscar_to_axis(band_path+'/POSCAR')
	rec_lat = reciprocal_lattice(axis)
	with open(band_path+'/Band_gap.log','r') as inp:
		lines = inp.readlines()
		VBM_k = lines[2].split()[1:4]
		VBM_E = float(lines[2].split()[5])
		VBM_n = int(lines[5].split()[1])
		VBM_s = int(lines[5].split()[3])

		CBM_k = lines[3].split()[1:4]
		CBM_E = float(lines[3].split()[5])
		CBM_n = int(lines[6].split()[1])
		CBM_s = int(lines[6].split()[3])

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

	with open(target+'/pocket.log','w') as out:
		for i in range(len(pocket)) :
			out.write(' '.join([str(x) for x in pocket[i][0]])+' :\t'+' '.join([str(x) for x in pocket[i][1:4]])+'\n')

def make_kpts_for_searching_space(target,search_grid_size):
	axis = poscar_to_axis(target+'/POSCAR')
	rec_lat = reciprocal_lattice(axis)
	pocket_kpt = []
	pocket_inform = []
	with open(target+'/pocket.log','r') as inp:
		for line in inp.readlines():
			pocket_kpt.append([float(x) for x in line.split()[0:3]])
			pocket_inform.append([int(x) for x in line.split()[4:6]]+[float(line.split()[6])])
	KPTD = [] # kpoints test direction
	KPTD.append([0.5, 0, 0])
	KPTD.append([0, 0.5, 0])
	KPTD.append([0, 0, 0.5])
	KPTD.append([0.5, 0.5, 0])
	KPTD.append([0.5, -0.5, 0])
	KPTD.append([0.5, 0, 0.5])
	KPTD.append([0.5, 0, -0.5])
	KPTD.append([0, 0.5, 0.5])
	KPTD.append([0, 0.5, -0.5])

	KPTDC = [] # cartessian KPTD
	for i in range(len(KPTD)):
		KPTDC.append(dir_to_cart(KPTD[i],rec_lat))

	KPT = []
	with open(target+'/NKPT','w') as out:
		out.write(str(len(pocket_kpt))+'\n')
		num_kp_for_test = []
		for i in range(len(KPTDC)):
			length_vec = dist_point(KPTDC[i],[0,0,0]) # length of KPTDC[i]
			num_kp_for_test.append(int(np.ceil(length_vec/search_grid_size)))

		for k in range(len(pocket_kpt)):
			pocket_kpt_cart = dir_to_cart(pocket_kpt[k],rec_lat)
			out.write('\t'.join([str(x) for x in pocket_kpt_cart])+'\t'+' '.join([str(x) for x in pocket_inform[k]])+'\n')
			for i in range(len(KPTDC)):
				for j in range(2*num_kp_for_test[i]+1):
					KPT.append([-KPTDC[i][x]+KPTDC[i][x]*float(j)/float(num_kp_for_test[i])+pocket_kpt_cart[x] for x in range(3)])

		out.write(' '.join([str(2*num_kp_for_test[x]+1) for x in range(len(num_kp_for_test))]))

	with open(target+'/KPOINTS','w') as kpt_out:
		kpt_out.write('K-pts for searching space\n')
		kpt_out.write('\t\t'+str(len(KPT))+'\n')
		kpt_out.write('car\n')
		for i in range(len(KPT)):
			kpt_out.write('\t'.join([str(KPT[i][x]) for x in range(3)])+'\t 1 \n')

def make_kpts_for_calculation(target,grid_size,max_E_diff):
	with open(target+'/NKPT','r') as inp:
		nkpt = inp.readlines()
		number_of_pocket = int(nkpt[0])
		pocket_kpt = []
		pocket_inform = []
		for k in range(number_of_pocket):
			pocket_kpt.append([float(x) for x in nkpt[k+1].split()[0:3]])
			pocket_inform.append([int(x)-1 for x in nkpt[k+1].split()[3:5]]+[float(nkpt[k+1].split()[5])]) # [band_idx,spin_idx,energy]
		num_kp_for_test = [int(x) for x in nkpt[-1].split()]

	Band = PROCAR_to_array(target+'/PROCAR') #Band[band_idx][kpt_idx][spin_idx]
	KPT = KPOINTS_to_array(target+'/KPOINTS')
	extreme_E = pocket_inform[0][2]


	shift = 0 
	KPT_for_effm = []
	num_kp_for_effm = []
	for k in range(number_of_pocket):
		kp_max = [x for x in pocket_kpt[k]]
		kp_min = [x for x in pocket_kpt[k]]
		for i in range(len(num_kp_for_test)):
			for j in range(num_kp_for_test[i]):
				if abs(extreme_E - Band[pocket_inform[k][0]][shift+j][pocket_inform[k][1]]) < max_E_diff:
					for k_id in range(3):
						if KPT[shift+j][k_id] > KPT[shift+(num_kp_for_test[i]-1)/2][k_id]:
							if kp_max[k_id] < KPT[shift+j][k_id]:
								kp_max[k_id] = KPT[shift+j][k_id]
						else:
							if kp_min[k_id] > KPT[shift+j][k_id]:
								kp_min[k_id] = KPT[shift+j][k_id]
			shift = shift + num_kp_for_test[i]

		num_kp_for_effm.append([]) # num_kp_for_effm[pocket_idx][kpt_idx][min or max]
		for k_id in range(3):
			num_kp_for_effm[k].append([int(np.ceil((pocket_kpt[k][k_id]-kp_min[k_id])/grid_size))+1,int(np.ceil((kp_max[k_id]-pocket_kpt[k][k_id])/grid_size))+1])
		for x_id in range(num_kp_for_effm[k][0][0]+num_kp_for_effm[k][0][1]+1):
			for y_id in range(num_kp_for_effm[k][1][0]+num_kp_for_effm[k][1][1]+1):
				for z_id in range(num_kp_for_effm[k][2][0]+num_kp_for_effm[k][2][1]+1):
					id_list = [x_id,y_id,z_id]
					KPT_for_effm.append([pocket_kpt[k][x]+float(id_list[x]-num_kp_for_effm[k][x][0])*grid_size for x in range(3)])

	with open(target+'/KPOINTS','w') as kpt_out:
		kpt_out.write('K-pts to calculate effective mass\n')
		kpt_out.write('\t\t'+str(len(KPT_for_effm))+'\n')
		kpt_out.write('car\n')
		for i in range(len(KPT_for_effm)):
			kpt_out.write('\t'.join([str(KPT_for_effm[i][x]) for x in range(3)])+'\t 1 \n')

	with open(target+'/inp_grid','w') as inp_effm:
		inp_effm.write(str(number_of_pocket)+'\n')
		inp_effm.write(str(grid_size)+'\n')
		inp_effm.write(str(max_E_diff)+'\n')
		for k in range(number_of_pocket):
			inp_effm.write('\t'.join([str(pocket_kpt[k][x]) for x in range(3)])+'\n')
			for i in range(3):
				inp_effm.write('\t'.join([str(num_kp_for_effm[k][i][x]) for x in range(2)])+'\n')
			inp_effm.write('\t'.join([str(pocket_inform[k][x]) for x in range(2)])+'\n')
		inp_effm.write(str(extreme_E))

def calc_effm(target,carrier_type,Temp):
	import scipy.constants as sc
	scp = sc.physical_constants
	with open(target+'/inp_grid') as inp:
		number_of_pocket = int(inp.readline())
		grid_size = float(inp.readline())
		max_E_diff = float(inp.readline())
		pocket_kpt = []
		pocket_inform = []
		num_kp = []
		for k in range(number_of_pocket):
			pocket_kpt.append([float(x) for x in inp.readline().split()])
			num_kp.append([])
			for i in range(3):
				num_kp[k].append(sum([int(x) for x in inp.readline().split()],1)) # num_negative_direction + num_positive_direction + 1
			pocket_inform.append([int(x) for x in inp.readline().split()])
		extreme_E = float(inp.readline())

	Band = PROCAR_to_array(target+'/PROCAR') #Band[band_idx][kpt_idx][spin_idx]
	diff_K = grid_size * np.pi * 2.0

	EN_tot = []
	Edk2_tot = []

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
			for x_id in range(1,num_kp[k][0]-1):
				for y_id in range(1,num_kp[k][1]-1):
					for z_id in range(1,num_kp[k][2]-1):
						energy = Band[n][idx_change(x_id,y_id,z_id,num_kp[k],shift)][spin_idx]
						if abs(energy-extreme_E) < max_E_diff:
							EN_tot.append(energy)
							Edk2_tot.append(cal_dk2(Band,n,x_id,y_id,z_id,num_kp[k],shift,spin_idx,diff_K))
		shift = shift + num_kp[k][0] * num_kp[k][1] * num_kp[k][2]

	sum_fermi = 0
	sum_deriv = [0,0,0,0,0,0]
	for i in range(len(EN_tot)):
		sum_fermi = sum_fermi + Fermi_dist(EN_tot[i],extreme_E,carrier_type,Temp)
		for j in range(6):
			sum_deriv[j] = sum_deriv[j] + Edk2_tot[i][j] * Fermi_dist(EN_tot[i],extreme_E,carrier_type,Temp)

	deriv_tensor = np.array([[sum_deriv[0],sum_deriv[3],sum_deriv[4]],[sum_deriv[3],sum_deriv[1],sum_deriv[5]],[sum_deriv[4],sum_deriv[5],sum_deriv[2]]])/sum_fermi
	deriv_mat = np.linalg.inv(deriv_tensor)
	deriv_dia = np.linalg.eigvals(deriv_tensor)

	effm_dia = []
	effm = []
	for i in range(3):
		effm_dia.append(scp['natural unit of action in eV s'][0]**2.0*scp['atomic unit of charge'][0]*1.0e20/(deriv_dia[i])/sc.m_e)
		effm.append([scp['natural unit of action in eV s'][0]**2.0*scp['atomic unit of charge'][0]*1.0e20*(deriv_mat[i][x])/sc.m_e for x in range(3)])

	return [effm_dia,effm]

def write_effm(effm_dia,effm,target,carrier_type):
	with open(target+'/effective_mass.log','w') as out:
		out.write(carrier_type+'\n')
		for i in range(3):
			out.write('\t'.join([str(effm[i][x]) for x in range(3)])+'\n')
		out.write('\t'.join([str(effm_dia[x]) for x in range(3)]))

def cal_dk2(Band,n,x_id,y_id,z_id,num,shift,spin_idx,diff_K):
	dxx = (Band[n][idx_change(x_id+1,y_id,z_id,num,shift)][spin_idx] + Band[n][idx_change(x_id-1,y_id,z_id,num,shift)][spin_idx] \
		- 2.0*Band[n][idx_change(x_id,y_id,z_id,num,shift)][spin_idx])/diff_K**2.0
	dyy = (Band[n][idx_change(x_id,y_id+1,z_id,num,shift)][spin_idx] + Band[n][idx_change(x_id,y_id-1,z_id,num,shift)][spin_idx] \
		- 2.0*Band[n][idx_change(x_id,y_id,z_id,num,shift)][spin_idx])/diff_K**2.0
	dzz = (Band[n][idx_change(x_id,y_id,z_id+1,num,shift)][spin_idx] + Band[n][idx_change(x_id,y_id,z_id-1,num,shift)][spin_idx] \
		- 2.0*Band[n][idx_change(x_id,y_id,z_id,num,shift)][spin_idx])/diff_K**2.0
	dxdy = (Band[n][idx_change(x_id+1,y_id+1,z_id,num,shift)][spin_idx] + Band[n][idx_change(x_id-1,y_id-1,z_id,num,shift)][spin_idx] \
		- Band[n][idx_change(x_id+1,y_id-1,z_id,num,shift)][spin_idx] - Band[n][idx_change(x_id-1,y_id+1,z_id,num,shift)][spin_idx])/diff_K**2.0/4.0
	dxdz = (Band[n][idx_change(x_id+1,y_id,z_id+1,num,shift)][spin_idx] + Band[n][idx_change(x_id-1,y_id,z_id-1,num,shift)][spin_idx] \
		- Band[n][idx_change(x_id+1,y_id,z_id-1,num,shift)][spin_idx] - Band[n][idx_change(x_id-1,y_id,z_id+1,num,shift)][spin_idx])/diff_K**2.0/4.0
	dydz = (Band[n][idx_change(x_id,y_id+1,z_id+1,num,shift)][spin_idx] + Band[n][idx_change(x_id,y_id-1,z_id-1,num,shift)][spin_idx] \
		- Band[n][idx_change(x_id,y_id+1,z_id-1,num,shift)][spin_idx] - Band[n][idx_change(x_id,y_id-1,z_id+1,num,shift)][spin_idx])/diff_K**2.0/4.0
	return [dxx,dyy,dzz,dxdy,dxdz,dydz]

def idx_change(x_id,y_id,z_id,num,shift): # idx change for Band
	new_idx = shift+x_id*num[2]*num[1]+y_id*num[2]+z_id
	return new_idx

def Fermi_dist(x,x0,carrier_type,Temp):
	import scipy.constants as sc
	scp = sc.physical_constants
	if carrier_type == 'hole':
		return 1.0-1.0/(np.exp((x-x0)/(scp['Boltzmann constant in eV/K'][0]*Temp))+1.0)
	else:
		if (x-x0) > 10.0: # if diff is too large, exp term become infinity.
			return 0.0
		else:
			return 1.0/(np.exp((x-x0)/(scp['Boltzmann constant in eV/K'][0]*Temp))+1.0)

def PROCAR_to_array(procar_file):
	# read head line
	lines = subprocess.check_output(['grep','# of',procar_file]).splitlines()
	spin = len(lines)
	line = lines[0].replace(':',': ').split() # due to printing format
	nband = int(line[7])
	nkpt = int(line[3])
	# read energy line
	lines = subprocess.check_output(['grep','# energy',procar_file]).splitlines()
	Band = []
	for n in range(nband):
		Band.append([])
		for k in range(nkpt):
			Band[n].append([])
			for i in range(spin):
				idx = i*nband*nkpt + k*nband + n
				Band[n][k].append(float(lines[idx].split()[4]))
	return Band

def KPOINTS_to_array(kpoints_file):
	with open(kpoints_file,'r') as inp:
		kpts = inp.readlines()
	nkpt = int(kpts[1])
	KPT = []	
	for i in range(nkpt):
		KPT.append([float(x) for x in kpts[i+3].split()[0:3]])

	return KPT

def check_vasp_done(target):
	# it can be used when KPOINTS is written by cart scale.
	if not os.path.isfile(target+'/OUTCAR'):
		return 0
	outcar_last_line = subprocess.check_output(['tail','-1',target+'/OUTCAR'])
	if not 'Voluntary' in outcar_last_line:
		return 0
	else:
		kpt_for_kpoints = subprocess.check_output(['head','-8',target+'/KPOINTS']).splitlines()[3:8]
		kpt_for_outcar =  subprocess.check_output(['grep','-A5','k-points in units of 2pi',target+'/OUTCAR']).splitlines()[1:6]
		for i in range(5):
			for j in range(3):
				if not round(float(kpt_for_kpoints[i].split()[j]),7) == round(float(kpt_for_outcar[i].split()[j]),7):
					return 2
		return 1

			
############################################################################################
def pocket_old(target,carrier_type,E_width,kspacing):
	spin = subprocess.check_output(['grep','ISPIN',target+'/OUTCAR']).split()[2]
	ncl = subprocess.check_output(['grep','NONCOL',target+'/OUTCAR']).split()[2]
	[KPT,Band,nelect] = EIGEN_to_array(target+'/EIGENVAL',spin)
	with open(target+'/Band_gap.log','r') as inp:
		lines = inp.readlines()
		VBM_E = float(lines[2].split()[5])
		VBM_n = int(lines[5].split()[1])
		VBM_s = int(lines[5].split()[3])

		CBM_E = float(lines[3].split()[5])
		CBM_n = int(lines[6].split()[1])
		CBM_s = int(lines[6].split()[3])

	targ_n = [0,0]
	if carrier_type == 'hole':
		sign = 1
		targ_E = VBM_E
		targ_n[VBM_s-1] = VBM_n-1
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
		if not CBM_s == VBM_s:
			targ_n[VBM_s-1] = VBM_n
		elif spin == '2' or not nelect%2 == 0 or ncl == 'T': # spin exists.
			for i in range(len(Band)):
				if Band[i][0][2-CBM_s] > targ_E:
					targ_n[2-CBM_s] = i
					break
	KPT_num = [[KPT[0],1,[0]]] # kpt, # of overlap kpt, index of overlap kpt
	for k in range(1,len(KPT)):
		swit = 0
		for j in range(len(KPT_num)):
			if KPT[k] == KPT_num[j][0]:
				KPT_num[j][1] = KPT_num[j][1] + 1
				KPT_num[j][2].append(k)
				swit = 1
				break
		if swit == 0:
			KPT_num.append([KPT[k],1,[k]])

	high_sym = []
	band_kpt = []
	for k in range(len(KPT_num)):
		if KPT_num[k][1] > 1:
			high_sym.append(KPT_num[k])
		else:
			band_kpt.append(KPT_num[k][2][0])

	pocket = []
	# kpt is high symmetry points
	for i in range(int(spin)):
		for k in range(len(high_sym)):
			if abs(Band[targ_n[i]][band_kpt[high_sym[k][2][0]]][i]-targ_E) < E_width:
				poc_on = 1 # pocket : 1, not : 0
				for idx_k in high_sym[k][2]:
					if idx_k == 0: # initial kpt, compare to next point
						if (Band[targ_n[i]][idx_k][i] - Band[targ_n[i]][idx_k+1][i]) * sign < 0:
							poc_on = 0
							break
					elif idx_k == len(KPT)-1: # last kpt, compare to previous point
						if (Band[targ_n[i]][idx_k][i] - Band[targ_n[i]][idx_k-1][i]) * sign < 0:
							poc_on = 0
							break
					else:
						if dist_point(KPT[idx_k],KPT[idx_k-1]) < kspacing*1.5: # check continuity with previous point
							if (Band[targ_n[i]][idx_k][i] - Band[targ_n[i]][idx_k-1][i]) * sign < 0:
								poc_on = 0
								break
						if dist_point(KPT[idx_k],KPT[idx_k+1]) < kspacing*1.5: # check continuity with next point
							if (Band[targ_n[i]][idx_k][i] - Band[targ_n[i]][idx_k+1][i]) * sign < 0:
								poc_on = 0
								break
				if poc_on == 1:
					pocket.append([KPT[high_sym[k][2][0]],targ_n[i]+1,i+1,Band[targ_n[i]][high_sym[k][2][0]][i]])

	# kpt is not high symmetry points
	for i in range(int(spin)):
		for k in range(len(band_kpt)):
			if abs(Band[targ_n[i]][band_kpt[k]][i]-targ_E) < E_width:
				if (Band[targ_n[i]][band_kpt[k]][i]-Band[targ_n[i]][band_kpt[k]-1][i])*(Band[targ_n[i]][band_kpt[k]+1][i]-Band[targ_n[i]][band_kpt[k]][i]) < 0:
					pocket.append([KPT[band_kpt[k]],targ_n[i]+1,i+1,Band[targ_n[i]][band_kpt[k]][i]])

	with open(target+'/pocket.log','w') as out:
		for i in range(len(pocket)) :
			out.write(' '.join([str(x) for x in pocket[i][0]])+' :\t'+' '.join([str(x) for x in pocket[i][1:4]])+'\n')


