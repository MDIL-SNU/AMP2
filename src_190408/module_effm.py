###########################################
### Date: 2019-01-22			###
### yybbyb@snu.ac.kr			###
###########################################
import subprocess,os
from module_band import EIGEN_to_array
from module_vector import *
from module_vasprun import *
from module_amp2_input import *
import numpy as np

def calc_max_E_diff(Fermi_cut,temp):
	import scipy.constants as sc
	scp = sc.physical_constants
	max_E_diff = -scp['Boltzmann constant in eV/K'][0]*temp*np.log(1./Fermi_cut-1.)
	return max_E_diff

def pocket(target,band_path,carrier_type,E_width,search_space):
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
#	KPTD.append([0.5, 0.5, 0.5])
#	KPTD.append([0.5, 0.5, -0.5])
#	KPTD.append([0.5, -0.5, 0.5])
#	KPTD.append([-0.5, 0.5, 0.5])

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

	with open(target+'/KPOINTS_searching','w') as kpt_out:
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
	kp_max_tot = []
	kp_min_tot = []
	final_pocket_kpt = []
	final_pocket_inform = []
	for k in range(number_of_pocket):
		sw_for_overlap = 0
		if not k == 0:
			for poc_id in range(len(kp_max_tot)):
				sw_for_overlap_single = 0
				for k_id in range(3):
					if pocket_kpt[k][k_id] > kp_max_tot[poc_id][k_id] or pocket_kpt[k][k_id] < kp_min_tot[poc_id][k_id]:
						sw_for_overlap_single = 1
						break
				if sw_for_overlap_single == 0:
					sw_for_overlap = 1
			
		if sw_for_overlap == 0:
			kp_max = [x for x in pocket_kpt[k]]
			kp_min = [x for x in pocket_kpt[k]]
			for i in range(len(num_kp_for_test)):
				## extreme to -direction
				for j in range((num_kp_for_test[i]-1)/2):
					kpj_idx = shift-j+(num_kp_for_test[i]-1)/2-1 # (shift + exterme point index - j), Ex> if num kp = 15; 6,5,4,...,0
					if abs(extreme_E - Band[pocket_inform[k][0]][kpj_idx][pocket_inform[k][1]]) < max_E_diff:
						for k_id in range(3):
							if KPT[kpj_idx][k_id] > KPT[shift+(num_kp_for_test[i]-1)/2][k_id]:
								if kp_max[k_id] < KPT[kpj_idx][k_id]:
									kp_max[k_id] = KPT[kpj_idx][k_id]
							else:
								if kp_min[k_id] > KPT[kpj_idx][k_id]:
									kp_min[k_id] = KPT[kpj_idx][k_id]
					else:
						break
				## extreme to +direction
				for j in range((num_kp_for_test[i]-1)/2):
					kpj_idx = shift+j+(num_kp_for_test[i]-1)/2+1 # (shift + exterme point index - j), Ex> if num kp = 15; 8,9,10,11
					if abs(extreme_E - Band[pocket_inform[k][0]][kpj_idx][pocket_inform[k][1]]) < max_E_diff:
						for k_id in range(3):
							if KPT[kpj_idx][k_id] > KPT[shift+(num_kp_for_test[i]-1)/2][k_id]:
								if kp_max[k_id] < KPT[kpj_idx][k_id]:
									kp_max[k_id] = KPT[kpj_idx][k_id]
							else:
								if kp_min[k_id] > KPT[kpj_idx][k_id]:
									kp_min[k_id] = KPT[kpj_idx][k_id]
					else:
						break

				shift = shift + num_kp_for_test[i]

			kp_max_tot.append(kp_max)
			kp_min_tot.append(kp_min)

			num_kp_for_effm.append([]) # num_kp_for_effm[pocket_idx][kpt_idx][min or max]
			final_pocket_kpt.append(pocket_kpt[k])
			final_pocket_inform.append(pocket_inform[k])
			for k_id in range(3):
				num_kp_for_effm[-1].append([int(np.ceil((pocket_kpt[k][k_id]-kp_min[k_id])/grid_size))+1,int(np.ceil((kp_max[k_id]-pocket_kpt[k][k_id])/grid_size))+1])
			for x_id in range(num_kp_for_effm[-1][0][0]+num_kp_for_effm[-1][0][1]+1):
				for y_id in range(num_kp_for_effm[-1][1][0]+num_kp_for_effm[-1][1][1]+1):
					for z_id in range(num_kp_for_effm[-1][2][0]+num_kp_for_effm[-1][2][1]+1):
						id_list = [x_id,y_id,z_id]
						KPT_for_effm.append([pocket_kpt[k][x]+float(id_list[x]-num_kp_for_effm[-1][x][0])*grid_size for x in range(3)])

	number_of_pocket = len(num_kp_for_effm)

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
			inp_effm.write('\t'.join([str(final_pocket_kpt[k][x]) for x in range(3)])+'\n')
			for i in range(3):
				inp_effm.write('\t'.join([str(num_kp_for_effm[k][i][x]) for x in range(2)])+'\n')
			inp_effm.write('\t'.join([str(final_pocket_inform[k][x]) for x in range(2)])+'\n')
		inp_effm.write(str(extreme_E))

def kpt_frac_in_cell(x_id_mod,y_id_mod,z_id_mod,grid_size,axis):
	new_kpt_cart = [float(x_id_mod)*grid_size,float(y_id_mod)*grid_size,float(z_id_mod)*grid_size]
	new_kpt_frac = cart_to_dir(new_kpt_cart,axis)
	sw = 0
	for i in range(3):
		if new_kpt_frac[i] < -0.5 or new_kpt_frac[i] >= 0.5:
			sw = 1
	if sw == 1:
		return 0 # out
	else:
		return 1 # in,do calculation

def calc_effm(target,carrier_type,Temp,oper):
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
	for k in range(number_of_pocket): # loop for pockets
		if carrier_type == 'hole':
			min_band_idx = 0
			max_band_idx = pocket_inform[k][0]+1
		else:
			min_band_idx = pocket_inform[k][0]
			max_band_idx = len(Band)

		# check pocket overlap
		sw_overlap = 0
		if not k == 0: # if the point is not first
			sym_kpts = symmetrical_kpts(frac_pocket_kpt[k],oper) # find duplicated kpts applying symmetry operator.
			for new_kpt in sym_kpts:
				for i in unique_kpt_idx:
					if valid_pocket(dir_to_cart(new_kpt,rec_axis),pocket_kpt[i],num_kp_sep[i],grid_size,rec_axis) == 0:
						sw_overlap = 1
						break
				if sw_overlap == 1:
					break
		# if sw_overlap is 0, the kpoint is unduplicated point.
		if sw_overlap == 0:
			unique_kpt_idx.append(k)
			sym_kpts = symmetrical_kpts(frac_pocket_kpt[k],oper)
			num_overlap_sym_kpt = 0 # number of pockets that are placed in the mesh grid from another pocket.
			for new_kpt in sym_kpts:
				if valid_pocket(dir_to_cart(new_kpt,rec_axis),pocket_kpt[k],num_kp_sep[k],grid_size,rec_axis) == 0:
					num_overlap_sym_kpt = num_overlap_sym_kpt +1
			weight = float(len(sym_kpts))/float(num_overlap_sym_kpt)

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
	# symmetry operation
	deriv_tensor_oper = np.zeros_like(deriv_tensor,dtype=float)
	for oper_line in oper:
		oper_xyz = np.matmul(np.linalg.inv(rec_axis),np.matmul(oper_line,rec_axis))
		deriv_tensor_oper = deriv_tensor_oper+np.matmul(np.transpose(oper_xyz),np.matmul(deriv_tensor,oper_xyz))
	deriv_tensor = deriv_tensor_oper/float(len(oper))

	deriv_mat = np.linalg.inv(deriv_tensor)
	deriv_dia = np.linalg.eigvals(deriv_mat)

	effm_dia = []
	effm = []
	for i in range(3):
		effm_dia.append(scp['natural unit of action in eV s'][0]**2.0*scp['atomic unit of charge'][0]*1.0e20*(deriv_dia[i])/sc.m_e)
		effm.append([scp['natural unit of action in eV s'][0]**2.0*scp['atomic unit of charge'][0]*1.0e20*(deriv_mat[i][x])/sc.m_e for x in range(3)])

	return [effm_dia,effm]

def write_effm(effm_dia,effm,target,carrier_type):
	with open(target+'/effective_mass.log','w') as out:
		out.write(carrier_type+'\n')
		for i in range(3):
			out.write(' '.join(['{:10.3f}'.format(effm[i][x]) for x in range(3)])+'\n')
#			out.write('\t'.join([str(effm[i][x]) for x in range(3)])+'\n')
		out.write('Diagonalized effective mass: '+' '.join(['{:10.3f}'.format(effm_dia[x]) for x in range(3)]))
#		out.write('\t'.join([str(effm_dia[x]) for x in range(3)]))

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
	lines = pygrep('# of',procar_file,0,0).splitlines()
#	lines = subprocess.check_output(['grep','# of',procar_file]).splitlines()
	spin = len(lines)
	line = lines[0].replace(':',': ').split() # due to printing format
	nband = int(line[7])
	nkpt = int(line[3])
	# read energy line
	lines = pygrep('# energy',procar_file,0,0).splitlines()
#	lines = subprocess.check_output(['grep','# energy',procar_file]).splitlines()
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
	outcar_last_line = pytail(target+'/OUTCAR')
#	outcar_last_line = subprocess.check_output(['tail','-1',target+'/OUTCAR'])
	if not 'Voluntary' in outcar_last_line:
		return 0
	else:
		kpt_for_kpoints = pyhead(target+'/KPOINTS',8).splitlines()[3:8]
#		kpt_for_kpoints = subprocess.check_output(['head','-8',target+'/KPOINTS']).splitlines()[3:8]
		kpt_for_outcar =  pygrep('k-points in units of 2pi',target+'/OUTCAR',0,5).splitlines()[1:6]
#		kpt_for_outcar =  subprocess.check_output(['grep','-A5','k-points in units of 2pi',target+'/OUTCAR']).splitlines()[1:6]
		for i in range(5):
			for j in range(3):
				if not round(float(kpt_for_kpoints[i].split()[j]),7) == round(float(kpt_for_outcar[i].split()[j]),7):
					return 2
		return 1

def read_operation(poscar):
	import spglib
	[axis,atom_pos] = read_poscar(poscar)

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

	rot_oper = np.ndarray.tolist(spglib.get_symmetry(Cell,symprec=1e-5)['rotations'])
	trans_oper = np.ndarray.tolist(spglib.get_symmetry(Cell,symprec=1e-5)['translations'])
	oper = []
	for i in range(len(rot_oper)):
		if sum([abs(x-0.5) for x in trans_oper[i]]) > 1.5-0.0000001: # considering the PBC (1.0 = 0.0)
			oper.append(rot_oper[i])
	# symmetry opperator according to a,b,c axis (not x,y,z axis)
	return oper

def valid_pocket(pocket_kpt,ref_kpt,num_kp_sep,grid_size,axis):
	frac_poc_tmp = cart_to_dir([pocket_kpt[x]-ref_kpt[x] for x in range(3)],axis)
	frac_poc = [x-np.floor(x+0.5) for x in frac_poc_tmp]
	cart_poc = dir_to_cart(frac_poc,axis)
	sw_for_overlap = 0 # 0 indicate overlapped k-points.
	for i in range(3):
		if cart_poc[i] > grid_size*float(num_kp_sep[i][1]) or cart_poc[i] < -grid_size*float(num_kp_sep[i][0]) :
			sw_for_overlap = 1
			break
	if sw_for_overlap == 1:
		return 1 # out,do calculation
	else :
		return 0 # in

def symmetrical_kpts(pocket_kpt,oper):
	poc = []
	poc.append([x for x in pocket_kpt])
	prev_poc_len = len(poc)
	while 1:
		prev_poc_len = len(poc)
		poc = dup_pos(poc,oper)
		if len(poc) == prev_poc_len:
			break
		else:
			prev_poc_len = len(poc)
	return poc

def vec_mod(vec):
	new_vec = []
	for xx in vec:
		new_vec.append(round(xx-np.floor(xx+0.5),7))
	return new_vec

def dup_pos(pos,oper):
	irre_pos = []
	made_pos = []
	idx = 0
	dic = {}

	for i in range(len(pos)):
		for j in range(len(oper)):
			made_pos.append(vec_mod(apply_operator(pos[i],oper[j])))
	for i in range(len(made_pos)):
		if len(irre_pos) == 0:
			irre_pos.append(made_pos[i])
		else:
			sw = 0
			for j in range(len(irre_pos)):
				if dist_point(made_pos[i],irre_pos[j]) < 0.000001:
					sw = 1
					break
			if sw == 0:
				irre_pos.append(made_pos[i])
	return irre_pos

def apply_operator(pos,operator):
	oper_pos = []
	for i in range(3):
		oper_pos.append(pos[0]*operator[0][i]+pos[1]*operator[1][i]+pos[2]*operator[2][i])
	return oper_pos

############################################################################################
def pocket_old(target,carrier_type,E_width,kspacing):
	spin = pygrep('ISPIN',target+'/OUTCAR',0,0).split()[2]
#	spin = subprocess.check_output(['grep','ISPIN',target+'/OUTCAR']).split()[2]
	ncl = pygrep('NONCOL',target+'/OUTCAR',0,0).split()[2]
#	ncl = subprocess.check_output(['grep','NONCOL',target+'/OUTCAR']).split()[2]
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

def mod_kpt_from_rec_to_car(kpt_file,pos_file):
	with open(kpt_file,'r') as kpt_inp:
		lines = kpt_inp.readlines()
