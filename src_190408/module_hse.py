###########################################
### Date: 2018-12-05			###
### yybbyb@snu.ac.kr			###
###########################################
from module_band import *
from module_vector import *
from module_vasprun import *

def make_kpts_for_oneshot(kpt_rlx_file,kpt_band_file,target):
	with open(kpt_rlx_file,'r') as kpt_rlx_inp:
		lines = kpt_rlx_inp.readlines()
	nkpt = int(lines[1].split()[0])
	kpt_total = []
	for line in lines[3:3+nkpt]:
		kpt_total.append([float(x) for x in line.split()[0:4]])

	with open(kpt_band_file,'r') as kpt_band_inp:
		lines = kpt_band_inp.readlines()
	kpt_band = []
	for line in lines[3:]:
		if line.split() > 3:
			kpt_band.append([float(x) for x in line.split()[0:4]])

	for i in range(len(kpt_band)):
		overlap_sw = 0
		for j in range(len(kpt_total)):
			if abs(kpt_band[i][0]-kpt_total[j][0]) + abs(kpt_band[i][1]-kpt_total[j][1]) + abs(kpt_band[i][2]-kpt_total[j][2]) < 0.0000001:
				overlap_sw = 1
		if overlap_sw == 0:
			kpt_total.append(kpt_band[i])
	with open(target+'/KPOINTS','w') as out:
		out.write('kpt for hse oneshot\n')
		out.write('             '+str(len(kpt_total))+'\n')
		out.write('Reciprocal\n')
		for i in range(len(kpt_total)):
			out.write('    '+'    '.join([str(x) for x in kpt_total[i]])+'\n')

def gap_estimation_hse(target,fermi,spin,ncl,KPT,Band,nelect):
	VBM = []
	CBM = []
	nVBM = []
	nCBM = []
	VBM_k = []
	CBM_k = []
	for i in range(int(spin)):
		VBM.append(-10.0**6.0)
		CBM.append(10.0**6.0)
	metal = 0
	# VBM & CBM search for spin-polarized calculation
	nVB = -1
	if spin == '1' and nelect%2 == 0 and ncl != 'T' :
		nVBM.append(nelect/2)
		nCBM.append(nelect/2+1)
		single_band = [Band[nVBM[0]-1][x][0] for x in range(len(KPT))]
		VBM[0] = max(single_band)
		VBM_k.append(KPT[single_band.index(VBM[0])])
		single_band = [Band[nCBM[0]-1][x][0] for x in range(len(KPT))]
		CBM[0] = min(single_band)
		CBM_k.append(KPT[single_band.index(CBM[0])])
		if CBM[0]-VBM[0] < 0.01:
			metal = 1
	else :
		for i in range(int(spin)) :
			for n in range(len(Band)) :
				if n == 0 :
					nVBM.append(n)
					nCBM.append(n)
					VBM_k.append(KPT[n][0][i])
					CBM_k.append(KPT[n][0][i])

				single_band = [Band[n][x][i] for x in range(len(KPT))]
				band_max = max(single_band)
				band_min = min(single_band)

				# nVBM and vCBM are determined at first kpt.
				# Fermi energy could be lower than VBM due to occupancy.
				if single_band[0] < fermi :
					if band_max >= VBM[i] :
						VBM[i] = band_max
						VBM_k[i] = KPT[single_band.index(band_max)]
						nVBM[i] = n+1
				elif single_band[1] > fermi :
					if band_min < CBM[i] :
						CBM[i] = band_min
						CBM_k[i] = KPT[single_band.index(band_min)]
						nCBM[i] = n+1
				else :
					metal = 1
					VBM[i] = band_max
					VBM_k[i] = KPT[single_band.index(band_max)]
					nVBM[i] = n+1
					CBM[i] = band_min
					CBM_k[i] = KPT[single_band.index(band_min)]
					nCBM[i] = n+1
					break

	if metal == 1:
		VB_max = []
		CB_min = []
		spin_index = []
		for i in range(len(KPT)):
			VB_max.append(Band[nVBM[0]-1][i][0])
			CB_min.append(Band[nVBM[0]][i][0])
#			spin_index.append([0,0])
			if not spin == '1':
				if VB_max[i] < Band[nVBM[1]-1][i][1]:
					VB_max[i] = Band[nVBM[1]-1][i][1]
#					spin_index[i][0] = 1
				if CB_min[i] > Band[nVBM[1]][i][1]:
					CB_min[i] = Band[nVBM[1]][i][1]
#					spin_index[i][1] = 1
			if CB_min[i] - VB_max[i] < 0.01:
				metal == 2
				break

	if metal == 0:
		total_VBM = max(VBM)
		total_CBM = min(CBM)
		total_VBM_spin = VBM.index(total_VBM)
		total_CBM_spin = CBM.index(total_CBM)
		total_VBM_index = nVBM[total_VBM_spin]
		total_CBM_index = nCBM[total_CBM_spin]
		total_VBM_kpt = VBM_k[total_VBM_spin]
		total_CBM_kpt = CBM_k[total_CBM_spin]
		gap = total_CBM-total_VBM

	elif metal == 1:
		total_VBM_kpt = KPT[VB_max.index(max(VB_max))]
		total_CBM_kpt = KPT[CB_min.index(min(CB_min))]
		gap = 0
	else:
		gap = 0
	gap_log = open(target+'/Band_gap.log', 'w')
	gap_simple = open(target+'/Band_gap', 'w')	# For auto_bin
	if metal != 0:
		gap_final = 'metal'
		gap_log.write('This system is metallic.\n')
		gap_simple.write('     metal\n')
		if metal == 2 :
			gap_log.write('! Additional search is required for HSE0\n')
			gap_simple.write('! Additional search is required for HSE0\n')
		else :
			gap_log.write("Use the 'KPT' file for HSE0 calculation\n")
	else :
		gap_final = str(gap)
		direct = 0
		for i in range(3) :
			total_VBM_kpt[i] = str(float(total_VBM_kpt[i]))
			total_CBM_kpt[i] = str(float(total_CBM_kpt[i]))
			if not total_VBM_kpt[i] == total_CBM_kpt[i]:
				direct = 1
		if direct == 0 :
			gap_log.write('Band gap: '+str(gap)+' eV (Direct)\n\n')
		else :
			gap_log.write('Band gap: '+str(gap)+' eV (Indirect)\n\n')
		gap_simple.write('     '+str(gap))
		gap_log.write('VBM: '+'  '.join(total_VBM_kpt)+'   : '+str(total_VBM)+' eV\n')
		gap_log.write('CBM: '+'  '.join(total_CBM_kpt)+'   : '+str(total_CBM)+' eV\n')
		gap_log.write('\nnVBM: '+str(total_VBM_index)+'  spin: '+str(total_VBM_spin+1)+'\n')
		gap_log.write('nCBM: '+str(total_CBM_index)+'  spin: '+str(total_CBM_spin+1)+'\n')
	gap_log.close()
	gap_simple.close()
	return gap_final

def DOS_ratio_fermi_to_vb(dos_file,fermi_width,vb_range):
	DOS = []
	with open(dos_file,'r') as dos_inp:
		for i in range(5):
			dos_inp.readline()
		line = dos_inp.readline().split()
		npoint = int(line[2])
		fermi = float(line[3])
		for i in range(npoint):
			line = dos_inp.readline().split()
			if len(line) == 3:
				DOS.append([float(line[0]),[float(line[1])]])
			else:
				DOS.append([float(line[0]),[float(line[1]),float(line[2])]])
	vb = [0,0]
	vb_ind = [0,0]
	mid = [0,0]
	mid_ind = [0,0]

	if vb_range[0] < 0:
		vb_range[0] = -vb_range[0]
		vb_range[1] = -vb_range[1]
	vb_range = sorted(vb_range)

	for i in range(npoint):
		if DOS[i][0] > fermi - vb_range[1] and DOS[i][0] < fermi - vb_range[0]:
			if vb_ind[0] == 0 :
				vb_ind[0] = i
			vb_ind[1] = i
			for j in range(len(DOS[0][1])):
				vb[j] = vb[j] + DOS[i][1][j]

		elif DOS[i][0] > fermi - fermi_width and DOS[i][0] < fermi + fermi_width:
			if mid_ind[0] == 0 :
				mid_ind[0] = i
			mid_ind[1] = i
			for j in range(len(DOS[0][1])):
				mid[j] = mid[j] + DOS[i][1][j]
 
	avg_vb = (vb[0]+vb[len(DOS[0][1])-1])/2./(DOS[vb_ind[1]][0]-DOS[vb_ind[0]][0])
	avg_mid = (mid[0]+mid[len(DOS[0][1])-1])/2./(DOS[mid_ind[1]][0]-DOS[mid_ind[0]][0])
	for i in range(1,len(DOS[0][1])-1):
		avg_vb = avg_vb + vb[i]/(DOS[vb_ind[1]][0]-DOS[vb_ind[0]][0])
		avg_mid = avg_mid + mid[i]/(DOS[mid_ind[1]][0]-DOS[mid_ind[0]][0])
	
	DF_DVB = avg_mid/(avg_vb+0.00001) # avoid to divide zero

	return DF_DVB

def find_extreme_kpt_for_hse(dir_band,e_width,search_space):
	fermi = float(subprocess.check_output(['head',dir_band+'/DOSCAR','-n','6']).splitlines()[-1].split()[3])
	spin = subprocess.check_output(['grep','ISPIN',dir_band+'/OUTCAR']).split()[2]
	ncl = subprocess.check_output(['grep','NONCOL',dir_band+'/OUTCAR']).split()[2]
	[KPT,Band,nelect] = EIGEN_to_array(dir_band+'/EIGENVAL',spin)

	num_new_kpt = 0
	new_kpt = []

	extreme_log = []

	axis = poscar_to_axis(dir_band+'/POSCAR')
	rec_lat = reciprocal_lattice(axis)

	for n in range(len(Band)):
		for i in range(len(Band[0][0])):
			for k in range(len(KPT)):
				if Band[n][k][i] > fermi-e_width and Band[n][k][i] < fermi+e_width :
					sign = 0
					ex_sw = 0
					for k2 in range(len(KPT)):
						if not k == k2 and dist_vec(KPT[k],KPT[k2],rec_lat) < search_space:
							if sign == 0:
								if Band[n][k][i] > Band[n][k2][i]: # top
									sign = 1 
								else: # bottom
									sign = -1
						
							if sign*(Band[n][k][i]-Band[n][k2][i]) < 0 : # not extreme point
								ex_sw = 1
								break
					if ex_sw == 0:
						extreme_log.append([n,k,i])

	reduced_log = []
	fin_kpt = []
	for i in range(len(extreme_log)):
		if not extreme_log[i][1] in [reduced_log[x][1] for x in range(len(reduced_log))]:
			reduced_log.append(extreme_log[i])
			overlap_sw = 0
			kpt_tmp = KPT[reduced_log[-1][1]]
			for j in range(len(fin_kpt)):
				print i,j,dist_vec(kpt_tmp,fin_kpt[j],rec_lat)
				if dist_vec(kpt_tmp,fin_kpt[j],rec_lat) < 0.0000001:
					overlap_sw = 1
					break
			if overlap_sw == 0:
				fin_kpt.append(KPT[reduced_log[-1][1]])
	with open(dir_band+'/KPT', 'w') as kpt_out:	# File for VBM & CBM k-point position
		kpt_out.write("Example file\n               0\nReciprocal\n")
		for i in range(len(fin_kpt)):
			kpt_out.write('    '+'    '.join(fin_kpt[i])+'\t0\n')
