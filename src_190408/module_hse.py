###########################################
### Date: 2018-12-05			###
### yybbyb@snu.ac.kr			###
###########################################
def make_kpts_for_oneshot(kpt_rlx_file,kpt_band_file,target):
	with open(kpt_rlx_file,'r') as kpt_rlx_inp:
		kpt_rlx = kpt_rlx_inp.readlines()
	with open(kpt_band_file,'r') as kpt_band_inp:
		kpt_band = kpt_band_inp.readlines()
	with open(target+'/KPOINTS','w') as out:
		out.write('kpt for hse oneshot\n')
		nkpt = int(kpt_rlx[1].split()[0])
		out.write('             '+str(nkpt+2)+'\n')
		out.write('Reciprocal\n')
		for i in range(nkpt):
			out.write(kpt_rlx[3+i])
		out.write(kpt_band[3])
		out.write(kpt_band[4])

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

