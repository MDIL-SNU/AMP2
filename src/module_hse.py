###########################################
### Date: 2018-12-05			###
### yybbyb@snu.ac.kr			###
###########################################
from module_band import *
from module_vector import *
from module_vasprun import *
import numpy as np
import sys,binascii,os,subprocess, glob

def make_kpts_for_hse(kpt_rlx_file,kpt_band_file,target,calc_type):
	# calc_type is band or oneshot
	with open(kpt_rlx_file,'r') as kpt_rlx_inp:
		lines = kpt_rlx_inp.readlines()
	nkpt = int(lines[1].split()[0])
	kpt_type = lines[2]
	kpt_total = []
	for line in lines[3:3+nkpt]:
		kpt_total.append([float(x) for x in line.split()[0:4]])

	with open(kpt_band_file,'r') as kpt_band_inp:
		lines = kpt_band_inp.readlines()
	kpt_band = []
	for line in lines[3:]:
		if line.split() > 3:
			kpt_band.append([float(x) for x in line.split()[0:3]]+[0.0])

	for i in range(len(kpt_band)):
		if calc_type == 'band':
			kpt_total.append(kpt_band[i])
		else:
			overlap_sw = 0
			for j in range(len(kpt_total)):
				if abs(kpt_band[i][0]-kpt_total[j][0]) + abs(kpt_band[i][1]-kpt_total[j][1]) + abs(kpt_band[i][2]-kpt_total[j][2]) < 0.0000001:
					overlap_sw = 1
			if overlap_sw == 0:
				kpt_total.append(kpt_band[i])

	with open(target+'/KPOINTS','w') as out:
		out.write('kpt for hse '+calc_type+'\n')
		out.write('             '+str(len(kpt_total))+'\n')
		out.write(kpt_type)
		for i in range(len(kpt_total)):
			out.write('    '+'    '.join([str(x) for x in kpt_total[i]])+'\n')

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
	spin = pygrep('ISPIN',dir_band+'/OUTCAR',0,0).split()[2]
#	spin = subprocess.check_output(['grep','ISPIN',dir_band+'/OUTCAR']).split()[2]
	ncl = pygrep('NONCOL',dir_band+'/OUTCAR',0,0).split()[2]
#	ncl = subprocess.check_output(['grep','NONCOL',dir_band+'/OUTCAR']).split()[2]
	[KPT,Band,nelect] = EIGEN_to_array(dir_band+'/EIGENVAL',spin)
	fermi = get_fermi_level(Band,nelect,ncl)

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
				if dist_vec(kpt_tmp,fin_kpt[j],rec_lat) < 0.0000001:
					overlap_sw = 1
					break
			if overlap_sw == 0:
				fin_kpt.append(KPT[reduced_log[-1][1]])
	with open(dir_band+'/KPT', 'w') as kpt_out:	# File for VBM & CBM k-point position
		kpt_out.write("Example file\n               0\nReciprocal\n")
		for i in range(len(fin_kpt)):
			kpt_out.write('    '+'    '.join(fin_kpt[i])+'\t0\n')

def calc_alpha_auto(diel_path):
	# alpha is defined as 1/e_inf
#	import numpy as np
	with open(diel_path,'r') as f:
		lines = f.readlines()[1:4]
	diel_e  = []
	for line in lines:
		diel_e.append([float(x) for x in line.split()])
	avg_diel = (diel_e[0][0]+diel_e[1][1]+diel_e[2][2])/3.0
	return 1.0/avg_diel

def read_wavecar(wave_file,min_band,max_band):
	AU2Ang = 0.529177249
	RY2Ev = 13.605826
	with open(wave_file,'rb') as f:
		f.seek(0)
		recl, nspin, rtag = np.array(np.fromfile(f,dtype=np.float, count=3),dtype = int)
		if rtag == 45200:
			wave_prec = np.complex64
		elif rtag == 45210:
			wave_prec = np.complex128

		f.seek(recl)
		dump = np.fromfile(f, dtype=np.float, count=12)
		nkpt  = int(dump[0])                     # No. of k-points
		nband = int(dump[1])                     # No. of bands
		encut  = dump[2]                          # Energy cutoff
		Acell  = dump[3:].reshape((3,3))          # real space supercell basis
		Omega  = np.linalg.det(Acell)       # real space supercell volume
		Bcell  = np.linalg.inv(Acell).T     # reciprocal space supercell volume

		# Minimum FFT grid size
		Anorm = np.linalg.norm(Acell, axis=1)
		CUTOF = np.ceil(np.sqrt(encut / RY2Ev) / (2.*np.pi / (Anorm / AU2Ang)))
		ngrid = np.array(2 * CUTOF + 1, dtype=int)

		num_of_plane_waves = np.zeros(nkpt, dtype = int)
		kvectors = np.zeros((nkpt,3), dtype = float)
		for k in range(nkpt):
			new_idx = wave_idx(nkpt,nband,0,k,0) - 1
			f.seek(new_idx * recl)
			tmp = np.fromfile(f,dtype = np.float, count = 4)
			num_of_plane_waves[k] = int(tmp[0])
			kvectors[k] = tmp[1:4]
		Unk = np.zeros((nspin,nkpt,nband,ngrid[0],ngrid[1],ngrid[2]),dtype = complex)
		for i in range(nspin):
			for k in range(nkpt):
				G = get_Gvector(kvectors[k],Bcell,ngrid,encut)
				for n in range(min_band,max_band):
					new_idx = wave_idx(nkpt,nband,i,k,n)
					f.seek(new_idx*recl)
					tmp = np.fromfile(f,dtype = wave_prec, count = num_of_plane_waves[k])
					tmp_cg = np.asarray(tmp,dtype = np.complex128)
					Unk[i,k,n,G[:,0],G[:,1],G[:,2]] = tmp_cg/np.linalg.norm(tmp_cg)

	return Unk

def get_Gvector(kvector,Bcell,ngrid,encut):
	AU2Ang = 0.529177249
	RY2Ev = 13.605826
	fx = [ii if ii < ngrid[0] / 2 + 1 else ii - ngrid[0] for ii in range(ngrid[0])]
	fy = [jj if jj < ngrid[1] / 2 + 1 else jj - ngrid[1] for jj in range(ngrid[1])]
	fz = [kk if kk < ngrid[2] / 2 + 1 else kk - ngrid[2] for kk in range(ngrid[2])]
	kgrid = np.array([(fx[ii], fy[jj], fz[kk]) for kk in range(ngrid[2]) for jj in range(ngrid[1]) for ii in range(ngrid[0])], dtype=float)
	kinetic_E = RY2Ev * AU2Ang ** 2.0 * np.linalg.norm(np.dot(kgrid + kvector[np.newaxis,:],2.0*np.pi*Bcell),axis = 1) ** 2.0
	Gvector = kgrid[np.where(kinetic_E < encut)[0]]
	return np.asarray(Gvector, dtype = int)

def wave_idx(nkpt,nband,spin_idx,kpt_idx,band_idx):
	new_idx = 2 + spin_idx*nkpt*(nband+1) + kpt_idx*(nband+1) + band_idx+1
	return new_idx

def read_kpt_file(kpt_file):
	kpt = []
	with open(kpt_file,'r') as f:
		lines = f.readlines()
	for line in lines[3:]:
		if len(line.split()) > 3:
			kpt.append([float(x) for x in line.split()[0:3]])
	return kpt

def read_xtic_file(xtic_file):
	kpt = []
	with open(xtic_file,'r') as f:
		lines = f.readlines()
	overlap_idx = []
	over_rm_idx = []
	xtic = []
	for i,line in enumerate(lines):
		if len(line.split()) == 4:
			xtic.append(float(line.split()[0]))
			kpt.append([float(x) for x in line.split()[1:4]])
			if sum([abs(kpt[0][x]-kpt[i][x]) for x in range(3)]) < 0.0000001:
				overlap_idx.append([i,0])
	band_line = [] 
	for i in range(len(overlap_idx)):
		check_dup = 0
		for j in range(len(band_line)):
			if overlap_idx[i][0] in band_line[j]:
				check_dup = 1
				over_rm_idx.append(overlap_idx[i])
				break
		if check_dup == 0:
			if not abs(xtic[overlap_idx[i][0]]-xtic[overlap_idx[i][0]+1]) < 0.0000001:
				for j in range(overlap_idx[i][0],len(xtic)-1):
					if abs(xtic[j]-xtic[j+1]) < 0.0000001:
						band_line.append([1,0]+range(overlap_idx[i][0]+1,j+1))
						break
					elif j == len(xtic)-2:
						band_line.append([1,0]+range(overlap_idx[i][0]+1,j+1))
						break
			if not abs(xtic[overlap_idx[i][0]]-xtic[overlap_idx[i][0]-1]) < 0.0000001:
				for j in range(overlap_idx[i][0],0,-1):
					if abs(xtic[j]-xtic[j-1]) < 0.0000001:
						band_line.append([1,0]+range(overlap_idx[i][0]-1,j-1,-1))
						break
	for item in over_rm_idx:
		overlap_idx.remove(item)
	for i in range(1,len(xtic)):
		if abs(xtic[i]-xtic[i-1]) < 0.0000001:
			# check that the point is alaeady included in the other line
			check_dup = 0
			for j in range(len(band_line)):
				if i in band_line[j] or i in [x[0] for x in overlap_idx]:
					check_dup = 1
					break
			# the point is not included in the other line
			if check_dup == 0:
				for k in range(len(band_line)):
					for j,idx in enumerate(band_line[k]):
						if sum([abs(kpt[i][x]-kpt[idx][x]) for x in range(3)]) < 0.0000001:
							overlap_idx.append([i,idx])
							if band_line[k][0] == j:
								address = [band_line[k][1],idx]
							else:
								address = [band_line[k][j-1],idx]
							break

				for j in range(i,len(xtic)-1):
					if abs(xtic[j]-xtic[j+1]) < 0.0000001:
						band_line.append(address+range(i+1,j+1))
						break
					elif j == len(xtic)-2:
						band_line.append(address+range(i+1,j+2))
						break
	return [band_line,overlap_idx]

def make_band_corrected_dat(xtic_file,Band,spin,target):
	xtic = []
	with open(xtic_file,'r') as inp:
		for line in inp:
			xtic.append(line.split()[0])
	with open(target+'/band_corrected.dat','w') as out:
		for k in range(len(xtic)):			# plot k-pts
			out.write(xtic[k])
			for i in range(int(spin)):		# spin
				for n in range(len(Band)):	# band idx
					out.write('\t'+str(Band[n][k][i]))
			out.write('\n')

def make_band_corrected_in(title,xlabel_file,fermi,gap,nband,spin,plot_range,target):
	with open(target+'/band_corrected.in','w') as out:
		out.write("set terminal pdfcairo enhanced color font 'Arial, 14' size 7.2,5.4\n")
#		out.write("set output '"+target+'/'+title+".eps'\n")
		out.write("set output '"+title+"_corrected.pdf'\n")
		title2 = '-'.join(title.split('_'))
		out.write("set title '"+title2+"' font 'Arial, 16'\n")
		if not spin == '2':
			out.write("set nokey\n")
		out.write('e = ' +fermi+'\n')
		out.write('set termoption dash\n')
		out.write('set yr[-'+str(plot_range[0])+':'+str(gap+plot_range[1])+']\n')
		with open(xlabel_file,'r') as inp:
			out.write(inp.read())
		out.write('f(x)=0\n')
		out.write('plot \\\n')
		if spin == '2':
			for i in range(nband-1):
				out.write("\t'band_corrected.dat' u 1:($"+str(i+2)+"-e) w l lt 1 lc rgb 'red' notitle,\\\n")
			out.write("\t'band_corrected.dat' u 1:($"+str(nband+1)+"-e) w l lt 1 lc rgb 'red' title 'up',\\\n")
			for i in range(nband-1):
				out.write("\t'band_corrected.dat' u 1:($"+str(nband+i+2)+"-e) w l lt 1 lc rgb 'blue' notitle,\\\n")
			out.write("\t'band_corrected.dat' u 1:($"+str(nband*2+1)+"-e) w l lt 1 lc rgb 'blue' title 'down',\\\n")
		else:
			for i in range(nband):
				out.write("\t'band_corrected.dat' u 1:($"+str(i+2)+"-e) w l lt 1 lc rgb 'black' notitle,\\\n")
		out.write("\tf(x) w l lt 0 notitle\n\n")
		out.write("set terminal pngcairo enhanced font 'Arial, 14' size 960,720\n")
		out.write("set output '"+title+"_corrected.png'\n")
		out.write("replot")

def plot_band_corrected_structure(spin,Band,fermi,xtic_file,xlabel_file,plot_range,target):
	# make band.dat
	make_band_corrected_dat(xtic_file,Band,spin,target)

	# make band.in
	nband = len(Band)
	if pyhead(target+'/Band_gap.log',1).split()[2] == 'is' :
#	if subprocess.check_output(['head','-n','1',target+'/Band_gap.log']).split()[2] == 'is' :
		gap = 0
		fermi = str(fermi)
	else :
		gap = round(float(subprocess.check_output(['head','-n','1',target+'/Band_gap.log']).split()[2]))
#		gap = round(float(subprocess.check_output(['head','-n','1',target+'/Band_gap.log']).split()[2]))
		fermi = pygrep('VBM',target+'/Band_gap.log',0,0).splitlines()[0].split()[-2]
#		fermi = subprocess.check_output(['grep','VBM',target+'/Band_gap.log']).splitlines()[0].split()[-2]

	title = target.split('/')[-2]
	make_band_corrected_in(title,xlabel_file,fermi,gap,nband,spin,plot_range,target)

def read_band_log(band_log):
	with open(band_log,'r') as f:
		lines = f.readlines()
	if 'metal' in lines[0]:
		sys.exit()
	else:
		gap = float(lines[0].split()[2])
		vbm_kp = [float(x) for x in lines[2].split()[1:4]]
		cbm_kp = [float(x) for x in lines[3].split()[1:4]]
		return [gap,vbm_kp,cbm_kp]

def make_xtic_hse(symk,order,rec,KSPACING,target):
	# make path from the combination of high symmetry kpts
	path = []
	for i in range(len(symk)):
		for j in range(i+1,len(symk)):
			path.append([999,i+1,j+1]) # 999 for sort [[index,kp_start,kp_end],...]
	# change path index for band structure
	for i in range(len(order)):
		path[order[i][0]-1][0] = abs(order[i][1])
		if order[i][1] < 0:
			path[order[i][0]-1][1],path[order[i][0]-1][2]=path[order[i][0]-1][2],path[order[i][0]-1][1]
	# sorting the path
	path=sorted(path)
	nkpt = []
	KPT = []
	section = []
	xtic = []
	for i in range(len(path)):
		kp_st = path[i][1]-1
		kp_end = path[i][2]-1
		dist_rec = [symk[kp_end][x]-symk[kp_st][x] for x in range(3)] # distance vector
		dist_val = dist_vec(symk[kp_end],symk[kp_st],rec) # distance value
		nsplit = round(dist_val/KSPACING,0)
		if nsplit == 0: # At the minimum, one kpt is included.
			nsplit = 1
		begin = 0
		if not i == 0:
			if kp_st == path[i-1][2]-1:
				begin = 1

		nkpt.append(nsplit+1-begin)
		for j in range(begin,int(nsplit)) :
			if j == 0 and i == 0:
				xtic.append(0)
			elif j == 0 :
				xtic.append(xtic[-1])
			else :
				xtic.append(xtic[-1]+dist_val/nsplit)
			KPT.append([symk[kp_st][x]+dist_rec[x]/nsplit*j for x in range(3)])

		xtic.append(xtic[-1]+dist_val/nsplit)
		KPT.append([symk[kp_end][x] for x in range(3)])

	with open(target+'/xtic_hse.dat','w') as out_xtic:
		for i in range(len(xtic)):
			out_xtic.write(str(xtic[i])+'\t'+'\t'.join([str(x) for x in KPT[i]])+'\n')

def get_band_reorder(Band,KPT,fermi,spin,target):
	overlap_cut = 0.7
	min_overlap_cut = 0.3
	min_band = -1
	max_band = 0
	ordering_window = 3

	for n in range(len(Band)):
		if min_band == -1 and Band[n][0][0] < fermi:
			sw = 0
		elif Band[n][0][0] > fermi:
			sw = 2
		else:
			sw = 1
		for k in range(len(Band[0])):
			if sw == 1:
				break
			for i in range(len(Band[0][0])):
				if sw == 0 and Band[n][k][i] > fermi-ordering_window:
					sw = 1
					min_band = n
					break
				elif sw == 2 and Band[n][k][i] < fermi+ordering_window:
					sw = 1
					max_band = n
					break
			if sw == 2:
				break
		if sw == 2:
			break

	[band_line,over_idx] = read_xtic_file(target+'/xtic_hse.dat')
	Unk_GGA = read_wavecar(target+'/WAVECAR',min_band,max_band)
	nband = len(Unk_GGA[0,0])

	band_idx = -1 * np.ones((len(Unk_GGA),len(KPT),nband),dtype=int)
	for i in range(len(Unk_GGA)):
		band_idx[i,0,:] = range(nband)
		band_idx[i,1,:] = range(nband)

	for i in range(len(Unk_GGA)):
		for idx_GGA in band_line:
			band_idx_prev = np.copy(band_idx[i,idx_GGA[0]])
			band_idx_prev[band_idx[i,idx_GGA[0]]] = band_idx[i,idx_GGA[1]]
			for k in range(2,len(idx_GGA)):
				M = np.zeros((nband,nband), dtype=float)
				M2 = np.zeros((nband,nband), dtype=float)
				for n1 in range(min_band,max_band):
					for n2 in range(min_band,max_band):
						M[n1,n2] = np.abs(np.sum(Unk_GGA[i,idx_GGA[k-1],n1, ...].conj() * Unk_GGA[i,idx_GGA[k],n2, ...])) ** 2.0
						M2[n1,n2] = np.abs(np.sum(Unk_GGA[i,idx_GGA[k-2],n1, ...].conj() * Unk_GGA[i,idx_GGA[k],n2, ...])) ** 2.0
				R = np.zeros_like(M, dtype=int)
				R[M >= overlap_cut] = 1
				R2 = np.zeros_like(M, dtype=int)
				R2[M2 >= overlap_cut] = 1
				flag = np.zeros(nband, dtype=int)
				M_left = np.copy(M)
	
				# R[idx_target, idx_from]
				# Ex. R = [[0,0,1],...] indicates band[3] is now band_reorder[1].
				for n1 in range(min_band,max_band):
					flag = flag + R[n1]

				for n1 in range(min_band,max_band):
					if not R[band_idx_prev[n1]].sum() == 1 and R2[n1].sum() == 1:
						for n2 in range(min_band,max_band):
							if R2[n1,n2] == 1 and flag[n2] == 0:
								R[band_idx_prev[n1]] = R2[n1]
		     						break
						flag = flag + R[band_idx_prev[n1]]

				for n1 in range(min_band,max_band):
					for n2 in range(min_band,max_band):
						if R[n1].sum() == 1 and not flag[n2] == 0:
							M_left[n1] = 0.0
							M_left[:,n2] = 0.0
				while 1 :
					ind = np.unravel_index(np.argmax(M_left),M_left.shape)
					if M_left[ind] < min_overlap_cut:
						break
					else:
						R[ind] = 1
						M_left[ind[0]] = 0.0
						M_left[:,ind[1]] = 0.0
						flag = flag + R[ind[0]]

				for n1 in range(min_band,max_band):
					if not R[n1].sum() == 1:
						for n2 in range(min_band,max_band):
							if flag[n2] == 0:
								R[n1,n2] = 1
								flag[n2] = 1
								break

				for n1 in range(min_band):
					R[n1,n1] = 1
				for n1 in range(max_band,nband):
					R[n1,n1] = 1

				band_idx_prev = np.dot(R,band_idx[i,0])
				band_idx[i,idx_GGA[k]] = band_idx_prev[band_idx[i,idx_GGA[k-1]]]

			for j in range(len(over_idx)):
				band_idx[i,over_idx[j][0]] = band_idx[i,over_idx[j][1]]

	Band_reorder = []
	for n in range(nband):
		Band_reorder.append([])
		for k in range(len(KPT)):
			Band_reorder[n].append([])
			for i in range(len(Band[0][0])):
				Band_reorder[n][k].append(0)

	for n in range(nband):
		for k in range(len(KPT)):
			for i in range(len(Band[0][0])):
				Band_reorder[n][k][i] = Band[band_idx[i,k,n]][k][i]

	return Band_reorder

def find_cb_gap(Band,fermi,dir_band):
	cb_idx = []
	vb_idx = []
	for i in range(len(Band[0][0])):
		cb_idx.append([])
		vb_idx.append([])
		for n in range(len(Band)):
			if Band[n][0][i] > fermi:
				cb_idx[i].append(n)
			else:
				vb_idx[i].append(n)

	eVBM = max([max([Band[vb_idx[i][-1]][x][y] for x in range(len(Band[0]))]) for y in range(len(Band[0][0]))])
	eCBM = min([min([Band[cb_idx[i][0]][x][y] for x in range(len(Band[0]))]) for y in range(len(Band[0][0]))])

	return [vb_idx,cb_idx,eVBM,eCBM]
	

def find_cb(Band,Band_re,KPT,fermi,hse_path,target):
	[gap_hse,vbm_kp,cbm_kp] = read_band_log(hse_path+'/Band_gap.log')
	ref_kp_idx = []
	ref_kp_search = []
	cb_idx = []
	vb_idx = []
	vb_kp_idx = -1
	cb_kp_idx = -1

	for i in range(len(Band[0][0])):
		ref_kp_idx.append(-1)
		ref_kp_search.append([])
		cb_idx.append([])
		vb_idx.append([])
		for k in range(len(Band[0])):
			if vb_kp_idx == -1 and sum([abs(vbm_kp[x]-float(KPT[k][x])) for x in range(3)]) < 0.00000001:
				vb_kp_idx = k
			if cb_kp_idx == -1 and sum([abs(cbm_kp[x]-float(KPT[k][x])) for x in range(3)]) < 0.00000001:
				cb_kp_idx = k
			for n in range(1,len(Band)):
				if Band[n-1][k][i] < fermi-0.5 and Band[n][k][i] > fermi+0.5:
					if n in [ref_kp_search[i][x][1] for x in range(len(ref_kp_search[i]))]:
						for x in range(len(ref_kp_search[i])):
							if ref_kp_search[i][x][1] == n:
								ref_kp_search[i][x][2] = ref_kp_search[i][x][2]+1
					else:					
						ref_kp_search[i].append([k,n,1])
					break
		if len(ref_kp_search[i]) == 0:
			# If band correction is difficult,
			return [-1,-1,0,0]
		ref_kp_search[i] = sorted(ref_kp_search[i],key=lambda x:x[2],reverse=True)
		ref_kp_idx[i] = ref_kp_search[i][0][0]

	for n in range(len(Band_re)):
		for i in range(len(Band_re[0][0])):
			k = ref_kp_idx[i]
			if Band_re[n][k][i] > fermi:
				cb_idx[i].append(n)
			else:
				vb_idx[i].append(n)

	eVBM = max([max([Band_re[x][vb_kp_idx][y] for x in vb_idx[y]]) for y in range(len(Band_re[0][0]))])
	eCBM = min([min([Band_re[x][cb_kp_idx][y] for x in cb_idx[y]]) for y in range(len(Band_re[0][0]))])

	return [vb_idx,cb_idx,eVBM,eCBM]


# check convergence for hse calculaton (only check energy)
def convergence_check_E(target):
    path_list = glob.glob(target+'/kptest/KP*')
    if os.path.isfile(target+'/kptest/KPOINTS_converged'):
        path_list.remove(target+'/kptest/KPOINTS_converged')
    path_list.sort()
    if os.path.isfile(target+'/kptest/kpoint.log'):
        ENCONV = float(pygrep('E/atom',target+'/kptest/kpoint.log',0,0).split()[-2])
    else:
        ENCONV = 0.01
    nion = int(pygrep('NION',path_list[0]+'/OUTCAR',0,0).splitlines()[-1].split()[11])
    check = 0
    for i in range(len(path_list[:-2])):
        ENERGY = []
        if ENCONV > 0:
            for j in range(3):
                ENERGY.append(float(pygrep('free  ',path_list[i+j]+'/OUTCAR',0,0).splitlines()[-1].split()[4]))
        converge = 0
        # Convergence check for energy/atom
        if (abs(ENERGY[0]-ENERGY[1])/nion < ENCONV and abs(ENERGY[0]-ENERGY[2])/nion < ENCONV) and check == 0:
#            print path_list[i]
            E_conv_kp = path_list[i]
            check = 1
    return E_conv_kp
