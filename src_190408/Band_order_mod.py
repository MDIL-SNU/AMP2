import sys,binascii
import numpy as np
import os,subprocess
from module_band import *
from module_vector import *
from module_vasprun import pygrep,pyhead,pytail

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

def make_band_new_dat(xtic_file,Band,spin,target):
	xtic = []
	with open(xtic_file,'r') as inp:
		for line in inp:
			xtic.append(line.split()[0])
	with open(target+'/band_new.dat','w') as out:
		for k in range(len(xtic)):			# plot k-pts
			out.write(xtic[k])
			for i in range(int(spin)):		# spin
				for n in range(len(Band)):	# band idx
					out.write('\t'+str(Band[n][k][i]))
			out.write('\n')


def make_band_new_in(title,xlabel_file,fermi,gap,nband,spin,plot_range,target):
	with open(target+'/band_new.in','w') as out:
		out.write("set terminal postscript enhanced color font 'Arial, 14' size 7.2,5.4\n")
#		out.write("set output '"+target+'/'+title+".eps'\n")
		out.write("set output '"+title+"_new.eps'\n")
		title2 = title.split('_')
		out.write("set title '"+title2[0]+', ICSD#'+title2[1]+"' font 'Arial, 16'\n")
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
				out.write("\t'band_new.dat' u 1:($"+str(i+2)+"-e) w l lt 1 lc rgb 'red' notitle,\\\n")
			out.write("\t'band_new.dat' u 1:($"+str(nband+1)+"-e) w l lt 1 lc rgb 'red' title 'up',\\\n")
			for i in range(nband-1):
				out.write("\t'band_new.dat' u 1:($"+str(nband+i+2)+"-e) w l lt 1 lc rgb 'blue' notitle,\\\n")
			out.write("\t'band_new.dat' u 1:($"+str(nband*2+1)+"-e) w l lt 1 lc rgb 'blue' title 'down',\\\n")
		else:
			for i in range(nband):
				out.write("\t'band_new.dat' u 1:($"+str(i+2)+"-e) w l lt 1 lc rgb 'black' notitle,\\\n")
		out.write("\tf(x) w l lt 0 notitle\n\n")
		out.write("set terminal pngcairo enhanced font 'Arial, 14' size 960,720\n")
		out.write("set output '"+title+"_new.png'\n")
		out.write("replot")

def plot_band_new_structure(spin,Band,fermi,xtic_file,xlabel_file,plot_range,target):
	# make band.dat
	make_band_new_dat(xtic_file,Band,spin,target)

	# make band.in
	nband = len(Band)
	if pyhead(target+'/Band_gap.log',1).split()[2] == 'is' :
#	if subprocess.check_output(['head','-n','1',target+'/Band_gap.log']).split()[2] == 'is' :
		gap = 0
		fermi = str(fermi)
	else :
		gap = round(float(pyhead(target+'/Band_gap.log',1).split()[2]))
#		gap = round(float(subprocess.check_output(['head','-n','1',target+'/Band_gap.log']).split()[2]))
		fermi = pygrep('VBM',target+'/Band_gap.log',0,0).splitlines()[0].split()[-2]
#		fermi = subprocess.check_output(['grep','VBM',target+'/Band_gap.log']).splitlines()[0].split()[-2]

	title = 'Band_corrected'
	make_band_new_in(title,xlabel_file,fermi,gap,nband,spin,plot_range,target)

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

def make_xtic_all(symk,order,rec,KSPACING):
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

	with open('xtic_all.dat','w') as out_xtic:
		for i in range(len(xtic)):
			out_xtic.write(str(xtic[i])+'\t'+'\t'.join([str(x) for x in KPT[i]])+'\n')


gga_path = sys.argv[1]

overlap_cut = 0.7
min_overlap_cut = 0.3
min_band = -1
max_band = 0
ordering_window = 3

fermi = float(pyhead(gga_path+'/DOSCAR',6).splitlines()[-1].split()[3])
#fermi = float(subprocess.check_output(['head',gga_path+'/DOSCAR','-n','6']).splitlines()[-1].split()[3])
spin = pygrep('ISPIN',gga_path+'/OUTCAR',0,0).split()[2]
#spin = subprocess.check_output(['grep','ISPIN',gga_path+'/OUTCAR']).split()[2]
[KPT,Band,nelect] = EIGEN_to_array(gga_path+'/EIGENVAL',spin)

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
#print min_band,max_band
#kpt = read_kpt_file(gga_path+'/KPOINTS')
[symk,order,xticlabel,rec] = make_symk(gga_path+'/sym')
#make_xtic_all(symk,order,rec,0.01)
make_xtic_all(symk,order,rec,0.02)

[band_line,over_idx] = read_xtic_file(gga_path+'/xtic_all.dat')

Unk_GGA = read_wavecar(gga_path+'/WAVECAR',min_band,max_band)

nband = len(Unk_GGA[0,0])

#idx_GGA = [x for x in range(len(KPT))]

band_idx = -1 * np.ones((len(Unk_GGA),len(KPT),nband),dtype=int)
for i in range(len(Unk_GGA)):
	band_idx[i,0,:] = range(nband)

for i in range(len(Unk_GGA)):
	for idx_GGA in band_line:
		for k in range(2,len(idx_GGA)):
			M = np.zeros((nband,nband), dtype=float)
			M2 = np.zeros((nband,nband), dtype=float)
			for n1 in range(min_band,max_band):
				for n2 in range(min_band,max_band):
					M[n1,n2] = np.abs(np.sum(Unk_GGA[i,idx_GGA[k],n1, ...].conj() * Unk_GGA[i,idx_GGA[k-1],n2, ...])) ** 2.0
					M2[n1,n2] = np.abs(np.sum(Unk_GGA[i,idx_GGA[k],n1, ...].conj() * Unk_GGA[i,idx_GGA[k-2],n2, ...])) ** 2.0
#			np.savetxt('Comp_Mat_{:d}_{:03d}.dat'.format(i+1,k+1),M,fmt='%4.2f')
#			if k < 40 and k > 31:
#				print [Band[x][k][i] for x in [25,26,27,28,29,30]]
#				np.savetxt('Mat1_{:d}_{:03d}.dat'.format(i+1,k+1),M[25:31,25:31],fmt='%4.2f')
#				np.savetxt('Mat2_{:d}_{:03d}.dat'.format(i+1,k+1),M2[25:31,25:31],fmt='%4.2f')
#				M3 = np.zeros((nband,nband), dtype=float)
#				M4 = np.zeros((nband,nband), dtype=float)
#				for n1 in range(nband):
#					for n2 in range(nband):
#						M3[n1,n2] = np.abs(np.sum(Unk_GGA[i,idx_GGA[k],n1, ...].conj() * Unk_GGA[i,idx_GGA[k-4],n2, ...])) ** 2.0
#						M4[n1,n2] = np.abs(np.sum(Unk_GGA[i,idx_GGA[k],n1, ...].conj() * Unk_GGA[i,idx_GGA[k-6],n2, ...])) ** 2.0
#				np.savetxt('Mat3_{:d}_{:03d}.dat'.format(i+1,k+1),M3[58:66,58:66],fmt='%4.2f')
#				np.savetxt('Mat4_{:d}_{:03d}.dat'.format(i+1,k+1),M4[58:66,58:66],fmt='%4.2f')
#			if k== 31:
#				sys.exit()

			R = np.zeros_like(M, dtype=int)
			R[M >= overlap_cut] = 1
			R2 = np.zeros_like(M, dtype=int)
			R2[M2 >= overlap_cut] = 1
			flag = np.zeros(nband, dtype=int)

			for n1 in range(min_band,max_band):
				flag = flag + R[n1]

			for n1 in range(min_band,max_band):
				if not R[n1].sum() == 1 and R2[n1].sum() == 1:
					for n2 in range(min_band,max_band):
						if R2[n1,n2] == 1 and flag[n2] == 0:
							R[n1] = R2[n1]
							break
					flag = flag + R[n1]


			for n1 in range(min_band,max_band):
				if not R[n1].sum() == 1:
					R[n1] = 0
					for n2 in np.argsort(M[n1])[::-1]:
						if flag[n2] == 0 and M[n1,n2] > min_overlap_cut:
							R[n1,n2] = 1
							break
				flag = flag + R[n1]

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

			band_idx[i,idx_GGA[k]] = np.dot(R, band_idx[i,idx_GGA[k-1]])
	for j in range(len(over_idx)):
		band_idx[i,over_idx[j][0]] = band_idx[i,over_idx[j][1]]

np.savetxt('total_idx.dat',band_idx[0] +1,fmt = '%3d')

[gap_hse,vbm_kp,cbm_kp] = read_band_log(sys.argv[2])

Band_reorder = []

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
		for n in range(1,nband):
			if Band[n-1][k][i] < fermi-0.5 and Band[n][k][i] > fermi+0.5:
				if n in [ref_kp_search[i][x][1] for x in range(len(ref_kp_search[i]))]:
					for x in range(len(ref_kp_search[i])):
						if ref_kp_search[i][x][1] == n:
							ref_kp_search[i][x][2] = ref_kp_search[i][x][2]+1
				else:					
					ref_kp_search[i].append([k,n,1])
				break
	ref_kp_search[i] = sorted(ref_kp_search[i],key=lambda x:x[2],reverse=True)
	ref_kp_idx[i] = ref_kp_search[i][0][0]

#print ref_kp_idx
#print KPT[ref_kp_idx[0]]
#print KPT[ref_kp_idx[1]]
for n in range(nband):
	Band_reorder.append([])
	for k in range(len(KPT)):
		Band_reorder[n].append([])
		for i in range(len(Band[0][0])):
			Band_reorder[n][k].append(0)

for n in range(nband):
	for k in range(len(KPT)):
		for i in range(len(Band[0][0])):
#			Band_reorder[n][k].append(Band[band_idx[i,k,n]][k][i])
			Band_reorder[band_idx[i,k,n]][k][i] = Band[n][k][i]
			if k == ref_kp_idx[i]:
				if Band_reorder[band_idx[i,k,n]][k][i] > fermi:
					cb_idx[i].append(band_idx[i,k,n])
				else:
					vb_idx[i].append(band_idx[i,k,n])

#for k in [32,33,34,35,36,37,38,39]:
#	print [Band_reorder[x][k][0] for x in [25,26,27,28,29,30]]

eVBM = max([max([Band_reorder[x][vb_kp_idx][y] for x in vb_idx[y]]) for y in range(len(Band[0][0]))])
eCBM = min([min([Band_reorder[x][cb_kp_idx][y] for x in cb_idx[y]]) for y in range(len(Band[0][0]))])

E_shift = gap_hse+eVBM-eCBM

for i in range(len(Band[0][0])):
	for n in cb_idx[i]:
		for k in range(len(KPT)):
			Band_reorder[n][k][i] = Band_reorder[n][k][i] + E_shift

plot_band_new_structure(spin,Band_reorder,eVBM,gga_path+'/xtic.dat',gga_path+'/xlabel.dat',[5,5],gga_path)
