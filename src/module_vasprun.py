###########################################
### Date: 2018-12-05			###
### yybbyb@snu.ac.kr			###
###########################################
# This is a package of modules to run VASP.
import subprocess
from module_log import *
from module_vector import *

# copy vasp input file (POSCAR)
def copy_input(source,target,pot_type):
	subprocess.call(['cp',source+'/INCAR',target+'/.'])
	subprocess.call(['cp',source+'/POSCAR',target+'/.'])
	subprocess.call(['cp',source+'/KPOINTS',target+'/.'])
	subprocess.call(['cp',source+'/POTCAR_'+pot_type,target+'/POTCAR'])

# copy vasp input file (CONTCAR)
def copy_input_cont(source,target):
	subprocess.call(['cp',source+'/INCAR',target+'/.'])
	subprocess.call(['cp',source+'/CONTCAR',target+'/POSCAR'])
	subprocess.call(['cp',source+'/KPOINTS',target+'/.'])
	subprocess.call(['cp',source+'/POTCAR',target+'/.'])

# copy vasp input file without KPOINTS (for kptest)
def copy_input_no_kp(source,target,pot_type):
	subprocess.call(['cp',source+'/INCAR',target+'/.'])
	subprocess.call(['cp',source+'/POSCAR',target+'/.'])
	subprocess.call(['cp',source+'/POTCAR_'+pot_type,target+'/POTCAR'])

# generate kpoints file from kpl and symmetry
def kpt_generation_for_relax(target,KPL,sym):
	kpoint = open(target+'/KPOINTS','w')
	# Gamma-centred mesh for hexagoanl and rhombohedral symmetry
	if sym==12 or sym==13 or sym==14 or sym=='gamma':
#	if sym==6 or sym==12 or sym==13 or sym==14 or sym==10 or sym==15:
		KPset = 'Gamma-centered'
	else :
		KPset = 'Monk-horst'
	axis = poscar_to_axis(target+'/POSCAR')
	recipro_latt = reciprocal_lattice(axis)
	l = []
	for i in range(3):
		l.append((recipro_latt[i][0]**2.+recipro_latt[i][1]**2.+recipro_latt[i][2]**2.)**0.5)
	KP = []
	for i in range(3) :
		if sym in [5,6,10]:  # symmetry BCT and ORCI
			KP.append(str(KPL))
#			KP.append(str(int(round(l[i]/(max(l)/KPL)))))
		else:
			KP.append(str(int(round(l[i]/(max(l)/KPL)))))
		if KP[-1] == '0' :
			KP[-1] = '1'
	# Write KPOINTS
	kpoint.write("Auto k-point\n 0\n"+KPset+"\n  "+KP[0]+"  "+KP[1]+"  "+KP[2]+"\n  0  0  0\n")
	kpoint.close()

def incar_from_yaml(target,yaml_incar):
	if bool(yaml_incar):
		for incar_set_key in list(yaml_incar.keys()):
			wincar(target+'/INCAR',target+'/INCAR',[[incar_set_key.upper(),str(yaml_incar[incar_set_key])]],[])

# set npar and kpar
def set_parallel(kpoints,incar,npar,kpar):
	with open(kpoints,'r') as kp_file:
		kp_line = kp_file.readlines()
	if kp_line[2][0].lower() == 'm' or kp_line[2][0].lower() == 'g': # Monk-horst of gamma
		KP = kp_line[3].split()
		if int(KP[0])*int(KP[1])*int(KP[2]) == 1 :
			nparf = str(npar*kpar)
			kparf = '#'
			wincar(incar,incar,[['NPAR',nparf],['KPAR',kparf]],[])
			return 1
		else :
			nparf = str(npar); kparf = str(kpar)
			wincar(incar,incar,[['NPAR',nparf],['KPAR',kparf]],[])
			return 0
	else:	# reciprocal or cart
		if int(kp_line[1]) == 1:
			nparf = str(npar*kpar)
			kparf = '#'
			wincar(incar,incar,[['NPAR',nparf],['KPAR',kparf]],[])
			if all([float(x) == 0.0 for x in kp_line[3].split()[0:3]]):
				return 1
			else:
				return 0
		else:
			nparf = str(npar); kparf = str(kpar)
			wincar(incar,incar,[['NPAR',nparf],['KPAR',kparf]],[])
			return 0

# run vasp	
def run_vasp(target,nproc,vasprun,mpi):
	# 0: well done, 1: error in vasp calculation
	mpi_core = {'mpirun':'-np','jsrun':'--np','srun':'-np','mpiexec':'-np','mpiexec.hydra':'-np','mpich':'-np'}
	os.chdir(target)
	if mpi in list(mpi_core.keys()):
		mpi_command = mpi+' '+mpi_core[mpi]
	else:
		mpi_command = mpi
	out = subprocess.call([mpi_command+' '+nproc+' '+vasprun+' >& stdout.x'], stdout=subprocess.PIPE, shell=True)
	if out == 0:
		out_res = pytail(target+'/OUTCAR').split()
		if len(out_res) > 0 and out_res[0] == 'Voluntary':
			make_amp2_log(target,'VASP calculation is performed successfully.')
			subprocess.call(['rm',target+'/vasprun.xml'])
			write_log_in_outcar(target+'/OUTCAR',target+'/amp2.log')
			return 0
		else:
			make_amp2_log(target,'ERROR occurs during vasp calculation. Check the calculation.')
			return 1
	else:
		make_amp2_log(target,'ERROR occurs. vasp calculation may be not started. Check the calculation.')
		return 1

# run vasp for relax (+check stdout.x)	
def run_vasp_rlx(target,nproc,vasprun,mpi):
	# 0: well done, 1: error in vasp calculation
	mpi_core = {'mpirun':'-np','jsrun':'--np','srun':'-np','mpiexec':'-np','mpiexec.hydra':'-np','mpich':'-np'}
	os.chdir(target)
	if mpi in list(mpi_core.keys()):
		mpi_command = mpi+' '+mpi_core[mpi]
	else:
		mpi_command = mpi
	out = subprocess.call([mpi_command+' '+nproc+' '+vasprun+' >& stdout.x'], stdout=subprocess.PIPE, shell=True)
	if out == 0:
		out_res = pytail(target+'/OUTCAR').split()
		std_res = pytail(target+'/stdout.x').split()
		if len(out_res) > 0 and out_res[0] == 'Voluntary':
			make_amp2_log(target,'VASP calculation is performed successfully.')
			subprocess.call(['rm',target+'/vasprun.xml'])
			write_log_in_outcar(target+'/OUTCAR',target+'/amp2.log')
			return 0
		elif (std_res[1] == 'POSCAR') and (std_res[3] == 'continue'):
			make_amp2_log(target,'VASP calculation is performed successfully.')
			subprocess.call(['rm',target+'/vasprun.xml'])
			write_log_in_outcar(target+'/OUTCAR',target+'/amp2.log')
			return 0
		else:
			make_amp2_log(target,'ERROR occurs during vasp calculation. Check the calculation.')
			return 1
	else:
		make_amp2_log(target,'ERROR occurs. vasp calculation may be not started. Check the calculation.')
		return 1

# check electronic step convergence
def electronic_step_convergence_check(target):
	# 0: well converged, 1: can try to calculation with other setting (retry), 2: cannot be converged (stop)
	with open(target+'/OSZICAR','r') as inp:
		fr_log = inp.readlines()[1:]
	spin = pygrep('ISPIN',target+'/OUTCAR',0,0).split()[2]
	elec_step = 0
	for ll in fr_log:
		if ll.split()[0] == '1':
			break
		elec_step = elec_step+1
	if elec_step == int(pygrep('NELM',target+'/OUTCAR',0,0).split(';')[0].split()[2]):
		make_amp2_log(target,'Electronic step is not converged.')
		algo = pygrep('ALGO',target+'/INCAR',0,0).split()[2]
		if algo == 'Normal' or algo == 'All':
			make_amp2_log(target,'Current ALGO is '+algo+' but it is not converged.')
			if spin == '2':
				if pygrep('BMIX_MAG',target+'/OUTCAR',0,0).split()[-1] == '1.00':
					make_amp2_log(target,'Change mixing parameter.')
					wincar(target+'/INCAR',target+'/INCAR',[['AMIX','0.2'],['BMIX','0.0001'],['AMIX_MAG','0.8'],['BMIX_MAG','0.0001']],[])
					return 1
				else:
					make_amp2_log(target,'We changed mixing parameters but it is not converged.')
					return 2
			else:
				return 2
		elif algo == 'Damped':
			wincar(target+'/INCAR',target+'/INCAR',[['ALGO','All']],[])
			make_amp2_log(target,'ALGO changes from '+algo+' to All.')
			return 1
		else:
			wincar(target+'/INCAR',target+'/INCAR',[['ALGO','Normal']],[])
			make_amp2_log(target,'ALGO changes from '+algo+' to Normal.')
			return 1
	else:
		return 0

# check electronic step convergence for CHGCAR
def electronic_step_convergence_check_CHGCAR(target):
	# 0: well converged, 1: can try to read CHGCAR
	with open(target+'/OSZICAR','r') as inp:
		fr_log = inp.readlines()[1:]
	spin = pygrep('ISPIN',target+'/OUTCAR',0,0).split()[2]
	elec_step = 0
	for ll in fr_log:
		if ll.split()[0] == '1':
			break
		elec_step = elec_step+1
	if elec_step == int(pygrep('NELM',target+'/OUTCAR',0,0).split(';')[0].split()[2]):
		wincar(target+'/INCAR',target+'/INCAR',[['NSW','0'],['ISYM','0'],['NELM','50'],['LCHARG','T'],['ICHARG','1']],[])
		return 1
	else:
		return 0

# check magnetism from relaxation
def check_magnet(dir_in,min_mom):
	nion = int(pygrep('NION',dir_in+'/OUTCAR',0,0).split()[-1])
	spin = pygrep('ISPIN',dir_in+'/OUTCAR',0,0).split()[2]
	mag_on = 0
	if spin == '2' :
		mag = [float(x.split()[-1]) for x in pygrep('magnetization (x)',dir_in+'/OUTCAR',0,nion+3).splitlines()[-nion:]]
		for i in range(len(mag)) :
			if abs(mag[i]) > min_mom :
				mag_on = 1
				break
	return mag_on

# make incar file including ncl calculation
def make_incar_for_ncl(dir_target,mag_on,kpar,npar,vasp_std,vasp_gam,vasp_ncl):
	# set parallel option
	gam = set_parallel(dir_target+'/KPOINTS',dir_target+'/INCAR',npar,kpar)
	if gam == 1:
		vasprun = vasp_gam
	else:
		vasprun = vasp_std

	# Fix INCAR
	if mag_on == 1 :
		spin = '2'
		magmom = '='
	elif mag_on == 0 :
		spin = '1'
		magmom = ''

	# Check LS coupling option
	with open(dir_target+'/INCAR','r') as inp:
		tmp = inp.read()
	if 'LSORBIT' in tmp and not vasp_ncl == vasp_std:
		vasprun = vasp_ncl
		SOC = '='
		maxmix = ''
		if 'T' in tmp.split('=')[1] :
			SOC = '.True.'
			if 'LMAXMIX' in tmp:
				maxmix = '!#'
			with open('SOC_note','w') as soc_note:
				soc_note.write('Spin-orbit coupling is considered!')
		# Fix INCAR
		wincar(dir_target+'/INCAR',dir_target+'/INCAR',[['MAGMOM',''],['LSORBIT',SOC],['LMAXMIX',maxmix],['ISPIN','2']],[])
	elif not mag_on == 2:
		wincar(dir_target+'/INCAR',dir_target+'/INCAR',[['ISPIN',spin],['MAGMOM',magmom]],[])

	return vasprun

# incar re-write 
def wincar(SOURCE,OUT,option,add) :
	incar = open(SOURCE,'r').readlines()
	out = open(OUT,'w')
	for i in range(len(option)) :
		check = 0
		option[i] = [str(x) for x in option[i]]
		for j in range(len(incar)) :
			if option[i][0].upper() in [x.upper() for x in incar[j].split()] :
				check = 1
				# Off the line
				if option[i][1] == '#' :
					incar[j] = '#'+incar[j].replace('#','')
				# Recover the line
				elif option[i][1] == '!#' :
					incar[j] = incar[j].replace('#','')
				# Remove the line
				elif option[i][1] == '' :
					incar[j] = ''
				# Leave the line
				elif option[i][1] == '=' :
					incar[j] = incar[j]
				else :
					incar[j] = '   '+option[i][0]+' = '+option[i][1]+'\n'
		if check == 0 and option[i][1] != '' :
			if option[i][1] == '#' :
				add.append('#   '+option[i][0]+' = ')
			elif option[i][1] != '=' :
				add.append('   '+option[i][0]+' = '+option[i][1])
	for j in range(len(incar)) :
		out.write(incar[j])
	for i in range(len(add)) :
		out.write(add[i]+'\n')
	out.close()
	return 0

def poscar_to_axis(POSCAR):
	with open(POSCAR,'r') as poscar:
		lines = poscar.readlines()
	axis_scale = float(lines[1].split()[0])
	axis = []
	for i in range(3):
		axis.append([float(x)*axis_scale for x in lines[2+i].split()])
	return axis

def count_line(filename):
	with open(filename,'r') as f:
		lines = 0
		for line in f:
			lines = lines + 1
	return lines

# make kpoints for dos and dielectric constant
def make_multiple_kpts(kp_log,kpt_file,pos_file,kp_multi,sym,gam_option):
	if os.path.isfile(kp_log):
		with open(kp_log,'r') as inp:
			KPL = kp_multi*int(inp.readlines()[-1].split()[-1])
		with open(kpt_file,'r') as kpt:
			khead = kpt.readlines()

		axis = poscar_to_axis(pos_file)
		recipro_latt = reciprocal_lattice(axis)
		l = []
		for i in range(3):
			l.append((recipro_latt[i][0]**2.+recipro_latt[i][1]**2.+recipro_latt[i][2]**2.)**0.5)
		idx = l.index(max(l))
		KP=[]
		for i in range(3) :
			if sym in [5,6,10]: # for BCT and ORCI
				KP.append(str(KPL))
			else:
				KP[i:] = [str(int(round(l[i]/(l[idx]/KPL))))]
				if KP[i] == '0' :
					KP[i] = '1'
	else:
		with open(kpt_file,'r') as kpt:
			khead = kpt.readlines()
		KP_ori = khead[3].split()
		KP=[]
		for i in range(3) :
			if sym == 5 or sym == 6:
				KP.append(str(max([int(x) for x in KP_ori[i]])*kp_multi))
			else:
				KP.append(str(int(KP_ori[i])*kp_multi))
	if gam_option == 'gamma':
		khead[2] = 'Gamma-centered\n'
	with open(kpt_file,'w') as kpt:
		kpt.write(khead[0]+khead[1]+khead[2]+"  "+KP[0]+"  "+KP[1]+"  "+KP[2]+"\n  0  0  0\n")

def incar_for_hse(incar_file):
	wincar(incar_file,incar_file,[['LDA+U',''],['LDAU',''],['LDAUTYPE',''],['LDAUL',''],['LDAUU',''],['LDAUJ',''],['LDAUPRINT','']],[])
	with open(incar_file,'r') as inc:
		if '\nHybrid calculation' in inc.read():
			wincar(incar_file,incar_file,[['LMAXMIX',''],['ISYM','3'],['LHFCALC','.T.'],['HFSCREEN','0.2'],['PRECFOCK','Normal'],['ALGO','Damped'],['AEXX','0.25']],[])
		else:
			wincar(incar_file,incar_file,[['ALGO',''],['LMAXMIX',''],['ISYM','3']],['\n\nHybrid calculation:\n   LHFCALC = .T.\n   HFSCREEN = 0.2\n   PRECFOCK = Normal\n   ALGO = Damped\n   AEXX = 0.25\n'])

def write_relaxed_poscar(target,pot_type):
	from module_amp2_input import read_poscar,write_poscar
	[axis_ori,atom_pos_ori] = read_poscar(target+'/INPUT0/POSCAR')
	[axis_rlx,atom_pos_rlx] = read_poscar(target+'/relax_'+pot_type+'/CONTCAR')
	new_atom_pos = []
	for i in range(len(atom_pos_rlx)):
		new_atom_pos.append(atom_pos_rlx[i][0:3]+atom_pos_ori[i][3:])
	write_poscar(axis_rlx,new_atom_pos,target+'/INPUT0/POSCAR_rlx_'+pot_type,'relaxed poscar')

# This function works like shell command in python script 
def pygrep(pattern,filename,prev_line,after_line):
	results = ''
	with open(filename,'r') as f:
		lines = f.readlines()
		for i,line in enumerate(lines):
			if pattern in line:
				for j in range(i-int(prev_line),i+1+int(after_line)):
					results = results+lines[j]
	return results

def pyhead(filename,num_line):
	results = ''
	with open(filename,'r') as f:
		for i in range(num_line):
			results = results + f.readline()
	return results

def pytail(filename):
	import os
	try:
		with open(filename,'rb') as f:
			f.seek(-2,os.SEEK_END)
			while f.read(1) != b'\n':
				f.seek(-2,os.SEEK_CUR)
			return f.readline().decode()
	except:
		with open(filename,'r') as f:
			return f.readline()
