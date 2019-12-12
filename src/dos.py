####################################
# Modifier : yybbyb@snu.ac.kr      #
# data : 2018-12-05                #
####################################
import shutil, os, sys, subprocess, yaml
from module_log import *
from module_vasprun import *
from module_dos import *
from input_conf import set_on_off
code_data = 'Version 0.9.4. Modified at 2019-11-28'

dir = sys.argv[1]

inp_file = sys.argv[2]
with open(inp_file,'r') as f:
	inp_yaml = yaml.safe_load(f)
ERROR_path = inp_yaml['directory']['error']
src_path = inp_yaml['directory']['src_path']
vasp_std = inp_yaml['program']['vasp_std']
vasp_gam = inp_yaml['program']['vasp_gam']
vasp_ncl = inp_yaml['program']['vasp_ncl']
mpi = inp_yaml['program']['mpi_command']
gnuplot = inp_yaml['program']['gnuplot']
npar = inp_yaml['vasp_parallel']['npar']
kpar = inp_yaml['vasp_parallel']['kpar']
inp_dos = inp_yaml['density_of_states']

node = node_simple(sys.argv[3])
nproc = sys.argv[4]

pot_type = sys.argv[5]
if pot_type == 'LDA':
	POT = 'LDA'
else:
	POT = 'GGA'

# Set directory for input structure and INCAR
dir_dos = dir+'/dos_'+pot_type

# Check existing data
if os.path.isdir(dir_dos) and os.path.isfile(dir_dos+'/Pdos_dat/Tot_dos.dat') :
	make_amp2_log_default(dir,src_path,'DOS calculation with '+pot_type+' potential.',node,code_data)
	make_amp2_log(dir,'DOS calculation is already done.\n')
	print 1
	sys.exit()

if not os.path.isdir(dir_dos):
	os.mkdir(dir_dos,0755)

os.chdir(dir_dos)

make_amp2_log_default(dir_dos,src_path,'DOS calculation with '+pot_type+' potential.',node,code_data)

# check band or relax calculation
no_rlx = 0
if not os.path.isdir(dir+'/relax_'+pot_type):
	make_amp2_log(dir_dos,'Relax directory does not exist.')
	no_rlx = 1
else:
	if not os.path.isfile(dir+'/relax_'+pot_type+'/CONTCAR'):
		make_amp2_log(dir_dos,'CONTCAR file in relaxation does not exist.')
		no_rlx = 1
	elif count_line(dir+'/relax_'+pot_type+'/CONTCAR') < 9:
		make_amp2_log(dir_dos,'CONTCAR file in relaxation is invalid.')
		no_rlx = 1

if no_rlx == 1 and set_on_off(inp_dos['relax_check']) == 1:
	print 0
	sys.exit()

dir_band = dir+'/band_'+pot_type

# Potential type
if pot_type == 'HSE':
	if no_rlx == 1:
		copy_input(dir+'/INPUT0',dir_dos,POT)
	else:
		copy_input_cont(dir+'/relax_'+pot_type,dir_dos)
	# make INCAR for CHGCAR
	incar_from_yaml(dir_dos,inp_dos['incar'])
	if no_rlx == 1:
		mag_on = 2
	else:
		mag_on = check_magnet(dir+'/relax_'+pot_type,inp_yaml['magnetic_ordering']['minimum_moment'])
	vasprun = make_incar_for_ncl(dir_dos,mag_on,kpar,npar,vasp_std,vasp_gam,vasp_ncl)
	with open(dir+'/INPUT0/sym','r') as symf:
		sym = int(symf.readline().split()[0])
	make_multiple_kpts(dir+'/kptest/kpoint.log',dir_dos+'/KPOINTS',dir_dos+'/POSCAR',inp_dos['kp_multiplier'],sym)
	incar_for_hse(dir_dos+'/INCAR')
	hse_algo = pygrep('ALGO',dir+'/relax_'+pot_type+'/INCAR',0,0).split()[2]
#	hse_algo = subprocess.check_output(['grep','ALGO',dir+'/relax_'+pot_type+'/INCAR']).split()[2]
	wincar(dir_dos+'/INCAR',dir_dos+'/INCAR',[['NSW','0'],['NEDOS','3001'],['ALGO',hse_algo]],[])
	incar_from_yaml(dir_dos,inp_dos['incar'])
	vasprun = vasp_std
	# VASP calculation
	out = run_vasp(dir_dos,nproc,vasprun,mpi)
	if out == 1:  # error in vasp calculation
		print 0
		sys.exit() 
	out = electronic_step_convergence_check(dir_dos)
	while out == 1:
		make_amp2_log(dir_dos,'Calculation options are changed. New calculation starts.')
		out = run_vasp(dir_dos,nproc,vasprun,mpi)
		if out == 1:  # error in vasp calculation
			print 0
			sys.exit()
		out = electronic_step_convergence_check(dir_dos)

	if out == 2:  # electronic step is not converged. (algo = normal)
		make_amp2_log(dir_dos,'The calculation stops but electronic step is not converged.')
		print 0
		sys.exit()

else:
	# check CHGCAR
	if os.path.isfile(dir_dos+'/CHGCAR') and os.path.getsize(dir_dos+'/CHGCAR') > 0 :
		make_amp2_log(dir_dos,'DOS calculation is performed by using existing CHGCAR file.')
	elif os.path.isdir(dir_band) and os.path.isfile(dir_band+'/CHGCAR') and os.path.getsize(dir_band+'/CHGCAR') > 0 :
		make_amp2_log(dir_dos,'DOS calculation is performed by using existing CHGCAR file for band calculation.')
		# Copy input data and write CHGCAR
		if no_rlx == 1:
			copy_input(dir+'/INPUT0',dir_dos,POT)
		else:
			copy_input_cont(dir+'/relax_'+pot_type,dir_dos)
			# make INCAR for dos
			incar_from_yaml(dir_dos,inp_dos['incar'])
			mag_on = check_magnet(dir+'/relax_'+pot_type,inp_yaml['magnetic_ordering']['minimum_moment'])
			vasprun = make_incar_for_ncl(dir_dos,mag_on,kpar,npar,vasp_std,vasp_gam,vasp_ncl)

		subprocess.call(['cp',dir_band+'/CHGCAR',dir_dos+'/.'])
	else:
		make_amp2_log(dir_dos,'VASP calculation for CHGCAR file.')
		# Copy input data and write CHGCAR
		if no_rlx == 1:
			copy_input(dir+'/INPUT0',dir_dos,POT)
		else:
			copy_input_cont(dir+'/relax_'+pot_type,dir_dos)
		# make INCAR for CHGCAR
		incar_from_yaml(dir_dos,inp_dos['incar'])
		if no_rlx == 1:
			mag_on = 2
		else:
			mag_on = check_magnet(dir+'/relax_'+pot_type,inp_yaml['magnetic_ordering']['minimum_moment'])
		vasprun = make_incar_for_ncl(dir_dos,mag_on,kpar,npar,vasp_std,vasp_gam,vasp_ncl)
		wincar(dir_dos+'/INCAR',dir_dos+'/INCAR',[['NSW','0'],['LCHARG','.T.']],[])

		# VASP calculation for CHGCAR
		out = run_vasp(dir_dos,nproc,vasprun,mpi)
		if out == 1:  # error in vasp calculation
			print 0
			sys.exit() 
		make_amp2_log(dir_dos,'CHGCAR file is generated successfully.')

	with open(dir+'/INPUT0/sym','r') as symf:
		sym = int(symf.readline().split()[0])
	make_multiple_kpts(dir+'/kptest/kpoint.log',dir_dos+'/KPOINTS',dir_dos+'/POSCAR',inp_dos['kp_multiplier'],sym)
	incar_from_yaml(dir_dos,inp_dos['incar'])

	wincar(dir_dos+'/INCAR',dir_dos+'/INCAR',[['NSW','0'],['ISTART','1'],['ICHARG','11'],['LCHARG','.F.'],['NEDOS','3001']],[])
	# make INCAR for CHGCAR
	if no_rlx == 1:
		mag_on = 2
	else:
		mag_on = check_magnet(dir+'/relax_'+pot_type,inp_yaml['magnetic_ordering']['minimum_moment'])
	vasprun = make_incar_for_ncl(dir_dos,mag_on,kpar,npar,vasp_std,vasp_gam,vasp_ncl)
	# VASP calculation
	out = run_vasp(dir_dos,nproc,vasprun,mpi)
	if out == 1:  # error in vasp calculation
		print 0
		sys.exit() 
	out = electronic_step_convergence_check(dir_dos)
	while out == 1:
		make_amp2_log(dir_dos,'Calculation options are changed. New calculation starts.')
		out = run_vasp(dir_dos,nproc,vasprun,mpi)
		if out == 1:  # error in vasp calculation
			print 0
			sys.exit()
		out = electronic_step_convergence_check(dir_dos)

	if out == 2:  # electronic step is not converged. (algo = normal)
		make_amp2_log(dir_dos,'The calculation stops but electronic step is not converged.')
		print 0
		sys.exit()

# set fermi level
fermi = float(pyhead(dir_dos+'/DOSCAR',6).splitlines()[-1].split()[3])
#fermi = float(subprocess.check_output(['head',dir_dos+'/DOSCAR','-n','6']).splitlines()[-1].split()[3])
gap = 0
if os.path.isdir(dir_band) and os.path.isfile(dir_band+'/Band_gap.log'):
	if not pyhead(dir_band+'/Band_gap.log',1,).split()[2] == 'is' :
#	if not subprocess.check_output(['head','-n','1',dir_band+'/Band_gap.log']).split()[2] == 'is' :
		fermi = float(pygrep('VBM',dir_band+'/Band_gap.log',0,0).splitlines()[0].split()[-2])
#		fermi = float(subprocess.check_output(['grep','VBM',dir_band+'/Band_gap.log']).splitlines()[0].split()[-2])
		gap = round(float(pyhead(dir_band+'/Band_gap.log',1).split()[2]))
#		gap = round(float(subprocess.check_output(['head','-n','1',dir_band+'/Band_gap.log']).split()[2]))
spin = pygrep('ISPIN',dir_dos+'/OUTCAR',0,0).split()[2]
#spin = subprocess.check_output(['grep','ISPIN',dir_dos+'/OUTCAR']).split()[2]
ncl = pygrep('NONCOL',dir_dos+'/OUTCAR',0,0).split()[2]
#ncl = subprocess.check_output(['grep','NONCOL',dir_dos+'/OUTCAR']).split()[2]
# read atom information from poscar
[atom_name,atom_num] = poscar_to_atom_inform(dir_dos+'/POSCAR')
# read dos information from doscar
[Ene,Tot_dos,par_dos] = make_dos_dat(dir_dos+'/DOSCAR',spin,atom_num,ncl)

# write dos dat files
write_tot_dos(Ene,Tot_dos,fermi,dir_dos)
write_par_dos(Ene,par_dos,atom_name,fermi,dir_dos)

make_dos_in(dir_dos,atom_name,spin,len(par_dos[0][0]),[inp_dos['y_min'],inp_dos['y_max']+gap])

if inp_yaml['calculation']['plot'] == 1:
	os.chdir(dir_dos+'/Pdos_dat')
	subprocess.call([gnuplot,dir_dos+'/Pdos_dat/dos.in'])

make_amp2_log(dir_dos,'DOS calculation is done.')

with open(dir_dos+'/amp2.log','r') as amp2_log:
	with open(dir+'/amp2.log','a') as amp2_log_tot:
		amp2_log_tot.write(amp2_log.read())
print 1
