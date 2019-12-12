###########################################
### Date: 2018-12-06			###
### yybbyb@snu.ac.kr			###
###########################################
import shutil, os, sys, subprocess, yaml
from module_log import *
from module_vasprun import *
from module_relax import *
code_data = 'Version 0.9.4. Modified at 2019-11-28'

dir = sys.argv[1]

inp_file = sys.argv[2]
with open(inp_file,'r') as f:
	inp_yaml = yaml.safe_load(f)
ERROR_path = inp_yaml['directory']['error']
src_path = inp_yaml['directory']['src_path']
vasp_std = inp_yaml['program']['vasp_std']
vasp_gam = inp_yaml['program']['vasp_gam']
mpi = inp_yaml['program']['mpi_command']
gnuplot = inp_yaml['program']['gnuplot']
npar = inp_yaml['vasp_parallel']['npar']
kpar = inp_yaml['vasp_parallel']['kpar']
inp_rlx = inp_yaml['relaxation']

node = node_simple(sys.argv[3])
nproc = sys.argv[4]

pot_type = sys.argv[5]
if pot_type == 'LDA':
	POT = 'LDA'
else:
	POT = 'GGA'

dir_relax = dir+'/relax_'+pot_type
# Check existing data
run_cont = 0
if os.path.isdir(dir_relax) and os.path.isfile(dir_relax+'/CONTCAR') and os.path.isfile(dir_relax+'/free') :
	if len(pygrep('free  ',dir_relax+'/free',0,0).splitlines()) > 0 and len(pygrep('free  ',dir_relax+'/free',0,0).splitlines()) <= inp_rlx['converged_ionic_step'] :
#	if len(subprocess.check_output(['grep','free  ',dir_relax+'/free']).splitlines()) > 0 :
		make_amp2_log_default(dir,src_path,'Relaxation with '+pot_type+' potential',node,code_data)
#		print('Success!')
		make_amp2_log(dir,'Already done')
		print 1
		sys.exit()

if os.path.isdir(dir_relax):
	if os.path.isfile(dir_relax+'/CONTCAR') and count_line(dir_relax+'/CONTCAR') > 9:
		make_amp2_log_default(dir_relax,src_path,'Relaxation with '+pot_type+' potential',node,code_data)
		make_amp2_log(dir_relax,'Calculation continue.')
		subprocess.call(['cp',dir_relax+'/CONTCAR',dir_relax+'/POSCAR'])
		run_cont = 1
		os.chdir(dir_relax)
	elif os.path.isfile(dir_relax+'/POSCAR') and count_line(dir_relax+'/POSCAR') > 9:
		make_amp2_log_default(dir_relax,src_path,'Relaxation with '+pot_type+' potential',node,code_data)
		make_amp2_log(dir_relax,'Calculation continue.')
		run_cont = 1
		os.chdir(dir_relax)
	else:
		subprocess.call(['rm','-r',dir_relax])

if run_cont == 0:
	os.mkdir(dir_relax,0755)
	os.chdir(dir_relax)
	make_amp2_log_default(dir_relax,src_path,'Relaxation with '+pot_type+' potential',node,code_data)
	make_amp2_log(dir_relax,'New calculation')
	copy_input(dir+'/INPUT0',dir_relax,POT)
	nsw = set_nsw(dir_relax+'/POSCAR',dir_relax+'/INCAR')
	wincar(dir_relax+'/INCAR',dir_relax+'/INCAR',[['LCHARG','.F.']],[])
	converge_condition = set_ediffg(dir_relax+'/POSCAR',inp_rlx['force'],inp_rlx['pressure'],inp_rlx['energy'])
	if not converge_condition == 0:
		wincar(dir_relax+'/INCAR',dir_relax+'/INCAR',[['EDIFFG',str(converge_condition)]],[])
	if pot_type == 'HSE':
		incar_for_hse(dir_relax+'/INCAR')

	incar_from_yaml(dir_relax,inp_rlx['incar'])

gam = set_parallel(dir_relax+'/KPOINTS',dir_relax+'/INCAR',npar,kpar)
if gam == 1:
	vasprun = vasp_gam
else:
	vasprun = vasp_std

out = run_vasp_rlx(dir_relax,nproc,vasprun,mpi)
if out == 1:  # error in vasp calculation
	print 0
	sys.exit() 

out = electronic_step_convergence_check(dir_relax)
while out == 1:
	make_amp2_log(dir_relax,'Calculation options are changed. New calculation starts.')
	out = run_vasp_rlx(dir_relax,nproc,vasprun,mpi)
	if out == 1:  # error in vasp calculation
		print 0
		sys.exit()
	out = electronic_step_convergence_check(dir_relax)

if out == 2:  # electronic step is not converged. (algo = normal)
	make_amp2_log(dir_relax,'The calculation stops but electronic step is not converged.')
	print 0
	sys.exit()

ionic_converge = 0
iteration = 1
energy = pygrep('free  ','OUTCAR',0,0).splitlines()
#energy = subprocess.check_output(['grep','free  ','OUTCAR']).splitlines()
while iteration < inp_rlx['max_iteration']:
#	with open(dir_relax+'/OUT_TOT','a') as out_log:
#		out_log.write(open(dir_relax+'/OUTCAR','r').read())
#	if len(energy) <= inp_rlx['converged_ionic_step']:
#		ionic_converge = 1
#		make_amp2_log(dir_relax,'Relaxation is done.')
#		break
	shutil.copyfile(dir_relax+'/CONTCAR',dir_relax+'/POSCAR')
	if bool(inp_rlx['incar']) and not 'EDIFFG' in inp_rlx['incar'].keys():
		converge_condition = set_ediffg(dir_relax+'/POSCAR',inp_rlx['force'],inp_rlx['pressure'],inp_rlx['energy'])
		if not converge_condition == 0:
			wincar(dir_relax+'/INCAR',dir_relax+'/INCAR',[['EDIFFG',str(converge_condition)]],[])

	out = run_vasp_rlx(dir_relax,nproc,vasprun,mpi)
	if out == 1:  # error in vasp calculation
		print 0
		sys.exit() 

	out = electronic_step_convergence_check(dir_relax)
	while out == 1:
		make_amp2_log(dir_relax,'Calculation options are changed. New calculation starts.')
		out = run_vasp_rlx(dir_relax,nproc,vasprun,mpi)
		if out == 1:  # error in vasp calculation
			print 0
			sys.exit()
		out = electronic_step_convergence_check(dir_relax)

	if out == 2:  # electronic step is not converged. (algo = normal)
		make_amp2_log(dir_relax,'The calculation stops but electronic step is not converged.')
		print 0
		sys.exit()

	energy = pygrep('free  ','OUTCAR',0,0).splitlines()
#	energy = subprocess.check_output(['grep','free  ','OUTCAR']).splitlines()
	iteration = iteration+1
	make_amp2_log(dir_relax,'Iteration number is '+str(iteration))
	with open(dir_relax+'/OUT_TOT','a') as out_log:
		out_log.write(open(dir_relax+'/OUTCAR','r').read())
	if len(energy) <= inp_rlx['converged_ionic_step']:
		ionic_converge = 1
		make_amp2_log(dir_relax,'Relaxation is done.')
		break

if ionic_converge == 0 and len(energy) <= inp_rlx['converged_ionic_step']:
	ionic_converge = 1
	make_amp2_log(dir_relax,'Relaxation is done.')

if ionic_converge == 0:
	make_amp2_log(dir_relax,'The iteration reach the max iteration number. You need to increase the max iteration number')
else:
	write_relaxed_poscar(dir,pot_type)

## Check the magnetic moment in relaxed cell
mag_on = check_magnet(dir_relax,inp_yaml['magnetic_ordering']['minimum_moment'])
if mag_on == 0 :
	wincar(dir+'/INPUT0/INCAR',dir+'/INPUT0/INCAR',[['MAGMOM',''],['ISPIN','1']],[])

with open('free','w') as fr_file:
	fr_file.write(pygrep('free  ','OUTCAR',0,0).splitlines()[-1])
#	fr_file.write(subprocess.check_output(['grep','free  ','OUTCAR']).splitlines()[-1])
warn = set_pos_compare(dir_relax+'/CONTCAR',dir+'/INPUT0/POSCAR',dir_relax,inp_rlx['pos_warning_percent'])

with open(dir_relax+'/amp2.log','r') as amp2_log:
	with open(dir+'/amp2.log','a') as amp2_log_tot:
		amp2_log_tot.write(amp2_log.read())

print 1
