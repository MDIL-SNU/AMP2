###########################################
### Date: 2018-12-05			###
### yybbyb@snu.ac.kr			###
###########################################
# This is for estimating dielectric tensor
import shutil, os, sys, subprocess, yaml
from module_log import *
from module_vasprun import *
from module_dielectric import *
from input_conf import set_on_off
from _version import __version__
code_data = 'Version '+__version__+'. Modified at 2019-12-17'

dir = sys.argv[1]

inp_file = sys.argv[2]
with open(inp_file,'r') as f:
	inp_yaml = yaml.safe_load(f)
ERROR_path = inp_yaml['directory']['error']
src_path = inp_yaml['directory']['src_path']
vasp_std = inp_yaml['program']['vasp_std']
mpi = inp_yaml['program']['mpi_command']
gnuplot = inp_yaml['program']['gnuplot']
kpar = inp_yaml['vasp_parallel']['kpar']
inp_diel = inp_yaml['dielectric']

node = node_simple(sys.argv[3])
nproc = sys.argv[4]

POT = sys.argv[5]

# Set directory for input structure and INCAR
dir_diel = dir+'/dielectric_'+POT

# Check existing data
if os.path.isdir(dir_diel) and os.path.isfile(dir_diel+'/dielectric.log') :
	make_amp2_log_default(dir,src_path,'Dielectric calculation with '+POT,node,code_data)
	make_amp2_log(dir,'Dielectric constant calculation with '+POT+' is already done.')
	print(1)
	sys.exit()

if not os.path.isdir(dir_diel):
	os.mkdir(dir_diel,0o755)

os.chdir(dir_diel)

make_amp2_log_default(dir_diel,src_path,'Dielectric calculation with '+POT,node,code_data)

no_rlx = 0
# check relax calculation
if not os.path.isdir(dir+'/relax_'+POT):
	make_amp2_log(dir_diel,'Relax directory does not exist.')
	no_rlx = 1
else:
	if not os.path.isfile(dir+'/relax_'+POT+'/CONTCAR'):
		make_amp2_log(dir_diel,'CONTCAR file in relaxation does not exist.')
		no_rlx = 1
	elif count_line(dir+'/relax_'+POT+'/CONTCAR') < 9:
		make_amp2_log(dir_diel,'CONTCAR file in relaxation is invalid.')
		no_rlx = 1

if no_rlx == 1 and set_on_off(inp_diel['relax_check']) == 1:
	print(0)
	sys.exit()

# check band gap calculation
if set_on_off(inp_diel['metal_check']) == 1:
	if os.path.isfile(dir+'/band_'+POT+'/Band_gap.log'):
		with open(dir+'/band_'+POT+'/Band_gap.log','r') as inp:
			gap_log = inp.readline()
	elif os.path.isfile(dir+'/band_GGA/Band_gap.log'):
		with open(dir+'/band_GGA/Band_gap.log','r') as inp:
			gap_log = inp.readline()
	else:
		make_amp2_log(dir_diel,'Band gap calculation should be performed.')
		print(0)
		sys.exit()
else:
	make_amp2_log(dir_diel,'Do not check the gap condition. (Metal or not)')
	gap_log = '0.0'

# start dielectric caclculation
if 'etal' in gap_log:
	make_amp2_log(dir_diel,'It is metallic band structure.\nDielectric calculation is unreasonable for a metallic system.')
	print(1)
	sys.exit()
else:
	if no_rlx == 1:
		copy_input(dir+'/INPUT0',dir_diel,POT)
	else:
		copy_input_cont(dir+'/relax_'+POT,dir_diel)
	with open(dir+'/INPUT0/sym','r') as symf:
		sym = int(symf.readline().split()[0])
	make_multiple_kpts(dir+'/kptest/kpoint.log',dir_diel+'/KPOINTS',dir_diel+'/POSCAR',inp_diel['kp_multiplier'],sym,'')
	incar_from_yaml(dir_diel,inp_diel['incar'])
	if not no_rlx == 1:
		mag_on = check_magnet(dir+'/relax_'+POT,inp_yaml['magnetic_ordering']['minimum_moment'])
		if mag_on == 0:
			wincar(dir_diel+'/INCAR',dir_diel+'/INCAR',[['ISPIN','1'],['MAGMOM','']],[])
	wincar(dir_diel+'/INCAR',dir_diel+'/INCAR',[['NSW','#'],['IBRION','8'],['NPAR','#'],['KPAR',kpar],['LEPSILON','.T.']],[])
	vasprun = vasp_std
	# VASP calculation
	out = run_vasp(dir_diel,nproc,vasprun,mpi)
	if out == 1:  # error in vasp calculation
		print(0)
		sys.exit() 
	out = electronic_step_convergence_check(dir_diel)
	while out == 1:
		make_amp2_log(dir_diel,'Calculation options are changed. New calculation starts.')
		out = run_vasp(dir_diel,nproc,vasprun,mpi)
		if out == 1:  # error in vasp calculation
			print(0)
			sys.exit()
		out = electronic_step_convergence_check(dir_diel)

	if out == 2:  # electronic step is not converged. (algo = normal)
		make_amp2_log(dir_diel,'The calculation stops but electronic step is not converged.')
		print(0)
		sys.exit()
	write_diel_log(dir_diel+'/OUTCAR',dir_diel)
	make_amp2_log(dir_diel,'Dielectric calculation is done.')

with open(dir_diel+'/amp2.log','r') as amp2_log:
	with open(dir+'/amp2.log','a') as amp2_log_tot:
		amp2_log_tot.write(amp2_log.read())

print(1)
