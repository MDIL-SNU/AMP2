###########################################
### Date: 2018-12-05			###
### yybbyb@snu.ac.kr			###
###########################################
import shutil, os, sys, subprocess, yaml
from module_log import *
from module_vasprun import *
code_data = 'Date: 2018-12-05'

dir = sys.argv[1]

inp_file = sys.argv[2]
with open(inp_file,'r') as f:
	inp_yaml = yaml.load(f)
ERROR_path = inp_yaml['directory']['error']
src_path = inp_yaml['directory']['src_path']
vasp_std = inp_yaml['program']['vasp_std']
vasp_gam = inp_yaml['program']['vasp_gam']
vasp_ncl = inp_yaml['program']['vasp_ncl']
mpi = inp_yaml['program']['mpi_command']
gnuplot = inp_yaml['program']['gnuplot']
npar = inp_yaml['vasp_parallel']['npar']
kpar = inp_yaml['vasp_parallel']['kpar']
inp_relax = inp_yaml['relaxation']
inp_band = inp_yaml['band_calculation']

node = node_simple(sys.argv[3])
nproc = sys.argv[4]

POT = sys.argv[5]

# Set directory for input structure and INCAR
dir_chg = dir+'/Optimized_'+POT

if not os.path.isdir(dir_chg):
	os.mkdir(dir_chg,0755)

os.chdir(dir_chg)

make_amp2_log_default(dir_chg,src_path,'Charge density calculation from optimized setup',node,code_data)

# check relax calculation
no_rlx = 0
if not os.path.isdir(dir+'/relax_'+POT):
	make_amp2_log(dir_band,'Relax directory does not exist.')
	no_rlx = 1
else:
	if not os.path.isfile(dir+'/relax_'+POT+'/CONTCAR'):
		make_amp2_log(dir_band,'CONTCAR file in relaxation does not exist.')
		no_rlx = 1
	elif count_line(dir+'/relax_'+POT+'/CONTCAR') < 9:
		make_amp2_log(dir_band,'CONTCAR file in relaxation is invalid.')
		no_rlx = 1

if no_rlx == 1 and inp_band['relax_check'] == 1:
	print 0
	sys.exit()

# check CHGCAR
if os.path.isfile(dir_chg+'/CHGCAR') and os.path.getsize(dir_chg+'/CHGCAR') > 0 :
	make_amp2_log(dir_chg,'The calculation was already done.')
else:
	make_amp2_log(dir_chg,'New calculation starts.')
	# Copy input data and write CHGCAR
	if no_rlx == 1:
		copy_input(dir+'/INPUT0',dir_chg,POT)
	else:
		copy_input_cont(dir+'/relax_'+POT,dir_chg)
	subprocess.call(['cp',dir+'/INPUT0/sym',dir_chg+'/.'])

	# make INCAR for CHGCAR
	incar_from_yaml(dir_chg,inp_relax['incar'])
	if no_rlx == 1:
		mag_on = 2
	else:
		mag_on = check_magnet(dir+'/relax_'+POT,inp_yaml['magnetic_ordering']['minimum_moment'])
	vasprun = make_incar_for_ncl(dir_chg,mag_on,kpar,npar,vasp_std,vasp_gam,vasp_ncl)
	wincar(dir_chg+'/INCAR',dir_chg+'/INCAR',[['NSW','0'],['LCHARG','.T.']],[])

	# VASP calculation for CHGCAR
	out = run_vasp(dir_chg,nproc,vasprun,mpi)
	if out == 1:  # error in vasp calculation
		print 0
		sys.exit() 
	make_amp2_log(dir_chg,'CHGCAR file is generated successfully.')


