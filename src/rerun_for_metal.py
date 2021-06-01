###########################################
### Date: 2019-01-23			###
### yybbyb@snu.ac.kr			###
###########################################
# This is a code to restart the all calculations without the on-site U term
# if the material was found to be metallic and U was applied.
import os, sys, subprocess, yaml, shutil, glob
from input_conf import input_conf
from module_amp2_input import *
from module_log import *
from module_band import check_half_metal
from input_conf import set_on_off
from _version import __version__
code_data = 'Version '+__version__+'. Modified at 2020-01-15'

# input from shell
inp_file = sys.argv[1]
node = sys.argv[2]
nproc = sys.argv[3]
target = sys.argv[4]
pypath = sys.executable

with open(inp_file,'r') as f:
	inp_yaml = yaml.safe_load(f)
cal_dic = inp_yaml['calculation']
src_path = inp_yaml['directory']['src_path']
ERROR_path = inp_yaml['directory']['error']
Done_path = inp_yaml['directory']['done']

# Check metal in HSE
for pot_type in inp_yaml['hybrid_oneshot']['potential_type']:
	if isinstance(pot_type,list):
		if len(pot_type) == 1:
			pot_cell = pot_type[0]
			pot_point = pot_type[0]
		else:
			pot_cell = pot_type[0]
			pot_point = pot_type[1]
	else:
		pot_cell = pot_type
		pot_point = pot_type
	if os.path.isfile(target+'/hybrid'+pot_cell+'_'+pot_point+'/Band_gap.log'):
		with open(target+'/hybrid'+pot_cell+'_'+pot_point+'/Band_gap.log','r') as inp:
			gap_log = inp.readline()
		if not 'etal' in gap_log:
			sys.exit()
# Check metal in GGA or LDA
if os.path.isfile(target+'/band_GGA/Band_gap.log'):
	with open(target+'/band_GGA/Band_gap.log','r') as inp:
		gap_log = inp.readline()
elif os.path.isfile(target+'/band_LDA/Band_gap.log'):
	with open(target+'/band_LDA/Band_gap.log','r') as inp:
		gap_log = inp.readline()
else:
	sys.exit()

if not 'etal' in gap_log:
	sys.exit()

# Check U
U_on = 0
with open(target+'/INPUT0/INCAR') as inp:
	lines = inp.readlines()
	for line in lines:
		if len(line.split()) > 0 and line.split()[0] == 'LDAU':
			if 'T' in line.split()[2]:
				U_on = 1
if U_on == 0 :
	sys.exit()

# Check half metal
if os.path.isfile(target+'/band_GGA/Band_gap.log'):
	if check_half_metal(target+'/band_GGA') == 1:
		sys.exit()
elif os.path.isfile(target+'/band_LDA/Band_gap.log'):
	if check_half_metal(target+'/band_LDA') == 1:
		sys.exit()
else:
	sys.exit()

# Backup +U calculation results
backup_path = target+'/'+target.split('/')[-1]+'_with_U'
if not os.path.isdir(backup_path):
	os.mkdir(backup_path,0o755)
path_list = glob.glob(target+'/*')
for path1 in path_list:
	if os.path.isdir(path1) and not path1 == backup_path:
		shutil.move(path1,backup_path+'/'+path1.split('/')[-1])

if os.path.isdir(backup_path+'/INPUT0_old'):
	shutil.copytree(backup_path+'/INPUT0_old',target+'/INPUT0')
else:
	shutil.copytree(backup_path+'/INPUT0',target+'/INPUT0')

make_amp2_log_default(target,src_path,'Checking metal with U.',node,code_data)
make_amp2_log(target,'It is metallic system. We rerun the calculation without U.')
# Modify U_note
with open(target+'/INPUT0/U_note','w') as out:
	out.write(' ')
# Modify INCAR
make_incar(target+'/INPUT0/POSCAR',target,src_path,inp_yaml['cif2vasp']['max_nelm'])

if os.path.isfile(target+'/INPUT0/CHGCAR_conv'): #Calculation with U would not have convergence problem
	os.remove(target+'/INPUT0/CHGCAR_conv')

calc_out = 0

if set_on_off(cal_dic['kp_test']) == 1:
	try:
		notice = subprocess.check_output([pypath,src_path+'/kpoint.py',target,inp_file,node,nproc],universal_newlines=True)
	except:
		notice = '0'
	if not notice.splitlines()[-1][0] == '1':
		shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
		sys.exit()

# check existance of follow calculation and K-pts file
if 1 in list(cal_dic.values()) and not os.path.isfile(target+'/INPUT0/KPOINTS'):
	make_amp2_log(target,'Warning!!! KPOINTS file should be located in INPUT0.')
	shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
	sys.exit()
# cutoff test
if set_on_off(cal_dic['encut_test']) == 1:
	try:
		notice = subprocess.check_output([pypath,src_path+'/cutoff.py',target,inp_file,node,nproc],universal_newlines=True)
	except:
		notice = '0'
	if not notice.splitlines()[-1][0] == '1':
		shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
		sys.exit()
# relaxation
if set_on_off(cal_dic['relaxation']) == 1:
	pot_type = inp_yaml['magnetic_ordering']['potential_type']
	if pot_type in inp_yaml['relaxation']['potential_type']:
		try:
			notice = subprocess.check_output([pypath,src_path+'/relax.py',target,inp_file,node,nproc,pot_type],universal_newlines=True)
		except:
			notice = '0'
		if not notice.splitlines()[-1][0] == '1':
			shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
			sys.exit()
# magnetic ordering
if set_on_off(cal_dic['magnetic_ordering']) == 1:
	if not os.path.isfile(target+'/INPUT0/KPOINTS'):
		sys.exit()
	try:
		notice = subprocess.check_output([pypath,src_path+'/magnetic_ordering.py',target,inp_file,node,nproc],universal_newlines=True)
	except:
		notice = '0'
	if notice.splitlines()[-1][0] == '2':
		shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
		sys.exit()
	elif not notice.splitlines()[-1][0] == '1':
		with open(target+'/magnetic_ordering/amp2.log','r') as amp2_log:
			with open(target+'/amp2.log','a') as amp2_log_tot:
				amp2_log_tot.write(amp2_log.read())
		make_amp2_log(target,'AMP2 failed to identify the most stable magnetic ordering. Ferromagnetic ordering is used.')
# relaxation if pot_type is different from pot_type of magnetic ordering
if set_on_off(cal_dic['relaxation']) == 1:
	for pot_type in inp_yaml['relaxation']['potential_type']:
		try:
			notice = subprocess.check_output([pypath,src_path+'/relax.py',target,inp_file,node,nproc,pot_type],universal_newlines=True)
		except:
			notice = '0'
		if not notice.splitlines()[-1][0] == '1':
			shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
			sys.exit()
# band calculation
if set_on_off(cal_dic['band']) == 1:
	for pot_type in inp_yaml['band_calculation']['potential_type']:
		try:
			notice = subprocess.check_output([pypath,src_path+'/band.py',target,inp_file,node,nproc,pot_type],universal_newlines=True)
		except:
			notice = '0'
		if not notice.splitlines()[-1][0] == '1':
			calc_out = 1
			break
	if calc_out == 1:
		calc_out = 0
		shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
		sys.exit()
# density of states
if set_on_off(cal_dic['density_of_states']) == 1:
	for pot_type in inp_yaml['density_of_states']['potential_type']:
		try:
			notice = subprocess.check_output([pypath,src_path+'/dos.py',target,inp_file,node,nproc,pot_type],universal_newlines=True)
		except:
			notice = '0'
		if not notice.splitlines()[-1][0] == '1':
			calc_out = 1
			break
	if calc_out == 1:
		calc_out = 0
		shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
		sys.exit()
# dielectric constant
if set_on_off(cal_dic['dielectric']) == 1:
	for pot_type in inp_yaml['dielectric']['potential_type']:
		try:
			notice = subprocess.check_output([pypath,src_path+'/dielectric.py',target,inp_file,node,nproc,pot_type],universal_newlines=True)
		except:
			notice = '0'
		if not notice.splitlines()[-1][0] == '1':
			calc_out = 1
			break
	if calc_out == 1:
		calc_out = 0
		shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
		sys.exit()
# HSE oneshot
if set_on_off(cal_dic['hse_oneshot']) == 1:
	for pot_type in inp_yaml['hybrid_oneshot']['potential_type']:
		if isinstance(pot_type,list):
			if len(pot_type) == 1:
				pot_cell = pot_type[0]
				pot_point = pot_type[0]
			else:
				pot_cell = pot_type[0]
				pot_point = pot_type[1]
		else:
			pot_cell = pot_type
			pot_point = pot_type
		try:
			notice = subprocess.check_output([pypath,src_path+'/hse_gap.py',target,inp_file,node,nproc,pot_cell,pot_point],universal_newlines=True)
		except:
			notice = '0'
		if not notice.splitlines()[-1][0] == '1':
			shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
			sys.exit()
