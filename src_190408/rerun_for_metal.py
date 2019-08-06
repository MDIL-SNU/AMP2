###########################################
### Date: 2019-01-23			###
### yybbyb@snu.ac.kr			###
###########################################
import os, sys, subprocess, yaml, shutil, glob
from input_conf import input_conf
from module_amp2_input import *
from module_log import *
code_data = 'Version xx, modified at 2019-07-16'

# input from shell
inp_file = sys.argv[1]
node = sys.argv[2]
nproc = sys.argv[3]
target = sys.argv[4]

with open(inp_file,'r') as f:
	inp_yaml = yaml.load(f)
cal_dic = inp_yaml['calculation']
src_path = inp_yaml['directory']['src_path']
ERROR_path = inp_yaml['directory']['error']
Done_path = inp_yaml['directory']['done']
large_off = inp_yaml['calculation']['large_off']

# Check metal in HSE
if os.path.isfile(target+'/HSE/Band_gap.log'):
	with open(target+'/HSE/Band_gap.log','r') as inp:
		gap_log = inp.readline()
# Check metal in GGA or LDA
elif os.path.isfile(target+'/band_GGA/Band_gap.log'):
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

# Backup +U calculation results
backup_path = target+'/'+target.split('/')[-1]+'_with_U'
if not os.path.isdir(backup_path):
	os.mkdir(backup_path,0755)
shutil.copytree(target+'/INPUT0',backup_path+'/INPUT0')
path_list = glob.glob(target+'/*')
for path1 in path_list:
	if os.path.isdir(path1) and not path1.split('/')[-1] == 'INPUT0' and not path1 == backup_path:
		shutil.move(path1,backup_path+'/'+path1.split('/')[-1])

make_amp2_log_default(target,src_path,'Checking metal with U.',node,code_data)
make_amp2_log(target,'It is metallic system. We rerun the calculation without U.')
# Modify U_note
with open(target+'/INPUT0/U_note','w') as out:
	out.write(' ')
# Modify INCAR
make_incar(target+'/INPUT0/POSCAR',target,src_path,inp_yaml['cif2vasp']['max_nelm'])

calc_out = 0

if cal_dic['kp_test'] == 1:
#	subprocess.call(['python',src_path+'/kpoint.py',target,inp_file,node,nproc])
	notice = subprocess.check_output(['python',src_path+'/kpoint.py',target,inp_file,node,nproc])
	if not notice.splitlines()[-1][0] == '1':
		shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
		sys.exit()
# check existance of follow calculation and K-pts file
if 1 in cal_dic.values() and not os.path.isfile(target+'/INPUT0/KPOINTS'):
	make_amp2_log(target,'Warning!!! KPOINTS file should be located in INPUT0.')
	shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
	sys.exit()
if cal_dic['encut_test'] == 1:
#	subprocess.call(['python',src_path+'/cutoff.py',target,inp_file,node,nproc])
	notice = subprocess.check_output(['python',src_path+'/cutoff.py',target,inp_file,node,nproc])
	if not notice.splitlines()[-1][0] == '1':
		shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
		sys.exit()
if cal_dic['relaxation'] == 1:
	pot_type = inp_yaml['magnetic_ordering']['potential_type']
	if pot_type in inp_yaml['relaxation']['potential_type']:
#		subprocess.call(['python',src_path+'/relax.py',target,inp_file,node,nproc,pot_type])
		if not os.path.isfile(target+'/INPUT0/POTCAR_'+pot_type):
			make_amp2_log(target,'POTCAR_'+pot_type+' file is missing.')
			sys.exit()
		notice = subprocess.check_output(['python',src_path+'/relax.py',target,inp_file,node,nproc,pot_type])
		if not notice.splitlines()[-1][0] == '1':
			shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
			sys.exit()
if cal_dic['magnetic_ordering'] == 1:
#	subprocess.call(['python',src_path+'/relax.py',target,inp_file,node,nproc,pot_type])
	if not os.path.isfile(target+'/INPUT0/KPOINTS'):
		sys.exit()
	notice = subprocess.check_output(['python',src_path+'/magnetic_ordering.py',target,inp_file,node,nproc])
	if not notice.splitlines()[-1][0] == '1':
		shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
		sys.exit()
if cal_dic['relaxation'] == 1:
	for pot_type in inp_yaml['relaxation']['potential_type']:
#		subprocess.call(['python',src_path+'/relax.py',target,inp_file,node,nproc,pot_type])
		if not pot_type == inp_yaml['magnetic_ordering']['potential_type']:
			if not os.path.isfile(target+'/INPUT0/POTCAR_'+pot_type):
				make_amp2_log(target,'POTCAR_'+pot_type+' file is missing.')
				sys.exit()
			notice = subprocess.check_output(['python',src_path+'/relax.py',target,inp_file,node,nproc,pot_type])
			if not notice.splitlines()[-1][0] == '1':
				shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
				sys.exit()
if cal_dic['band'] == 1:
	for pot_type in inp_yaml['band_calculation']['potential_type']:
#		subprocess.call(['python',src_path+'/band.py',target,inp_file,node,nproc,pot_type])
		if not os.path.isfile(target+'/INPUT0/POTCAR_'+pot_type):
			make_amp2_log(target,'POTCAR_'+pot_type+' file is missing.')
			continue
		notice = subprocess.check_output(['python',src_path+'/band.py',target,inp_file,node,nproc,pot_type])
		if not notice.splitlines()[-1][0] == '1':
			calc_out = 1
			break
	if calc_out == 1:
		calc_out = 0
		shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
		sys.exit()
if cal_dic['density_of_states'] == 1:
	for pot_type in inp_yaml['density_of_states']['potential_type']:
		if not os.path.isfile(target+'/INPUT0/POTCAR_'+pot_type):
			make_amp2_log(target,'POTCAR_'+pot_type+' file is missing.')
			continue
		notice = subprocess.check_output(['python',src_path+'/dos.py',target,inp_file,node,nproc,pot_type])
		if not notice.splitlines()[-1][0] == '1':
			calc_out = 1
			break
	if calc_out == 1:
		calc_out = 0
		shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
		sys.exit()
if cal_dic['dielectric'] == 1:
	for pot_type in inp_yaml['dielectric']['potential_type']:
		if not os.path.isfile(target+'/INPUT0/POTCAR_'+pot_type):
			make_amp2_log(target,'POTCAR_'+pot_type+' file is missing.')
			continue
#		subprocess.call(['python',src_path+'/dielectric.py',target,inp_file,node,nproc,pot_type])
		notice = subprocess.check_output(['python',src_path+'/dielectric.py',target,inp_file,node,nproc,pot_type])
		if not notice.splitlines()[-1][0] == '1':
			calc_out = 1
			break
	if calc_out == 1:
		calc_out = 0
		shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
		sys.exit()
if cal_dic['hse_oneshot'] == 1:
	for pot_type in inp_yaml['hybrid_oneshot']['potential_type']:
		notice = subprocess.check_output(['python',src_path+'/hse_gap.py',target,inp_file,node,nproc,pot_type])
		if not notice.splitlines()[-1][0] == '1':
			shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
			sys.exit()

