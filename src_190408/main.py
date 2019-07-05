###########################################
### Date: 2019-01-23			###
### yybbyb@snu.ac.kr			###
###########################################
import os, sys, subprocess, yaml, shutil
from input_conf import input_conf
from module_amp2_input import *
from module_log import *
# input from shell
conf_file = sys.argv[1]
node = sys.argv[2]
nproc = sys.argv[3]

home = os.getcwd()

inp_file = input_conf(conf_file)
with open(inp_file,'r') as f:
	inp_yaml = yaml.load(f)
cal_dic = inp_yaml['calculation']
src_path = inp_yaml['directory']['src_path']
ERROR_path = inp_yaml['directory']['error']
Done_path = inp_yaml['directory']['done']
large_off = inp_yaml['calculation']['large_off']

calc_list = make_list(inp_file)
calc_out = 0

while len(make_list(inp_file)) > 0:
	target_mat = make_list(inp_file)[0]
	target = subprocess.check_output(['python',src_path+'/amp2_input.py',inp_file,node,target_mat[0],target_mat[1]]).splitlines()[0]
	if target == '0':
		shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
		continue
	if large_off > 1:
		num_pos = len(read_poscar(target+'/INPUT0/POSCAR')[1])
		if num_pos > large_off:
			shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
			continue
	if cal_dic['kp_test'] == 1:
#		subprocess.call(['python',src_path+'/kpoint.py',target,inp_file,node,nproc])
		notice = subprocess.check_output(['python',src_path+'/kpoint.py',target,inp_file,node,nproc])
		if not notice.splitlines()[-1][0] == '1':
			shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
			continue
	# check existance of follow calculation and K-pts file
	if 1 in cal_dic.values() and not os.path.isfile(target+'/INPUT0/KPOINTS'):
		make_amp2_log(target,'Warning!!! KPOINTS file should be located in INPUT0.')
		shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
		continue
	if cal_dic['encut_test'] == 1:
#		subprocess.call(['python',src_path+'/cutoff.py',target,inp_file,node,nproc])
		notice = subprocess.check_output(['python',src_path+'/cutoff.py',target,inp_file,node,nproc])
		if not notice.splitlines()[-1][0] == '1':
			shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
			continue
	if cal_dic['relaxation'] == 1:
		if 'GGA' in inp_yaml['relaxation']['potential_type']:
#			subprocess.call(['python',src_path+'/relax.py',target,inp_file,node,nproc,pot_type])
			pot_type = 'GGA'
			if not os.path.isfile(target+'/INPUT0/POTCAR_'+pot_type):
				make_amp2_log(target,'POTCAR_'+pot_type+' file is missing.')
				continue
			notice = subprocess.check_output(['python',src_path+'/relax.py',target,inp_file,node,nproc,pot_type])
			if not notice.splitlines()[-1][0] == '1':
				shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
				continue
	if cal_dic['magnetic_ordering'] == 1:
#		subprocess.call(['python',src_path+'/relax.py',target,inp_file,node,nproc,pot_type])
		if not os.path.isfile(target+'/INPUT0/KPOINTS'):
			continue
		notice = subprocess.check_output(['python',src_path+'/magnetic_ordering.py',target,inp_file,node,nproc])
		if not notice.splitlines()[-1][0] == '1':
			shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
			continue
	if cal_dic['relaxation'] == 1:
		if 'LDA' in inp_yaml['relaxation']['potential_type']:
#			subprocess.call(['python',src_path+'/relax.py',target,inp_file,node,nproc,pot_type])
			pot_type = 'LDA'
			if not os.path.isfile(target+'/INPUT0/POTCAR_'+pot_type):
				make_amp2_log(target,'POTCAR_'+pot_type+' file is missing.')
				continue
			notice = subprocess.check_output(['python',src_path+'/relax.py',target,inp_file,node,nproc,pot_type])
			if not notice.splitlines()[-1][0] == '1':
				shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
				continue
	if cal_dic['band'] == 1:
		for pot_type in inp_yaml['band_calculation']['potential_type']:
#			subprocess.call(['python',src_path+'/band.py',target,inp_file,node,nproc,pot_type])
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
			continue
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
			continue
	if cal_dic['dielectric'] == 1:
		for pot_type in inp_yaml['dielectric']['potential_type']:
			if not os.path.isfile(target+'/INPUT0/POTCAR_'+pot_type):
				make_amp2_log(target,'POTCAR_'+pot_type+' file is missing.')
				continue
#			subprocess.call(['python',src_path+'/dielectric.py',target,inp_file,node,nproc,pot_type])
			notice = subprocess.check_output(['python',src_path+'/dielectric.py',target,inp_file,node,nproc,pot_type])
			if not notice.splitlines()[-1][0] == '1':
				calc_out = 1
				break
		if calc_out == 1:
			calc_out = 0
			shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
			continue
	if cal_dic['hse_oneshot'] == 1:
		notice = subprocess.check_output(['python',src_path+'/hse_gap.py',target,inp_file,node,nproc])
		if not notice.splitlines()[-1][0] == '1':
			shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
			continue
	if cal_dic['effective_mass'] == 1:
		for pot_type in inp_yaml['effective_mass']['potential_type']:
			if not os.path.isfile(target+'/INPUT0/POTCAR_'+pot_type):
				make_amp2_log(target,'POTCAR_'+pot_type+' file is missing.')
				continue
			for carrier_type in inp_yaml['effective_mass']['carrier_type']:
#				subprocess.call(['python',src_path+'/effm.py',target,inp_file,node,nproc,pot_type,carrier_type])
				notice = subprocess.check_output(['python',src_path+'/effm.py',target,inp_file,node,nproc,pot_type,carrier_type])
				if not notice.splitlines()[-1][0] == '1':
					calc_out = 1
					break
			if calc_out == 1:
				break
		if calc_out == 1:
			calc_out = 0
			shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
			continue

	shutil.move(target,Done_path+'/'+target.split('/')[-1])
