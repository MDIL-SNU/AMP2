###########################################
### Date: 2019-04-08			###
### yybbyb@snu.ac.kr			###
###########################################
import math, glob, os, shutil, sys, subprocess
import yaml
from module_vasprun import wincar,incar_from_yaml
from module_log import *
from module_amp2_input import *
code_data = 'Version 0.9.1. Modified at 2019-11-12'

# set file path from input file
inp_file = sys.argv[1]
with open(inp_file,'r') as f:
	inp_yaml = yaml.safe_load(f)
Output_path = inp_yaml['directory']['output']
ERROR_path = inp_yaml['directory']['error']
src_path = inp_yaml['directory']['src_path']
pot_path_gga = inp_yaml['directory']['pot_path_gga']
pot_path_lda = inp_yaml['directory']['pot_path_lda']

node = node_simple(sys.argv[2])

target_mat = sys.argv[3]
target_idx = sys.argv[4]

# In this step, if user give a specific material as Submit_path, Submit_path is not given directory.
# Instead, we set to be the directory where the file is being as Submit_path.
Submit_path = '/'.join(target_mat.split('/')[:-1])

if target_idx == '0':	# directories
	title = target_mat.split('/')[-1]
	target = Output_path+'/'+title
	shutil.move(target_mat,target)
	make_amp2_log_default(target,src_path,'AMP2 continue',node,code_data)

elif target_idx == '1': # cifs
	title = target_mat.split('/')[-1][:-4]
	target = Output_path+'/'+title
	subprocess.call(['mkdir',target])		#make directory	
	subprocess.call(['mkdir',target+'/INPUT0'])
	make_amp2_log_default(target+'/INPUT0',src_path,'AMP2 from cif',node,code_data)
	out_mk_poscar = make_poscar_from_cif(target_mat,target)
	if out_mk_poscar == 1:
		shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
		print 0
		sys.exit()
	[axis,atom_pos] = read_poscar(target+'/INPUT0/POSCAR_conv')
	[prim_axis,prim_atom_pos,sym] = get_primitive_cell(axis,atom_pos)
	write_file(str(sym),target+'/INPUT0/sym')
	write_poscar(prim_axis,prim_atom_pos,target+'/INPUT0/POSCAR','Primitive Cell')
	out_mk_potcar = make_potcar(target+'/INPUT0/POSCAR',pot_path_gga,pot_path_lda,target,src_path,inp_yaml['cif2vasp']['pot_name'])
	if out_mk_potcar == 1:
		shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
		print 0
		sys.exit()
	out_mk_in_note = make_incar_note(target+'/INPUT0/POSCAR',target,inp_yaml['cif2vasp']['soc_target'],inp_yaml['cif2vasp']['u_value'],inp_yaml['cif2vasp']['magmom'],src_path)
	make_incar(target+'/INPUT0/POSCAR',target,src_path,inp_yaml['cif2vasp']['max_nelm'])

	# if KPOINTS, POTCAR, INCAR are provided.
	submit_list = ['KPOINTS','POTCAR_GGA','POTCAR_LDA','INCAR']
	for sub_file in submit_list:
		if os.path.isfile(Submit_path+'/'+sub_file+'_'+title):
			shutil.move(Submit_path+'/'+sub_file+'_'+title, target+'/INPUT0/'+sub_file)

	incar_from_yaml(target+'/INPUT0',inp_yaml['cif2vasp']['incar'])
	shutil.move(target_mat, target+'/'+target_mat.split('/')[-1])

else:	## For POSCAR type input
	poscar_error = poscar_error_check(target_mat,ERROR_path)
	if poscar_error == 1:
		print 0
		sys.exit()
	poscar = target_mat.split('/')[-1].split('_')
	title = '_'.join(poscar[1:])
	target = Output_path+'/'+title
	subprocess.call(['mkdir',target])		#make directory	
	subprocess.call(['mkdir',target+'/INPUT0'])
	make_amp2_log_default(target+'/INPUT0',src_path,'AMP2 from POSCAR',node,code_data)
	shutil.copyfile(target_mat, target+'/INPUT0/POSCAR_conv')

	[axis,atom_pos] = read_poscar(target+'/INPUT0/POSCAR_conv')
	[prim_axis,prim_atom_pos,sym] = get_primitive_cell(axis,atom_pos)
	write_file(str(sym),target+'/INPUT0/sym')
	write_poscar(prim_axis,prim_atom_pos,target+'/INPUT0/POSCAR','Primitive Cell')

	out_mk_potcar = make_potcar(target+'/INPUT0/POSCAR',pot_path_gga,pot_path_lda,target,src_path,inp_yaml['cif2vasp']['pot_name'])
	if out_mk_potcar == 1:
		shutil.move(target,ERROR_path+'/'+target.split('/')[-1])
		print 0
		sys.exit()
	out_mk_in_note = make_incar_note(target+'/INPUT0/POSCAR',target,inp_yaml['cif2vasp']['soc_target'],inp_yaml['cif2vasp']['u_value'],inp_yaml['cif2vasp']['magmom'],src_path)
	make_incar(target+'/INPUT0/POSCAR',target,src_path,inp_yaml['cif2vasp']['max_nelm'])

	# if KPOINTS, POTCAR, INCAR are provided.
	submit_list = ['KPOINTS','POTCAR_GGA','POTCAR_LDA','INCAR']
	for sub_file in submit_list:
		if os.path.isfile(Submit_path+'/'+sub_file+'_'+title):
			shutil.move(Submit_path+'/'+sub_file+'_'+title, target+'/INPUT0/'+sub_file)

	incar_from_yaml(target+'/INPUT0',inp_yaml['cif2vasp']['incar'])
	shutil.move(target_mat, target+'/'+target_mat.split('/')[-1])

if not target_idx == '0':
	with open(target+'/INPUT0/amp2.log','r') as amp2_log:
		with open(target+'/amp2.log','a') as amp2_log_tot:
			amp2_log_tot.write(amp2_log.read())

print target

