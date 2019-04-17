import os,sys,shutil,subprocess,glob,math
from module_amp2_input import *

cifs = glob.glob(sys.argv[1]+'/*.cif')
Output_path = sys.argv[2]
for target_mat in cifs:
	title = target_mat.split('/')[-1][:-4]
	target = Output_path+'/'+title
	subprocess.call(['mkdir',target])		#make directory	
	subprocess.call(['mkdir',target+'/INPUT0'])
	out_mk_poscar = make_poscar_from_cif(target_mat,target)
	[axis,atom_pos] = read_poscar(target+'/INPUT0/POSCAR_conv')
	[prim_axis,prim_atom_pos,sym] = get_primitive_cell(axis,atom_pos)
	with open(target+'/INPUT0/sym','w') as out:
		out.write(str(sym))
	write_poscar(prim_axis,prim_atom_pos,target+'/INPUT0/POSCAR','Primitive Cell')
	shutil.copy2(target_mat,target+'/INPUT0/'+title+'.cif')
