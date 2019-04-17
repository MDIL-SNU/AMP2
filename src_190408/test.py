import os, sys, subprocess, yaml,shutil,math
from module_log import *
from module_vasprun import *
from module_converge import *
from module_relax import *
from module_band import *
from module_effm import *
from module_AF import *
from module_amp2_input import *

#target_mat = sys.argv[1]
#Output_path = sys.argv[2]
#title = target_mat.split('/')[-1][:-4]
#target = Output_path+'/'+title
#subprocess.call(['mkdir',target])		#make directory	
#subprocess.call(['mkdir',target+'/INPUT0'])
#out_mk_poscar = make_poscar_from_cif(target_mat,target)
#[axis,atom_pos] = read_poscar(target+'/INPUT0/POSCAR_conv')
#[prim_axis,prim_atom_pos,sym] = get_primitive_cell(axis,atom_pos)
#write_poscar(prim_axis,prim_atom_pos,target+'/INPUT0/POSCAR_prim_pp','Primitive Cell')


[A,B,C] = find_pair(sys.argv[1],5,['V1','V2','V3','V4','V5','V6'],{'V1':3,'V2':3,'V3':3,'V4':3,'V5':3,'V6':3},0.01)
for line in A:
	print line
