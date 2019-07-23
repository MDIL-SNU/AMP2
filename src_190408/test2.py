import os, sys, subprocess, yaml
from module_log import *
from module_vasprun import *
from module_converge import *
from module_relax import *
from module_band import *
from module_effm import *
from module_AF import *
import spglib
import numpy as np

def sym_rot(pos,rot):
	pos_new = []
	for rot_tar in rot:
		pos_new.append(sum([pos[x]*rot_tar[x] for x in range(3)]))
	return pos_new


[axis,atom_pos] = read_poscar(sys.argv[1])

L = np.mat(axis)
pos = []
atom_type = []
atom_dic = {}
type_dic = {}
index = 0
for line in atom_pos:
	pos.append(line[0:3])
	if not line[4] in atom_dic.keys():
		atom_dic[line[4]] = index
		type_dic[index] = [line[3],line[4]]
		index = index+1
	atom_type.append(atom_dic[line[4]])
D = np.mat(pos)
Cell = (L,D,atom_type)

info = spglib.get_symmetry(Cell)


pos = [0.2,0.2,0]
pos_cart = dir_to_cart(pos,axis)

sym_pos = []
sym_pos_cart = []
for sym in info['rotations']:
	sym_pos.append(sym_rot(pos,sym)+atom_pos[0][3:])
	sym_pos_cart.append(cart_to_dir(sym_rot(pos_cart,sym),axis)+atom_pos[0][3:])
	print sym_pos[-1],sym_pos_cart[-1]

write_poscar(axis,sym_pos,'POSCAR_sym','test1')
write_poscar(axis,sym_pos_cart,'POSCAR_sym_cart','test2')

#	print sym
#print info['rotations']
