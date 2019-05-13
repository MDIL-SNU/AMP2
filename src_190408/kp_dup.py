import os, sys, subprocess, yaml,shutil,math
from module_log import *
from module_vasprun import *
from module_converge import *
from module_relax import *
from module_band import *
from module_effm import *
from module_AF import *
from module_amp2_input import *
import spglib
import numpy as np

def vec_mod(vec):
	new_vec = []
	for xx in vec:
		new_vec.append(round(xx-np.floor(xx+0.5),7))
	return new_vec

def vec_same(vec1,vec2):
	summ = sum([abs(vec1[x]-vec2[x]) for x in range(len(vec1))])
	if summ < 0.00001:
		return 1
	else:
		return 0

def dir_to_cart(vec,axis):
	cart_vec = []
	for i in range(3):
		cart_vec.append(vec[0]*axis[0][i]+vec[1]*axis[1][i]+vec[2]*axis[2][i])
	return cart_vec


def dup_pos(pos):
	irre_pos = []
	made_pos = []
	idx = 0
	dic = {}

	for i in range(len(pos)):
		for j in range(len(oper)):
			made_pos.append(vec_mod(dir_to_cart(pos[i],oper[j])))
	for i in range(len(made_pos)):
		if len(irre_pos) == 0:
			irre_pos.append(made_pos[i])
		else:
			sw = 0
			for j in range(len(irre_pos)):
				if vec_same(made_pos[i],irre_pos[j]):
					sw = 1
					break
			if sw == 0:
				irre_pos.append(made_pos[i])
	return irre_pos


pocket_path = sys.argv[2]


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

rot_oper = spglib.get_symmetry(Cell)['rotations']
oper = np.ndarray.tolist(rot_oper)

with open(pocket_path,'r') as f:
	lines = f.readlines()

poc_kp = []
for line in lines:
	poc_kp.append([[float(x) for x in line.split()[0:3]]])

weight = []
for poc in poc_kp:
	prev_poc_len = len(poc)
	while 1:
		prev_poc_len = len(poc)
		poc = dup_pos(poc)
		if len(poc) == prev_poc_len:
			break
		else:
			prev_poc_len = len(poc)
	print poc
	print '-------------------'

	weight.append(len(poc))

print len(poc)
