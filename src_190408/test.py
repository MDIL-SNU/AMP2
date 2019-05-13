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

print info['rotations']
