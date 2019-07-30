import os, sys, subprocess, yaml,shutil,math
from module_log import *
from module_vasprun import *
from module_converge import *
from module_relax import *
from module_band import *
from module_effm import *
from module_AF import *
from module_amp2_input import *
from module_vector import *
import spglib
import numpy as np
import math
# input from shell
dir = sys.argv[1]
POT = 'GGA'
if os.path.isfile(dir+'/INPUT0/POSCAR_rlx_'+POT):
	oper = read_operation(dir+'/INPUT0/POSCAR_rlx_'+POT)
else:
	oper = read_operation(dir+'/INPUT0/POSCAR')
dir_effm = dir+'/effm_GGA/hole'
[effm_dia,effm] = calc_effm(dir_effm,'hole',300,oper)
print effm_dia,effm

