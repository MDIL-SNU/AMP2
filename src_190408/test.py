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
from input_conf import *
import spglib
import numpy as np
import math
# input from shell
conf = sys.argv[1]
with open(conf,'r') as inp:
	inp_yaml = yaml.load(inp)
print inp_yaml

inp_yaml_low = dic_to_lowercase(inp_yaml)

print inp_yaml_low

