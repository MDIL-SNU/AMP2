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

print pygrep(sys.argv[1],sys.argv[2],0,int(sys.argv[3]))

