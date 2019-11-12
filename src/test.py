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

#target = '/data/opusg514/AMP2_snumat/ERROR/Au2Nb1_58558/relax_GGA'

target = '/data/opusg514/AMP2_snumat/ERROR/rlx_pos2cont/Ag1Ce1_57368/relax_GGA'
target2 = '/data/opusg514/AMP2_snumat/ERROR/rlx_pos2cont/Ag1Ce1_57368/kptest/KP1'
a = 1
def test_elec (target):
    # 0: well converged, 1: can try to calculation with other setting (retry), 2: cannot be converged (stop)
    with open(target+'/OSZICAR','r') as inp:
        fr_log = inp.readlines()[1:]
    spin = pygrep('ISPIN',target+'/OUTCAR',0,0).split()[2]
#   spin = subprocess.check_output(['grep','ISPIN',target+'/OUTCAR']).split()[2]
    elec_step = 0
    for ll in fr_log:
        if ll.split()[0] == '1':
            break
        elec_step = elec_step+1
    print elec_step
    if elec_step == int(pygrep('NELM',target+'/OUTCAR',0,0).split(';')[0].split()[2]):
#   if elec_step == int(subprocess.check_output(['grep','NELM',target+'/OUTCAR']).split(';')[0].split()[2]):
#        make_amp2_log(target,'Electronic step is not converged.')
        algo = pygrep('ALGO',target+'/INCAR',0,0).split()[2]
#       algo = subprocess.check_output(['grep','ALGO',target+'/INCAR']).split()[2]
        if algo == 'Normal' or algo == 'All':
#            make_amp2_log(target,'Current ALGO is '+algo+' but it is not converged.')
            if spin == '2':
                if pygrep('BMIX_MAG',target+'/OUTCAR',0,0).split()[-1] == '1.00':
#               if subprocess.check_output(['grep','BMIX_MAG',target+'/OUTCAR']).split()[-1] == '1.00':
#                    make_amp2_log(target,'Change mixing parameter.')
#                    wincar(target+'/INCAR',target+'/INCAR',[['AMIX','0.2'],['BMIX','0.0001'],['AMIX_MAG','0.8'],['BMIX_MAG','0.0001']],[])
                    return 1
                else:
#                    make_amp2_log(target,'We changed mixing parameters but it is not converged.')
                    return 2
            else:
                return 2
        elif algo == 'Damped':
#            wincar(target+'/INCAR',target+'/INCAR',[['ALGO','All']],[])
#            make_amp2_log(target,'ALGO changes from '+algo+' to All.')
            return 1
        else:
#            wincar(target+'/INCAR',target+'/INCAR',[['ALGO','Normal']],[])
#            make_amp2_log(target,'ALGO changes from '+algo+' to Normal.')
            return 1
    else:
        return 0


print test_elec (target)
