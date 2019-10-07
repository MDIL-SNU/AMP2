###########################################
### Date: 2018-12-05			###
### yybbyb@snu.ac.kr			###
###########################################
from module_log import *
from module_vasprun import *
from module_amp2_input import read_poscar
from module_vector import *

# set nsw for relax
def set_nsw(poscar,incar):
	pos = open(poscar,'r').readlines()[6].split()
	natom = 0
	for i in range(len(pos)) :
		natom = natom + int(pos[i])
	nsw = str(natom*3)
	if nsw < 20:
		nsw = 20
	nsw = 100
	wincar(incar,incar,[['NSW',nsw]],[])
	return int(nsw)

# Lattice mismatch check between relaxed poscar and original poscar
def set_pos_compare(poscar_rlx,poscar_ori,target,err_percent):
	axis_rlx = poscar_to_axis(poscar_rlx)
	axis_ori = poscar_to_axis(poscar_ori)

	min_diff = 0.01

	warning = 0
	for i in range(3):
		for j in range(3):
			if abs(axis_rlx[i][j]-axis_ori[i][j]) > min_diff and abs((axis_rlx[i][j]-axis_ori[i][j])/(axis_ori[i][j]+0.001)) > float(err_percent/100):
				warning = 1
	if warning == 1:
		make_amp2_log(target,'Warning!! More than One of lattice parameter changes more than 1 %.')
	return warning

def set_ediffg(poscar,force,pressure,energy):
	if pressure < 0 and force < 0 and energy < 0:
		return 0
	elif pressure < 0 and force < 0 and enegy > 0:
		return energy
	else:
		[axis,atom_pos] = read_poscar(poscar)
		volume = calc_volume(axis)
		nion = len(atom_pos)
		conv_press = pressure*volume/nion/1602.17733 # it is identical to VASP precision.
		if conv_press < force:
			return -1.0*conv_press
		else:
			return -1.0*force
	
