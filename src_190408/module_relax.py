###########################################
### Date: 2018-12-05			###
### yybbyb@snu.ac.kr			###
###########################################
from module_log import *
from module_vasprun import *

# set nsw for relax
def set_nsw(poscar,incar):
	pos = open(poscar,'r').readlines()[6].split()
	natom = 0
	for i in range(len(pos)) :
		natom = natom + int(pos[i])
	nsw = str(natom*3)
	if nsw < 20:
		nsw = 20
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
