import os, sys, subprocess, yaml
from module_log import *
from module_vasprun import *
from module_converge import *
from module_relax import *
from module_band import *
from module_effm import *
from module_AF import *

dir_band = sys.argv[1]

# Extract data from VASP output files
fermi = float(subprocess.check_output(['head',dir_band+'/DOSCAR','-n','6']).splitlines()[-1].split()[3])
spin = subprocess.check_output(['grep','ISPIN',dir_band+'/OUTCAR']).split()[2]
ncl = subprocess.check_output(['grep','NONCOL',dir_band+'/OUTCAR']).split()[2]
[KPT,Band,nelect] = EIGEN_to_array(dir_band+'/EIGENVAL',spin)

# Error checking for deviated eigen value
warn = band_warning(Band,dir_band)

# Band gap calculation
gap = gap_estimation(dir_band,fermi,spin,ncl,KPT,Band,nelect) # gap is string

make_amp2_log(dir_band,'Band calculaton is done.\nBand gap is '+gap)

# Drawing band structure
plot_band_structure(spin,Band,fermi,dir_band+'/xtic.dat',dir_band+'/xlabel.dat',-3,3,dir_band)

os.chdir(dir_band)
subprocess.call([gnuplot,dir_band+'/band.in'])


#if mat1 == mat2:
#	print 'equal'
#else:
#	print 'non equal'
#for i in range(len(mat1)):
#	for j in range(len(mat1[0])):
#		print mat1[i][j],'     ',mat2[i][j]

#[axis,atom_pos] = read_poscar(sys.argv[1])

#[prim_axis,prim_atom_pos] = get_primitive_cell(axis,atom_pos)

#write_poscar(prim_axis,prim_atom_pos,sys.argv[2],'Primitive Cell')

#print prim_axis
#print prim_atom_pos

#write_poscar(prim_axis,prim_atom_pos,sys.argv[2],'Primitive cell')
#[mag_atom_list,mag_val] = check_spin(sys.argv[1],sys.argv[2])
#print mag_atom_list
#mag_atom_list = ['Ni1']
#poscar = sys.argv[1]
#list1 = set_supercell_list(poscar,mag_atom_list)
#print list1
