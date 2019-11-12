from module_dos import *
import sys, subprocess

dir_dos = sys.argv[1]
spin = subprocess.check_output(['grep','ISPIN',dir_dos+'/OUTCAR']).split()[2]
ncl = subprocess.check_output(['grep','NONCOL',dir_dos+'/OUTCAR']).split()[2]
# read atom information from poscar
[atom_name,atom_num] = poscar_to_atom_inform(dir_dos+'/POSCAR')
# read dos information from doscar
[Ene,Tot_dos,par_dos] = make_dos_dat(dir_dos+'/DOSCAR',spin,atom_num,ncl)
print ncl
fermi = 5.82327940 
#fermi = 0
# write dos dat files
write_tot_dos(Ene,Tot_dos,fermi,dir_dos)
write_par_dos(Ene,par_dos,atom_name,fermi,dir_dos)
