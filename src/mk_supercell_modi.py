# This is a code to build supercell for the Ising coefficient
import os, math, sys, yaml
import numpy as np
import module_subr
from itertools import permutations
from itertools import combinations
import operator
from time import strftime
t_start = strftime("%y%m%d-%H%M%S")
try:
   pos_file = str(sys.argv[1])
except:
   sys.exit()

len_min = 5.0
out_pos_file = sys.argv[2]
len_iter = float(sys.argv[3])
latt_diff = 1.0

len_crit = len_min+len_iter+latt_diff
min_natom = 16    # number of atoms in  primitive cell
max_angle = 100
min_angle = 60
# Using relaxed position file 
poslines = module_subr.read_file(pos_file)

# Read lattice vector & atom position
lat_prim=module_subr.read_lat(poslines)
atnames= poslines[5].split()
totnum= sum([int(x) for x in poslines[6].split()])
at_prim = module_subr.read_at_pos(poslines,totnum)
vol_prim = np.inner(lat_prim[0],np.cross(lat_prim[1],lat_prim[2]))

# Make combination of rotation vector ex) [-2,1,2]
m = 1
len_check = 1
tmp = []
bf_list = []
# only min_length check
while len_check:
    rot = module_subr.mk_combination(m)
    if m!=1:
        rot_small = module_subr.mk_combination(m-1)
        rot_new = [x for x in rot if x not in rot_small] 
        bf_list = new_list
        rot = rot_new
    
    # Remove duplicate or all element < 0
    rot_sort = rot
    # Check length criterion
    [new_list,min_sw] = module_subr.check_length(rot_sort, lat_prim, len_min, len_crit)
    if (len(bf_list) !=0) and (len(new_list) == 0):  #until do not append new list, it takes loop infinitely
        len_check = 0
    else:
        if min_sw == 0:
            print(1)
            sys.exit()
        m = m+1

    tmp = tmp+new_list
    

new_list = tmp
# make rotation matrix using basis
# format : [[rot_mat],[lattice_mat]]
basis_mat = module_subr.mk_rot_matrix(new_list)

# remove det = 0
basis_mat = module_subr.rm_det_0(basis_mat)
# remove too sharp shape structures
basis_mat = module_subr.check_shape(basis_mat, min_angle, max_angle)

# Check the number of atoms:
basis_mat = module_subr.check_natom(basis_mat, lat_prim, totnum, min_natom)

# Check defect distance:
basis_mat = module_subr.check_defect_length(basis_mat, len_min)    

# Check pair:
#cri_natom = 100
#best_set = []
#lattice diff check and tot number check
#for i in range(len(basis_mat)):
#    lat_mat = basis_mat[i][1]
#    length = module_subr.cal_lat_len(lat_mat)
#    tot_atom = module_subr.cal_volume(lat_mat)/vol_prim*totnum

#    if (max(length)-min(length))/max(length) < latt_diff:
#        if tot_atom < cri_natom :
#            best_set.append([basis_mat[i][0],basis_mat[i][1]])


#find small natoms
#basis_mat = best_set
target = sys.argv[2].split('POSCAR_param')[0]
with open (target+'tmp_mag_list.dat','r') as f1:
    mag_atom_list = f1.readlines()[0].split()

check_set = module_subr.check_pair(basis_mat,at_prim,mag_atom_list) 

best_natom = 10e5
best = 'Null'
best_set = []


for i in range(len(check_set)):
    lat_mat = check_set[i][1]
    length = module_subr.cal_lat_len(lat_mat)
    tot_atom = module_subr.cal_volume(lat_mat)/vol_prim*totnum

    if tot_atom < best_natom:
        best_natom = tot_atom
        best_set = [[check_set[i][0],check_set[i][1]]]
    elif tot_atom == best_natom:
        best_set.append([check_set[i][0],check_set[i][1]])
#best_set = check_set


#best_set = module_subr.check_pair(best_set,at_prim,mag_atom_list) 
best_latt_sum = -10e5

for i in range(len(best_set)):
    latt_sum = 0
    for j in range(3):
        for k in range(3):
            latt_sum = latt_sum + best_set[i][1][j][k]

    if latt_sum > best_latt_sum:
        best_latt_sum = latt_sum
        best = best_set[i]
det = int(round(best_natom))//int(totnum)
# Set atomic position
if best == 'Null':
    print(1)
    sys.exit()
else:
    at = module_subr.set_atom(at_prim, best[0])
lat_mat = best[1]
# Create POSCAR files (GGA/HSE)
with open(out_pos_file,'w') as f1:
    f1.write('Clean supercell of '+str(atnames)+'\n')
    f1.write('   1.00000000000000\n')
    f1.write("{:21.13f}{:19.13f}{:19.13f}\n".format(lat_mat[0][0],lat_mat[0][1],lat_mat[0][2]))
    f1.write("{:21.13f}{:19.13f}{:19.13f}\n".format(lat_mat[1][0],lat_mat[1][1],lat_mat[1][2]))
    f1.write("{:21.13f}{:19.13f}{:19.13f}\n".format(lat_mat[2][0],lat_mat[2][1],lat_mat[2][2]))
    f1.write(poslines[5])
    f1.write('   '+'   '.join([str(int(x)*det) for x in poslines[6].split()])+'\n')
    f1.write('Selective dynamics\nDirect\n')
    for t in range(len(at)):
        f1.write("{:19.15f}{:19.15f}{:19.15f}   T   T   T\t! {}\n".format(at[t][0],at[t][1],at[t][2],at[t][3]))

t_end = strftime("%y%m%d-%H%M%S")
print(0)
