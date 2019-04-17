import os, math, sys, yaml
import numpy as np
import z_subr
from itertools import permutations
from itertools import combinations
import operator
from time import strftime
t_start = strftime("%y%m%d-%H%M%S")
#print(t_start)
try:
   pos_file = str(sys.argv[1])
except:
   sys.exit()

kp_file = sys.argv[2]
len_min = float(sys.argv[3])
out_pos_file = sys.argv[4]

gamma_on = 0

#len_min = 10
latt_diff = 0.5
len_crit = len_min+2+latt_diff
min_natom = 100
#max_natom = 300
max_angle = 115
min_angle = 65
#max_angle = np.pi*0.75
#min_angle = np.pi*0.25
# Using relaxed position file 
poslines = z_subr.read_file(pos_file)

# Read lattice vector & atom position
lat_prim=z_subr.read_lat(poslines)
atnames= poslines[5].split()
totnum= sum([int(x) for x in poslines[6].split()])
#print "# of atoms in unit: "+str(totnum)
at_prim = z_subr.read_at_pos(poslines,totnum)
vol_prim = np.inner(lat_prim[0],np.cross(lat_prim[1],lat_prim[2]))
#print vol_prim

# Make combination of rotation vector ex) [-2,1,2]
m = 1
len_check = 1
tmp = []
bf_list = []
while len_check:
#    print m
    rot = z_subr.mk_combination(m)
    if m!=1:
        rot_small = z_subr.mk_combination(m-1)
        rot_new = [x for x in rot if x not in rot_small] 
        bf_list = new_list
        rot = rot_new
#    print len(rot)
    
    # Remove duplicate or all element < 0
#    rot_sort = z_subr.rm_neg_rot(rot, lat_prim)
    rot_sort = rot
#    print "Remove duplicate or all element < 0"
#    print len(rot), len(rot_sort)
    
    # Check length criterion
    new_list = z_subr.check_length(rot_sort, lat_prim, len_min, len_crit)
#    print "Check length criterion"
#    print len(rot_sort), len(new_list)
#    print len(bf_list), len(new_list)
#    print "-------"
    if (len(bf_list) !=0) and (len(new_list) == 0):
        len_check = 0
    else:
        m = m+1

    tmp = tmp+new_list


#print "-------"
new_list = tmp
#print len(new_list)
#print "-------"
#for i in range(len(new_list)):
#    print new_list[i]

# make rotation matrix using basis
# format : [[rot_mat],[lattice_mat]]
basis_mat = z_subr.mk_rot_matrix(new_list)
#print "Make rotation matrix using basis"
#print len(basis_mat)

#print basis_mat

# remove det = 0
basis_mat = z_subr.rm_det_0(basis_mat)
#print "Remove det = 0"
#print len(basis_mat)

# remove too sharp shape structures
basis_mat = z_subr.check_shape(basis_mat, min_angle, max_angle)
#print "Remove too sharp shape structures"
#print len(basis_mat)

# Check the number of atoms:
basis_mat = z_subr.check_natom(basis_mat, lat_prim, totnum, min_natom)
#print "Check the number of atoms"
#print len(basis_mat)

# Check defect distance:
basis_mat = z_subr.check_defect_length(basis_mat, len_min)    
#print "Check defect distance"
#print len(basis_mat)

# Check kpt
reci_prim = z_subr.cal_reci_vec(lat_prim)
KPT = [float(x) for x in open(kp_file,'r').readlines()[3].split()]
min_dk =  max([x/y for x,y in zip(reci_prim,KPT)])

if gamma_on == 1:
    tmp = []
    for i in range(len(basis_mat)):
        lat_mat = basis_mat[i][1]
        reci_mat = z_subr.cal_reci_vec(lat_mat)
        kpt_new = [x/min_dk for x in reci_mat]
        if max(kpt_new) < 1.01:
            tmp.append(basis_mat[i])
            lat_mat = basis_mat[i][1]
            natoms = z_subr.cal_volume(lat_mat)/vol_prim*totnum
    basis_mat = tmp

if len(basis_mat) == 0 :
    print 1
    sys.exit()

# Select mininum total atom number
best_natom = 10e5
best = 'Null'
best_set = []
for i in range(len(basis_mat)):
    lat_mat = basis_mat[i][1]
    length = z_subr.cal_lat_len(lat_mat)
    tot_atom = z_subr.cal_volume(lat_mat)/vol_prim*totnum

    if (max(length)-min(length))/max(length) < latt_diff:
        if tot_atom < best_natom:
            best_natom = tot_atom
            best_set = [[basis_mat[i][0],basis_mat[i][1]]]
        elif tot_atom == best_natom:
            best_set.append([basis_mat[i][0],basis_mat[i][1]])

#print best_set
best_latt_sum = -10e5
for i in range(len(best_set)):
    latt_sum = 0
    for j in range(3):
        for k in range(3):
            latt_sum = latt_sum + best_set[i][1][j][k]

#    print latt_sum
    if latt_sum > best_latt_sum:
        best_latt_sum = latt_sum
        best = best_set[i]
det = int(round(best_natom))/int(totnum)
#print "Atomic position"
# Set atomic position
at = z_subr.set_atom(at_prim, best[0])
lat_mat = best[1]

# Create POSCAR files (GGA/HSE)
#print poslines[6].split()
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
#print(t_start)
#print(t_end)
print 0
