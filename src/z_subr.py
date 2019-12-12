import numpy as np
import math
from itertools import product
from itertools import combinations

def dotproduct(v1, v2):
    return sum((a*b) for a, b in zip(v1, v2))

def length(v):
    return math.sqrt(dotproduct(v, v))

def angle(v1, v2):
    value =  dotproduct(v1, v2)/(length(v1)*length(v2))
    if value > 1.0: value = 1.0 
    if value < -1.0: value = -1.0 
    return (math.acos(value))*180/np.pi

def read_file(target):
    with open(target, 'r') as f1:
        lines=f1.readlines()
    return lines

def read_lat(poscar):
    lat_prim=[]
    for i in range(3):
        lat_prim.append([float(x) for x in poscar[2+i].split()])
    return lat_prim

def read_at_pos(poscar,totnum):
    at_prim=[]
    for i in range(totnum):
        tmp = poscar[9+i].split()
        at_prim.append([float(x) for x in tmp[0:3]]+[tmp[-1]])
    return at_prim

def read_tot(poscar,totnum):
    lat_prim=[]
    for i in range(3):
        lat_prim.append([float(x) for x in poscar[2+i].split()])
    at_prim=[]
    for i in range(totnum):
        tmp = poscar[9+i].split()
        at_prim.append([float(x) for x in tmp[0:3]]+[tmp[-1]])
    return lat_prim, at_prim

def cal_lat_len(lattice):
    lat_length = []
    for i in range(3):
        lat_length.append(np.sqrt(sum([x**2 for x in lattice[i]])))
    return lat_length

def cal_angle(lattice):
    alpha = (angle(lattice[1],lattice[2]))
    beta = (angle(lattice[0],lattice[2]))
    gamma = (angle(lattice[0],lattice[1]))
    return [alpha, beta, gamma]

def rotate_mat(lat_prim, rot_mat):
    lat_rot=np.array(np.matmul(np.array(lat_prim).transpose(), rot_mat)).transpose()
    return lat_rot 

def mk_combination(m):
    items= [range(-m,m+1)]*3
    rot = list(product(*items))
    return rot

def rm_neg_rot(rot,lat_prim):
    rot_sort = []
    for i in range(len(rot)):
        count_neg = sum(1 for x in rot[i] if x <= 0)
        if not count_neg == len(rot[i]):
            if [-x for x in rot[i]] not in rot_sort:
                rot_sort.append([x for x in rot[i]])
    return rot_sort

def mk_new_basis(lat_prim, mat):
    new = []
    for i in range(3):
        tmp = 0
        for j in range(3):
            tmp = tmp+lat_prim[j][i]*mat[j]
        new.append(tmp)
    return new

def check_length(rot_sort, lat_prim, len_min, len_crit):
    new_list = []
    min_sw = 0
    for i in range(len(rot_sort)):
        # produce new basis using linear combination
        new_basis = mk_new_basis(lat_prim,rot_sort[i])
        # Check length criteria
        lat_length=np.sqrt(sum([x**2 for x in new_basis]))
        if lat_length > len_min and lat_length < len_crit:
            new_list.append([lat_length, rot_sort[i], new_basis])
        if min_sw == 0 and lat_length < len_crit:
            min_sw = 1
    return [new_list,min_sw]

def mk_rot_matrix(new_list):
    # Sort new_list
    new_list = sorted(new_list, reverse=True)
    basis_list = [x[1:] for x in new_list]
    # basis_mat format: [[rot1, lat1],[rot2, lat2],[rot3,lat3]]
    basis_mat = (list(combinations(basis_list,3)))
    # Change format of basis_mat: [[rot_mat],[lattice_mat]]
    tmp_mat = []
    for i in range(len(basis_mat)):
        lattice = []
        rot_vec = []
        for j in range(3):
            rot_vec.append(basis_mat[i][j][0])
            lattice.append(basis_mat[i][j][1])
        lat_sum = 0
        for j in range(3):
            lat_sum = lat_sum+sum(lattice[j])
	if lat_sum > 0:
            tmp_mat.append([rot_vec, lattice])
    return tmp_mat

def rm_det_0(basis_mat):
    tmp = []
    for i in range(len(basis_mat)):
        rot_mat = basis_mat[i][0]
        lat_mat = basis_mat[i][1]
        det = np.linalg.det(rot_mat)
        if abs(det) > 1e-10:
            tmp.append(basis_mat[i])
    basis_mat = tmp
    return basis_mat

def check_shape(basis_mat,min_ang,max_ang):
    tmp = []
    for i in range(len(basis_mat)):
        rot_mat = basis_mat[i][0]
        lat_mat = basis_mat[i][1]
        angle = cal_angle(lat_mat)
        lat_len = cal_lat_len(lat_mat)
        if (min(angle) > min_ang) and (max(angle) < max_ang):
            vol_mat = np.inner(lat_mat[0],np.cross(lat_mat[1],lat_mat[2]))
            # from left-hand to right-hand rule
            if vol_mat < 0:
                basis_mat[i] = [[rot_mat[0],rot_mat[2],rot_mat[1]],[lat_mat[0],lat_mat[2],lat_mat[1]]]
            tmp.append(basis_mat[i])
    basis_mat = tmp
    return basis_mat

def sum_vector(list_of_lists):
#    tmp=[]
#    for i in range(len(a)):
#        tmp.append(a[i]+b[i])


    tmp = [sum(x) for x in zip(*list_of_lists)]
    return tmp

def check_defect_length(basis_mat, len_min):
    tmp=[]
    for i in range(len(basis_mat)):
        dis = []
        [a,b,c]=basis_mat[i][1]
        edge = [a, b, c, sum_vector([a,b]), sum_vector([b,c]), sum_vector([a,c]), sum_vector([a,b,c])]
        for j in range(len(edge)):
            dis.append(np.sqrt(sum([x**2 for x in edge[j]])))
        if min(dis) > len_min:
            tmp.append(basis_mat[i])
#        else:
#            test =  basis_mat[i]
#            print "-----"
#            print [round(x,5) for x in cal_lat_len(test[1])]
#            print [round(x*180/np.pi,4) for x in cal_angle(test[1])]
#            print dis

    basis_mat = tmp
    return basis_mat

def cal_volume(lat):
    return np.inner(lat[0],np.cross(lat[1],lat[2]))

def check_natom(basis_mat, lat_prim, natom, min_atom):
    tmp = []
    vol_prim = cal_volume(lat_prim)
    for i in range(len(basis_mat)):
        vol_mat = cal_volume(basis_mat[i][1])
        natom_mat = vol_mat/vol_prim*natom
        if natom_mat >= min_atom:
            tmp.append(basis_mat[i])

    return tmp

def cal_reci_vec(lat):
    [aa, bb, cc] = lat
    reci = []
    reci.append(np.cross(bb,cc)/(np.inner(aa,np.cross(bb,cc)))) 
    reci.append(np.cross(cc,aa)/(np.inner(bb,np.cross(cc,aa))))
    reci.append(np.cross(aa,bb)/(np.inner(cc,np.cross(aa,bb))))
#    print aa, bb, cc

    len_reci = []
    for i in range(3):
        vector = math.sqrt(sum(map(lambda x: x**2, reci[i][0:3])))
#        print vector
        len_reci.append(vector)
    return len_reci

def set_atom(prim_at, M):
    inv_M = np.linalg.inv(M)
    det_M = np.linalg.det(M)
#    print inv_M
#    print det_M
    MM = M[0]+M[1]+M[2]
    
    rr =  max([abs(x) for x in MM])*4
#    rr = 20
#    print 'rr = '+str(rr)
    sup_at_lines=[]
    sup_at_fin=[]
    for i in range(len(prim_at)):
        for j in range(-rr,rr):
            for k in range(-rr,rr):
                for l in range(-rr,rr): 
                    sup_at_lines.append([prim_at[i][0]+j,prim_at[i][1]+k, prim_at[i][2]+l,prim_at[i][3]])
    count = 0
    for i in range(len(sup_at_lines)):
        [xx,yy,zz] = [round(x,6) for x in np.matmul(sup_at_lines[i][:3],inv_M)]
        if xx >= 0. and xx < 1. and yy >= 0. and yy < 1. and zz >= 0. and zz < 1. :
            sup_at_fin.append([xx,yy,zz,sup_at_lines[i][3]])
            count = count+1
#            print count, sup_at_lines[i], xx, yy, zz
#    print count,len(sup_at_fin)
#    print np.linalg.det(M)
    if int(round(np.linalg.det(M))) * len(prim_at) != count:
        sup_at_fin = []
#        print "rr values is not enough"
#    print sup_at_fin
    return sup_at_fin



