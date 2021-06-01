####################################
# date : 2019-04-08                #
# Author : yybbyb@snu.ac.kr        #
####################################
# This is a package of modules for identifying the most stable magnetic spin ordering.
import os,sys,yaml,subprocess, re
import numpy as np
from module_vector import *
from module_vasprun import *
from module_amp2_input import *
from module_GA import *

def check_spin(ref_dir,inp_pos,min_mom):
	[axis,atom_pos] = read_poscar(inp_pos)
	magnet_out = pygrep('magnetization (x)',ref_dir+'/OUTCAR',0,len(atom_pos)+3).splitlines()[-len(atom_pos):]
	magnet_list = [float(x.split()[-1]) for x in magnet_out]
	type_name = []
	mag_val = {}
	for i in range(len(atom_pos)):
		over = 0
		for j in range(len(type_name)):
			if type_name[j] == atom_pos[i][4]:
				over = 1
				break
		if over == 0:
			type_name.append(atom_pos[i][4])
			mag_val[atom_pos[i][4]] = magnet_list[i]
	new_type_name = []
	for typ in type_name:
		if abs(mag_val[typ]) > min_mom:
			new_type_name.append(typ)
		else:
			mag_val[typ] = 0
	return new_type_name,mag_val

# This function is for finding the possible magnetic pairs
def find_pair(poscar,cutoff,mag_atom_list,mag_val,tol):
	[axis,atom_pos] = read_poscar(poscar)
	sole_list = []
	pair_list = []
	mag_list = []
	for i in range(len(atom_pos)):
		if not atom_pos[i][4] in mag_atom_list:
			mag_list.append(0)
		else:
			mag_list.append(mag_val[atom_pos[i][4]])
			if len(sole_list) == 0:
				sole_list.append([atom_pos[i][4],i])
			else:
				if not atom_pos[i][4] in [x[0] for x in sole_list]:
					sole_list.append([atom_pos[i][4],i])
			for j in range(i,len(atom_pos)):
				if atom_pos[j][4] in mag_atom_list:
					length = short_dist(atom_pos[i][0:3],atom_pos[j][0:3],axis)
					if length <= cutoff and length > 0.1:
						if len(pair_list) == 0:
							pair_list.append([length,atom_pos[i][4],atom_pos[j][4],i,j])
						else:
							dup = 0
							for k in range(len(pair_list)):
								if abs(length-pair_list[k][0]) < tol:
									if atom_pos[i][4] == pair_list[k][1] and atom_pos[j][4] == pair_list[k][2] :
										dup = 1
										break
									elif atom_pos[i][4] == pair_list[k][2] and atom_pos[j][4] == pair_list[k][1] :
										dup = 1
										break
							if dup == 0:
								pair_list.append([length,atom_pos[i][4],atom_pos[j][4],i,j])
	return [pair_list,sole_list,mag_list]

# This function is for writing given magnetic pair list
def write_pair_list(pair_list,sole_list,mag_list,target):
	with open(target+'/pair_list_table.dat','w') as plw:
		for line in pair_list:
			plw.write('\t'.join([str(x) for x in line])+'\n')
	tot_mag_list = []
	[spin,num_spin] = reduce_mag_list(mag_list)
	tot_mag_list.append(' '.join([str(num_spin[x])+'*'+str(spin[x]) for x in range(len(spin))]))
	for line in sole_list:
		new_mag_list = mag_list[0:line[1]]+[-mag_list[line[1]]]+mag_list[line[1]+1:]
		[spin,num_spin] = reduce_mag_list(new_mag_list)
		tot_mag_list.append(' '.join([str(num_spin[x])+'*'+str(spin[x]) for x in range(len(spin))]))
	for line in pair_list:
		new_mag_list = mag_list[0:line[3]]+[-mag_list[line[3]]]+mag_list[line[3]+1:line[4]]+[-mag_list[line[4]]]+mag_list[line[4]+1:]
		[spin,num_spin] = reduce_mag_list(new_mag_list)
		tot_mag_list.append(' '.join([str(num_spin[x])+'*'+str(spin[x]) for x in range(len(spin))]))
	with open(target+'/mag_list.dat','w') as mlw:
		for line in tot_mag_list:
			mlw.write(line+'\n')
	return tot_mag_list

def reduce_mag_list(mag_list):
	spin = []
	num_spin = []
	for i in mag_list:
		if len(spin) == 0:
			spin.append(i)
			num_spin.append(1)
		else:
			if i == spin[-1]:
				num_spin[-1] = num_spin[-1]+1
			else:
				spin.append(i)
				num_spin.append(1)
	return [spin,num_spin]
## for modified Ising model ##
def mk_p_list(target):
    with open(target+'/pair_list_table.dat','r') as f1:
        plt_lines = f1.readlines()
    p_list = []
    for i in range(len(plt_lines)):
        tmp_list = [0 for j in range(3)]
        tmp = plt_lines[i].split()
        for j in range(1,3):
            tmp_list[j-1] = float(re.findall('\d+',tmp[j])[0])
        tmp_list[2] = round(float(tmp[0]),2)
        p_list.append(tmp_list[0:3])
    return p_list
 
def mk_modi_ising_list(mag_lines,p_list,target,mag_atom_list,axis,atom_pos):
    pattern = re.compile(r'((-)?\d{1,3}(,\d{3})*(\.\d+)?)')
    m_atoms = mag_lines.split()
    m_atoms = list(map(float,[x.split('*')[0] for x in mag_lines.split()]))
    s_atoms = list(map(float,[x.split('*')[1] for x in mag_lines.split()]))
    atom_pos_cart = []
    for j in range(len(m_atoms)):
        for k in range(int(sum(m_atoms[0:j])),int(sum(m_atoms[0:j+1]))):
            if s_atoms[j] < 0:
                atom_pos_cart.append(dir_to_cart(atom_pos[k][0:3],axis)+atom_pos[k][3:])
                atom_pos_cart[k][-1] = atom_pos_cart[k][-1]+'_down'
            else:
                atom_pos_cart.append(dir_to_cart(atom_pos[k][0:3],axis)+atom_pos[k][3:])
                atom_pos_cart[k][-1] = atom_pos_cart[k][-1]+'_up'
    mag_true = []
    mag_sum = 0    
    for j in range(len(atom_pos_cart)):
        if atom_pos_cart[j][4] in mag_atom_list:
            mag_true.append(1)
            mag_sum = mag_sum +1
        else:
            mag_true.append(0)
    info_pair = pairCategorize(atom_pos_cart,axis,mag_true,5.0)
    count = {}
    for line in info_pair:
        if line[5] in list(count.keys()):
            count[line[5]] = count[line[5]] +1
        else:
            count[line[5]] = 1
    keys = list(count.keys())
    tot_list=[]
    for j in range(len(keys)): #info_pair to list, format ex) [[1,2,3.66(pair length),4(count),Fe1_upFe2_down3.66(info_pair key)]] 
        tmp_list = []
        idx_check = pattern.findall(keys[j])
        if idx_check[0][0] >= idx_check[1][0]:
            tmp_list.append(float(idx_check[1][0]))
            tmp_list.append(float(idx_check[0][0]))
            tmp_list.append(float(idx_check[2][0]))
        else:
            tmp_list.append(float(idx_check[0][0]))
            tmp_list.append(float(idx_check[1][0]))
            tmp_list.append(float(idx_check[2][0]))
        tmp_list.append(count[keys[j]])
        tmp_list.append(keys[j])
        tot_list.append(tmp_list)
    
    m_list = [[] for j in range(len(p_list))]
    for j in range(len(tot_list)):
        for k in range(len(p_list)):
            if tot_list[j][0:3] == p_list[k]:
                m_list[k].append(tot_list[j])
    return m_list

def cal_pinv(e_list,N_list,target):
    e_array = array(e_list)
    N_array = array(N_list)
    JJ = dot(linalg.pinv(N_array),e_array)
    return JJ

def mk_e_list(target):
    with open(target+'/energy','r') as f1:
        e_lines = f1.readlines()
    e_list = []
    for i in range(len(e_lines)):
        e_list.append(float(e_lines[i].split()[-1]))
    return e_list 

# This function is for calculating ising parameter from magnetic pairs
def calc_ising_param(sole_list,pair_list,ene_file):
	with open(ene_file,'r') as enef:
		lines = enef.readlines()
	ENE = {}
	ENE['Ferro'] = float(lines[0].split()[1])
	for i in range(len(sole_list)):
		ENE[sole_list[i][0]] = float(lines[1+i].split()[1])
	ene_pair = []
	JJ = []
	for i in range(len(pair_list)):
		pair_ene = float(lines[len(sole_list)+i+1].split()[1])
		JJ.append(pair_list[i][0:3]+[(ENE['Ferro']+pair_ene-ENE[pair_list[i][1]]-ENE[pair_list[i][2]])/4])
	return JJ

def write_ising_param(JJ,target):
	with open(target+'/Pair_coeff.dat','w') as pcw:
		for line in JJ:
			pcw.write(str(line[0])+'\t'+str(line[3])+'\t'+line[1]+'\t'+line[2]+'\n')

# This function is for calculating modified ising parameter from magnetic pairs
def calc_ising_param_modi(mag_atom_list,target):
    up_p = re.compile(r'up')
    down_p = re.compile(r'down')
    spin_atom_list = []
    for i in range(len(mag_atom_list)):
        spin_atom_list.append(mag_atom_list[i]+'_up')
        spin_atom_list.append(mag_atom_list[i]+'_down')
    p_list = mk_p_list(target)
    [axis, atom_pos] = read_poscar(target+'/POSCAR_param')
    with open(target+'/mag','r') as f1:
        mag_lines = f1.readlines()
    N_list = []
    for i in range(len(mag_lines)):
        # m_list format = [[1,1,3.82,4,Fe1_upFe1_up3.82][1,1,3.02,6,Fe1_upFe1_down3.02]]..[mag_atom_a,mag_atom_b,pair_length,count,spin_info]
        m_list = mk_modi_ising_list(mag_lines[i],p_list,target,spin_atom_list,axis,atom_pos)
        N_tmp = [1]
        for j in range(len(p_list)):
            ferro = 0
            af = 0
            for k in range(len(m_list[j])):
                a = up_p.findall(m_list[j][k][-1])
                b = down_p.findall(m_list[j][k][-1])
                if len(a) == 2 or len(b) == 2:
                    ferro += m_list[j][k][-2]
                if len(a) == 1 and len(b) == 1:
                    af += m_list[j][k][-2]
            tmp = (ferro-af)/2
            N_tmp.append(tmp)
        N_list.append(N_tmp)
    e_list = mk_e_list(target)
    JJ = cal_pinv(e_list,N_list,target)
    return JJ

def write_ising_param_modi(JJ,target):
    with open(target+'/pair_list_table.dat','r') as f1:
        pt_lines = f1.readlines()
    with open(target+'/Pair_coeff.dat','w') as f1:
        for i in range(len(pt_lines)):
            p_lines = pt_lines[i].split()
            f1.write('{}\t {} \t {} \t {} \n'.format(p_lines[0],JJ[i+1],p_lines[1],p_lines[2]))

# This function is for writing inputs for genetic algorithm
def write_inp_for_GA(mag_atom_list,inp_af,target):
	inp_yaml = {}
	for key in list(inp_af['genetic_algorithm'].keys()):
		inp_yaml[key] = inp_af['genetic_algorithm'][key]
	inp_yaml['mag_atom_list'] = mag_atom_list
	inp_yaml['cutoff'] = inp_af['cutoff_for_parameter']
	inp_yaml['pair_coeff'] = target+'/Pair_coeff.dat'
	with open(target+'/input_GA.yaml','w') as inp_GA:
		yaml.dump(inp_yaml,inp_GA,default_flow_style = False)

def count_mag_atom_num(poscar,mag_atom_list):
	[axis,atom_pos] = read_poscar(poscar)
	mag_atom_num = 0
	for line in atom_pos:
		if line[4] in mag_atom_list:
			mag_atom_num = mag_atom_num + 1
	return mag_atom_num

def set_supercell_list(poscar,mag_atom_list):
	sup_size_list = [1,2]
	supercell_list = []
	mag_atom_num = count_mag_atom_num(poscar,mag_atom_list)
	for i in sup_size_list:
		for j in sup_size_list:
			for k in sup_size_list:
				if i*j*k*mag_atom_num % 2 == 0:
					supercell_list.append([i,j,k,i*j*k])
	sort_supercell_list = sorted(supercell_list,key=lambda x : x[3])
	supercell_list = []
	for line in sort_supercell_list:
		supercell_list.append(line[0:3])
	return supercell_list

# This function is for making distance matrix from poscar
def get_bond_matrix(poscar):
	[axis,atom_pos] = read_poscar(poscar)
	bond_matrix = []
	for i in range(len(atom_pos)):
		tmp = []
		for j in range(len(atom_pos)):
			tmp.append([round(short_dist(atom_pos[i][0:3],atom_pos[j][0:3],axis),6),atom_pos[i][4],atom_pos[j][4]])
		bond_matrix.append(sorted(sorted(tmp,key = lambda x:x[0]),key = lambda x:x[2]))

	sorted_bond_matrix = sorted(sorted(bond_matrix,key = lambda x:x[0][0]),key = lambda x:x[0][1])
	
	return sorted_bond_matrix

# This function is for modifying kpoints corresponding target path
def resized_kpoints(ref_path,targ_path):
	sym = int(open(targ_path+'/sym','r').readline())
	KPT_ref = [float(x) for x in open(ref_path+'/KPOINTS','r').readlines()[3].split()]
	axis = poscar_to_axis(ref_path+'/POSCAR')
	recipro_latt = reciprocal_lattice(axis)
	l = []
	for i in range(3):
		l.append((recipro_latt[i][0]**2.+recipro_latt[i][1]**2.+recipro_latt[i][2]**2.)**0.5)
	min_dk = max([x/y for x,y in zip(l,KPT_ref)])
	KP = []

	axis = poscar_to_axis(targ_path+'/POSCAR')
	recipro_latt = reciprocal_lattice(axis)
	l = []
	for i in range(3):
		l.append((recipro_latt[i][0]**2.+recipro_latt[i][1]**2.+recipro_latt[i][2]**2.)**0.5)
	
	kpoint = open(targ_path+'/KPOINTS','w')
	# Gamma-centred mesh for hexagoanl symmetry
	if sym==12 or sym==13 or sym==14:
		KPset = 'Gamma-centered'
	else :
		KPset = 'Monk-horst'
	for i in range(3) :
		if sym in [5,6,10]: # symmetry BCT
			KP.append(str(max([int(round(x/min_dk)) for x in l])))
		else:
			KP.append(str(int(round(l[i]/min_dk))))
		if KP[-1] == '0' :
			KP[-1] = '1'
	# Write KPOINTS
	kpoint.write("Auto k-point\n 0\n"+KPset+"\n  "+KP[0]+"  "+KP[1]+"  "+KP[2]+"\n  0  0  0\n")
	kpoint.close()

def set_mag_info(atom_pos,target,mag_val):
	magmom = []
	for line in atom_pos:
		if len(line[4].split('_')) > 1:
			sline = line[4].split('_')
			mag_type = '_'.join(sline[:-1])
			if sline[-1] == 'up':
				magmom.append(mag_val[mag_type])
			elif sline[-1] == 'down': 
				magmom.append(-1.0*mag_val[mag_type])
			else:
				magmom.append(mag_val[line[4]])
		else:
			magmom.append(mag_val[line[4]])
	[spin,num_spin] = reduce_mag_list(magmom)
	with open(target+'/mag_info','w') as inp:
		inp.write(' '.join([str(num_spin[x])+'*'+str(spin[x]) for x in range(len(spin))]))

def mk_kpoints_modi(c_dir):
    axis_conv = poscar_to_axis(c_dir+'/relax_GGA/POSCAR')
    recipro_latt_conv = reciprocal_lattice(axis_conv)
    l_conv = []
    for i in range(3):
        l_conv.append((recipro_latt_conv[i][0]**2.+recipro_latt_conv[i][1]**2.+recipro_latt_conv[i][2]**2.)**0.5)
    axis_modi = poscar_to_axis(c_dir+'/magnetic_ordering/POSCAR_param')
    recipro_latt_modi = reciprocal_lattice(axis_modi)
    l_modi = []
    for i in range(3):
        l_modi.append((recipro_latt_modi[i][0]**2.+recipro_latt_modi[i][1]**2.+recipro_latt_modi[i][2]**2.)**0.5)
    l_ratio = [l_modi[x]/l_conv[x] for x in range(3) ]
    with open(c_dir+'/INPUT0/KPOINTS','r') as f1:
        kp_lines = f1.readlines()
    kp_conv = kp_lines[3].split()
    kp_conv = list(map(int,kp_conv))
    with open(c_dir+'/INPUT0/KPOINTS_modi','w') as kp_modi:
        for j in range(3):
            kp_modi.write(kp_lines[j])
        KP = []
        for j in range(3):
            KP.append(str(int(round(kp_conv[j]*l_ratio[j]))))
            if KP[-1] == '0':
                KP[-1] = '1'
        kp_modi.write("  "+KP[0]+"  "+ KP[1]+"  "+KP[2]+"\n  0  0  0\n")

def mk_mag_list_modi(targ_dir):
    datfile = open(targ_dir+'/OUTCAR','r') 
    tmp_mag = {}
    it_mag = 0
    for line in datfile:
        if 'NIONS = ' in line:
            NIONS=int(line.split()[-1])
        if 'magnetization (x)' in line:
            Null = datfile.readline()
            Null = datfile.readline()
            Null = datfile.readline()
            tmp_list = []
            for i in range(int(NIONS)):
                tmp_list.append(datfile.readline().split()[-1])
            tmp_mag[it_mag] = tmp_list
            it_mag +=1
    datfile.close()
    tot_mag = tmp_mag[it_mag-1]
    mag_lines = ''
    for i in range(len(tot_mag)):
        mag_lines += str(' 1*'+tot_mag[i])
    return mag_lines 
