####################################
# Modifier : yybbyb@snu.ac.kr      #
# data : 2019-08-08                #
####################################
import os, sys, yaml, shutil, glob, math, subprocess
from operator import itemgetter
from module_log import *
from module_vector import *
from module_vasprun import pygrep,pyhead,pytail

def make_list(inp_file):
	with open(inp_file,'r') as f:
		inp_yaml = yaml.load(f)
	Submit_path = inp_yaml['directory']['submit']
	ERROR_path = inp_yaml['directory']['error']
	# Submit_path is file (cif or POSCAR)
	cifs = []
	POSCARs = []
	dirs = []
	if os.path.isfile(Submit_path):
		if '.cif' in Submit_path[-4:]:
			cifs = [Submit_path]
		elif 'POSCAR_' in Submit_path.split('/')[-1][:7]:
			POSCARs = [Submit_path]
	# Submit_path is directory for a single material.
	elif os.path.isdir(Submit_path+'/INPUT0'):
		dirs = [Submit_path]
	
	elif os.path.isdir(Submit_path):
		cifs = glob.glob(Submit_path+'/*.cif')
		POSCARs = glob.glob(Submit_path+'/POSCAR_*')

		dir_list = glob.glob(Submit_path+'/*')
		dirs = []
		for ll in dir_list:
			title = ll.split('/')[-1]
			if os.path.isdir(ll) and len(title.split('_')) >= 2:
				if os.path.isdir(ll+'/INPUT0'):
					dirs.append(Submit_path+'/'+title)
				else:
					make_amp2_log(ll,'ERROR: There is no INPUT0.')
					shutil.move(ll,ERROR_path+'/'+title)
	calc_list = []
	for i in range(len(dirs)):
		calc_list.append([dirs[i],'0'])
	for i in range(len(cifs)):
		calc_list.append([cifs[i],'1'])
	for i in range(len(POSCARs)):
		calc_list.append([POSCARs[i],'2'])
	return calc_list

def make_sym(sp_group):
	# Symmetry k-points table number
	k_table = {'1' : '19', '2' : '19', '3' : '15', '4' : '15', '5' : '16', '6' : '15', '7' : '15', '8' : '16',
 '9' : '16', '10' : '15', '11' : '15', '12' : '16', '13' : '15', '14' : '15', '15' : '16', '16' : '7', '17' : '7',
 '18' : '7', '19' : '7', '20' : '11', '21' : '11', '22' : '8', '23' : '10', '24' : '10', '25' : '7', '26' : '7',
 '27' : '7', '28' : '7', '29' : '7', '30' : '7', '31' : '7', '32' : '7', '33' : '7', '34' : '7', '35' : '11',
 '36' : '11', '37' : '11', '38' : '11', '39' : '11', '40' : '11', '41' : '11', '42' : '9', '43' : '8', '44' : '10',
 '45' : '10', '46' : '10', '47' : '7', '48' : '7', '49' : '7', '50' : '7', '51' : '7', '52' : '7', '53' : '7',
 '54' : '7', '55' : '7', '56' : '7', '57' : '7', '58' : '7', '59' : '7', '60' : '7', '61' : '7', '62' : '7',
 '63' : '11', '64' : '11', '65' : '11', '66' : '11', '67' : '11', '68' : '11', '69' : '9', '70' : '8',
 '71' : '10', '72' : '10', '73' : '10', '74' : '10', '75' : '4', '76' : '4', '77' : '4', '78' : '4',
 '79' : '6', '80' : '6', '81' : '4', '82' : '5', '83' : '4', '84' : '4', '85' : '4', '86' : '4', '87' : '6',
 '88' : '6', '89' : '4', '90' : '4', '91' : '4', '92' : '4', '93' : '4', '94' : '4', '95' : '4', '96' : '4',
 '97' : '6', '98' : '6', '99' : '4', '100' : '4', '101' : '4', '102' : '4', '103' : '4', '104' : '4',
 '105' : '4', '106' : '4', '107' : '6', '108' : '6', '109' : '6', '110' : '6', '111' : '4', '112' : '4',
 '113' : '4', '114' : '4', '115' : '4', '116' : '4', '117' : '4', '118' : '4', '119' : '5', '120' : '5',
 '121' : '5', '122' : '5', '123' : '4', '124' : '4', '125' : '4', '126' : '4', '127' : '4', '128' : '4',
 '129' : '4', '130' : '4', '131' : '4', '132' : '4', '133' : '4', '134' : '4', '135' : '4', '136' : '4',
 '137' : '4', '138' : '4', '139' : '6', '140' : '6', '141' : '6', '142' : '6', '143' : '12', '144' : '12',
 '145' : '12', '146' : '13', '147' : '12', '148' : '13', '149' : '12', '150' : '12', '151' : '12',
 '152' : '12', '153' : '12', '154' : '12', '155' : '13', '156' : '12', '157' : '12', '158' : '12',
 '159' : '12', '160' : '13', '161' : '13', '162' : '12', '163' : '12', '164' : '12', '165' : '12',
 '166' : '13', '167' : '13', '168' : '12', '169' : '12', '170' : '12', '171' : '12', '172' : '12',
 '173' : '12', '174' : '12', '175' : '12', '176' : '12', '177' : '12', '178' : '12', '179' : '12',
 '180' : '12', '181' : '12', '182' : '12', '183' : '12', '184' : '12', '185' : '12', '186' : '12',
 '187' : '12', '188' : '12', '189' : '12', '190' : '12', '191' : '12', '192' : '12', '193' : '12',
 '194' : '12', '195' : '1', '196' : '2', '197' : '3', '198' : '1', '199' : '3', '200' : '1', '201' : '1',
 '202' : '2', '203' : '2', '204' : '3', '205' : '1', '206' : '3', '207' : '1', '208' : '1', '209' : '2',
 '210' : '2', '211' : '3', '212' : '1', '213' : '1', '214' : '3', '215' : '1', '216' : '2', '217' : '3',
 '218' : '1', '219' : '2', '220' : '3', '221' : '1', '222' : '1', '223' : '1', '224' : '1', '225' : '2',
 '226' : '2', '227' : '2', '228' : '2', '229' : '3', '230' : '3'}
	table = k_table[sp_group]
	return table

def make_poscar_from_cif(cif,target):
	# Read cif file
	f = open(cif,'r')
	line = f.readline()
	check_error = 0
	while line :
		### read lattice parameters
		if "_cell_length_a " in line :
			a = read_cif_lattice(line)
			if a == -1:
				make_amp2_log(target,'ERROR. Invalid cif file')
				sys.exit()
		if "_cell_length_b " in line :
			b = read_cif_lattice(line)
			if b == -1:
				make_amp2_log(target,'ERROR. Invalid cif file')
				sys.exit()
		if "_cell_length_c " in line :
			c = read_cif_lattice(line)
			if c == -1:
				make_amp2_log(target,'ERROR. Invalid cif file')
				sys.exit()
		if "_cell_angle_alpha " in line :
			alpha = read_cif_lattice(line)
			if alpha == -1:
				make_amp2_log(target,'ERROR. Invalid cif file')
				sys.exit()
		if "_cell_angle_beta " in line :
			beta = read_cif_lattice(line)
			if beta == -1:
				make_amp2_log(target,'ERROR. Invalid cif file')
				sys.exit()
		if "_cell_angle_gamma " in line :
			gamma = read_cif_lattice(line)
			if gamma == -1:
				make_amp2_log(target,'ERROR. Invalid cif file')
				sys.exit()

		### read symmetry information
		# from ICSD cif
		if "_symmetry_equiv_pos_as_xyz" in line :
			sym = []
			while True :
				tmp = f.readline()
				tmp = tmp.replace('\n','').replace('\r','')
				if not "'" in tmp:
					break
				tmp = tmp.split("'")[1]
				tmp = tmp.replace(" ","")
				tmp = tmp.replace("/",".0/")
				sym.append(tmp.split(','))

		# from VESTA cif
		if "_space_group_symop_operation_xyz" in line :
			sym = []
			while True :
				tmp = f.readline()
				tmp = tmp.replace('\n','').replace('\r','')
				if not len(tmp.split()) > 0:
					break
				if not tmp.split()[0][0] == "'":
					break
				tmp = ''.join(tmp.split())
				tmp = tmp.replace("'","")
				tmp = tmp.replace("/",".0/")
				sym.append(tmp.split(','))


		### read atom information
		if "loop_" in line:
			tmp = f.readline()
			if "_atom_site_" in tmp:
				atom_info_idx_cif = [-1,-1,-1,-1,-1,-1] # label, frac_pos_x, frac_pos_y, frac_pos_z, occupancy, type_symbol
				atom_info_idx_num = 0
				while len(tmp.split()) == 1:
					val = atom_inform_idx(tmp)
					if not val == -1:
						atom_info_idx_cif[val] = atom_info_idx_num
					atom_info_idx_num = atom_info_idx_num+1

					tmp = f.readline()

				if sum(atom_info_idx_cif) == -6:
					break
				elif -1 in atom_info_idx_cif[0:5]:
					make_amp2_log(target,'ERROR. There is not enough information if cif file.')
					sys.exit()

				atoms = []	# Atom postition lines
				while len(tmp.split()) == atom_info_idx_num:
					# atom_info = [name,index,charge,occupancy,pos_vector]
					atom_info = read_cif_position(tmp,atom_info_idx_cif)

					# Partial occupation case (not used) 
					check = 0
					if atom_info[3] != 1.0 :
						make_amp2_log(target,'Warning. Partially occufied structure. Check '+atom_info[0]+atom_info[1])
						for atom_list in atoms:
							if atom_list[4] == atom_info[4]:
								if not atom_list[0] == atom_info[0]:
									make_amp2_log(target,'ERROR. Two different type atoms cannot be placed at the same position.')
									sys.exit()
								else:
									make_amp2_log(target,'Warning. Same type but different index atoms are placed at the same position.')
									make_amp2_log(target,'One is '+atom_list[0]+atom_list[1]+'. The other is '+atom_info[0]+atom_info[1]+'.')
									check = 1
					if not check == 1:
						atoms.append(atom_info)
					tmp = f.readline()
			else:
				line = tmp
		else:
			line = f.readline()

	### make poscar from atomic information
	# lattice parameter	
        alpha = math.radians(alpha)
        beta = math.radians(beta)
        gamma = math.radians(gamma)
	
	aa = [a,0.0,0.0]
	bb = [b*math.cos(gamma),b*math.sin(gamma),0.0]
	c1 = c*math.cos(beta)
	c2 = (c*math.cos(alpha)-math.cos(gamma)*c*math.cos(beta))/math.sin(gamma)
	c3 = (c**2.0-c1**2.0-c2**2.0)**0.5
	cc = [c1,c2,c3]
	# atomic postion from symmetry information
	with open(target+'/INPUT0/final.dat','w') as out_final:
		for i in range(len(atoms)):
			out_final.write('\t'.join([str(x) for x in atoms[i][:3]])+'\n')

	atom_z = []		# Final atom names
	atom_cnt = []	# Final atom counts
	atom_pos = []
	k = 0
	chem = atoms[0][0]
	atom_z.append(chem)
	atom_cnt.append(0)
	for index in range(len(atoms)) :
		if chem != atoms[index][0] :
			chem = atoms[index][0]
			atom_z.append(chem)
			atom_cnt.append(0)
			k = k+1
		# atom_pos[i] = [x,y,z,atom_type,atom_name]
		for j in range(len(sym)) :
			[x, y, z] = atoms[index][4]
			exec 'xx='+sym[j][0]
			exec 'yy='+sym[j][1]
			exec 'zz='+sym[j][2]
			[xx,yy,zz] = setBasis(xx,yy,zz)
			check = 0
			for l in range(len(atom_pos)) :
				if abs(xx-atom_pos[l][0]) < 0.001 and abs(yy-atom_pos[l][1]) < 0.001 and abs(zz-atom_pos[l][2]) < 0.001 :
					check = 1
			if check == 0 :
				atom_pos.append([xx, yy, zz, atoms[index][0], atoms[index][1]])
				atom_cnt[k] = atom_cnt[k]+1

	# Write POSCAR & POTCAR
	write_poscar([aa,bb,cc],atom_pos,target+'/INPUT0/POSCAR_conv','cif to poscar')
	return 0

# poscar file error check
def poscar_error_check(poscar_full,target):
	poscar_file_name = poscar_full.split('/')[-1]
	poscar = poscar_file_name.split('_')
	if len(poscar) < 2 :
		make_amp2_log(target,'Error for '+poscar_full)
		make_amp2_log(target,'Name of the POSCAR file must be POSCAR_title.')
		shutil.move(poscar_full, target+'/'+poscar_file_name)
		return 1
		sys.exit()
	pos = open(poscar_full,'r').readlines()
	atom_z = pos[6].split()
	atom_cnt = pos[6].split()
	if not len(atom_z) == len(atom_cnt):
		make_amp2_log(target,'Error for '+poscar_full)
		make_amp2_log(target,'Not enough composition information. (Name of atom type should be included.)')
		shutil.move(poscar_full, target+'/'+poscar_file_name)
		return 1
		sys.exit()
	return 0

# make potcar
def make_potcar(poscar,pot_path_gga,pot_path_lda,target,src_path,pot_name):
	POT_LDA = yaml.load(open(src_path+'/pot_table.yaml','r'))['LDA']
	POT_GGA = yaml.load(open(src_path+'/pot_table.yaml','r'))['GGA']
	if not pot_name['GGA'] is None:
		for pot_key in pot_name['GGA'].keys() :
			POT_GGA[pot_key] = pot_name['GGA'][pot_key]
	if not pot_name['LDA'] is None:
		for pot_key in pot_name['LDA'].keys() :
			POT_GGA[pot_key] = pot_name['LDA'][pot_key]

	POT_dir_GGA = pot_path_gga+'/'
	POT_dir_LDA = pot_path_lda+'/'

	pos = open(poscar,'r').readlines()
	atom_z = pos[5].split()
	pot_gga_none = 0
	pot_lda_none = 0
	for atom_name in atom_z :
		if atom_name in POT_LDA and os.path.isfile(POT_dir_LDA+POT_LDA[atom_name]+'/POTCAR'):
			with open(target+'/INPUT0/POTCAR_LDA', 'a') as pot_out_L:
				pot_out_L.write(open(POT_dir_LDA+POT_LDA[atom_name]+'/POTCAR','r').read())
		else:
			pot_lda_none = 1
		if atom_name in POT_GGA and os.path.isfile(POT_dir_GGA+POT_GGA[atom_name]+'/POTCAR'):
			with open(target+'/INPUT0/POTCAR_GGA', 'a') as pot_out_G:
				pot_out_G.write(open(POT_dir_GGA+POT_GGA[atom_name]+'/POTCAR','r').read())
		else:
			pot_gga_none = 1

	if pot_lda_none == 1:
		make_amp2_log(target+'/INPUT0','Cannot make POTCAR for LDA.')

	if pot_gga_none == 1:
		make_amp2_log(target+'/INPUT0','Cannot make POTCAR for GGA.')
		return 1
		sys.exit()

	if pot_gga_none == 1 and pot_lda_none == 1:
		make_amp2_log(target+'/INPUT0','There is no POTCAR file.')
		return 1
		sys.exit()
	return 0


# make incar note file including spin, +u and soc information.
def make_incar_note(poscar,target,soc_target,u_value,magmom_def,src_path):
	# U-J values for LDA+U method
	TMU = yaml.load(open(src_path+'/U_table.yaml','r'))
	if not u_value is None:
		if 'All' in u_value.keys() :
			for u_key in TMU.keys():
				TMU[u_key] = u_value['All']
		elif 'all' in u_value.keys() :
			for u_key in TMU.keys():
				TMU[u_key] = u_value['all']
		for u_key in u_value.keys() :
			TMU[u_key] = u_value[u_key]
		for u_key in TMU.keys():
			if TMU[u_key] == 0:
				TMU.pop(u_key)
	
	if soc_target is None:
		SOC = []
	else:
		SOC = soc_target

	TMd = ['Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg']
	TMf = ['La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu']

	NVAL = {'H':1, 'He':0, 'Li':1, 'Be':2, 'B':3, 'C':4, 'N':5, 'O':6, 'F':7, 'Ne':0, 'Na':1, 'Mg':2, 'Al':3,
        'Si':4, 'P':5, 'S':6, 'Cl':7, 'Ar':0, 'K':1, 'Ca':2, 'Sc':3, 'Ti':4, 'V':5, 'Cr':6, 'Mn':7,
        'Fe':8, 'Co':9, 'Ni':10, 'Cu':11, 'Zn':2, 'Ga':3, 'Ge':4, 'As':5, 'Se':6, 'Br':7, 'Kr':0,
        'Rb':1, 'Sr':2, 'Y':3, 'Zr':4, 'Nb':5, 'Mo':6, 'Tc':7, 'Ru':8, 'Rh':9, 'Pd':10, 'Ag':11,
        'Cd':2, 'In':3, 'Sn':4, 'Sb':5, 'Te':6, 'I':7, 'Xe':0, 'Cs':1, 'Ba':2, 'La':3, 'Ce':4,
        'Pr':5, 'Nd':6, 'Pm':7, 'Sm':8, 'Eu':9, 'Gd':10, 'Tb':11, 'Dy':12, 'Ho':13, 'Er':14, 'Tm':15,
        'Yb':16, 'Lu':17, 'Hf':4, 'Ta':5, 'W':6, 'Re':7, 'Os':8, 'Ir':9, 'Pt':10, 'Au':11, 'Hg':12,
        'Tl':3, 'Pb':4, 'Bi':5, 'Po':6, 'At':7, 'Rn':0, 'Ac':3, 'Th':4, 'Pa':5, 'U':6, 'Pu':8, 'No':16}

	pos = open(poscar,'r').readlines()
	atom_z = pos[5].split()
	atom_cnt = pos[6].split()

	# Check U_note for 3d transition metals	
	u_note = ""
	u_val1 = ""
	u_val2 = ""
	u_on = 0
	mag_on = 0
	soc_on = 0
	nelect = 0

	for j in range(len(atom_z)) :
		pot_line = pygrep('ZVAL',target+'/INPUT0/POTCAR_GGA',0,0).splitlines()[j].split()[5]
#		pot_line = subprocess.check_output(['grep','ZVAL',target+'/INPUT0/POTCAR_GGA']).splitlines()[j].split()[5]
		nelect = nelect + int(atom_cnt[j])*int(float(pot_line))
		if atom_z[j] in TMd :
			if atom_z[j] in TMU :
				if u_on == 0:
					u_on = 1
				u_note = u_note + ' 2'
				u_val1 = u_val1 + ' ' + str(TMU[atom_z[j]]+1.0)
				u_val2 = u_val2 + ' 1.0'
			else:
				u_note = u_note + ' -1'
				u_val1 = u_val1 + ' 0.0'
				u_val2 = u_val2 + ' 0.0'
			mag_on = 1
		elif atom_z[j] in TMf :
			if atom_z[j] in TMU :
				u_on = 2
				u_note = u_note + ' 3'
				u_val1 = u_val1 + ' ' + str(TMU[atom_z[j]]+1.0)
				u_val2 = u_val2 + ' 1.0'
			else :
				u_note = u_note + ' -1'
				u_val1 = u_val1 + ' 0.0'
				u_val2 = u_val2 + ' 0.0'
			mag_on = 1
		else :
			u_note = u_note + ' -1'
			u_val1 = u_val1 + ' 0.0'
			u_val2 = u_val2 + ' 0.0'
		if atom_z[j] in SOC :
			soc_on = 1

	# Checking metallic compounds (If all components are metallic elements, we do not impose +U.)
	non_metal = ['H','He','B','C','N','O','F','Ne','Si','P','S','Cl','Ar','Ge','As','Se','Br','Kr','Sb','Te','I','Xe','At','Rn']
	if not any([x in non_metal for x in atom_z]):
		u_on = 0

	maxmix = str(2+2*u_on)
	with open(target+'/INPUT0/U_note', 'w') as out:
		if not u_on == 0 :
			out.write('\n LDA+U calculation:\n   LDAU      = .TRUE. !Treat L(S)DA+U\n   LDAUTYPE  = 2\n')
			out.write('   LDAUL   ='+u_note+'\n')
			out.write('   LDAUU   ='+u_val1+'\n')
			out.write('   LDAUJ   ='+u_val2+'\n')
			out.write('   LDAUPRINT = 0\n   LMAXMIX = '+maxmix+'\n\n')


	oxi_check = 0	# 0: no spin-polarizaton, 1: use MAGMOM, 2: no MAGMOM
	tot_atom = sum([int(x) for x in atom_cnt])

	if mag_on == 1 and os.path.isfile(target+'/INPUT0/final.dat'):
		make_amp2_log(target+'/INPUT0','MAGMOM is determined by cif file\n')
		with open(target+'/INPUT0/final.dat','r') as fin_file:
			final = []
			for line in fin_file.readlines():
				final.append(line.split())
		if len(final[0]) == 3:
			mag = []	# [atom_oxi, atom, mag, n_mag]
			fin_idx = 0
			for pos_idx in range(9,tot_atom+9):
				if pos_idx != 9 and pos[pos_idx].split()[-1] == mag[-1][0] :
					mag[-1][3] = mag[-1][3] + 1
				else:
#					if len(final[fin_idx][2]) < 2 :
					if len(final[fin_idx][2]) < 2 or float(final[fin_idx][2][:-1]) == 0.0 :
						oxi_check = 2	# No oxi information
						break
					mag.append([final[fin_idx][1],final[fin_idx][0],float(final[fin_idx][2][-1]+final[fin_idx][2][:-1]),1])	# float part: number sign --> sign number
					fin_idx = fin_idx+1


			if oxi_check != 2 :
				for j in range(len(mag)) :
					if mag[j][2] > 0 :
						mag[j][2] = NVAL[mag[j][1]] - mag[j][2]
						if mag[j][1] in TMd and mag[j][2] > 5 and mag[j][2] <= 10 :
							mag[j][2] = 10 - mag[j][2]
						elif mag[j][1] in TMf and mag[j][2] > 7 and mag[j][2] <= 14 :
							mag[j][2] = 14 - mag[j][2]
						elif (not mag[j][1] in TMd) and (not mag[j][1] in TMf) :
							if mag[j][2] == 2 :
								mag[j][2] = 0.0
						if mag[j][2] > 0 :
							mag[j][2] = mag[j][2]+0.5
							oxi_check = 1
					else :
						mag[j][2] = 0.0

	# write spin note for unpaired spin electrons
	if mag_on == 1:
		if not os.path.isfile(target+'/INPUT0/final.dat') or not len(final[0]) == 3 or oxi_check == 2:
			make_amp2_log(target+'/INPUT0','MAGMOM is determined by magmom default value\n')
			oxi_check = 1
			mag = []
			for j in range(len(atom_z)) :
				mag.append([atom_z[j],'',0,atom_cnt[j]])
				if atom_z[j] in TMd :
					mag[-1][2] = magmom_def
				elif atom_z[j] in TMf :
					mag[-1][2] = magmom_def

	with open(target+'/INPUT0/spin_note','w') as out2:
		if oxi_check == 1 :
			out2.write(" Spin polarized caculation w/ initial magnetic moment:\n")
			out2.write("   ISPIN = 2\n   MAGMOM = ")
			for j in range(len(mag)) :
				out2.write(str(mag[j][3])+"*"+str(round(mag[j][2],2))+" ")
			out2.write('\n')
		elif oxi_check == 2 :	# Manual setting part
			out2.write("   ISPIN = 2\n ")
			with open(target+'/valence_notice!', 'w') as out3:
				out3.write("Unpaired spin atoms can not be clearly identified!\n")
		elif nelect%2 != 0 :
			out2.write(" Spin polarized caculation:\n")
			out2.write("   ISPIN = 2\n ")
		if soc_on == 1 :
			out2.write("#   LSORBIT = .True.\n")
			if u_on == 1 :
				out2.write("#   LMAXMIX = "+maxmix+"\n")

# make incar file in INCAR0 directory. (reference incar for target meterial)
def make_incar(poscar,target,src_path,max_nelm):
	from module_vasprun import wincar
	pos = open(poscar,'r').readlines()
	atom_z = pos[5].split()
	atom_cnt = pos[6].split()
	nelect = 0

	for j in range(len(atom_z)) :
		pot_line = pygrep('ZVAL',target+'/INPUT0/POTCAR_GGA',0,0).splitlines()[j].split()[5]
#		pot_line = subprocess.check_output(['grep','ZVAL',target+'/INPUT0/POTCAR_GGA']).splitlines()[j].split()[5]
		nelect = nelect + int(atom_cnt[j])*int(float(pot_line))

	# Set maximum number of electronic steps
	nelm = max_nelm
#	if nelect > int(max_nelm/3) :
#		nelm = max_nelm
#	else :
#		nelm = nelect*3
#	if nelm < 30:
#		nelm = 30
	wincar(src_path+'/INCAR0',target+'/INPUT0/INCAR0',[['SYSTEM',target.split('/')[-1]],['NELM',str(nelm)]],[])
	with open(target+'/INPUT0/INCAR','w') as out_inc:
		with open(target+'/INPUT0/INCAR0','r') as f:
			out_inc.write(f.read())
		with open(target+'/INPUT0/U_note','r') as f:
			out_inc.write(f.read())
		with open(target+'/INPUT0/spin_note','r') as f:
			out_inc.write(f.read())
#		subprocess.call(['cat',target+'/INPUT0/INCAR0',target+'/INPUT0/U_note',target+'/INPUT0/spin_note'], stdout=out_inc)

def read_cif_lattice(line):
	a = line.replace('\n','').replace('\r','')
	a = a.split()[1]
	if '(' in line :
		a = a[:a.index('(')]
	if a.replace('.','').isdigit() :
		a = float(a)
		return a
	else :
		return -1

def atom_inform_idx(line):
	if '_atom_site_label' in line:
		return 0
	elif '_atom_site_fract_x' in line:
		return 1
	elif '_atom_site_fract_y' in line:
		return 2
	elif '_atom_site_fract_z' in line:
		return 3
	elif '_atom_site_occupancy' in line:
		return 4
	elif '_atom_site_type_symbol' in line:
		return 5
	else:
		return -1

def read_cif_position(line,inform):
	tmp = line.split()
	if tmp[inform[0]][1].isdigit():
		atom_name_length = 1
	else:
		atom_name_length = 2
	atom_name = tmp[inform[0]][0:atom_name_length]
	atom_index = tmp[inform[0]][atom_name_length:]

	pos = []
	for i in range(3):
		if '(' in tmp[inform[1+i]] :
			tmp[inform[1+i]] = tmp[inform[1+i]][:tmp[inform[1+i]].index('(')]
		if tmp[inform[1+i]] == '.':
			tmp[inform[1+i]] = '0.0'
		pos.append(float(tmp[inform[1+i]]))

	if '(' in tmp[inform[4]] :
		tmp[inform[4]] = tmp[inform[4]][:tmp[inform[4]].index('(')]
	occupancy = float(tmp[inform[4]])

	if not inform[5] == -1:
		atom_charge = tmp[inform[5]][atom_name_length:]
	else:
		atom_charge = 0

	out = [atom_name,atom_name+atom_index,atom_charge,occupancy,pos]
	return out

def read_poscar(poscar):
	with open(poscar,'r') as pos:
		lines = pos.readlines()
	axis_scale = float(lines[1].split()[0])
	axis = []
	for i in range(3):
		axis.append([float(x)*axis_scale for x in lines[2+i].split()])
	type_name = lines[5].split()
	type_num = [int(x) for x in lines[6].split()]
	if lines[7][0] == 'S' or lines[7][0] == 's':
		cart_type = lines[8][0]
		cart_idx = 8
	else:
		cart_type = lines[7][0]
		cart_idx = 7
	atom_pos = []
	atom_idx = cart_idx
	sym_info = 0
	for i in range(len(type_num)):
		for j in range(type_num[i]):
			atom_idx = atom_idx+1
			if len(lines[atom_idx].split('!')) > 1: 
				typ_idx = lines[atom_idx].split('!')[-1].split()[0]
			else:
				# symmetry check should be included
				typ_idx = type_name[i]+'1'
				sym_info = 1
			atom_pos_line = [float(x) for x in lines[atom_idx].split()[0:3]]
			if cart_type == 'C' or cart_type == 'c':
				atom_pos_line = [float(x)*axis_scale for x in lines[atom_idx].split()[0:3]]
				atom_pos_line = cart_to_dir(atom_pos_line,axis)
			atom_pos.append(atom_pos_line+[type_name[i],typ_idx])
	# atom_pos = [[x,y,z,type_name,detail_type_name]*# of atoms]
	if sym_info == 1:
		atom_pos = impose_atom_type_index(axis,atom_pos)
	return [axis,atom_pos]

def impose_atom_type_index(axis,atom_pos):
	import spglib
	import numpy as np
	L = np.mat(axis)
	pos = []
	atom_type = []
	atom_dic = {}
	type_dic = {}
	index = 0
	for line in atom_pos:
		pos.append(line[0:3])
		if not line[4] in atom_dic.keys():
			atom_dic[line[4]] = index
			type_dic[index] = [line[3],line[4]]
			index = index+1
		atom_type.append(atom_dic[line[4]])
	D = np.mat(pos)
	Cell = (L,D,atom_type)
	equ_atoms = spglib.get_symmetry_dataset(Cell,symprec=2e-3)['equivalent_atoms']
	new_type_dic = {}
	new_index = [1 for x in range(index)]
	new_atom_pos = []
	for i in range(len(atom_pos)):
		if not str(atom_type[i])+'_'+str(equ_atoms[i]) in new_type_dic.keys():
			new_type_dic[str(atom_type[i])+'_'+str(equ_atoms[i])] = [atom_pos[i][3],atom_pos[i][3]+str(new_index[atom_type[i]])]
			new_index[atom_type[i]] = new_index[atom_type[i]]+1
		new_atom_pos.append(atom_pos[i][0:3]+new_type_dic[str(atom_type[i])+'_'+str(equ_atoms[i])])
	return new_atom_pos

def write_poscar(axis,atom_pos,out_pos,title):
	atom_type = []
	atom_num = []
	atom_type.append(atom_pos[0][3])
	atom_num.append(1)
	for i in range(1,len(atom_pos)):
		if atom_type[-1] == atom_pos[i][3]:
			atom_num[-1] = atom_num[-1] + 1
		else:
			atom_type.append(atom_pos[i][3])
			atom_num.append(1)

	with open(out_pos,'w') as pos:
		pos.write(title+'\n')
		pos.write('   1.000000000\n')
		for i in range(3):
			pos.write('      '+'    '.join([str(x) for x in axis[i]])+'\n')
		pos.write('    '+'    '.join(atom_type)+'\n')
		pos.write('    '+'    '.join([str(x) for x in atom_num])+'\n')
		pos.write('Selective dynamics\n')
		pos.write('Direct\n')
		for line in atom_pos:
			pos.write('    '+'    '.join([str(x) for x in line[0:3]])+'  T  T  T ! '+line[4]+'\n')

def get_primitive_cell(axis,atom_pos):
	import spglib
	import numpy as np
	L = np.mat(axis)
	pos = []
	atom_type = []
	atom_dic = {}
	type_dic = {}
	index = 0
	for line in atom_pos:
		pos.append(line[0:3])
		if not line[4] in atom_dic.keys():
			atom_dic[line[4]] = index
			type_dic[index] = [line[3],line[4]]
			index = index+1
		atom_type.append(atom_dic[line[4]])
	D = np.mat(pos)
	Cell = (L,D,atom_type)
	prim_cell = spglib.find_primitive(Cell,symprec=2e-3)
	equ_atoms = spglib.get_symmetry_dataset(prim_cell,symprec=2e-3)['equivalent_atoms']
	prim_axis = prim_cell[0].tolist()
	prim_pos = prim_cell[1].tolist()
	prim_type = prim_cell[2].tolist()

	prim_atom_pos = []
	for i in range(len(prim_pos)):
		prim_atom_pos.append(prim_pos[i]+type_dic[prim_type[i]])

	sym = spglib.get_spacegroup(Cell,symprec=2e-3).split('(')[1].split(')')[0]

	[new_axis,new_pos,sym_typ] = axis_rotation(prim_axis,prim_atom_pos,sym)
	return [new_axis,new_pos,sym_typ]

def axis_rotation(axis,pos,sym):
	brava_sym = int(make_sym(sym))
	sym = int(sym)
	pos_new = []
	axis_new = []
	axis_order = [0,1,2]
	# Determine barava sym of BCT 
	if brava_sym in [5,6]:
		if abs(axis[0][0]) > abs(axis[0][2]):
			brava_sym = 5
		else:
			brava_sym = 6

	# sort axis order a < b < c for orthorhombic
	elif brava_sym == 7:
		dic = [[axis[0][0],0],[axis[1][1],1],[axis[2][2],2]]
		order = sorted(dic)
		axis_order = [order[x][1] for x in range(3)]

	# Determine brava sym of ORCF
	elif brava_sym in [8,9]:
		dic = [[axis[1][0],0],[axis[0][1],1],[axis[0][2],2]]
		order = sorted(dic)
		axis_order = [order[x][1] for x in range(3)]
		if 1./(sum([axis[x][axis_order[0]] for x in range(3)])**2.0) < 1./(sum([axis[x][axis_order[1]] for x in range(3)])**2.0) + 1./(sum([axis[x][axis_order[2]] for x in range(3)])**2.0):
			brava_sym = 9
		else:
			brava_sym = 8

	elif brava_sym == 10:
		dic = [[abs(axis[0][0]),0],[abs(axis[1][1]),1],[abs(axis[2][2]),2]]
		order = sorted(dic)
		axis_order = [order[x][1] for x in range(3)]

	# sort axis order a < b for c-centered orthorhombic
	elif brava_sym == 11:
		# a-centered orthorhombic --> c-centered orthorhombic
		if sym in [38,39,40,41]:
			axis_order = [1,2,0]
		# sort axis order a < b 
		if abs(axis[axis_order[0]][axis_order[0]]) > abs(axis[axis_order[0]][axis_order[1]]):
			tmp = [x for x in axis_order]
			axis_order = [tmp[1],tmp[0],tmp[2]]

	# Determine brava sym of RHL
	elif brava_sym in [13,14]:
		angle = calc_angle(axis[0][:],axis[1][:])
		if math.degrees(angle) < 90.:
			brava_sym = 13
		else:
			brava_sym = 14

	# change a and b axis for simple monoclinic phase
	elif brava_sym in [15,16,17,18]:
		axis_order = [1,0,2]

	# Determine brava sym of TRI
	elif brava_sym in [19,20]:
		reci_lat = reciprocal_lattice(axis)
		k_alpha = calc_angle(reci_lat[2],reci_lat[1])
		if math.degrees(k_alpha) > 90.:
			brava_sym = 19
		else:
			brava_sym = 20


	# redefine axis and position
	for i in range(3):
		axis_new.append([axis[axis_order[i]][x] for x in axis_order])
	for line in pos:
		pos_new.append([line[axis_order[x]] for x in range(3)]+line[3:])

	# rotate axis for alpha < 90
	if brava_sym in [15,16,17,18]:
		axis_new[2][1] = -axis_new[2][1]
		for i in range(len(pos)):
			pos_new[i][2] = 1.-pos_new[i][2]

	# Determine brava sym of MCLC
	if brava_sym in [16,17,18]:
		conv_lat = [[abs(axis_new[0][0])*2.0,0,0],[0,abs(axis_new[1][1])*2.0,0],axis_new[2]]
		len_c = dist_point(axis_new[2],[0,0,0])
		reci_lat = reciprocal_lattice(conv_lat)
		k_gamma = calc_angle(reci_lat[0],reci_lat[1])
		if math.degrees(k_gamma) >= 90.:
			brava_sym = 16
		elif conv_lat[1][1]/len_c*conv_lat[2][1]/len_c+(conv_lat[1][1]/conv_lat[0][0]*conv_lat[2][2]/len_c)**2.0 <=1.:
			brava_sym = 17
		else:
			brava_sym = 18

	return [axis_new,pos_new,brava_sym]

def setBasis(x,y,z) :
	x = x - int(x)
	y = y - int(y)
	z = z - int(z)
	if x < 0 :
		x = x+1
	if y < 0 :
		y = y+1
	if z < 0 :
		z = z+1
	if abs(x-1) < 0.0001 :
		x = 0
	if abs(y-1) < 0.0001 :
		y = 0
	if abs(z-1) < 0.0001 :
		z = 0
	return [x,y,z]

def write_file(line,filename):
	with open(filename,'w') as out:
		out.write(line)
