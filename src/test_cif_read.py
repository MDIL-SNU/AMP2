import sys
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

def make_poscar_from_cif(cif):
	# Read cif file
	f = open(cif,'r')
	line = f.readline()
	check_error = 0
	while line :
		line = f.readline()
		### read lattice parameters
		if "_cell_length_a " in line :
			a = read_cif_lattice(line)
		if "_cell_length_b " in line :
			b = read_cif_lattice(line)
		if "_cell_length_c " in line :
			c = read_cif_lattice(line)
		if "_cell_angle_alpha " in line :
			alpha = read_cif_lattice(line)
		if "_cell_angle_beta " in line :
			beta = read_cif_lattice(line)
		if "_cell_angle_gamma " in line :
			gamma = read_cif_lattice(line)

		### read symmetry information
		# from ICSD cif
		if "_symmetry_equiv_pos_as_xyz" in line :
			sym = []
			while True :
				tmp = f.readline()
				tmp = tmp.replace('\n','').replace('\r','')
				if not len(tmp.split()) > 1:
					break
				if not tmp.split()[1][0] == "'":
					break
				tmp = ''.join(tmp.split()[1:])
				tmp = tmp.replace("'","")
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

				print atom_info_idx_cif
				if sum(atom_info_idx_cif) == -6:
					break
				elif -1 in atom_info_idx_cif[0:5]:
					print 'not enough information'
					sys.exit()

				atoms = []	# Atom postition lines
				while len(tmp.split()) == atom_info_idx_num:
					# atom_info = [name,index,charge,occupancy,pos_vector]
					atom_info = read_cif_position(tmp,atom_info_idx_cif)

					# Partial occupation case (not used) 
					check = 0
					if atom_info[3] != 1.0 :
						for atom_list in atoms:
							if atom_list[4] == atom_info[4]:
								if not atom_list[0] == atom_info[0]:
									sys.exit()
								else:
									check = 1
					if not check == 1:
						atoms.append(atom_info)
					tmp = f.readline()
	print atoms

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

make_poscar_from_cif(sys.argv[1])
