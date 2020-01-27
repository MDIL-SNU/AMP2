import os,sys,json,shutil
from module_vector import *
from module_vasprun import poscar_to_axis
from _version import __version__

def write_formatted_band_log(band_log,band_log_form):
	with open(band_log,'r') as f:
		lines = f.readlines()
	with open(band_log_form,'w') as g:
		if 'etal' in lines[0]:
			for line in lines:
				g.write(line)
		else:
			gap = float(lines[0].split()[2])
			E_vb = float(lines[2].split()[5])
			E_cb = float(lines[3].split()[5])
			g.write('Band gap: '+'{:10.3f}'.format(gap)+' eV '+lines[0].split()[4]+'\n\n')
			g.write('VBM:'+lines[2].split(':')[1]+': '+'{:10.3f}'.format(E_vb)+' eV\n')
			g.write('CBM:'+lines[3].split(':')[1]+': '+'{:10.3f}'.format(E_cb)+' eV\n\n')
			g.write(lines[5]+lines[6])

target = os.path.abspath(sys.argv[1])

if not os.path.isdir(target+'/Results'):
	os.mkdir(target+'/Results')

res_path = target+'/Results'
DB = {}
name =  target.split('/')[-1]
DB['Version'] = __version__
DB['Material_name'] = name

pot_list = ['GGA','LDA','HSE']
for POT in pot_list:
	# Crystal information
	if os.path.isdir(target+'/relax_'+POT):
		axis = poscar_to_axis(target+'/relax_'+POT+'/CONTCAR')
		lat_pos = []
		for i in range(3):
			lat_pos.append(dist_point(axis[i],[0,0,0]))
		ang_pos = []
		ang_pos.append(calc_angle(axis[1],axis[2]))
		ang_pos.append(calc_angle(axis[0],axis[2]))
		ang_pos.append(calc_angle(axis[1],axis[0]))
		DB['Lattice_parameter_'+POT] = lat_pos
		DB['Angle_parameter_'+POT] = ang_pos
		if os.path.isfile(target+'/INPUT0/POSCAR_rlx_'+POT):
			shutil.copy(target+'/INPUT0/POSCAR_rlx_'+POT,res_path+'/POSCAR_'+POT)
		else:
			shutil.copy(target+'/relax_'+POT+'/CONTCAR',res_path+'/POSCAR_'+POT)

	# Band_gap
	if os.path.isfile(target+'/band_'+POT+'/Band_gap.log'):
		with open(target+'/band_'+POT+'/Band_gap.log','r') as inp:
			gap_log = inp.readline()
		if 'etal' in gap_log:
			gga_gap = 0.0
			gga_eg_type = 'Null'
		else:
			gga_gap = float(gap_log.split()[2])
			gga_eg_type = str(gap_log.split()[4])
		DB['Band_gap_'+POT] = gga_gap
		DB['Band_gap_D/I_'+POT] = gga_eg_type
		write_formatted_band_log(target+'/band_'+POT+'/Band_gap.log',res_path+'/Band_gap_'+POT+'.log')
#		shutil.copy(target+'/band_'+POT+'/Band_gap.log',res_path+'/Band_gap_'+POT+'.log')

	# Band_structure
	if os.path.isfile(target+'/band_'+POT+'/'+name+'.png'):
		shutil.copy(target+'/band_'+POT+'/'+name+'.png',res_path+'/band_'+POT+'.png')
		shutil.copy(target+'/band_'+POT+'/'+name+'.pdf',res_path+'/band_'+POT+'.pdf')

	# HSE band gap
	if os.path.isfile(target+'/hybrid_'+POT+'/Band_gap.log'):
		with open(target+'/hybrid_'+POT+'/Band_gap.log','r') as inp:
			gap_log = inp.readline()
		if 'etal' in gap_log:
			hse_gap = 0.0
			hse_eg_type = 'Null'
		else:
			hse_gap = float(gap_log.split()[2])
			hse_eg_type = str(gap_log.split()[4])
		DB['Band_gap_hybrid_'+POT] = hse_gap
		DB['Band_gap_D/I_hybrid_'+POT] = hse_eg_type
		write_formatted_band_log(target+'/hybrid_'+POT+'/Band_gap.log',res_path+'/Band_gap_hybrid_'+POT+'.log')

	for POT2 in pot_list:
		if os.path.isfile(target+'/hybrid_'+POT+'_'+POT2+'/Band_gap.log'):
			with open(target+'/hybrid_'+POT+'_'+POT2+'/Band_gap.log','r') as inp:
				gap_log = inp.readline()
			if 'etal' in gap_log:
				hse_gap = 0.0
				hse_eg_type = 'Null'
			else:
				hse_gap = float(gap_log.split()[2])
				hse_eg_type = str(gap_log.split()[4])
			DB['Band_gap_hybrid_'+POT+'_'+POT2] = hse_gap
			DB['Band_gap_D/I_hybrid_'+POT+'_'+POT2] = hse_eg_type
			write_formatted_band_log(target+'/hybrid_'+POT+'_'+POT2+'/Band_gap.log',res_path+'/Band_gap_hybrid_'+POT+'_'+POT2+'.log')
#			shutil.copy(target+'/hybrid_'+POT+'_'+POT2+'/Band_gap.log',res_path+'/Band_gap_hybrid_'+POT+'_'+POT2+'.log')

	# Corrected band_structure
	if os.path.isfile(target+'/band_'+POT+'/'+name+'_corrected.png'):
		shutil.copy(target+'/band_'+POT+'/'+name+'_corrected.png',res_path+'/band_'+POT+'_corrected.png')
		shutil.copy(target+'/band_'+POT+'/'+name+'_corrected.pdf',res_path+'/band_'+POT+'_corrected.pdf')

	# Density of states
	if os.path.isfile(target+'/dos_'+POT+'/Pdos_dat/dos_'+name+'.png'):
		shutil.copy(target+'/dos_'+POT+'/Pdos_dat/dos_'+name+'.png',res_path+'/dos_'+POT+'.png')
		shutil.copy(target+'/dos_'+POT+'/Pdos_dat/dos_'+name+'.pdf',res_path+'/dos_'+POT+'.pdf')

	# Dielectric constant
	if os.path.isfile(target+'/dielectric_'+POT+'/dielectric.log'):
		with open(target+'/dielectric_'+POT+'/dielectric.log','r') as inp:
			lines = inp.readlines()
		diel_ele = []
		diel_ion = []
		diel_avg = float(lines[9].split()[-1])
		for i in range(3):
			ll1 = lines[1+i].split()
			ll2 = lines[5+i].split()
			diel_ele.append([float(ll1[0]),float(ll1[1]),float(ll1[2])])
			diel_ion.append([float(ll2[0]),float(ll2[1]),float(ll2[2])])
		DB['Dielectric_tensor_'+POT] = diel_avg
		DB['Dielectric_tensor_electronic_'+POT] = diel_ele
		DB['Dielectric_tensor_ionic_'+POT] = diel_ion
		shutil.copy(target+'/dielectric_'+POT+'/dielectric.log',res_path+'/dielectric_'+POT+'.log')

	# Effective mass
	effm_typ = ['hole','electron']
	for typ in effm_typ:
		if os.path.isfile(target+'/effm_'+POT+'/'+typ+'/effective_mass.log'):
			effm_ten = []
			effm_avg = 0.0
			with open(target+'/effm_'+POT+'/'+typ+'/effective_mass.log','r') as inp:
				lines = inp.readlines()

			if lines[4].split() == 3:
				effm = [float(x) for x in lines[4].split()]
			else:
				effm = [float(x) for x in lines[4].split()[-3:]]
			for i in range(3):
				effm_ten.append([float(x) for x in lines[i+1].split()])
				effm_avg = effm_avg+1.0/effm[i]
			effm_avg = 3.0/effm_avg
			DB[typ+'_effective_mass_averaged_'+POT] = effm_avg
			DB[typ+'_effective_mass_diagonalized_'+POT] = effm
			DB[typ+'_effective_mass_tensor_'+POT] = effm_ten
			shutil.copy(target+'/effm_'+POT+'/'+typ+'/effective_mass.log',res_path+'/effective_mass_'+typ+'_'+POT+'.log')

with open(res_path+'/Properties.json','w') as out:
#	out.write(json.dump(DB,indent='\t'))
	out.write(json.dumps(DB,sort_keys=True,indent=4))
