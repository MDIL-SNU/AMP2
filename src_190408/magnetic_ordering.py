####################################
# date : 2019-04-08                #
# Author : yybbyb@snu.ac.kr        #
####################################
import os, sys, subprocess, yaml, glob, shutil
from module_log import *
from module_vasprun import *
from module_converge import *
from module_AF import *
from module_amp2_input import *
from module_relax import *
from input_conf import set_on_off
code_data = 'Version xx. Modified at 2019-07-17'

# Set input
dir = sys.argv[1]

inp_file = sys.argv[2]
with open(inp_file,'r') as f:
	inp_yaml = yaml.safe_load(f)
ERROR_path = inp_yaml['directory']['error']
src_path = inp_yaml['directory']['src_path']
vasp_std = inp_yaml['program']['vasp_std']
vasp_gam = inp_yaml['program']['vasp_gam']
mpi = inp_yaml['program']['mpi_command']
gnuplot = inp_yaml['program']['gnuplot']
npar = inp_yaml['vasp_parallel']['npar']
kpar = inp_yaml['vasp_parallel']['kpar']

max_nelm = str(inp_yaml['cif2vasp']['max_nelm'] * 2)

inp_af = inp_yaml['magnetic_ordering']
inp_rlx = inp_yaml['relaxation']
cutoff_length = inp_af['cutoff_for_parameter']
pot_type = inp_af['potential_type']
if pot_type == 'LDA':
	POT = 'LDA'
else:
	POT = 'GGA'

node = node_simple(sys.argv[3])
nproc = sys.argv[4]

# Check existing data
if os.path.isdir(dir+'/magnetic_ordering') and os.path.isfile(dir+'/magnetic_ordering/POSCAR_spin') :
	make_amp2_log_default(dir,src_path,'Magnetic ordering',node,code_data)
#	print('Success!')
	make_amp2_log(dir,'Already done')
	print 1
	sys.exit()

if not os.path.isdir(dir+'/magnetic_ordering') :
	os.mkdir(dir+'/magnetic_ordering', 0755)
	make_amp2_log_default(dir+'/magnetic_ordering',src_path,'Magnetic ordering',node,code_data)
	make_amp2_log(dir+'/magnetic_ordering','New calculation.')
else:
	make_amp2_log_default(dir+'/magnetic_ordering',src_path,'Magnetic ordering',node,code_data)
	make_amp2_log(dir+'/magnetic_ordering','Calculation continue.')


os.chdir(dir+'/magnetic_ordering')
target = dir+'/magnetic_ordering'

### Check magnetic elements
# define path for OUTCAR and POSCAR (atom type should be including (! xx1))
if set_on_off(inp_af['from_relax']) == 1:
	if os.path.isfile(dir+'/INPUT0/POSCAR_rlx_'+pot_type):
		make_amp2_log(target,'Magnetic elements are determined by OUTCAR in relax_'+pot_type+'.')
		ref_dir = dir+'/relax_'+pot_type
		inp_pos = dir+'/INPUT0/POSCAR_rlx_'+pot_type
	else:
		make_amp2_log(target,'Cannot find relaxed POSCAR.')
		print 0
		sys.exit()
else:
	ref_dir = ''
	inp_pos = dir+'/INPUT0/POSCAR'
	if os.path.isdir(dir+'/cutoff') and os.path.isfile(dir+'/cutoff/cutoff.log'):
		with open(dir+'/cutoff/cutoff.log','r') as cutoff_log:
			line = cutoff_log.readlines()[-1].split()
		if line[0] == 'Converged':
			ref_dir = dir+'/cutoff/EN'+str(line[2])
			make_amp2_log(target,'Magnetic elements are determined by OUTCAR in encut test.')
	if ref_dir == '':
		with open(dir+'/kptest/kpoint.log','r') as kpoint_log:
			line = kpoint_log.readlines()[-1].split()
		if line[0] == 'Converged':
			ref_dir = dir+'/kptest/KP'+str(line[2])
			make_amp2_log(target,'Magnetic elements are determined by OUTCAR in kptest.')

if ref_dir == '':
	make_amp2_log(target,'Cannot find reference OUTCAR.')
	print 0
	sys.exit()

spin = pygrep('ISPIN',ref_dir+'/OUTCAR',0,0).split()[2]
#spin = subprocess.check_output(['grep','ISPIN',ref_dir+'/OUTCAR']).split()[2]
if spin == '1':
	make_amp2_log(target,'ISPIN was turned off.')
	with open(target+'/amp2.log','r') as amp2_log:
		with open(dir+'/amp2.log','a') as amp2_log_tot:
			amp2_log_tot.write(amp2_log.read())
	print 1
	sys.exit()

[mag_atom_list,mag_val] = check_spin(ref_dir,inp_pos,inp_af['minimum_moment'])
if len(mag_atom_list) == 0:
	make_amp2_log(target,'There is no magnetic atom.')
	with open(target+'/amp2.log','r') as amp2_log:
		with open(dir+'/amp2.log','a') as amp2_log_tot:
			amp2_log_tot.write(amp2_log.read())
	print 1
	sys.exit()
else:
	make_amp2_log(dir+'/magnetic_ordering','Magnetic atoms are '+', '.join(mag_atom_list)+'.')


# Make supercell 
if not os.path.isfile(target+'/POSCAR_param'):
	min_cell_length = cutoff_length*2
	while 1:
		out = int(subprocess.check_output(['python',src_path+'/mk_supercell.py',inp_pos,dir+'/INPUT0/KPOINTS',str(min_cell_length),target+'/POSCAR_param']).split()[0])
		if out == 0:
			break
		else:
			min_cell_length = min_cell_length + 2

[axis_param,pos_param] = read_poscar(target+'/POSCAR_param')
if isinstance(inp_af['maximum_atoms_number'],int) and inp_af['maximum_atoms_number'] > 0 and len(pos_param) > inp_af['maximum_atoms_number']:
	make_amp2_log(target,'The number of atoms in supercell is larger than the criteria.')
	print 0
	sys.exit()

[pair_list,sole_list,mag_list] = find_pair(target+'/POSCAR_param',cutoff_length,mag_atom_list,mag_val,inp_af['tolerance'])
if len(pair_list) == 0:
	make_amp2_log(target,'There is no magnetic atom pair within cutoff radius.')
	with open(target+'/amp2.log','r') as amp2_log:
		with open(dir+'/amp2.log','a') as amp2_log_tot:
			amp2_log_tot.write(amp2_log.read())
	print 1
	sys.exit()

if isinstance(inp_af['maximum_pair_type'],int) and inp_af['maximum_pair_type'] > 0 and len(pair_list) > inp_af['maximum_pair_type']:
	make_amp2_log(target,'The number of pairs is larger than the criteria.')
	print 0
	sys.exit()

tot_mag_list = write_pair_list(pair_list,sole_list,mag_list,target)

for i in range(len(tot_mag_list)):
	targ_dir = target+'/Spin_'+str(i)
	if not os.path.isdir(targ_dir):   # check calculation is done or path exists
		os.mkdir(targ_dir)
	if not os.path.isfile(targ_dir+'/energy') or not os.path.getsize(targ_dir+'/energy') > 1:
		os.chdir(targ_dir)
		copy_input_no_kp(dir+'/INPUT0',targ_dir,POT)
		subprocess.call(['cp',target+'/POSCAR_param',targ_dir+'/POSCAR'])
		subprocess.call(['cp',src_path+'/KPOINTS_gamma',targ_dir+'/KPOINTS'])
		wincar(targ_dir+'/INCAR',targ_dir+'/INCAR',[['LCHARG','F'],['LWAVE','F'],['ALGO','Normal'],['NELM',max_nelm]],[])
		if pot_type == 'HSE':
			incar_for_hse(targ_dir+'/INCAR')
			wincar(targ_dir+'/INCAR',targ_dir+'/INCAR',[['ALGO','All']],[])
		incar_from_yaml(targ_dir,inp_af['incar'])

		wincar(targ_dir+'/INCAR',targ_dir+'/INCAR',[['MAGMOM',tot_mag_list[i]]],[])
		gam = set_parallel(targ_dir+'/KPOINTS',targ_dir+'/INCAR',npar,kpar)
		if gam == 1:
			vasprun = vasp_gam
		else:
			vasprun = vasp_std
		wincar(targ_dir+'/INCAR',targ_dir+'/INCAR',[['NSW','0'],['ISYM','0']],[])

		out = run_vasp(targ_dir,nproc,vasprun,mpi)
		if out == 1:  # error in vasp calculation
			print 0
			sys.exit() 

		out = electronic_step_convergence_check(targ_dir)
		while out == 1:
			make_amp2_log(targ_dir,'Calculation options are changed. New calculation starts.')
			out = run_vasp(targ_dir,nproc,vasprun,mpi)
			if out == 1:  # error in vasp calculation
				print 0
				sys.exit()
			out = electronic_step_convergence_check(targ_dir)

		if out == 2:  # electronic step is not converged. (algo = normal)
			make_amp2_log(targ_dir,'The calculation stops but electronic step is not converged.')
			print 0
			sys.exit()

		energy = pygrep('free  ','OUTCAR',0,0).splitlines()[-1].split()[4]
#		energy = subprocess.check_output(['grep','free  ','OUTCAR']).splitlines()[-1].split()[4]
		with open(targ_dir+'/energy','w') as enef:
			enef.write('Spin_'+str(i)+'\t'+energy+'\n')

	if i == 0:
		with open(target+'/energy','w') as enef:
			with open(targ_dir+'/energy','r') as ener:
				enef.write(ener.read())
	else:
		with open(target+'/energy','a') as enef:
			with open(targ_dir+'/energy','r') as ener:
				enef.write(ener.read())

JJ = calc_ising_param(sole_list,pair_list,target+'/energy')
write_ising_param(JJ,target)
write_inp_for_GA(mag_atom_list,inp_af,target)

### genetic algorithm in supercell
supercell_list = set_supercell_list(inp_pos,mag_atom_list)

min_energy_in_ga = 10.0e5
energy_tolerance = str(inp_af['energy_tolerance']) # eV/atom
optimized_supercell = []
for cell in supercell_list:
	ga_path = target+'/GA_'+str(cell[0])+str(cell[1])+str(cell[2])
	if not os.path.isdir(ga_path):
		os.mkdir(ga_path)
	os.chdir(ga_path)

	# shuld be checked
	if not os.path.isfile(ga_path+'/POSCAR_spin0'):
		subprocess.call(['python',src_path+'/make_supercell.py',inp_pos,str(cell[0]),str(cell[1]),str(cell[2]),ga_path+'/POSCAR_ref'])
		subprocess.call(['python',src_path+'/genetic_algorithm.py',target+'/input_GA.yaml',ga_path+'/POSCAR_ref',energy_tolerance])

	fin_energy = float(pytail(ga_path+'/energy.dat').split()[0])
#	fin_energy = float(subprocess.check_output(['tail','-1',ga_path+'/energy.dat']).split()[0])
	fin_energy = fin_energy/float(cell[0]*cell[1]*cell[2])
	if round(fin_energy-min_energy_in_ga,8) < 0 :
		min_energy_in_ga = fin_energy
		optimized_supercell = cell

## run vasp for stable spin configuration in optimized supercell
make_amp2_log(target,'Smallest cell is in GA_'+''.join([str(x) for x in optimized_supercell])+'.')
ground_path = target+'/GA_'+''.join([str(x) for x in optimized_supercell])
os.chdir(ground_path)
POSCARs = glob.glob(ground_path+'/POSCAR_spin*')
pos_atom_num = []
for i in range(len(POSCARs)):
	if not os.path.isdir(ground_path+'/Stable'+str(i)):
		os.mkdir(ground_path+'/Stable'+str(i))
	[axis,atom_pos] = read_poscar(POSCARs[i])
	[prim_axis,prim_atom_pos,sym] = get_primitive_cell(axis,atom_pos)
	write_file(str(sym),ground_path+'/Stable'+str(i)+'/sym')
	set_mag_info(prim_atom_pos,ground_path+'/Stable'+str(i),mag_val)
	if not os.path.isfile(ground_path+'/Stable'+str(i)+'/POSCAR'):
		write_poscar(prim_axis,prim_atom_pos,ground_path+'/Stable'+str(i)+'/POSCAR','Primitive cell')
	if not os.path.isfile(ground_path+'/Stable'+str(i)+'/POSCAR0'):
		write_poscar(prim_axis,prim_atom_pos,ground_path+'/Stable'+str(i)+'/POSCAR0','Primitive cell')
	pos_atom_num.append(len(prim_atom_pos))


stable_ene = [999,9999999.9]
for i in range(len(POSCARs)):
	calc_path = ground_path+'/Stable'+str(i)
	if not pos_atom_num[i] == min(pos_atom_num):
		make_amp2_log(target,'There is smaller primitive cell than POSCAR_spin'+str(i)+'.')
		shutil.rmtree(calc_path)
	else:
		os.chdir(calc_path)
		if os.path.isfile(calc_path+'/CONTCAR') and os.path.isfile(calc_path+'/free') and len(pygrep('free  ',calc_path+'/free',0,0).splitlines()) > 0 and len(pygrep('free  ',calc_path+'/free',0,0).splitlines()) <= inp_rlx['converged_ionic_step'] :
			energy = pygrep('free  ','OUTCAR',0,0).splitlines()
#			energy = subprocess.check_output(['grep','free  ','OUTCAR']).splitlines()
		else:
			run_cont = 0
			if os.path.isfile(calc_path+'/CONTCAR') and count_line(calc_path+'/CONTCAR') > 9:
				make_amp2_log_default(calc_path,src_path,'Relaxation with '+pot_type+' potential',node,code_data)
				make_amp2_log(calc_path,'Calculation continue.')
				subprocess.call(['cp',calc_path+'/CONTCAR',calc_path+'/POSCAR'])
				run_cont = 1

			if run_cont == 0:
				make_amp2_log_default(calc_path,src_path,'Relaxation with '+pot_type+' potential',node,code_data)
				make_amp2_log(calc_path,'New calculation')
				subprocess.call(['cp',dir+'/INPUT0/INCAR',calc_path+'/INCAR'])
				subprocess.call(['cp',dir+'/INPUT0/POTCAR_'+POT,calc_path+'/POTCAR'])
				resized_kpoints(dir+'/INPUT0',calc_path)
				nsw = set_nsw(calc_path+'/POSCAR',calc_path+'/INCAR')
				converge_condition = set_ediffg(calc_path+'/POSCAR',inp_rlx['force'],inp_rlx['pressure'],inp_rlx['energy'])
				if not converge_condition == 0:
					wincar(calc_path+'/INCAR',calc_path+'/INCAR',[['EDIFFG',str(converge_condition)]],[])
				wincar(calc_path+'/INCAR',calc_path+'/INCAR',[['LCHARG','.F.'],['ALGO','Normal'],['NELM',max_nelm]],[])
				if pot_type == 'HSE':
					incar_for_hse(calc_path+'/INCAR')
					wincar(calc_path+'/INCAR',calc_path+'/INCAR',[['ALGO','All']],[])
				incar_from_yaml(calc_path,inp_rlx['incar'])

				with open(calc_path+'/mag_info','r') as inp:
					mag_info = inp.readline()
				wincar(calc_path+'/INCAR',calc_path+'/INCAR',[['MAGMOM',mag_info]],[])

			gam = set_parallel(calc_path+'/KPOINTS',calc_path+'/INCAR',npar,kpar)
			if gam == 1:
				vasprun = vasp_gam
			else:
				vasprun = vasp_std

			out = run_vasp(calc_path,nproc,vasprun,mpi)
			if out == 1:  # error in vasp calculation
				print 0
				sys.exit() 
			out = electronic_step_convergence_check(calc_path)
			while out == 1:
				make_amp2_log(calc_path,'Calculation options are changed. New calculation starts.')
				out = run_vasp(calc_path,nproc,vasprun,mpi)
				if out == 1:  # error in vasp calculation
					print 0
					sys.exit()
				out = electronic_step_convergence_check(calc_path)

			if out == 2:  # electronic step is not converged. (algo = normal)
				make_amp2_log(calc_path,'The calculation stops but electronic step is not converged.')
				print 0
				sys.exit()

			ionic_converge = 0
			iteration = 1
			energy = pygrep('free  ','OUTCAR',0,0).splitlines()
#			energy = subprocess.check_output(['grep','free  ','OUTCAR']).splitlines()
			while iteration < inp_rlx['max_iteration']:
				with open(calc_path+'/OUT_TOT','a') as out_log:
					out_log.write(open(calc_path+'/OUTCAR','r').read())
				if len(energy) <= inp_rlx['converged_ionic_step']:
					ionic_converge = 1
					make_amp2_log(calc_path,'Relaxation is done.')
					break
				shutil.copyfile(calc_path+'/CONTCAR',calc_path+'/POSCAR')
				if bool(inp_rlx['incar']) and not 'EDIFFG' in inp_rlx['incar'].keys():
					converge_condition = set_ediffg(calc_path+'/POSCAR',inp_rlx['force'],inp_rlx['pressure'],inp_rlx['energy'])
					if not converge_condition == 0:
						wincar(calc_path+'/INCAR',calc_path+'/INCAR',[['EDIFFG',str(converge_condition)]],[])

				out = run_vasp(calc_path,nproc,vasprun,mpi)
				if out == 1:  # error in vasp calculation
					print 0
					sys.exit() 
				out = electronic_step_convergence_check(calc_path)
				while out == 1:
					make_amp2_log(calc_path,'Calculation options are changed. New calculation starts.')
					out = run_vasp(calc_path,nproc,vasprun,mpi)
					if out == 1:  # error in vasp calculation
						print 0
						sys.exit()
					out = electronic_step_convergence_check(calc_path)
				if out == 2:  # electronic step is not converged. (algo = normal)
					make_amp2_log(calc_path,'The calculation stops but electronic step is not converged.')
					print 0
					sys.exit()

				energy = pygrep('free  ','OUTCAR',0,0).splitlines()
#				energy = subprocess.check_output(['grep','free  ','OUTCAR']).splitlines()
				iteration = iteration+1
				make_amp2_log(calc_path,'Iteration number is '+str(iteration))

		if stable_ene[1] > float(energy[-1].split()[4]):
			stable_ene = [i,float(energy[-1].split()[4])]
		with open('free','w') as fr_file:
			fr_file.write(pygrep('free  ','OUTCAR',0,0).splitlines()[-1])
#			fr_file.write(subprocess.check_output(['grep','free  ','OUTCAR']).splitlines()[-1])

with open(target+'/amp2.log','r') as amp2_log:
	with open(dir+'/amp2.log','a') as amp2_log_tot:
		amp2_log_tot.write(amp2_log.read())

## check magnetic ordering (ferro or not)
stable_path = ground_path+'/Stable'+str(stable_ene[0])
shutil.copy(stable_path+'/POSCAR',target+'/POSCAR_spin')
[new_mag_atom_list,new_mag_val] = check_spin(stable_path,stable_path+'/POSCAR0',inp_af['minimum_moment'])
ferro = 0 # 0 indicates Ferromagnetic ordering
if min([new_mag_val[x] for x in new_mag_val])*max([new_mag_val[x] for x in new_mag_val]) < -0.01:
		ferro = 1 # 1 indicates non-ferromagnetic ordering.
if ferro == 0:
	make_amp2_log(dir,'The ground state spin ordering is ferromagnetic ordering.')
	if os.path.isdir(dir+'/relax_'+pot_type):
		if os.path.isfile(dir+'/relax_'+pot_type+'/free'):
			make_amp2_log(dir,'We did not change the files in relax_'+pot_type+' for magnetic calculation.')
		else:
			make_amp2_log(dir,'We changed the files in relax_'+pot_type+' from magnetic calculation.')
			shutil.move(dir+'/relax_'+pot_type,dir+'/relax_'+pot_type+'_old')
			shutil.copytree(stable_path,dir+'/relax_'+pot_type)
	else:
		make_amp2_log(dir,'We changed the files in relax_'+pot_type+' from magnetic calculation.')
		shutil.copytree(stable_path,dir+'/relax_'+pot_type)
else:
	if os.path.isdir(dir+'/relax_'+pot_type):
		make_amp2_log(dir,'We changed the files in relax_'+pot_type+' from magnetic calculation.')
		shutil.move(dir+'/relax_'+pot_type,dir+'/relax_'+pot_type+'_old')
		shutil.copytree(stable_path,dir+'/relax_'+pot_type)
	else:
		make_amp2_log(dir,'We changed the files in relax_'+pot_type+' from magnetic calculation.')
		shutil.copytree(stable_path,dir+'/relax_'+pot_type)
	shutil.copytree(dir+'/INPUT0',dir+'/INPUT0_old')
	shutil.copy(stable_path+'/POSCAR0',dir+'/INPUT0/POSCAR')
	write_relaxed_poscar(dir,pot_type)
	shutil.copy(stable_path+'/INCAR',dir+'/INPUT0/INCAR')
	shutil.copy(stable_path+'/KPOINTS',dir+'/INPUT0/KPOINTS')
	shutil.copy(stable_path+'/sym',dir+'/INPUT0/sym')
	with open(stable_path+'/mag_info','r') as inp:
		mag_info = inp.readline()
	wincar(dir+'/INPUT0/spin_note',dir+'/INPUT0/spin_note',[['MAGMOM',mag_info]],[])

print 1
