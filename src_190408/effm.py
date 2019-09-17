###########################################
### Date: 2018-12-05			###
### yybbyb@snu.ac.kr			###
###########################################
import shutil, os, sys, subprocess, yaml
from module_log import *
from module_vasprun import *
from module_effm import *
code_data = 'Version xx. Modified at 2019-08-07.'

dir = sys.argv[1]

inp_file = sys.argv[2]
with open(inp_file,'r') as f:
	inp_yaml = yaml.safe_load(f)
ERROR_path = inp_yaml['directory']['error']
src_path = inp_yaml['directory']['src_path']
vasp_std = inp_yaml['program']['vasp_std']
vasp_gam = inp_yaml['program']['vasp_gam']
vasp_ncl = inp_yaml['program']['vasp_ncl']
mpi = inp_yaml['program']['mpi_command']
gnuplot = inp_yaml['program']['gnuplot']
npar = inp_yaml['vasp_parallel']['npar']
kpar = inp_yaml['vasp_parallel']['kpar']
inp_effm = inp_yaml['effective_mass']

node = node_simple(sys.argv[3])
nproc = sys.argv[4]
pot_type = sys.argv[5]
if pot_type == 'LDA':
	POT = 'LDA'
else:
	POT = 'GGA'
carrier_type = sys.argv[6]
max_E_diff = calc_max_E_diff(inp_effm['fermi_for_cutoff'],inp_effm['temperature_for_fermi'])

# Set directory for input structure and INCAR
dir_effm = dir+'/effm_'+pot_type+'/'+carrier_type

# Check existing data
if os.path.isdir(dir_effm) and os.path.isfile(dir_effm+'/effective_mass.log') :
	make_amp2_log_default(dir,src_path,'Effective mass for '+carrier_type+' calculation with '+pot_type,node,code_data)
	make_amp2_log(dir,'Effective mass for '+carrier_type+' calculation with '+pot_type+' is already done.')
	print 1
	sys.exit()

if not os.path.isdir(dir+'/effm_'+pot_type):
	os.mkdir(dir+'/effm_'+pot_type,0755)
if not os.path.isdir(dir_effm):
	os.mkdir(dir_effm,0755)

os.chdir(dir_effm)

make_amp2_log_default(dir_effm,src_path,'Effective mass for '+carrier_type+' calculation with '+pot_type,node,code_data)
########
if pot_type == 'HSE':
	make_amp2_log(dir_effm,'Effective mass calculation with HSE is not supported yet.')
	print 1
	sys.exit()

########
with open(dir+'/INPUT0/INCAR','r') as inp:
	tmp = inp.read()
if 'LSORBIT' in tmp:
	make_amp2_log(dir_effm,'Effective mass calculation with SOC is not supported yet.')
	print 1
	sys.exit()

# check band calculation
if not os.path.isdir(dir+'/band_'+pot_type):
	make_amp2_log(dir_effm,'Band directory does not exist.')
	print 0
	sys.exit()
else:
	if not os.path.isfile(dir+'/band_'+pot_type+'/CONTCAR'):
		make_amp2_log(dir_effm,'CONTCAR file in band does not exist.')
		print 0
		sys.exit()
	elif count_line(dir+'/band_'+pot_type+'/CONTCAR') < 9:
		make_amp2_log(dir_effm,'CONTCAR file in band is invalid.')
		print 0
		sys.exit()

	if os.path.isfile(dir+'/band_'+pot_type+'/Band_gap.log'):
		with open(dir+'/band_'+pot_type+'/Band_gap.log','r') as inp:
			gap_log = inp.readline()
	else:
		make_amp2_log(dir_effm,'Band gap calculation should be performed.')
		print 0
		sys.exit()



if 'etal' in gap_log:
	make_amp2_log(dir_effm,'It is metallic band structure.\nEffective mass cannot be estimated in metallic system.')
	print 1
	sys.exit()
else:
	if not pot_type == 'HSE':
		# check CHGCAR
		if os.path.isfile(dir_effm+'/CHGCAR') and os.path.getsize(dir_effm+'/CHGCAR') > 0 :
			make_amp2_log(dir_effm,'Effective mass calculation is performed by using existing CHGCAR file.')
		elif os.path.isfile(dir+'/band_'+pot_type+'/CHGCAR') and os.path.getsize(dir+'/band_'+pot_type+'/CHGCAR') > 0 :
			copy_input_cont(dir+'/band_'+pot_type,dir_effm)
			shutil.copy(dir+'/INPUT0/INCAR',dir_effm+'/INCAR')
			subprocess.call(['cp',dir+'/band_'+pot_type+'/CHGCAR',dir_effm+'/.'])
			make_amp2_log(dir_effm,'Effective mass calculation is performed by using existing CHGCAR file in band calculation.')
		else:
			make_amp2_log(dir_effm,'VASP calculation is conducted for CHGCAR file.')
			# Copy input data and write CHGCAR
			copy_input_cont(dir+'/band_'+pot_type,dir_effm)
			shutil.copy(dir+'/INPUT0/INCAR',dir_effm+'/INCAR')
			subprocess.call(['cp',dir+'/INPUT0/KPOINTS',dir_effm+'/.'])
			# make INCAR for CHGCAR
			incar_from_yaml(dir_effm,inp_effm['incar'])
			if os.path.isfile(dir+'/band_'+pot_type+'/OUTCAR'):
				mag_on = check_magnet(dir+'/band_'+pot_type,inp_yaml['magnetic_ordering']['minimum_moment'])
			else:
				mag_on = 2
			vasprun = make_incar_for_ncl(dir_effm,mag_on,kpar,npar,vasp_std,vasp_gam,vasp_ncl)
			incar_from_yaml(dir_effm,inp_effm['incar'])
			wincar(dir_effm+'/INCAR',dir_effm+'/INCAR',[['NSW','0'],['LCHARG','.T.']],[])

			# VASP calculation for CHGCAR
			out = run_vasp(dir_effm,nproc,vasprun,mpi)
			if out == 1:  # error in vasp calculation
				print 0
				sys.exit() 
			make_amp2_log(dir_effm,'CHGCAR file is generated successfully.')

	# make pocket.log
	if os.path.isfile(dir_effm+'/pocket.log') and os.path.getsize(dir_effm+'/pocket.log') > 0 :
		make_amp2_log(dir_effm,'Existing pocket.log is used.')
	else:
		pocket(dir_effm,dir+'/band_'+pot_type,carrier_type,max_E_diff,inp_effm['pocket_search_dist'])
		make_amp2_log(dir_effm,'New pocket.log is generated.')

	# make k-points to determine searching space
	if os.path.isfile(dir_effm+'/NKPT') and os.path.getsize(dir_effm+'/NKPT') > 0 :
		make_amp2_log(dir_effm,'Existing NKPT is used.')
	else:
		make_kpts_for_searching_space(dir_effm,inp_effm['grid_size_for_searching'])
		if pot_type == 'HSE':
			make_kpts_for_hse(dir+'/relax_'+pot_type+'/IBZKPT',dir_effm+'/KPOINTS_searching',dir_effm,'band')
		else:
			shutil.copy(dir_effm+'/KPOINTS_searching',dir_effm+'/KPOINTS')
		make_amp2_log(dir_effm,'New NKPT is generated.')

	# make k-points to calculate effective mass	
	if os.path.isfile(dir_effm+'/inp_grid') and os.path.getsize(dir_effm+'/inp_grid') > 0 :
		make_amp2_log(dir_effm,'Existing inp_grid is used.')
	else:
		# check the vasp calculation to determine searching space
		if not check_vasp_done(dir_effm) == 1:
			if os.path.isfile(dir+'/band_'+pot_type+'/OUTCAR'):
				mag_on = check_magnet(dir+'/band_'+pot_type,inp_yaml['magnetic_ordering']['minimum_moment'])
			else:
				mag_on = 2
			vasprun = make_incar_for_ncl(dir_effm,mag_on,kpar,npar,vasp_std,vasp_gam,vasp_ncl)
			if pot_type == 'HSE':
				incar_for_hse(dir_effm+'/INCAR')
				wincar(dir_effm+'/INCAR',dir_effm+'/INCAR',[['NSW','0'],['ALGO','All']],[])
				incar_from_yaml(dir_effm,inp_effm['incar'])
			else:
				incar_from_yaml(dir_effm,inp_effm['incar'])
				wincar(dir_effm+'/INCAR',dir_effm+'/INCAR',[['NSW','0'],['ISTART','1'],['ICHARG','11'],['LCHARG','.F.']],[])
			out = run_vasp(dir_effm,nproc,vasprun,mpi)
			if out == 1:  # error in vasp calculation
				print 0
				sys.exit() 
		else:
			make_amp2_log(dir_effm,'VASP calculation is already done.')
		make_kpts_for_calculation(dir_effm,inp_effm['grid_size_for_calculation'],max_E_diff)
		make_amp2_log(dir_effm,'New inp_grid is generated.')

	# check the vasp calculation to calculate effective mass
	if not check_vasp_done(dir_effm) == 1:
		if os.path.isfile(dir+'/band_'+pot_type+'/OUTCAR'):
			mag_on = check_magnet(dir+'/band_'+pot_type,inp_yaml['magnetic_ordering']['minimum_moment'])
		else:
			mag_on = 2
		vasprun = make_incar_for_ncl(dir_effm,mag_on,kpar,npar,vasp_std,vasp_gam,vasp_ncl)
		incar_from_yaml(dir_effm,inp_effm['incar'])
		wincar(dir_effm+'/INCAR',dir_effm+'/INCAR',[['NSW','0'],['ISTART','1'],['ICHARG','11'],['LCHARG','.F.']],[])
		out = run_vasp(dir_effm,nproc,vasprun,mpi)
		if out == 1:  # error in vasp calculation
			print 0
			sys.exit() 
	else:
		make_amp2_log(dir_effm,'VASP calculation is already done.')

	if os.path.isfile(dir+'/INPUT0/POSCAR_rlx_'+POT):
		oper = read_operation(dir+'/INPUT0/POSCAR_rlx_'+POT)
	else:
		oper = read_operation(dir+'/INPUT0/POSCAR')
	[effm_dia,effm] = calc_effm(dir_effm,carrier_type,inp_effm['temperature_for_fermi'],oper)
	write_effm(effm_dia,effm,dir_effm,carrier_type)

	make_amp2_log(dir_effm,'Effective mass for '+carrier_type+' calculation is done.')

with open(dir_effm+'/amp2.log','r') as amp2_log:
	with open(dir+'/amp2.log','a') as amp2_log_tot:
		amp2_log_tot.write(amp2_log.read())

print 1
