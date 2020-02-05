###########################################
### Date: 2018-12-05			###
### yybbyb@snu.ac.kr			###
###########################################
import shutil, os, sys, subprocess, yaml
import numpy as np
from module_log import *
from module_vasprun import *
from module_band import *
from module_hse import *
from input_conf import set_on_off
from module_relax import *
from _version import __version__
code_data = 'Version '+__version__+'. Modified at 2019-12-17'

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
inp_hse = inp_yaml['hybrid_oneshot']
inp_band = inp_yaml['band_calculation']
node = node_simple(sys.argv[3])
nproc = sys.argv[4]
pot_cell = sys.argv[5]
pot_point = sys.argv[6]

if pot_cell == 'HSE' and pot_point == 'HSE':
	print(1)
	sys.exit()

# Set directory for input structure and INCAR
if pot_cell == pot_point:
	dir_hse = dir+'/hybrid_'+pot_cell
else:
	dir_hse = dir+'/hybrid_'+pot_cell+'_'+pot_point

# Check existing data
if os.path.isdir(dir_hse) and os.path.isfile(dir_hse+'/Band_gap.log') :
	make_amp2_log_default(dir,src_path,'Hybrid oneshot calculation',node,code_data)
	make_amp2_log(dir,'Hybrid oneshot calculation is already done.')
	print(1)
	sys.exit()

if not os.path.isdir(dir_hse):
	os.mkdir(dir_hse,0o755)

# Alpha auto setting (From dielectric constant)
PBE0_on = 0
if inp_hse['alpha'] in ['Auto','auto','AUTO','A','.A.'] or inp_hse['alpha'] == 0:
	if os.path.isdir(dir+'/dielectric_'+pot_point):
		PBE0_on = 1
		diel_path = dir+'/dielectric_'+pot_point
	else:
		inp_hse['alpha'] = 0.25

os.chdir(dir_hse)

make_amp2_log_default(dir_hse,src_path,'Hybrid oneshot calculation',node,code_data)

# check relax calculation
if not os.path.isdir(dir+'/relax_'+pot_cell):
	make_amp2_log(dir_hse,'Relax directory does not exist.')
	print(0)
	sys.exit()
else:
	if not os.path.isfile(dir+'/relax_'+pot_cell+'/CONTCAR'):
		make_amp2_log(dir_hse,'CONTCAR file in relaxation does not exist.')
		print(0)
		sys.exit()
	elif count_line(dir+'/relax_'+pot_cell+'/CONTCAR') < 9:
		make_amp2_log(dir_hse,'CONTCAR file in relaxation is invalid.')
		print(0)
		sys.exit()

if os.path.isfile(dir+'/band_'+pot_point+'/Band_gap.log'):
	with open(dir+'/band_'+pot_point+'/Band_gap.log','r') as inp:
		gap_log = inp.readline()
else:
	make_amp2_log(dir_hse,'Band gap calculation should be performed.')
	print(0)
	sys.exit()

# Odd number electrons and non-magnetic system --> always metal.
spin = pygrep('ISPIN',dir+'/band_'+pot_point+'/OUTCAR',0,0).split()[2]
ncl = pygrep('NONCOL',dir+'/band_'+pot_point+'/OUTCAR',0,0).split()[2]
nelect = pygrep('NELECT',dir+'/band_'+pot_point+'/OUTCAR',0,0).split()[2]
if spin == '1' and int(float(nelect))%2 == 1 and ncl == 'F':
	make_amp2_log(dir_hse,'The band gap cannot be opened in this system.')
	print(1)
	sys.exit()

# Identify the candidates of which band gap can open in hse calculation.
if 'etal' in gap_log:
	dir_dos = dir_hse+'/dos_indicator'
	if not os.path.isdir(dir_dos):
		os.mkdir(dir_dos,0o755)
	if not os.path.isfile(dir_dos+'/DOSCAR') or os.path.getsize(dir_dos+'/DOSCAR') < 5000:
		os.chdir(dir_dos)
		make_amp2_log(dir_hse,'The calculation for DOS indicator starts.')
		copy_input_cont(dir+'/band_'+pot_point,dir_dos)
		subprocess.call(['cp',dir+'/INPUT0/KPOINTS',dir_dos+'/.'])
		subprocess.call(['cp',dir+'/band_'+pot_point+'/CHGCAR',dir_dos+'/.'])
		mag_on = int(spin)-1
		vasprun = make_incar_for_ncl(dir_dos,mag_on,kpar,npar,vasp_std,vasp_gam,vasp_ncl)
		wincar(dir_dos+'/INCAR',dir_dos+'/INCAR',[['NEDOS','3001'],['ISMEAR','0'],['SIGMA','0.05'],['LWAVE','F']],[])
		with open(dir+'/INPUT0/sym','r') as symf:
			sym = int(symf.readline().split()[0])
		make_multiple_kpts(dir+'/kptest/kpoint.log',dir_dos+'/KPOINTS',dir_dos+'/POSCAR',2,sym,'')
		# VASP calculation
		out = run_vasp(dir_dos,nproc,vasprun,mpi)
		if out == 1:  # error in vasp calculation
			print(0)
			sys.exit() 
		out = electronic_step_convergence_check(dir_dos)
		while out == 1:
			make_amp2_log(dir_dos,'Calculation options are changed. New calculation starts.')
			out = run_vasp(dir_dos,nproc,vasprun,mpi)
			if out == 1:  # error in vasp calculation
				print(0)
				sys.exit()
			out = electronic_step_convergence_check(dir_dos)

		if out == 2:  # electronic step is not converged. (algo = normal)
			make_amp2_log(dir_dos,'The calculation stops but electronic step is not converged.')
			print(0)
			sys.exit()
		os.chdir(dir_hse)

	DF_DVB = round(DOS_ratio_fermi_to_vb(dir_dos+'/DOSCAR',inp_hse['fermi_width'],[inp_hse['vb_dos_min'],inp_hse['vb_dos_max']]),4)
	if DF_DVB < inp_hse['cutoff_df_dvb']:
		make_amp2_log(dir_hse,'DF/DVB is '+str(DF_DVB)+'. Band_gap can open.')
		find_extreme_kpt_for_hse(dir+'/band_'+pot_point,inp_hse['energy_width_for_extreme'],inp_hse['search_space_for_extreme'])
	else:
		make_amp2_log(dir_hse,'DF/DVB is '+str(DF_DVB)+'. It is difficult for band gap to open.')
		print(1)
		sys.exit()

make_amp2_log(dir_hse,'This hybrid calculation is performed in the cell with '+pot_cell+' potential at the VBM and CBM with '+pot_point+' potential.')

# we perform the hse calculation only for the materials of which band gap can open.
if os.path.isfile(dir+'/band_'+pot_point+'/KPT') and count_line(dir+'/band_'+pot_point+'/KPT') > 3 :
	if not os.path.isfile(dir_hse+'/POSCAR'):
		# copy VASP input
		copy_input_cont(dir+'/relax_'+pot_cell,dir_hse)
		subprocess.call(['cp',dir+'/INPUT0/POTCAR_GGA',dir_hse+'/POTCAR'])  #HSE always needs POTCAR_GGA
	# make reduce_KPOINTS
	reduce_kpt =  convergence_check_E(dir)
	make_kpts_for_hse(reduce_kpt+'/IBZKPT',dir+'/band_'+pot_point+'/KPT',dir_hse,'oneshot')
	# make INCAR
	# delete LDAU tag
	wincar(dir_hse+'/INCAR',dir_hse+'/INCAR',[['LDA+U',''],['LDAU',''],['LDAUTYPE',''],['LDAUL',''],['LDAUU',''],['LDAUJ',''],['LDAUPRINT','']],[])
	if PBE0_on == 1:
		inp_hse['alpha'] = calc_alpha_auto(diel_path+'/dielectric.log')
		wincar(dir_hse+'/INCAR',dir_hse+'/INCAR',[['NSW','0'],['ALGO',''],['LMAXMIX',''],['ISYM','3']],['\n\nHybrid calculation:\n   LHFCALC= .T.\n   HFSCREEN = 0.0\n   PRECFOCK = Normal\n   ALGO = ALL\n   AEXX = '+str(inp_hse['alpha'])+'\n'])
	else:
		wincar(dir_hse+'/INCAR',dir_hse+'/INCAR',[['NSW','0'],['ALGO',''],['LMAXMIX',''],['ISYM','3']],['\n\nHybrid calculation:\n   LHFCALC= .T.\n   HFSCREEN = 0.2\n   PRECFOCK = Normal\n   ALGO = ALL\n   AEXX = '+str(inp_hse['alpha'])+'\n'])
	incar_from_yaml(dir_hse,inp_hse['incar'])
	mag_on = check_magnet(dir+'/relax_'+pot_cell,inp_yaml['magnetic_ordering']['minimum_moment'])
	vasprun = make_incar_for_ncl(dir_hse,mag_on,kpar,npar,vasp_std,vasp_gam,vasp_ncl)
	make_amp2_log(dir_hse,'Run VASP calculation.')
	# VASP calculation for HSE
	out = run_vasp(dir_hse,nproc,vasprun,mpi)
	if out == 1:  # error in vasp calculation
		print(0)
		sys.exit() 
	out = electronic_step_convergence_check(dir_hse)
	if not out == 0:
		make_amp2_log(dir_hse,'Electronic step is not converged.')
		print(0)
		sys.exit()

	# Extract data from VASP output files
	spin = pygrep('ISPIN',dir_hse+'/OUTCAR',0,0).split()[2]
	ncl = pygrep('NONCOL',dir_hse+'/OUTCAR',0,0).split()[2]
	[KPT,Band,nelect] = EIGEN_to_array(dir_hse+'/EIGENVAL',spin)
	fermi = get_fermi_level(Band,nelect,ncl)

	# Band gap calculation
	gap = gap_estimation(dir_hse,fermi,spin,ncl,KPT,Band,nelect)
	if gap == 'metal':
		make_amp2_log(dir_hse,'HSE band gap calculation is done but it is metallic.')
		gap = '0.0'
	else:
		make_amp2_log(dir_hse,'HSE band gap calculation is done.\nBand gap is '+gap)

	# Small gap correction
	if set_on_off(inp_hse['band_structure_correction']) == 1 and float(gap) > 0.01:
		dir_band = dir+'/band_'+pot_point
		# Insulator in GGA. Band reodering is not required.
		if not 'etal' in gap_log:
			ncl = pygrep('NONCOL',dir_band+'/OUTCAR',0,0).split()[2]
			spin = pygrep('ISPIN',dir_band+'/OUTCAR',0,0).split()[2]
			[KPT,Band,nelect] = EIGEN_to_array(dir_band+'/EIGENVAL',spin)
			fermi = get_fermi_level(Band,nelect,ncl)
			[vb_idx,cb_idx,eVBM,eCBM] = find_cb_gap(Band,fermi,dir_band)
			E_shift = float(gap)+eVBM-eCBM
			for i in range(len(Band[0][0])):
				for n in cb_idx[i]:
					for k in range(len(KPT)):
						Band[n][k][i] = Band[n][k][i] + E_shift
			plot_band_corrected_structure(spin,Band,eVBM,dir_band+'/xtic.dat',dir_band+'/xlabel.dat',[inp_band['y_min'],inp_band['y_max']+E_shift],dir_band)
			if inp_yaml['calculation']['plot'] == 1:
				os.chdir(dir_band)
				subprocess.call([gnuplot,dir_band+'/band_corrected.in'])

		# Metal in GGA. Band reordering is required.
		else:
			# WAVECAR is required for band reordering.
			if os.path.isfile(dir_band+'/WAVECAR') and os.path.getsize(dir_band+'/WAVECAR') > 1000:
				# Make new xtic including all k-points in band caculation
				if int(pyhead(dir_band+'/KPOINTS_band',2).splitlines()[-1].split()[0]) == count_line(dir_band+'/xtic.dat'):
					shutil.copy(dir_band+'/xtic.dat',dir_band+'/xtic_hse.dat')
				else:
					[symk,order,xticlabel,rec] = make_symk(dir_band+'/sym')
					make_xtic_hse(symk,order,rec,inp_band['kspacing_for_band'],dir_band)

				ncl = pygrep('NONCOL',dir_band+'/OUTCAR',0,0).split()[2]
				spin = pygrep('ISPIN',dir_band+'/OUTCAR',0,0).split()[2]
				[KPT,Band,nelect] = EIGEN_to_array(dir_band+'/EIGENVAL',spin)
				fermi = get_fermi_level(Band,nelect,ncl)
				Band_reorder = get_band_reorder(Band,KPT,fermi,spin,dir_band)
				[vb_idx,cb_idx,eVBM,eCBM] = find_cb(Band,Band_reorder,KPT,fermi,dir_hse,dir_band)
				if isinstance(vb_idx,list):
					E_shift = float(gap)+eVBM-eCBM
					num_kpt_for_image = count_line(dir_band+'/xtic.dat')
					Band_reorder = np.array(Band_reorder)
					Band_reorder = Band_reorder[:,:num_kpt_for_image,:]
					for i in range(len(Band[0][0])):
						for n in cb_idx[i]:
							for k in range(len(KPT[:num_kpt_for_image])):
								Band_reorder[n][k][i] = Band_reorder[n][k][i] + E_shift
					fermi = get_fermi_level(Band_reorder,nelect,ncl)
					if calc_gap(fermi,spin,ncl,KPT[:num_kpt_for_image],Band_reorder,nelect) > 0.01:
						plot_band_corrected_structure(spin,Band_reorder,eVBM,dir_band+'/xtic.dat',dir_band+'/xlabel.dat',[inp_band['y_min'],inp_band['y_max']+float(gap)],dir_band)
						if inp_yaml['calculation']['plot'] == 1:
							os.chdir(dir_band)
							subprocess.call([gnuplot,dir_band+'/band_corrected.in'])
					else:
						make_amp2_log(dir_hse,'Warning. We cannot correct the band structure.')
				else:
					make_amp2_log(dir_hse,'Warning. We cannot correct the band structure.')

			else:
				make_amp2_log(dir_hse,'WAVECAR file is missing in band directory')
else:
	if 'etal' in gap_log:
		make_amp2_log(dir_hse,'It is metallic band structure.')
	else:		
		make_amp2_log(dir_hse,'KPT file is empty or do not exist.')
		print(0)
		sys.exit()

with open(dir_hse+'/amp2.log','r') as amp2_log:
	with open(dir+'/amp2.log','a') as amp2_log_tot:
		amp2_log_tot.write(amp2_log.read())

print(1)

