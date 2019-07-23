###########################################
### Date: 2018-12-05			###
### yybbyb@snu.ac.kr			###
###########################################
import shutil, os, sys, subprocess, yaml
from module_log import *
from module_vasprun import *
from module_band import *
from module_hse import *
from module_relax import *
code_data = 'Version xx. Modified at 2019-07-16'

dir = sys.argv[1]

inp_file = sys.argv[2]
with open(inp_file,'r') as f:
	inp_yaml = yaml.load(f)
ERROR_path = inp_yaml['directory']['error']
src_path = inp_yaml['directory']['src_path']
vasp_std = inp_yaml['program']['vasp_std']
vasp_gam = inp_yaml['program']['vasp_gam']
vasp_ncl = inp_yaml['program']['vasp_ncl']
gnuplot = inp_yaml['program']['gnuplot']
npar = inp_yaml['vasp_parallel']['npar']
kpar = inp_yaml['vasp_parallel']['kpar']
inp_hse = inp_yaml['Hybrid_oneshot']
inp_band = inp_yaml['band_calculation']
pot_type = 'GGA'
if isinstance(pot_type,list):
	if len(pot_type) == 1:
		pot_cell = pot_type[0]
		pot_point = pot_type[0]
	else:
		pot_cell = pot_type[0]
		pot_point = pot_type[1]
else:
	pot_cell = pot_type
	pot_point = pot_type

# Set directory for input structure and INCAR
dir_hse = dir+'/hybrid_'+pot_cell+'_'+pot_point

os.chdir(dir_hse)

# check relax calculation
if not os.path.isdir(dir+'/relax_'+pot_cell):
	make_amp2_log(dir_hse,'Relax directory does not exist.')
	print 0
	sys.exit()
else:
	if not os.path.isfile(dir+'/relax_'+pot_cell+'/CONTCAR'):
		make_amp2_log(dir_hse,'CONTCAR file in relaxation does not exist.')
		print 0
		sys.exit()
	elif count_line(dir+'/relax_'+pot_cell+'/CONTCAR') < 9:
		make_amp2_log(dir_hse,'CONTCAR file in relaxation is invalid.')
		print 0
		sys.exit()

if os.path.isfile(dir+'/band_'+pot_point+'/Band_gap.log'):
	with open(dir+'/band_'+pot_point+'/Band_gap.log','r') as inp:
		gap_log = inp.readline()
else:
	make_amp2_log(dir_hse,'Band gap calculation should be performed.')
	print 0
	sys.exit()

# Identify the candidates of which band gap can open in hse calculation.
if 'etal' in gap_log and os.path.isdir(dir+'/dos_'+pot_point) and os.path.isdir(dir+'/dos_'+pot_point+'/Pdos_dat'):
	DF_DVB = round(DOS_ratio_fermi_to_vb(dir+'/dos_'+pot_point+'/DOSCAR',inp_hse['fermi_width'],[inp_hse['vb_dos_min'],inp_hse['vb_dos_max']]),4)
	if DF_DVB < inp_hse['cutoff_DF_DVB']:
		make_amp2_log(dir_hse,'DF/DVB is '+str(DF_DVB)+'. Band_gap can open.')
		find_extreme_kpt_for_hse(dir+'/band_GGA',inp_hse['energy_width_for_extreme'],inp_hse['search_space_for_extreme'])
	else:
		make_amp2_log(dir_hse,'DF/DVB is '+str(DF_DVB)+'. It is difficult for band gap to open.')
		print 1
		sys.exit()

# we perform the hse calculation only for the materials of which band gap can open.
if os.path.isfile(dir+'/band_'+pot_point+'/KPT') and count_line(dir+'/band_'+pot_point+'/KPT') > 3  :
	# Extract data from VASP output files
	spin = subprocess.check_output(['grep','ISPIN',dir_hse+'/OUTCAR']).split()[2]
	ncl = subprocess.check_output(['grep','NONCOL',dir_hse+'/OUTCAR']).split()[2]
	[KPT,Band,nelect] = EIGEN_to_array(dir_hse+'/EIGENVAL',spin)
	fermi = get_fermi_level(Band,nelect,ncl)

	# Band gap calculation
	gap = gap_estimation(dir_hse,fermi,spin,ncl,KPT,Band,nelect)
	make_amp2_log(dir_hse,'HSE band gap calculaton is done.\nBand gap is '+gap)

	# Small gap correction
	if inp_hse['band_structure_correction'] == 1 and float(gap) > 0.01:
		dir_band = dir+'/band_'+pot_point
		# Insulator in GGA. Band reodering is not required.
		if not 'etal' in gap_log:
			ncl = subprocess.check_output(['grep','NONCOL',dir_band+'/OUTCAR']).split()[2]
			spin = subprocess.check_output(['grep','ISPIN',dir_band+'/OUTCAR']).split()[2]
			[KPT,Band,nelect] = EIGEN_to_array(dir_band+'/EIGENVAL',spin)
			fermi = get_fermi_level(Band,nelect,ncl)
			[vb_idx,cb_idx,eVBM,eCBM] = find_cb(Band,Band,KPT,fermi,dir_hse,dir_band)
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
				if int(subprocess.check_output(['head',dir_band+'/KPOINTS_band','-n','2']).splitlines()[-1].split()[0]) == count_line(dir_band+'/xtic.dat'):
					shutil.copy(dir_band+'/xtic.dat',dir_band+'/xtic_hse.dat')
				else:
					[symk,order,xticlabel,rec] = make_symk(dir_band+'/sym')
					make_xtic_hse(symk,order,rec,inp_band['kspacing_for_band'],dir_band)

				ncl = subprocess.check_output(['grep','NONCOL',dir_band+'/OUTCAR']).split()[2]
				spin = subprocess.check_output(['grep','ISPIN',dir_band+'/OUTCAR']).split()[2]
				[KPT,Band,nelect] = EIGEN_to_array(dir_band+'/EIGENVAL',spin)
				fermi = get_fermi_level(Band,nelect,ncl)
				Band_reorder = get_band_reorder(Band,KPT,fermi,spin,dir_band)
				[vb_idx,cb_idx,eVBM,eCBM] = find_cb(Band,Band_reorder,KPT,fermi,dir_hse,dir_band)

				E_shift = float(gap)+eVBM-eCBM

				for i in range(len(Band[0][0])):
					for n in cb_idx[i]:
						for k in range(len(KPT)):
							Band_reorder[n][k][i] = Band_reorder[n][k][i] + E_shift

				fermi = get_fermi_level(Band_reorder,nelect,ncl)
				if calc_gap(fermi,spin,ncl,KPT,Band_reorder,nelect) > 0.01:
					plot_band_corrected_structure(spin,Band_reorder,eVBM,dir_band+'/xtic.dat',dir_band+'/xlabel.dat',[inp_band['y_min'],inp_band['y_max']+float(gap)],dir_band)
					if inp_yaml['calculation']['plot'] == 1:
						os.chdir(dir_band)
						subprocess.call([gnuplot,dir_band+'/band_corrected.in'])
				else:
					make_amp2_log(dir_hse,'ERROR. We cannot correct the band structure.')
			else:
				make_amp2_log(dir_hse,'WAVECAR file is missing in band directory')
else:
	if 'etal' in gap_log:
		make_amp2_log(dir_hse,'It is metallic band structure.')
	else:		
		make_amp2_log(dir_hse,'KPT file is empty or do not exist.')
		print 0
		sys.exit()

with open(dir_hse+'/amp2.log','r') as amp2_log:
	with open(dir+'/amp2.log','a') as amp2_log_tot:
		amp2_log_tot.write(amp2_log.read())

print 1

