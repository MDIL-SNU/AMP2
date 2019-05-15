###########################################
### Date: 2018-12-05			###
### yybbyb@snu.ac.kr			###
###########################################
import shutil, os, sys, subprocess, yaml
from module_log import *
from module_vasprun import *
from module_band import *
from module_hse import *
code_data = 'Date: 2018-12-05'

dir = sys.argv[1]

inp_file = sys.argv[2]
with open(inp_file,'r') as f:
	inp_yaml = yaml.load(f)
ERROR_path = inp_yaml['directory']['error']
src_path = inp_yaml['directory']['src_path']
vasp_std = inp_yaml['program']['vasp_std']
vasp_gam = inp_yaml['program']['vasp_gam']
vasp_ncl = inp_yaml['program']['vasp_ncl']
npar = inp_yaml['vasp_parallel']['npar']
kpar = inp_yaml['vasp_parallel']['kpar']
inp_hse = inp_yaml['Hybrid_oneshot']

node = node_simple(sys.argv[3])
nproc = sys.argv[4]

# Set directory for input structure and INCAR
dir_hse = dir+'/HSE'

# Check existing data
if os.path.isdir(dir_hse) and os.path.isfile(dir_hse+'/Band_gap.log') :
	make_amp2_log_default(dir,src_path,'Hybrid oneshot calculation',node,code_data)
	make_amp2_log(dir,'Hybrid oneshot calculation is already done.')
	print 1
	sys.exit()

if not os.path.isdir(dir_hse):
	os.mkdir(dir_hse,0755)

os.chdir(dir_hse)

make_amp2_log_default(dir_hse,src_path,'Hybrid oneshot calculation',node,code_data)

# check relax calculation
if not os.path.isdir(dir+'/relax_GGA'):
	make_amp2_log(dir_hse,'Relax directory does not exist.')
	print 0
	sys.exit()
else:
	if not os.path.isfile(dir+'/relax_GGA/CONTCAR'):
		make_amp2_log(dir_hse,'CONTCAR file in relaxation does not exist.')
		print 0
		sys.exit()
	elif count_line(dir+'/relax_GGA/CONTCAR') < 9:
		make_amp2_log(dir_hse,'CONTCAR file in relaxation is invalid.')
		print 0
		sys.exit()

if os.path.isfile(dir+'/band_GGA/Band_gap.log'):
	with open(dir+'/band_GGA/Band_gap.log','r') as inp:
		gap_log = inp.readline()
else:
	make_amp2_log(dir_hse,'Band gap calculation should be performed.')
	print 0
	sys.exit()

# Identify the candidates of which band gap can open in hse calculation.
if os.path.isdir(dir+'/dos_GGA') and os.path.isdir(dir+'/dos_GGA/Pdos_dat'):
	DF_DVB = DOS_ratio_fermi_to_vb(dir+'/dos_GGA/DOSCAR',inp_hse['fermi_width'],[inp_hse['vb_dos_min'],inp_hse['vb_dos_max']])
	if DF_DVB < inp_hse['cutoff_DF_DVB']:
		make_amp2_log(dir_hse,'DF/DVB is '+str(DF_DVB)+'. Band_gap can open.')
		find_extreme_kpt_for_hse(dir+'/band_GGA',inp_hse['energy_width_for_extreme'],inp_hse['search_space_for_extreme'])
	else:
		make_amp2_log(dir_hse,'DF/DVB is '+str(DF_DVB)+'. It is difficult for band gap to open.')


if os.path.isfile(dir+'/band_GGA/KPT') and os.path.getsize(dir+'/band_GGA/KPT') > 0 :
	# copy VASP input
	copy_input_cont(dir+'/relax_GGA',dir_hse)
	# make KPOINTS
	make_kpts_for_oneshot(dir+'/relax_GGA/IBZKPT',dir+'/band_GGA/KPT',dir_hse)
	# make INCAR
	wincar(dir_hse+'/INCAR',dir_hse+'/INCAR',[['NSW','0'],['ALGO',''],['LDA',''],['LMAXMIX','']],['\n\nHybrid calculation:\n   LHFCALC= .T.\n   HFSCREEN = 0.2\n   PRECFOCK = Normal\n   ALGO = ALL\n   AEXX = '+str(inp_hse['alpha'])+'\n'])
	incar_from_yaml(dir_hse,inp_hse['INCAR'])
	mag_on = check_magnet(dir+'/relax_GGA')
	vasprun = make_incar_for_ncl(dir_hse,mag_on,kpar,npar,vasp_std,vasp_gam,vasp_ncl)
	make_amp2_log(dir_hse,'Run VASP calculation.')
	# VASP calculation for HSE
	out = run_vasp(dir_hse,nproc,vasprun)
	if out == 1:  # error in vasp calculation
		print 0
		sys.exit() 

	# Extract data from VASP output files
	fermi = float(subprocess.check_output(['head',dir_hse+'/DOSCAR','-n','6']).splitlines()[-1].split()[3])
	spin = subprocess.check_output(['grep','ISPIN',dir_hse+'/OUTCAR']).split()[2]
	ncl = subprocess.check_output(['grep','NONCOL',dir_hse+'/OUTCAR']).split()[2]
	[KPT,Band,nelect] = EIGEN_to_array(dir_hse+'/EIGENVAL',spin)

	# Band gap calculation
	gap = gap_estimation_hse(dir_hse,fermi,spin,ncl,KPT,Band,nelect)
	make_amp2_log(dir_hse,'HSE band gap calculaton is done.\nBand gap is '+gap)

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

