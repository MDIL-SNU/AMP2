####################################
# date : 2018-12-06                #
# Author : feihoon@snu.ac.kr       #
# Modifier : yybbyb@snu.ac.kr      #
####################################
import os, sys, subprocess, yaml
from module_log import *
from module_vasprun import *
from module_converge import *
import math
from _version import __version__
code_data = 'Version '+__version__+'. Modified at 2019-12-17'

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

inp_conv = inp_yaml['convergence_test']
ENCONV = inp_conv['enconv']
PRCONV = inp_conv['prconv']
FOCONV = inp_conv['foconv']
ENSTART = inp_conv['enstart']
ENSTEP = inp_conv['enstep']
ENMAX = inp_conv['enmax']
pot_type = inp_conv['potential_type']
if pot_type == 'LDA':
	POT = 'LDA'
else:
	POT = 'GGA'

node = node_simple(sys.argv[3])
nproc = sys.argv[4]

# Check existing data
if os.path.isdir(dir+'/cutoff') and os.path.isfile(dir+'/cutoff/cutoff.log') :
	if len(pygrep('Converged',dir+'/cutoff/cutoff.log',0,0).splitlines()) > 0 :
		make_amp2_log_default(dir,src_path,'cutoff energy test',node,code_data)
		make_amp2_log(dir,'Already done')
#		print('Success!')
		print(1)
		sys.exit()
if not os.path.isdir(dir+'/cutoff') :
	os.mkdir(dir+'/cutoff', 0o755)
	make_amp2_log_default(dir+'/cutoff',src_path,'cutoff energy test',node,code_data)
	make_amp2_log(dir+'/cutoff','New calculation.')
else:
	make_amp2_log_default(dir+'/cutoff',src_path,'cutoff energy test',node,code_data)
	make_amp2_log(dir+'/cutoff','Calculation continue.')

os.chdir(dir+'/cutoff')
enlog = dir+'/cutoff/cutoff.log'
if os.path.isfile(enlog):
	os.remove(enlog)

### cutoff energy test setting ###
ENPOT = pygrep('ENMAX',dir+'/INPUT0/POTCAR_'+POT,0,0).splitlines()
EN_recommend = max([float(ENPOT[x].split()[2].split(';')[0]) for x in range(len(ENPOT))]) # Maximum value in the ENMAX
if EN_recommend - 50 > ENSTART:
	ENSTART = int(math.ceil((EN_recommend-50)/ENSTEP))*ENSTEP

ENCUT = ENSTART
loopnum = 0
convergence = 1

while convergence == 1 :
	loopnum = loopnum+1
	if ENCUT > ENMAX :
		make_amp2_log(dir+'/cutoff','Too high cut-off energy is required!')
#		print ("ERROR: Too high cut-off energy is required!")
		print(0)
		sys.exit()

	if not os.path.isdir(dir+'/cutoff/EN'+str(ENCUT)) :
		os.mkdir(dir+'/cutoff/EN'+str(ENCUT), 0o755)

	make_amp2_log(dir+'/cutoff','EN'+str(ENCUT)+' calculation start')
	now_path = dir+'/cutoff/EN'+str(ENCUT)
	os.chdir(now_path)

	rerun = 0 
	if os.path.isfile(now_path+'/OUTCAR'):
		if 'Voluntary' in pytail(now_path+'/OUTCAR'):
			rerun = 1
			write_conv_result(now_path,enlog)
	if rerun == 0:
		copy_input(dir+'/INPUT0',now_path,POT)
		gam = set_parallel(now_path+'/KPOINTS',now_path+'/INCAR',npar,kpar)
		if gam == 1:
			vasprun = vasp_gam
		else:
			vasprun = vasp_std

		if pot_type == 'HSE':
			incar_for_hse(now_path+'/INCAR')
		incar_from_yaml(now_path,inp_conv['incar'])
		wincar(now_path+'/INCAR',now_path+'/INCAR',[['ENCUT',str(ENCUT)],['NSW','0'],['LWAVE','F'],['LCHARG','F']],[])
		# Running vasp
		out = run_vasp(now_path,nproc,vasprun,mpi)
		if ENCUT < EN_recommend:
			if out == 1:
				make_amp2_log(dir+'/cutoff','VASP error occurs, but pass because of lower encut than ENMAX in POTCAR')
				out = 3
				loopnum = 0
			else:
				out = electronic_step_convergence_check(now_path)

			while out == 1:
				make_amp2_log(dir+'/cutoff','Calculation options are changed. New calculation starts.')
				out = run_vasp(now_path,nproc,vasprun,mpi)
				if out == 1:  # error in vasp calculation
					make_amp2_log(dir+'/cutoff','VASP error occurs, but pass because of lower encut than ENMAX in POTCAR')
					out = 3
					loopnum = 0
				else:
					out = electronic_step_convergence_check(now_path)

			if out == 2:  # electronic step is not converged. (algo = normal)
				make_amp2_log(dir+'/cutoff','The calculation stops but electronic step is not converged, but pass because of lower encut than ENMAX in POTCAR')
				loopnum = 0
			elif out < 2:
				write_conv_result(now_path,enlog)

		else:
			if out == 1:  # error in vasp calculation
				print(0)
				sys.exit() 

			out = electronic_step_convergence_check(now_path)
			while out == 1:
				make_amp2_log(dir+'/cutoff','Calculation options are changed. New calculation starts.')
				out = run_vasp(now_path,nproc,vasprun,mpi)
				if out == 1:  # error in vasp calculation
					print(0)
					sys.exit()
				out = electronic_step_convergence_check(now_path)

			if out == 2:  # electronic step is not converged. (algo = normal)
				make_amp2_log(dir+'/cutoff','The calculation stops but electronic step is not converged.')
				print(0)
				sys.exit()
			write_conv_result(now_path,enlog)

	# electronic step is converged.
	if loopnum >= 3:
		convergence = convergence_check(now_path,dir+'/cutoff/EN'+str(ENCUT-ENSTEP),dir+'/cutoff/EN'+str(ENCUT-2*ENSTEP),ENCONV,PRCONV,FOCONV)
	ENCUT = ENCUT + ENSTEP
	os.chdir('../')

make_amp2_log(dir+'/cutoff','cut off energy test is done')

with open(enlog,'a') as result:
	result.write("\nConvergence criterion: E/atom < "+str(ENCONV)+" eV\n")
	if PRCONV > 0 :
		result.write("                       Pressure < "+str(PRCONV)+" kB\n")
	if FOCONV > 0 :
		result.write("                       Force < "+str(FOCONV)+" eV/Angst\n")
	result.write("Converged ENCUT: "+str(ENCUT-3*ENSTEP)+' eV\n')
	wincar(dir+'/INPUT0/INCAR',dir+'/INPUT0/INCAR',[['ENCUT',str(ENCUT-3*ENSTEP)]],[])

make_conv_dat(dir+'/cutoff','cutoff')
if inp_yaml['calculation']['plot'] == 1:
	os.chdir(dir+'/cutoff')
	subprocess.call([gnuplot,dir+'/cutoff/conv_plot.in'])

with open(dir+'/cutoff/amp2.log','r') as amp2_log:
	with open(dir+'/amp2.log','a') as amp2_log_tot:
		amp2_log_tot.write(amp2_log.read())

print(1)
