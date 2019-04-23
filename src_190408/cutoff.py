####################################
# date : 2018-12-06                #
# Author : feihoon@snu.ac.kr       #
# Modifier : yybbyb@snu.ac.kr      #
####################################
import os, sys, subprocess, yaml
from module_log import *
from module_vasprun import *
from module_converge import *
code_data = 'Date: 2018-12-06'

# Set input
dir = sys.argv[1]

inp_file = sys.argv[2]
with open(inp_file,'r') as f:
	inp_yaml = yaml.load(f)
ERROR_path = inp_yaml['directory']['error']
src_path = inp_yaml['directory']['src_path']
vasp_std = inp_yaml['program']['vasp_std']
vasp_gam = inp_yaml['program']['vasp_gam']
mpi_type = inp_yaml['program']['mpi_type']
gnuplot = inp_yaml['program']['gnuplot']
npar = inp_yaml['vasp_parallel']['npar']
kpar = inp_yaml['vasp_parallel']['kpar']

inp_conv = inp_yaml['converge_test']
ENCONV = inp_conv['ENCONV']
PRCONV = inp_conv['PRCONV']
FOCONV = inp_conv['FOCONV']
ENSTART = inp_conv['ENSTART']
ENSTEP = inp_conv['ENSTEP']
ENMAX = inp_conv['ENMAX']

node = node_simple(sys.argv[3])
nproc = sys.argv[4]

# Check existing data
if os.path.isdir(dir+'/cutoff') and os.path.isfile(dir+'/cutoff/cutoff.log') :
	if len(subprocess.check_output(['grep','Converged',dir+'/cutoff/cutoff.log']).splitlines()) > 0 :
		make_amp2_log_default(dir,src_path,'cutoff energy test',node,code_data)
		make_amp2_log(dir,'Already done')
#		print('Success!')
		print 1
		sys.exit()
if not os.path.isdir(dir+'/cutoff') :
	os.mkdir(dir+'/cutoff', 0755)
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
ENPOT = subprocess.check_output(['grep','ENMIN',dir+'/INPUT0/POTCAR_GGA']).splitlines()
for i in range(len(ENPOT)) :
	if float(ENPOT[i].split()[5]) > ENSTART :
		ENSTART = int(float(ENPOT[i].split()[5])/ENSTEP)*ENSTEP

ENCUT = ENSTART
loopnum = 0
convergence = 1

while convergence == 1 :
	loopnum = loopnum+1
	if ENCUT > ENMAX :
		make_amp2_log(dir+'/cutoff','Too high cut-off energy is required!')
#		print ("ERROR: Too high cut-off energy is required!")
		print 0
		sys.exit()

	if not os.path.isdir(dir+'/cutoff/EN'+str(ENCUT)) :
		os.mkdir(dir+'/cutoff/EN'+str(ENCUT), 0755)

	make_amp2_log(dir+'/cutoff','EN'+str(ENCUT)+' calculation start')
	now_path = dir+'/cutoff/EN'+str(ENCUT)
	copy_input(dir+'/INPUT0',now_path,'GGA')

	os.chdir(now_path)

	rerun = 0 
	if os.path.isfile(now_path+'/OUTCAR'):
		if 'Voluntary' in subprocess.check_output(['tail','-n','1',now_path+'/OUTCAR']):
			rerun = 1
	if rerun == 0:
		gam = set_parallel(now_path+'/KPOINTS',now_path+'/INCAR',npar,kpar)
		if gam == 1:
			vasprun = vasp_gam
		else:
			vasprun = vasp_std

		wincar(now_path+'/INCAR',now_path+'/INCAR',[['ENCUT',str(ENCUT)],['NSW','0'],['LWAVE','F'],['LCHARG','F']],[])
		# Running vasp
		out = run_vasp(now_path,nproc,vasprun,mpi_type)
		if out == 1:  # error in vasp calculation
			print 0
			sys.exit() 

		out = electronic_step_convergence_check(now_path)
		if out == 2:  # electronic step is not converged. (algo = normal)
			make_amp2_log(dir+'/cutoff','Electronic step is not converged. ALGO is already normal.')
			print 0
			sys.exit()
		elif out == 1:  # elctronic step is not converged. (algo = fast) Algo changes to normal and rerun.
			make_amp2_log(dir+'/cutoff','Electronic step is not converged. ALGO changes to normal.')
			out = run_vasp(now_path,nproc,vasprun,mpi_type)
			if out == 1:  # error in vasp calculation
				print 0
				sys.exit()
			out = electronic_step_convergence_check(now_path)
			if out == 2:  # electronic step is not converged. (algo = normal)
				make_amp2_log(dir+'/cutoff','Electronic step is not converged. ALGO is already normal.')
				print 0
				sys.exit()
	# electronic step is converged.
	write_conv_result(now_path,enlog)
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

print 1
