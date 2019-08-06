####################################
# date : 2018-12-06                #
# Author : feihoon@snu.ac.kr       #
# Modifier : yybbyb@snu.ac.kr      #
####################################
import os, sys, subprocess, yaml
from module_log import *
from module_vasprun import *
from module_converge import *
code_data = 'Version xx. Modified at 2019-07-17'

# Set input
dir = sys.argv[1]

inp_file = sys.argv[2]
with open(inp_file,'r') as f:
	inp_yaml = yaml.load(f)
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
pot_type = inp_conv['potential_type']
if pot_type == 'LDA':
	POT = 'LDA'
else:
	POT = 'GGA'

node = node_simple(sys.argv[3])
nproc = sys.argv[4]

# Check existing data
if os.path.isdir(dir+'/kptest') and os.path.isfile(dir+'/kptest/KPOINTS_converged') :
	make_amp2_log_default(dir,src_path,'KPOINTS test',node,code_data)
#	print('Success!')
	make_amp2_log(dir,'Already done')
	print 1
	sys.exit()

if not os.path.isdir(dir+'/kptest') :
	os.mkdir(dir+'/kptest', 0755)
	make_amp2_log_default(dir+'/kptest',src_path,'KPOINTS test',node,code_data)
	make_amp2_log(dir+'/kptest','New calculation.')
else:
	make_amp2_log_default(dir+'/kptest',src_path,'KPOINTS test',node,code_data)
	make_amp2_log(dir+'/kptest','Calculation continue.')

os.chdir(dir+'/kptest')
kplog = dir+'/kptest/kpoint.log'
if os.path.isfile(kplog):
	os.remove(kplog)

KPL = inp_conv['initial_kpl']
loopnum = 0
convergence = 1

while convergence == 1 :
	loopnum = loopnum+1
	if KPL > inp_conv['max_kpl'] :
		make_amp2_log(dir+'/kptest','Too many k-point mesh is required!')
#		print("Too many k-point mesh is required!")
		print 0
		sys.exit()

	if not os.path.isdir(dir+'/kptest/KP'+str(KPL)) :
		os.mkdir(dir+'/kptest/KP'+str(KPL), 0755)

	make_amp2_log(dir+'/kptest','KP'+str(KPL)+' calculation start')
	now_path = dir+'/kptest/KP'+str(KPL)

	os.chdir(now_path)
	
	rerun = 0 
	if os.path.isfile(now_path+'/OUTCAR'):
		if 'Voluntary' in subprocess.check_output(['tail','-n','1',now_path+'/OUTCAR']):
			rerun = 1
	if rerun == 0:
		copy_input_no_kp(dir+'/INPUT0',now_path,POT)
		sym = int(open(dir+'/INPUT0/sym','r').readline())
		kpt_generation_for_relax(now_path,KPL,sym)
		gam = set_parallel(now_path+'/KPOINTS',now_path+'/INCAR',npar,kpar)
		if gam == 1:
			vasprun = vasp_gam
		else:
			vasprun = vasp_std

		if pot_type == 'HSE':
			incar_for_hse(now_path+'/INCAR')
		incar_from_yaml(now_path,inp_conv['incar'])
		wincar(now_path+'/INCAR',now_path+'/INCAR',[['NSW','0'],['LCHARG','F'],['LWAVE','F']],[])
		out = run_vasp(now_path,nproc,vasprun,mpi)
		if out == 1:  # error in vasp calculation
			print 0
			sys.exit() 

		out = electronic_step_convergence_check(now_path)
		while out == 1:
			make_amp2_log(dir+'/kptest','Calculation options are changed. New calculation starts.')
			out = run_vasp(now_path,nproc,vasprun,mpi)
			if out == 1:  # error in vasp calculation
				print 0
				sys.exit()
			out = electronic_step_convergence_check(now_path)

		if out == 2:  # electronic step is not converged. (algo = normal)
			make_amp2_log(dir+'/kptest','The calculation stops but electronic step is not converged.')
			print 0
			sys.exit()

	# electronic step is converged.
	write_conv_result(now_path,kplog)
	if loopnum >= 3:
		convergence = convergence_check(now_path,dir+'/kptest/KP'+str(KPL-1),dir+'/kptest/KP'+str(KPL-2),ENCONV,PRCONV,FOCONV)
	KPL = KPL + 1 
	os.chdir('../')

make_amp2_log(dir+'/kptest','kpoint test is done')

with open(kplog,'a') as result:
	result.write("\nConvergence criterion: E/atom < "+str(ENCONV)+" eV\n")
	if PRCONV > 0 :
		result.write("                       Pressure < "+str(PRCONV)+" kB\n")
	if FOCONV > 0 :
		result.write("                       Force < "+str(FOCONV)+" eV/Angst\n")
	result.write("Converged KPL: "+str(KPL-3)+'\n')
	subprocess.call(['cp',dir+'/kptest/KP'+str(KPL-3)+'/KPOINTS',dir+'/kptest/KPOINTS_converged'])
	subprocess.call(['cp',dir+'/kptest/KPOINTS_converged',dir+'/INPUT0/KPOINTS'])

make_conv_dat(dir+'/kptest','kpoint')
if inp_yaml['calculation']['plot'] == 1:
	os.chdir(dir+'/kptest')
	subprocess.call([gnuplot,dir+'/kptest/conv_plot.in'])

with open(dir+'/kptest/amp2.log','r') as amp2_log:
	with open(dir+'/amp2.log','a') as amp2_log_tot:
		amp2_log_tot.write(amp2_log.read())

print 1
