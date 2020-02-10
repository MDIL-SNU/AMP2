###########################################
### Date: 2018-12-31			###
### yybbyb@snu.ac.kr			###
###########################################
# This is a package of modules for convergence test.
import subprocess,sys
from module_log import *
from module_vector import dist_point
from module_vasprun import pygrep,pyhead,pytail

# check convergence for kpoints and cutoff energy
def convergence_check(target,be1,be2,ENCONV,PRCONV,FOCONV):
	path_list = [target,be1,be2]
	nion = int(pygrep('NION',target+'/OUTCAR',0,0).splitlines()[-1].split()[11])
	# Read energy, pressure and force
	ENERGY = []
	PRESS = []
	FORCE = []
	for path in path_list:
		if ENCONV > 0:
			ENERGY.append(float(pygrep('free  ',path+'/OUTCAR',0,0).splitlines()[-1].split()[4]))
		if PRCONV > 0:
			line = pygrep('in kB',path+'/OUTCAR',0,0).splitlines()[-1].split()
			if len(line) == 8:
				PRESS.append([float(x) for x in line[2:]])
			else:
				PRESS.append([99999.,99999.,99999.,99999.,99999.,99999.])
		if FOCONV > 0:
			force_lines = pygrep('TOTAL-F',path+'/OUTCAR',0,nion+1).splitlines()[2:nion+2]
			FORCE_one_sample = []
			for line in force_lines:
				if len(line.split()) == 6:
					FORCE_one_sample.append([float(x) for x in line.split()[3:6]])
				else:
					FORCE_one_sample.append([999.,999.,999.])
			FORCE.append(FORCE_one_sample)
	converge = 0
	# Convergence check for energy/atom
	if ENCONV > 0 :
		if abs(ENERGY[1]-ENERGY[2])/nion > ENCONV or abs(ENERGY[0]-ENERGY[2])/nion > ENCONV or abs(ENERGY[0]-ENERGY[1])/nion > ENCONV:
			make_amp2_log(target,'Energy is not yet converged')
			return 1
			sys.exit()
	# Convergence check for pressure
	if PRCONV > 0 :
		for i in range(6) :
			if abs(PRESS[1][i]-PRESS[2][i]) > PRCONV or abs(PRESS[0][i]-PRESS[2][i]) > PRCONV or abs(PRESS[0][i]-PRESS[1][i]) > PRCONV:
#				if abs(PRESS[0][i]-PRESS[1][i])/abs(PRESS[0][i]) > 0.1 or abs(PRESS[0][i]-PRESS[2][i])/abs(PRESS[0][i]) > 0.1:
#					make_amp2_log(target,'Pressure is not yet converged')
#					return 1
#					sys.exit()
				make_amp2_log(target,'Pressure is not yet converged')
				return 1
				sys.exit()
		
	# Convergence check for force
	if FOCONV > 0 :
		for i in range(nion):
			for j in range(3) :
				if abs(FORCE[1][i][j]-FORCE[2][i][j]) > FOCONV or abs(FORCE[0][i][j]-FORCE[2][i][j]) > FOCONV or abs(FORCE[0][i][j]-FORCE[1][i][j]) > FOCONV:
					make_amp2_log(target,'Force is not yet converged')
					return 1
					sys.exit()
	return 0

# make a log file for convergence test
def write_conv_result(target,logfile):
	head = target.split('/')[-1]
	energy = float(pygrep('free  ',target+'/OUTCAR',0,0).splitlines()[-1].split()[4])
	line = pygrep('in kB',target+'/OUTCAR',0,0).splitlines()[-1].split()
	if len(line) == 8:
		press = [float(x) for x in line[2:]]
	else:
		press = [99999.,99999.,99999.,99999.,99999.,99999.]
	with open(logfile,'a') as log:
		log.write(head+': '+str(energy)+' eV')
		log.write('\t/ ')
		for ll in press:
			log.write(str(ll)+' ')
		log.write('kB\n')

# This is function for making data file for plotting convergence test
def make_conv_dat(target,typ):
	# typ is 'kpoint' or 'cutoff'
	conv = []
	with open(target+'/'+typ+'.log','r') as inp:
		for line in inp.readlines():
			tmp = line.split()
			if len(tmp) < 2:
				break
			else:
				conv.append(tmp)

	test_val = []
	with open(target+'/'+conv[0][0].split(':')[0]+'/POSCAR','r') as pos_f:
		lines = pos_f.readlines()
	natom = sum([int(x) for x in lines[6].split()])

	with open(target+'/conv_plot.dat','w') as out:
		for i in range(len(conv)):
			target_path = conv[i][0].split(':')[0]
			if typ =='kpoint':
				with open(target+'/'+target_path+'/KPOINTS','r') as kp_in:
					kpt = [float(x) for x in kp_in.readlines()[3].split()]
					test_val.append(dist_point(kpt,[0,0,0]))
			else:
				test_val.append(float(target_path[2:]))
			out.write(target_path+'\t'+str(test_val[i])+'\t')		# title
			out.write(str(float(conv[i][1])/float(natom))+'\t')				# energy
			max_press = max([float(x) for x in conv[i][4:7]])
			out.write(str(max_press)+'\t')				# pressure
			out.write(str(max_force(target+'/'+target_path+'/OUTCAR')[0])+'\n')	# force
	xr = [round(test_val[0]-(test_val[1]-test_val[0])/2,1),round(test_val[-1]+(test_val[1]-test_val[0])/2,1)]
	make_conv_plot(target,typ,xr)

# This is function for plotting convergence test
def make_conv_plot(target,typ,region):
	if typ == 'kpoint':
		xlabel = '(# of k-points)^{1/3}'
		xtic = '1'
	else:
		xlabel = 'cut off energy (eV)'
		xtic = '50'

	with open(target+'/conv_plot.in','w') as inp:
		inp.write("set term pngcairo size 720,720 font 'Arial, 20'\n")
		inp.write("set output 'Energy_conv.png'\n")
		inp.write("set xr["+str(region[0])+':'+str(region[1])+"]\n")
		inp.write("set tics nomirror out\n")
		inp.write("set xlabel '"+xlabel+"'\n")
		inp.write("set ylabel 'Energy (eV/atom)'\n")
		inp.write("set border 3 back\n")
		inp.write("set xtics "+xtic+" nomirror\n")
		inp.write("p 'conv_plot.dat' u 2:3 w linespoints pt 5 ps 2.5 lc rgb '#00CC00' title 'Energy'\n\n")

		inp.write("set term pngcairo size 720,720 font 'Arial, 20'\n")
		inp.write("set output 'Pressure_conv.png'\n")
		inp.write("set xr["+str(region[0])+':'+str(region[1])+"]\n")
		inp.write("set tics nomirror out\n")
		inp.write("set xlabel '"+xlabel+"'\n")
		inp.write("set ylabel 'Pressure (kB)'\n")
		inp.write("set border 3 back\n")
		inp.write("set xtics "+xtic+" nomirror\n")
		inp.write("p 'conv_plot.dat' u 2:4 w linespoints pt 7 ps 2.5 lc rgb 'red' title 'Pressure'\n\n")

		inp.write("set term pngcairo size 720,720 font 'Arial, 20'\n")
		inp.write("set output 'Force_conv.png'\n")
		inp.write("set xr["+str(region[0])+':'+str(region[1])+"]\n")
		inp.write("set tics nomirror out\n")
		inp.write("set xlabel '"+xlabel+"'\n")
		inp.write("set ylabel 'Force (eV/Angst)'\n")
		inp.write("set border 3 back\n")
		inp.write("set xtics "+xtic+" nomirror\n")
		inp.write("p 'conv_plot.dat' u 2:5 w linespoints pt 9 ps 2.5 lc rgb 'blue' title 'Force'\n")

# This function is finding maximum force in OUTCAR
def max_force(outcar_path):
	nion = int(pygrep('NION',outcar_path,0,0).split()[-1])
	force_list = pygrep('TOTAL-FORCE',outcar_path,0,nion+1).splitlines()
	max_f = [0,[0,0,0]]
	for line in force_list[2:]:
		if len(line.split()) == 6:
			force = [float(x) for x in line.split()[3:6]]
		else:
			force = [999.,999.,999.]
		if dist_point(force,[0,0,0]) > max_f[0]:
			max_f = [dist_point(force,[0,0,0]),force]
	return max_f
