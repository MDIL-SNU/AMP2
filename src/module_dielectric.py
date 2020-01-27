###########################################
### Date: 2018-12-05			###
### yybbyb@snu.ac.kr			###
###########################################
import os
from module_vasprun import *
from module_vector import *

def write_diel_log(outcar_file,target):
	import numpy as np
	diel_file_error = 0 # 0 is no error
	try:
		diel = pygrep('DIELECTRIC TENSOR (including',outcar_file,0,4).splitlines()[-3:]
		diel_e = []
		for i in range(3) :
			diel_e.append(diel[i].split())
	except:
		diel_e[[0,0,0],[0,0,0],[0,0,0]]
		diel_file_error = 1
	try:
		diel = pygrep('DIELECTRIC TENSOR IONIC',outcar_file,0,4).splitlines()[-3:]
		diel_i = []
		for i in range(3) :
			diel_i.append(diel[i].split())
	except:
		diel_i[[0,0,0],[0,0,0],[0,0,0]]
		diel_file_error = 1

	diel_sum = [[],[],[]]
	for i in range(3) :
		diel_e[i] = [float(diel_e[i][0]), float(diel_e[i][1]), float(diel_e[i][2])]
		diel_i[i] = [float(diel_i[i][0]), float(diel_i[i][1]), float(diel_i[i][2])]
		diel_sum[i] = [diel_e[i][0]+diel_i[i][0],diel_e[i][1]+diel_i[i][1],diel_e[i][2]+diel_i[i][2]]
	diel_e_array = np.real(np.asarray(diel_e))
	diel_i_array = np.real(np.asarray(diel_i))
	diel_sum_array = np.real(np.asarray(diel_sum))
	diel_e_dia = np.real(np.linalg.eigvals(diel_e))
	diel_i_dia = np.real(np.linalg.eigvals(diel_i))
	diel_sum_dia = np.real(np.linalg.eigvals(diel_sum))

	diel0 = sum(diel_sum_dia)/3.0
	mode_check = check_imaginary(outcar_file,target)
	with open(target+'/dielectric.log', 'w') as diel_log:
		if diel_file_error == 1:
			diel_log.write('One of dielectric tensor is missing. ')
		if mode_check == 2:
			diel_log.write('DO NOT BELIEVE IT!! Read amp2.log   ')
		elif mode_check == 1:
			diel_log.write('USE IT CAUTIOUSLY!! Read amp2.log   ')
		diel_log.write('Dielectric tensor (electronic contribution):\n')
		for i in range(3) :
			diel_log.write('  '.join(['{:10.3f}'.format(diel_e[i][x]) for x in range(3)])+'\n')
#			diel_log.write('  '+str(round(diel_e[i][0],3))+'\t'+str(round(diel_e[i][1],3))+'\t'+str(round(diel_e[i][2],3))+'\n')
		diel_log.write('Dielectric tensor (ionic contribution):\n')
		for i in range(3) :
			diel_log.write('  '.join(['{:10.3f}'.format(diel_i[i][x]) for x in range(3)])+'\n')
#			diel_log.write('  '+str(round(diel_i[i][0],3))+'\t'+str(round(diel_i[i][1],3))+'\t'+str(round(diel_i[i][2],3))+'\n')

		diel_log.write('\nDielectric constant diagonalization (electronic): '+' '.join(['{:10.3f}'.format(diel_e_dia[x]) for x in range(3)])+'\n')
		diel_log.write('Dielectric constant diagonalization (ionic): '+' '.join(['{:10.3f}'.format(diel_i_dia[x]) for x in range(3)])+'\n')
		diel_log.write('\nAveraged static dielectric constant: '+'{:10.3f}'.format(diel0)+'\n')

def check_imaginary(outcar_file,target):
	# return 0: no probelm, return 1: there are imaginary modes.
	lines = pygrep('f/i',outcar_file,0,0).splitlines()
	max_mode = 0
	if len(lines) > 3:
		make_amp2_log(target,'Warning!! Imaginary phonon mode is observed.')
		for line in lines:
			if float(line.split()[-2]) > max_mode:
				max_mode = float(line.split()[-2])
			if float(line.split()[-2]) > 0.1 :
				make_amp2_log(target,'Warning!! Maximum frequency of imaginary mode is over 0.1 meV.')
				return 2
	else:
		return 0
	make_amp2_log(target,'Maximum frequency of imaginary mode is '+str(max_mode)+'. The results should be used cautiously.')
	return 1

