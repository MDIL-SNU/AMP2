###########################################
### Date: 2018-12-05			###
### yybbyb@snu.ac.kr			###
###########################################
import os
from module_vasprun import *
from module_vector import *

def write_diel_log(outcar_file,target):
	diel = subprocess.check_output(['grep','-A4','DIELECTRIC TENSOR',outcar_file]).splitlines()[-9:]
	diel_e = []
	diel_i = []
	for i in range(3) :
		diel_e.append(diel[i].split())
		diel_i.append(diel[6+i].split())
	diel_sum = [[],[],[]]
	for i in range(3) :
		diel_e[i] = [float(diel_e[i][0]), float(diel_e[i][1]), float(diel_e[i][2])]
		diel_i[i] = [float(diel_i[i][0]), float(diel_i[i][1]), float(diel_i[i][2])]
		diel_sum[i] = [diel_e[i][0]+diel_i[i][0],diel_e[i][1]+diel_i[i][1],diel_e[i][2]+diel_i[i][2]]
	diel0 = (diel_sum[0][0]+diel_sum[1][1]+diel_sum[2][2])/3.0
	with open(target+'/dielectric.log', 'w') as diel_log:
		diel_log.write('Dielectric tensor (electronic contribution):\n')
		for i in range(3) :
			diel_log.write('  '+str(round(diel_e[i][0],3))+'\t'+str(round(diel_e[i][1],3))+'\t'+str(round(diel_e[i][2],3))+'\n')
		diel_log.write('Dielectric tensor (ionic contribution):\n')
		for i in range(3) :
			diel_log.write('  '+str(round(diel_i[i][0],3))+'\t'+str(round(diel_i[i][1],3))+'\t'+str(round(diel_i[i][2],3))+'\n')
		diel_log.write('\nAveraged static dielectric constant: '+str(round(diel0,3))+'\n')

