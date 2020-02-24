####################################
# Modifier : yybbyb@snu.ac.kr      #
# data : 2018-12-05                #
####################################
# This is a package of modules for drawing band structure and calculating band gap.
from module_vector import *
from module_vasprun import poscar_to_axis,pygrep,pyhead,pytail
from module_log import *
import os,sys,math,subprocess

def make_sym_for_band(poscar_file,sym_file,target):
	axis = poscar_to_axis(poscar_file)
	with open(sym_file,'r') as sym:
		table = int(sym.readline().split()[0])

### set axis length and angle of conventional lattice
	if table in [5,6,8,9,10]: 
		a = abs(axis[1][0])*2.
		b = abs(axis[0][1])*2.
		c = abs(axis[0][2])*2.
		alpha = math.radians(90)
	elif table == 11 :
		a = abs(axis[0][0])*2.
		b = abs(axis[1][1])*2.
		c = axis[2][2]
		alpha = math.radians(90)
	elif table in [16,17,18]:
		a = abs(axis[0][0])*2.
		b = abs(axis[1][1])*2.
		c = 0
		for i in range(3) :
			c = c + axis[2][i]**2.0
		c = c**(0.5)
		alpha = math.asin(axis[2][2]/c)
	else :
		a = 0
		for i in range(3) :
			a = a + axis[0][i]**2.0
		a = a**(0.5)
		b = 0
		for i in range(3) :
		        b = b + axis[1][i]**2.0
		b = b**(0.5)
		c = 0
		for i in range(3) :
		        c = c + axis[2][i]**2.0
		c = c**(0.5)
		alpha = math.acos((axis[1][0]*axis[2][0]+axis[1][1]*axis[2][1]+axis[1][2]*axis[2][2])/(b*c))
###############################

	[b1,b2,b3] = reciprocal_lattice(axis)

#### define constant for kpoint and subdevide the symmetry number
	kvar = []
	if table == 5 :
		kvar.append((1.+c**2./a**2.)/4.)
	elif table == 6 :
		kvar.append((1.+a**2./c**2.)/4.)
		kvar.append(a**2./(2.*c**2.))
	elif table == 8 :
		kvar.append((1.+a**2./b**2.-a**2./c**2.)/4.)
		kvar.append((1.+a**2./b**2.+a**2./c**2.)/4.)
	elif table == 9 :
		kvar.append((1.+a**2./b**2.-a**2./c**2.)/4.)
		kvar.append((1.+c**2./b**2.-c**2./a**2.)/4.)
		kvar.append((1.+b**2./a**2.-b**2./c**2.)/4.)
	elif table == 10 :
		kvar.append((1.+a**2./c**2.)/4.)
		kvar.append((1.+b**2./c**2.)/4.)
		kvar.append((b**2.-c**2.)/(4.*c**2.))
		kvar.append((a**2.+b**2.)/(4.*c**2.))
	elif table == 11 :
		kvar.append((1.+a**2./b**2.)/4.)
	elif table == 13:
		kvar.append((1.+4.*math.cos(alpha))/(2.+4.*math.cos(alpha)))
		kvar.append(3./4.0-kvar[0]/2.)
	elif table == 14 :
		kvar.append(1./(2.*(math.tan(alpha/2.))**2.))
		kvar.append(3./4.0-kvar[0]/2.)
	elif table == 15 :
		kvar.append((1.-b*math.cos(alpha)/c)/(2.*(math.sin(alpha))**2.))
		kvar.append(1./2.0-kvar[0]*c*math.cos(alpha)/b)
	elif table == 16 :
		kvar.append((2-b*math.cos(alpha)/c)/(4*(math.sin(alpha))**2))
		kvar.append(1/2.0+2*kvar[0]*c*math.cos(alpha)/b)
		kvar.append(3/4.0-a**2/(4*b**2*(math.sin(alpha))**2))
		kvar.append(kvar[2]+(3/4.0-kvar[2])*b*math.cos(alpha)/c)
	elif table == 17 :
		kvar.append((1.+b**2./a**2.)/4.)
		kvar.append(b*c*math.cos(alpha)/(2.*a**2.))
		kvar.append(kvar[0]-1./4.0+(1.-b*math.cos(alpha)/c)/(4.*(math.sin(alpha))**2.))
		kvar.append(1./2.0+2.*kvar[2]*c*math.cos(alpha)/b)
		kvar.append(1.+kvar[2]-2.*kvar[0])
		kvar.append(kvar[3]-2.*kvar[1])
	elif table == 18:
		kvar.append((b**2./a**2.+(1.-b*math.cos(alpha)/c)/(math.sin(alpha))**2.)/4.)
		kvar.append((1./2.0+2.*kvar[0]*c*math.cos(alpha)/b)/2.+b**2./(4.*a**2.)-b*c*math.cos(alpha)/(2.*a**2.))
		kvar.append((4.*(2.*kvar[1]-kvar[0])-1.-b**2.*(math.sin(alpha))**2./a**2.)*c/(2.*b*math.cos(alpha)))
		kvar.append(1./2.0+2.*kvar[0]*c*math.cos(alpha)/b)
		kvar.append(kvar[0]*c*math.cos(alpha)/b+kvar[2]/2.-1./4.0)
		kvar.append(2.*kvar[1]-kvar[0])
		kvar.append(1.-kvar[0]*a**2./b**2.)

##############################################################

	with open(target+'/sym','w') as out:
		out.write(str(table)+'\n')  # symmetry number
		if len(kvar) >= 1 : # constant for making kpoints for band
			for j in range(len(kvar)) :
				out.write(' '+str(kvar[j]))
#		reciprocal lattice
		out.write('\n'+str(round(b1[0],8))+' '+str(round(b1[1],8))+' '+str(round(b1[2],8)))
		out.write('\n'+str(round(b2[0],8))+' '+str(round(b2[1],8))+' '+str(round(b2[2],8)))
		out.write('\n'+str(round(b3[0],8))+' '+str(round(b3[1],8))+' '+str(round(b3[2],8)))

# This function is writing kpoints path corresponding symmetry
def make_symk(sym_file):
	sym = open(sym_file,'r').readlines()
	table = sym[0].split()[0]
	const = [float(x) for x in sym[1].split()]
	symk = []
	if table == '1' :	#CUB
		symk.append([0.,0.,0.]) #gamma
		symk.append([1./2,1./2,0.]) #M
		symk.append([1./2,1./2,1./2]) #R
		symk.append([0.,1./2,0.]) #X
		order=[[3,1],[5,-2],[1,-3],[2,4],[6,5],[4,6]]
		xticlabel = ['{/Symbol G}', 'X', 'M', '{/Symbol G}', 'R', 'X|M', 'R']
	elif table == '2' :	#FCC
		symk.append([0.,0.,0.]) #gamma
		symk.append([3./8,3./8,3./4]) #K
		symk.append([1./2,1./2,1./2]) #L
		symk.append([5./8,1./4,5./8]) #U
		symk.append([1./2,1./4,3./4]) #W
		symk.append([1./2,0.,1./2]) #X
		order=[[5,1],[15,-2],[8,-3],[1,-4],[2,5],[10,6],[13,7],[11,-8],[6,-9],[14,10]]
		xticlabel=['{/Symbol G}', 'X', 'W', 'K', '{/Symbol G}', 'L', 'U', 'W', 'L', 'K|U', 'X']
	elif table == '3' :	#BCC
		symk.append([0.,0.,0.]) #gamma
		symk.append([1./2,-1./2,1./2]) #H
		symk.append([1./4,1./4,1./4]) #P
		symk.append([0.,0.,1./2]) #N
		order=[[1,1],[5,2],[3,-3],[2,4],[4,-5],[6,6]]
		xticlabel=['{/Symbol G}', 'H', 'N', '{/Symbol G}', 'P', 'H|P', 'N']
	elif table == '4' :	#TET
		symk.append([0.,0.,0.]) #gamma
		symk.append([1./2,1./2,1./2]) #A
		symk.append([1./2,1./2,0.]) #M
		symk.append([0.,1./2,1./2]) #R
		symk.append([0.,1./2,0.]) #X
		symk.append([0.,0.,1./2]) #Z
		order=[[4,1],[11,-2],[2,-3],[5,4],[14,-5],[7,-6],[9,7],[13,-8],[6,-9]]
		xticlabel=['{/Symbol G}', 'X', 'M', '{/Symbol G}', 'Z', 'R', 'A', 'Z|X', 'R|M', 'A']
	elif table == '5' :	#BCT1
		symk.append([0.,0.,0.]) #gamma
		symk.append([-1./2,1./2,1./2]) #M
		symk.append([0.,1./2,0.]) #N
		symk.append([1./4,1./4,1./4]) #P
		symk.append([0.,0.,1./2]) #X
		symk.append([const[0],const[0],-const[0]]) #Z
		symk.append([-const[0],1.-const[0],const[0]]) #Z_1
		order=[[4,1],[9,-2],[1,-3],[5,4],[17,-5],[12,-6],[15,7],[11,-8],[16,-9]]
		xticlabel=['{/Symbol G}', 'X', 'M', '{/Symbol G}', 'Z', 'P', 'N', 'Z_1', 'M|X', 'P']
	elif table == '6' :	#BCT2
		symk.append([0.,0.,0.]) #gamma
		symk.append([0.,1./2,0.]) #N
		symk.append([1./4,1./4,1./4]) #P
		symk.append([-const[0],const[0],const[0]]) #sigma
		symk.append([const[0],1.-const[0],-const[0]]) #sigma_1
		symk.append([0.,0.,1./2]) #X
		symk.append([-const[1],const[1],1./2]) #Y
		symk.append([1./2,1./2,-const[1]]) #Y1
		symk.append([1./2,1./2,-1./2]) #Z
		order=[[5,1],[31,2],[24,-3],[3,-4],[8,5],[30,-6],[11,-7],[9,8],[20,9],[36,10],[18,-11]]
		xticlabel=['{/Symbol G}', 'X', 'Y', '{/Symbol S}', '{/Symbol G}', 'Z', '{/Symbol S}_1', 'N', 'P', 'Y_1', 'Z|X', 'P']
	elif table == '7' :	#ORC
		symk.append([0.,0.,0.]) #gamma
		symk.append([1./2,1./2,1./2]) #R
		symk.append([1./2,1./2,0.]) #S
		symk.append([0.,1./2,1./2]) #T
		symk.append([1./2,0.,1./2]) #U
		symk.append([1./2,0.,0.]) #X
		symk.append([0.,1./2,0.]) #Y
		symk.append([0.,0.,1./2]) #Z
		order=[[5,1],[16,-2],[17,3],[6,-4],[7,5],[25,-6],[10,-7],[9,8],[22,9],[21,-10],[23,11],[8,-12]]
		xticlabel=['{/Symbol G}', 'X', 'S', 'Y', '{/Symbol G}', 'Z', 'U', 'R', 'T', 'Z|Y', 'T|U', 'X|S', 'R']
	elif table == '8' :	#ORCF1
		symk.append([0.,0.,0.]) #gamma
		symk.append([1./2,1./2+const[0],const[0]]) #A
		symk.append([1./2,1./2-const[0],1.-const[0]]) #A_1
		symk.append([1./2,1./2,1./2]) #L
		symk.append([1.,1./2,1./2]) #T
		symk.append([0.,const[1],const[1]]) #X
		symk.append([1.,1.-const[1],1.-const[1]]) #X_1
		symk.append([1./2,0.,1./2]) #Y
		symk.append([1./2,1./2,0.]) #Z
		order=[[7,1],[29,-2],[30,3],[8,-4],[5,5],[18,-6],[20,7],[28,8],[12,-9],[15,10],[3,-11]]
		xticlabel=['{/Symbol G}', 'Y', 'T', 'Z', '{/Symbol G}', 'X', 'A_1', 'Y|T', 'X_1|X', 'A', 'Z|L', '{/Symbol G}']
	elif table == '9' :	#ORCF2
		symk.append([0.,0.,0.]) #gamma
		symk.append([1./2,1./2-const[0],1.-const[0]]) #C
		symk.append([1./2,1./2+const[0],const[0]]) #C_1
		symk.append([1./2-const[2],1./2,1.-const[2]]) #D
		symk.append([1./2+const[2],1./2,const[2]]) #D_1
		symk.append([1./2,1./2,1./2]) #L
		symk.append([1.-const[1],1./2-const[1],1./2]) #H
		symk.append([const[1],1./2+const[1],1./2]) #H_1
		symk.append([0.,1./2,1./2]) #X
		symk.append([1./2,0.,1./2]) #Y
		symk.append([1./2,1./2,0.]) #Z
		order=[[9,1],[18,-2],[12,3],[32,4],[8,-5],[10,6],[40,-7],[36,8],[15,-9],[27,10],[50,-11],[48,12],[5,-13]]
		xticlabel=['{/Symbol G}', 'Y', 'C', 'D', 'X', '{/Symbol G}', 'Z', 'D_!', 'H', 'C|C_1', 'Z', 'X', 'H_1|H', 'Y|L', '{/Symbol G}']
	elif table == '10' :	#ORCI
		symk.append([0.,0.,0.]) #gamma
		symk.append([-const[3],const[3],1./2-const[2]]) #L
		symk.append([const[3],-const[3],1./2+const[2]]) #L_1
		symk.append([1./2-const[2],1./2+const[2],-const[3]]) #L_2
		symk.append([0.,1./2,0.]) #R
		symk.append([1./2,0.,0.]) #S
		symk.append([0.,0.,1./2]) #T
		symk.append([1./4,1./4,1./4]) #W
		symk.append([-const[0],const[0],const[0]]) #X
		symk.append([const[0],1.-const[0],-const[0]]) #X_1
		symk.append([const[1],-const[1],const[1]]) #Y
		symk.append([1.-const[1],const[1],-const[1]]) #Y_1
		symk.append([1./2,1./2,-1./2]) #Z
		order=[[8,1],[19,-2],[17,3],[58,4],[45,-5],[47,6],[75,7],[12,-8],[10,9],[55,-10],[52,11],[31,12],[78,13]]
		xticlabel=['{/Symbol G}', 'X', 'L', 'T', 'W', 'R', 'X_1', 'Z', '{/Symbol G}', 'Y_1', 'S', 'W|L_1', 'Y|Y_1', 'Z']
	elif table == '11' :	#ORCC
		symk.append([0.,0.,0.]) #gamma
		symk.append([const[0],const[0],1./2]) #A
		symk.append([-const[0],1.-const[0],1./2]) #A_1
		symk.append([0.,1./2,1./2]) #R
		symk.append([0.,1./2,0.]) #S
		symk.append([-1./2,1./2,1./2]) #T
		symk.append([const[0],const[0],0.]) #X
		symk.append([-const[0],1.-const[0],0.]) #X_1
		symk.append([-1./2,1./2,0.]) #Y
		symk.append([0.,0.,1./2]) #Z
		order=[[6,1],[32,-2],[25,-3],[11,-4],[17,5],[9,-6],[8,7],[43,-8],[22,-9],[20,10],[38,11],[39,-12]]
		xticlabel=['{/Symbol G}', 'X', 'S', 'R', 'A', 'Z', '{/Symbol G}', 'Y', 'X_1', 'A_1', 'T', 'Y|Z', 'T']
	elif table == '12' :	#HEX
		symk.append([0.,0.,0.]) #gamma
		symk.append([0.,0.,1./2]) #A
		symk.append([1./3,1./3,1./2]) #H
		symk.append([1./3,1./3,0.]) #K
		symk.append([1./2,0.,1./2]) #L
		symk.append([1./2,0.,0.]) #M
		order=[[5,1],[14,-2],[3,-3],[1,4],[8,5],[11,-6],[6,-7],[15,8],[10,-9]]
		xticlabel=['{/Symbol G}', 'M', 'K', '{/Symbol G}', 'A', 'L', 'H', 'A|L', 'M|K', 'H']
	elif table == '13' :	#RHL1
		symk.append([0.,0.,0.]) #gamma
		symk.append([const[0],1./2,1.-const[0]]) #B
		symk.append([1./2,1.-const[0],const[0]-1.]) #B_1
		symk.append([1./2,1./2,0.]) #F
		symk.append([1./2,0.,0.]) #L
		symk.append([0.,0.,-1./2]) #L_1
		symk.append([const[0],const[1],const[1]]) #P
		symk.append([1.-const[1],1.-const[1],1.-const[0]]) #P_1
		symk.append([const[1],const[1],const[0]-1.]) #P_2
		symk.append([1.-const[1],const[1],0.]) #Q
		symk.append([const[1],0.,-const[1]]) #X
		symk.append([1./2,1./2,1./2]) #Z
		order=[[4,1],[23,-2],[21,3],[11,-4],[10,5],[36,-6],[34,7],[60,8],[40,9]]
		xticlabel=['{/Symbol G}', 'L', 'B_1|B', 'Z', '{/Symbol G}', 'X|Q', 'F', 'P_1', 'Z|L', 'P']
	elif table == '14' :	#RHL2
		symk.append([0.,0.,0.]) #gamma
		symk.append([1./2,-1./2,0.]) #F
		symk.append([1./2,0.,0.]) #L
		symk.append([1.-const[1],-const[1],1.-const[1]]) #P
		symk.append([const[1],const[1]-1.,const[1]-1.]) #P_1
		symk.append([const[0],const[0],const[0]]) #Q
		symk.append([1.-const[0],-const[0],-const[0]]) #Q_1
		symk.append([1./2,-1./2,1./2]) #Z
		order=[[3,1],[22,2],[27,-3],[5,-4],[1,5],[10,6],[24,7],[17,-8],[18,9]]
		xticlabel=['{/Symbol G}', 'P', 'Z', 'Q', 'R', 'F', 'P_1', 'Q_1', 'L', 'Z']
	elif table == '15' :	#MCL
		symk.append([0.,0.,0.]) #gamma
		symk.append([1./2,1./2,0.]) #A
		symk.append([0.,1./2,1./2]) #C
		symk.append([1./2,0.,1./2]) #D
		symk.append([1./2,0.,-1./2]) #D_1
		symk.append([1./2,1./2,1./2]) #E
		symk.append([0.,const[0],1.-const[1]]) #H
		symk.append([0.,1.-const[0],const[1]]) #H_1
		symk.append([0.,const[0],-const[1]]) #H_2
		symk.append([1./2,const[0],1.-const[1]]) #M
		symk.append([1./2,1.-const[0],const[1]]) #M_1
		symk.append([1./2,const[0],-const[1]]) #M_2
		symk.append([0.,1./2,0.]) #X
		symk.append([0.,0.,1./2]) #Y
		symk.append([0.,0.,-1./2]) #Y_1
		symk.append([1./2,0.,0.]) #Z
		order=[[13,1],[82,-2],[33,-3],[32,4],[70,5],[24,-6],[26,7],[89,-8],[48,-9],[54,10],[52,-11]]
		xticlabel=['{/Symbol G}', 'Y', 'H', 'C', 'E', 'M_1', 'A', 'X', 'H_1|M', 'D', 'Z|Y', 'D']
	elif table == '16' :	#MCLC1
		symk.append([0.,0.,0.]) #gamma
		symk.append([1./2,0.,0.]) #N
		symk.append([0.,-1./2,0.]) #N_1
		symk.append([1.-const[0],1.-const[0],1.-const[1]]) #F
		symk.append([const[0],const[0],const[1]]) #F_1
		symk.append([-const[0],-const[0],1.-const[1]]) #F_2
		symk.append([1.-const[0],-const[0],1.-const[1]]) #F_3
		symk.append([const[3],1.-const[3],1./2]) #I
		symk.append([1.-const[3],const[3]-1.,1./2]) #I_1
		symk.append([1./2,1./2,1./2]) #L
		symk.append([1./2,0.,1./2]) #M
		symk.append([1.-const[2],const[2]-1.,0.]) #X
		symk.append([const[2],1.-const[2],0.]) #X_1
		symk.append([const[2]-1.,-const[2],0.]) #X_2
		symk.append([1./2,1./2,0.]) #Y
		symk.append([-1./2,-1./2,0.]) #Y_1
		symk.append([0.,0.,1./2]) #Z
		order=[[14,1],[56,-2],[51,3],[93,-4],[108,5],[70,-6],[128,-7],[11,-8],[1,9],[10,-10]]
		xticlabel=['{/Symbol G}', 'Y', 'F', 'L', 'I|I_1', 'Z', 'F_1|Y', 'X_1|X', '{/Symbol G}', 'N|M', '{/Symbol G}']
	elif table == '17' :	#MCLC3
		symk.append([0.,0.,0.]) #gamma
		symk.append([1.-const[4],1.-const[4],1.-const[5]]) #F
		symk.append([const[4],const[4]-1.,const[5]]) #F_1
		symk.append([1.-const[4],-const[4],1.-const[5]]) #F_2
		symk.append([const[2],const[2],const[3]]) #H
		symk.append([1.-const[2],-const[2],1.-const[3]]) #H_1
		symk.append([-const[2],-const[2],1.-const[3]]) #H_2
		symk.append([1./2,-1./2,1./2]) #I
		symk.append([1./2,0.,1./2]) #M
		symk.append([1./2,0.,0.]) #N
		symk.append([0.,-1./2,0.]) #N_1
		symk.append([1./2,-1./2,0.]) #X
		symk.append([const[0],const[0],const[1]]) #Y
		symk.append([1.-const[0],-const[0],-const[1]]) #Y_1
		symk.append([-const[0],-const[0],-const[1]]) #Y_2
		symk.append([const[0],const[0]-1.,const[1]]) #Y_3
		symk.append([0.,0.,1./2]) #Z
		order=[[12,1],[27,-2],[19,3],[70,4],[100,-5],[36,-6],[78,7],[123,-8],[11,-9],[9,10],[8,-11]]
		xticlabel=['{/Symbol G}', 'Y', 'F', 'H', 'Z', 'I', 'F_1|H_1', 'Y_1', 'X', 'Ga', 'N}M', '{/Symbol G}']
	elif table == '18' :	#MCLC5
		symk.append([0.,0.,0.]) #gamma
		symk.append([const[5],const[5],const[2]]) #F
		symk.append([1.-const[5],1.-const[5],1.-const[2]]) #F_1
		symk.append([const[5],const[5]-1.,const[2]]) #F_2
		symk.append([const[0],const[0],const[3]]) #H
		symk.append([1.-const[0],-const[0],1.-const[3]]) #H_1
		symk.append([-const[0],-const[0],1.-const[3]]) #H_2
		symk.append([const[6],1.-const[6],1./2]) #I
		symk.append([1.-const[6],const[6]-1.,1./2]) #I_1
		symk.append([1./2,1./2,1./2]) #L
		symk.append([1./2,0.,1./2]) #M
		symk.append([1./2,0.,0.]) #N
		symk.append([0.,-1./2,0.]) #N_1
		symk.append([1./2,-1./2,0.]) #X
		symk.append([const[1],const[1],const[4]]) #Y
		symk.append([1.-const[1],-const[1],-const[4]]) #Y_1
		symk.append([-const[1],-const[1],-const[4]]) #Y_2
		symk.append([const[1],const[1]-1.,const[4]]) #Y_3
		symk.append([0.,0.,1./2]) #Z
		order=[[14,1],[31,-2],[26,3],[107,-4],[126,5],[80,-6],[37,-7],[90,8],[158,-9],[13,-10],[11,11],[10,-12]]
		xticlabel=['{/Symbol G}', 'Y', 'F', 'L', 'I|I_1', 'Z', 'H', 'F_1|H_1', 'Y_1', 'X', '{/Symbol G}', 'N|M', '{/Symbol G}']
	elif table == '19' :	#TRIa
		symk.append([0.,0.,0.]) #gamma
		symk.append([1./2,1./2,0.]) #L
		symk.append([0.,1./2,1./2]) #M
		symk.append([1./2,0.,1./2]) #N
		symk.append([1./2,1./2,1./2]) #R
		symk.append([1./2,0.,0.]) #X
		symk.append([0.,1./2,0.]) #Y
		symk.append([0.,0.,1./2]) #Z
		order=[[5,-1],[6,2],[1,-3],[7,4],[3,-5],[2,6],[4,-7]]
		xticlabel=['X', '{/Symbol G}', 'Y|L', '{/Symbol G}', 'Z|N', '{/Symbol G}', 'M|R', '{/Symbol G}']
	elif table == '20' :	#TRIb
		symk.append([0.,0.,0.]) #gamma
		symk.append([1./2,-1./2,0.]) #L
		symk.append([0.,0.,1./2]) #M
		symk.append([-1./2,-1./2,1./2]) #N
		symk.append([0.,-1./2,1./2]) #R
		symk.append([0.,-1./2,0.]) #X
		symk.append([1./2,0.,0.]) #Y
		symk.append([-1./2,0.,1./2]) #Z
		order=[[5,-1],[6,2],[1,-3],[7,4],[3,-5],[2,6],[4,-7]]
		xticlabel=['X', '{/Symbol G}', 'Y|L', '{/Symbol G}', 'Z|N', '{/Symbol G}', 'M|R', '{/Symbol G}']
	elif table == '21' :
		symk.append([0.,0.,0.]) #gamma
		symk.append([0.,1./2,0.]) #F
		symk.append([0.,1./2,1./2]) #Q
		symk.append([0.,0.,1./2]) #Z
		order=[[1,1],[4,2],[6,3],[3,-4]]
		xticlabel=['{/Symbol G}', 'F', 'Q', 'Z', '{/Symbol G}']

	# set the reciprocal axis
	if len(sym[2].split()) < 3 :
		idx = 3
	else : 
		idx = 2
	rec1=[float(x) for x in sym[idx].split()]
	rec2=[float(x) for x in sym[idx+1].split()]
	rec3=[float(x) for x in sym[idx+2].split()]
	rec=[rec1,rec2,rec3]

	return [symk,order,xticlabel,rec]


def make_kp_for_band(symk,order,xticlabel,rec,KSPACING,opt,target):
	# make path from the combination of high symmetry kpts
	path = []
	for i in range(len(symk)):
		for j in range(i+1,len(symk)):
			path.append([999,i+1,j+1]) # 999 for sort [[index,kp_start,kp_end],...]
	# change path index for band structure
	for i in range(len(order)):
		path[order[i][0]-1][0] = abs(order[i][1])
		if order[i][1] < 0:
			path[order[i][0]-1][1],path[order[i][0]-1][2]=path[order[i][0]-1][2],path[order[i][0]-1][1]
	# sorting the path
	path=sorted(path)
	# set the length of kp path. if band, includes kpts for band figure.
	if opt == 'band' :
		path = path[:len(order)]

	nkpt = []
	KPT = []
	section = []
	xtic = []
	xlabel = []
	for i in range(len(path)):
		kp_st = path[i][1]-1
		kp_end = path[i][2]-1
		dist_rec = [symk[kp_end][x]-symk[kp_st][x] for x in range(3)] # distance vector
		dist_val = dist_vec(symk[kp_end],symk[kp_st],rec) # distance value
		nsplit = round(dist_val/KSPACING,0)
		if nsplit == 0: # At the minimum, one kpt is included.
			nsplit = 1
		begin = 0
		if not i == 0:
			if kp_st == path[i-1][2]-1:
				begin = 1

		nkpt.append(nsplit+1-begin)
		for j in range(begin,int(nsplit)) :
			if j == 0 and i == 0:
				xtic.append(0)
				xlabel.append(0)
			elif j == 0 :
				xtic.append(xtic[-1])
			else :
				xtic.append(xtic[-1]+dist_val/nsplit)
			KPT.append([symk[kp_st][x]+dist_rec[x]/nsplit*j for x in range(3)])

		xtic.append(xtic[-1]+dist_val/nsplit)
		xlabel.append(xtic[-1])
		KPT.append([symk[kp_end][x] for x in range(3)])

	N_plot = int(sum(nkpt[:len(order)]))

	with open('xtic.dat','w') as out_xtic:
		for i in range(N_plot):
			out_xtic.write(str(xtic[i])+'\t'+'\t'.join([str(x) for x in KPT[i]])+'\n')

	with open('xlabel.dat','w') as out_xlabel:
		out_xlabel.write("set xtics('")
		for i in range(len(order)):
			out_xlabel.write(xticlabel[i]+"' "+str(xlabel[i])+", '")
		out_xlabel.write(xticlabel[len(order)]+"' "+str(xlabel[len(order)])+')\n')
		for i in range(1,len(order)):
			out_xlabel.write("set arrow from "+str(xlabel[i])+",graph(0,0) to "+str(xlabel[i])+",graph(1,1) nohead lt 0 lc 0\n")

	with open('KPOINTS_band', 'w') as out:
		out.write('k-points for band\n')
		out.write('  '+str(int(sum(nkpt)))+'\nReciprocal\n')
		for i in range(len(KPT)) :
			out.write('  '+'  '.join(["%0.8f" % (round(x,8)) for x in KPT[i]])+'   1\n')

	return N_plot

def EIGEN_to_array(eigen_file,spin):
	with open(eigen_file,'r') as eig:
		for i in range(5) :
			eig.readline()
		[nelect, nkpt, nband] = eig.readline().split()[0:3]
		nelect = int(nelect)
		nkpt = int(nkpt)
		nband = int(nband)
		KPT = [[] for i in range(nkpt)]
		Band = [[] for i in range(nband)]
		for k in range(nkpt) :
			eig.readline()
			KPT[k] = eig.readline().split()[0:3]
			for n in range(nband) :
				if spin == '2' :	# Read eigenvalues for spin-polarized calculation
					tmp = eig.readline().split()[1:3]
					Band[n].append([float(tmp[0]),float(tmp[1])])
				else :
					Band[n].append([float(eig.readline().split()[1])])
	return [KPT,Band,nelect] # KPT[kpts_idx,axis_idx],Band[band_idx,kpts_idx,spin_idx]

def get_fermi_level(Band,nelect,ncl):
	import numpy as np
	num_kpt = len(Band[0])
	spin = len(Band[0][0])
	Band_reshape = np.array(Band).reshape([-1])
	Band_reshape = np.sort(Band_reshape)
	if ncl == 'T':
		occupied_index = nelect*num_kpt
		fermi =  (Band_reshape[occupied_index-1] + Band_reshape[occupied_index])/2.0
	elif nelect%2 == 1 and spin == 1 and num_kpt%2 == 1:
		occupied_index = (nelect*num_kpt+1)//2
		fermi = Band_reshape[occupied_index-1]
	else:
		occupied_index = nelect*num_kpt//(3-spin)
		fermi =  (Band_reshape[occupied_index-1] + Band_reshape[occupied_index])/2.0
	return fermi

def calc_gap(fermi,spin,ncl,KPT,Band,nelect):
	VBM = []
	CBM = []
	nVBM = []
	nCBM = []
	VBM_k = []
	CBM_k = []
	for i in range(int(spin)):
		VBM.append(-10.0**6.0)
		CBM.append(10.0**6.0)
	metal = 0
	# VBM & CBM search for spin-polarized calculation
	nVB = -1
	for i in range(int(spin)) :
		for n in range(len(Band)) :
			if n == 0 :
				nVBM.append(n)
				nCBM.append(n)
				VBM_k.append(KPT[n][0][i])
				CBM_k.append(KPT[n][0][i])

			single_band = [Band[n][x][i] for x in range(len(KPT))]
			band_max = max(single_band)
			band_min = min(single_band)

			if band_max < fermi :
				if band_max >= VBM[i] :
					VBM[i] = band_max
					VBM_k[i] = KPT[single_band.index(band_max)]
					nVBM[i] = n+1
			elif band_min > fermi :
				if band_min < CBM[i] :
					CBM[i] = band_min
					CBM_k[i] = KPT[single_band.index(band_min)]
					nCBM[i] = n+1
			else :
				metal = 1
				VBM[i] = band_max
				VBM_k[i] = KPT[single_band.index(band_max)]
				nVBM[i] = n+1
				CBM[i] = band_min
				CBM_k[i] = KPT[single_band.index(band_min)]
				nCBM[i] = n+1
				break

	if metal == 1:
		VB_max = []
		CB_min = []
		spin_index = []
		for i in range(len(KPT)):
			VB_max.append(Band[nVBM[0]-1][i][0])
			CB_min.append(Band[nVBM[0]][i][0])
#			spin_index.append([0,0])
			if not spin == '1':
				if VB_max[i] < Band[nVBM[1]-1][i][1]:
					VB_max[i] = Band[nVBM[1]-1][i][1]
#					spin_index[i][0] = 1
				if CB_min[i] > Band[nVBM[1]][i][1]:
					CB_min[i] = Band[nVBM[1]][i][1]
#					spin_index[i][1] = 1
			if CB_min[i] - VB_max[i] < 0.01:
				metal = 2
				break

	if metal == 0:
		total_VBM = max(VBM)
		total_CBM = min(CBM)
		gap = total_CBM-total_VBM

	elif metal == 1:
		gap = 0
	else:
		gap = 0

	return gap

# This function is for calulating band gap by finding CBM and VBM of band structure
def gap_estimation(target,fermi,spin,ncl,KPT,Band,nelect):
	VBM = []
	CBM = []
	nVBM = []
	nCBM = []
	VBM_k = []
	CBM_k = []
	for i in range(int(spin)):
		VBM.append(-10.0**6.0)
		CBM.append(10.0**6.0)
	metal = 0
	# VBM & CBM search for spin-polarized calculation
	nVB = -1
	for i in range(int(spin)) :
		for n in range(len(Band)) :
			if n == 0 :
				nVBM.append(n)
				nCBM.append(n)
				VBM_k.append(KPT[n][0][i])
				CBM_k.append(KPT[n][0][i])

			single_band = [Band[n][x][i] for x in range(len(KPT))]
			band_max = max(single_band)
			band_min = min(single_band)

			if band_max < fermi :
				if band_max >= VBM[i] :
					VBM[i] = band_max
					VBM_k[i] = KPT[single_band.index(band_max)]
					nVBM[i] = n+1
			elif band_min > fermi :
				if band_min < CBM[i] :
					CBM[i] = band_min
					CBM_k[i] = KPT[single_band.index(band_min)]
					nCBM[i] = n+1
			else :
				metal = 1
				VBM[i] = band_max
				VBM_k[i] = KPT[single_band.index(band_max)]
				nVBM[i] = n+1
				CBM[i] = band_min
				CBM_k[i] = KPT[single_band.index(band_min)]
				nCBM[i] = n+1
				break

	if metal == 1:
		VB_max = []
		CB_min = []
		spin_index = []
		for i in range(len(KPT)):
			VB_max.append(Band[nVBM[0]-1][i][0])
			CB_min.append(Band[nVBM[0]][i][0])
#			spin_index.append([0,0])
			if not spin == '1':
				if VB_max[i] < Band[nVBM[1]-1][i][1]:
					VB_max[i] = Band[nVBM[1]-1][i][1]
#					spin_index[i][0] = 1
				if CB_min[i] > Band[nVBM[1]][i][1]:
					CB_min[i] = Band[nVBM[1]][i][1]
#					spin_index[i][1] = 1
			if CB_min[i] - VB_max[i] < 0.01:
				metal = 2
				break

	if metal == 0:
		total_VBM = max(VBM)
		total_CBM = min(CBM)
		total_VBM_spin = VBM.index(total_VBM)
		total_CBM_spin = CBM.index(total_CBM)
		total_VBM_index = nVBM[total_VBM_spin]
		total_CBM_index = nCBM[total_CBM_spin]
		total_VBM_kpt = VBM_k[total_VBM_spin]
		total_CBM_kpt = CBM_k[total_CBM_spin]
		gap = total_CBM-total_VBM
		if gap < 0.01:
			gap = 0.0
			metal = 1
	elif metal == 1:
		total_VBM_kpt = KPT[VB_max.index(max(VB_max))]
		total_CBM_kpt = KPT[CB_min.index(min(CB_min))]
		gap = 0
	else:
		gap = 0

	gap_log = open(target+'/Band_gap.log', 'w')
	gap_simple = open(target+'/Band_gap', 'w')	# For auto_bin
	kpt_out = open(target+'/KPT', 'w')	# File for VBM & CBM k-point position
	if metal != 0:
		gap_final = 'metal'
		gap_log.write('This system is metallic.\n')
		gap_simple.write('     metal\n')
		if metal == 2 :
			gap_log.write('! If it is not hybrid calculation, additional search is required for hybrid calculation.\n')
			gap_simple.write('! If it is not hybrid calculation, additional search is required for hybrid calculation.\n')
		else :
			gap_log.write("Use the 'KPT' file for hybrid calculation\n")
			kpt_out.write("Example file\n               0\nReciprocal\n")
			kpt_out.write('    '+'    '.join(total_VBM_kpt)+'\t0\n')
			kpt_out.write('    '+'    '.join(total_CBM_kpt)+'\t0\n')
	else :
		gap_final = str(gap)
		direct = 0
		for i in range(3) :
			total_VBM_kpt[i] = str(float(total_VBM_kpt[i]))
			total_CBM_kpt[i] = str(float(total_CBM_kpt[i]))
			if not total_VBM_kpt[i] == total_CBM_kpt[i]:
				direct = 1
		if direct == 0 :
			gap_log.write('Band gap: '+str(gap)+' eV (Direct)\n\n')
		else :
			gap_log.write('Band gap: '+str(gap)+' eV (Indirect)\n\n')
		gap_simple.write('     '+str(gap))
		gap_log.write('VBM: '+'  '.join(total_VBM_kpt)+'   : '+str(total_VBM)+' eV\n')
		gap_log.write('CBM: '+'  '.join(total_CBM_kpt)+'   : '+str(total_CBM)+' eV\n')
		gap_log.write('\nnVBM: '+str(total_VBM_index)+'  spin: '+str(total_VBM_spin+1)+'\n')
		gap_log.write('nCBM: '+str(total_CBM_index)+'  spin: '+str(total_CBM_spin+1)+'\n')
		kpt_out.write("Example file\n               0\nReciprocal\n")
		kpt_out.write('    '+'    '.join(total_VBM_kpt)+'\t0\n')
		kpt_out.write('    '+'    '.join(total_CBM_kpt)+'\t0\n')
	gap_log.close()
	gap_simple.close()
	kpt_out.close()
	return gap_final

# This function is for making band data file for plotting
def make_band_dat(xtic_file,Band,spin,target):
	xtic = []
	with open(xtic_file,'r') as inp:
		for line in inp:
			xtic.append(line.split()[0])
	with open(target+'/band.dat','w') as out:
		for k in range(len(xtic)):			# plot k-pts
			out.write(xtic[k])
			for i in range(int(spin)):		# spin
				for n in range(len(Band)):	# band idx
					out.write('\t'+str(Band[n][k][i]))
			out.write('\n')

# This function is for making input file for gnuplot
def make_band_in(title,xlabel_file,fermi,gap,nband,spin,plot_range,target):
	with open(target+'/band.in','w') as out:
		out.write("set terminal pdfcairo enhanced color font 'Arial, 14' size 7.2,5.4\n")
#		out.write("set output '"+target+'/'+title+".eps'\n")
		out.write("set output '"+title+".pdf'\n")
		title2 = '-'.join(title.split('_'))
		out.write("set title '"+title2+"' font 'Arial, 16'\n")
		if not spin == '2':
			out.write("set nokey\n")
		out.write('e = ' +fermi+'\n')
		out.write('set termoption dash\n')
		out.write('set yr[-'+str(plot_range[0])+':'+str(gap+plot_range[1])+']\n')
		with open(xlabel_file,'r') as inp:
			out.write(inp.read())
		out.write('f(x)=0\n')
		out.write('plot \\\n')
		if spin == '2':
			for i in range(nband-1):
				out.write("\t'band.dat' u 1:($"+str(i+2)+"-e) w l lt 1 lc rgb 'red' notitle,\\\n")
			out.write("\t'band.dat' u 1:($"+str(nband+1)+"-e) w l lt 1 lc rgb 'red' title 'up',\\\n")
			for i in range(nband-1):
				out.write("\t'band.dat' u 1:($"+str(nband+i+2)+"-e) w l lt 1 lc rgb 'blue' notitle,\\\n")
			out.write("\t'band.dat' u 1:($"+str(nband*2+1)+"-e) w l lt 1 lc rgb 'blue' title 'down',\\\n")
		else:
			for i in range(nband):
				out.write("\t'band.dat' u 1:($"+str(i+2)+"-e) w l lt 1 lc rgb 'black' notitle,\\\n")
		out.write("\tf(x) w l lt 0 notitle\n\n")
		out.write("set terminal pngcairo enhanced font 'Arial, 14' size 960,720\n")
		out.write("set output '"+title+".png'\n")
		out.write("replot")

def plot_band_structure(spin,Band,fermi,xtic_file,xlabel_file,plot_range,target):
	# make band.dat
	make_band_dat(xtic_file,Band,spin,target)

	# make band.in
	nband = len(Band)
	if pyhead(target+'/Band_gap.log',1).split()[2] == 'is' :
		gap = 0
		fermi = str(fermi)
	else :
		gap = round(float(pyhead(target+'/Band_gap.log',1).split()[2]))
		fermi = pygrep('VBM',target+'/Band_gap.log',0,0).splitlines()[0].split()[-2]

	title = target.split('/')[-2]
	make_band_in(title,xlabel_file,fermi,gap,nband,spin,plot_range,target)

def band_warning(Band,target):
	warn = []
	err_E_diff = 0.5	# if the energy differences between neighboring two points are larger than max_E_diff, return warning.
	for i in range(len(Band[0][0])):
		for n in range(len(Band)):
			for k in range(1,len(Band[0])-1):
				if abs(Band[n][k][i]-Band[n][k-1][i]) > err_E_diff and abs(Band[n][k][i]-Band[n][k+1][i]) > err_E_diff:
					if (Band[n][k][i]-Band[n][k-1][i])*(Band[n][k][i]-Band[n][k+1][i]) > 0:	# Error point should be extereme point.
						swit = 1
						for j in range(len(warn)):
							if k == warn[j][0]:
								warn[j][1].append(n)
								swit = 0
							break
						if swit == 1:
							warn.append([k,[n]])
	if len(warn) > 0:
		make_amp2_log(target,'Warning!! More than One of eigen value is deviated significantly from the line.\nThe list is followed:')
		for warn_list in warn:
			make_amp2_log(target,'\tIndex of kpts: '+str(warn_list[0]+1)+'\t from '+str(warn_list[1][0]+1)+' th band')

	return warn

def check_half_metal(target):
	spin = pygrep('ISPIN',target+'/OUTCAR',0,0).split()[2]
	ncl = pygrep('NONCOL',target+'/OUTCAR',0,0).split()[2]
	if not (spin == '2' and ncl =='F'):
		half_metal = 0
	else:
		[KPT,Band,nelect] = EIGEN_to_array(target+'/EIGENVAL',spin)
		fermi = get_fermi_level(Band,nelect,ncl)
		half_metal = 0
		for i in range(int(spin)) :
			metal = 0
			for n in range(len(Band)) :
				single_band = [Band[n][x][i] for x in range(len(KPT))]
				band_max = max(single_band)
				band_min = min(single_band)
				if band_max > fermi and band_min < fermi:
					metal = 1
			if metal = 0:
				half_metal = 1
	return half_metal
