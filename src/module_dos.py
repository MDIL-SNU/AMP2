###########################################
### Date: 2018-12-05			###
### yybbyb@snu.ac.kr			###
###########################################
# This is a package of modules for drawing density of states.
import os

# This function is for reading atomic information from poscar
def poscar_to_atom_inform(poscar_file):
	with open(poscar_file,'r') as pos_inp:
		pos_lines = pos_inp.readlines()
	atom_name = pos_lines[5].split()
	atom_num = [int(x) for x in pos_lines[6].split()]
	return [atom_name,atom_num]

# This function is for making total and partial dos data file from DOSCAR
def make_dos_dat(dos_file,spin,atom_num,ncl):
	with open(dos_file,'r') as dos_inp:
		dos = dos_inp.readlines()
	ndos = int(dos[5].split()[2]) + 1
	Ene = []
	Tot_dos = []
	# read total dos
	for i in range(6,5+ndos):
		Ene.append(float(dos[i].split()[0]))
		if spin == '1':
			Tot_dos.append([float(dos[i].split()[1])])
		else:
			Tot_dos.append([float(dos[i].split()[1]),float(dos[i].split()[2])])

	# read partial dos
	Obt_num = {'s':1,'p':3,'d':5,'f':7} 
	if ncl == 'F':
		Obt_dic = {3:[1,'sum','s','p','d'],6:[2,'sum','s','p','d'],9:[1,'no_sum','s','p','d'],18:[2,'no_sum','s','p','d'],4:[1,'sum','s','p','d','f'],8:[2,'sum','s','p','d','f'], 16:[1,'no_sum','s','p','d','f'], 32:[2,'no_sum','s','p','d','f']}
	else:
		Obt_dic = {12:[4,'sum','s','p','d'], 36:[4,'no_sum','s','p','d'], 16:[4,'sum','s','p','d','f'],64:[4,'no_sum','s','p','d','f']}
	par_dos = []
	idf = len(dos[6+ndos].split()) - 1
	for typ in range(len(atom_num)):
		atom_sum = sum(atom_num[:typ])+1
		par_dos.append([])
		for i in range(ndos-1):
			if ncl == 'T':
				dos_tmp = [[0] for x in range(len(Obt_dic[idf])-2)]
			else:
				dos_tmp = [[0 for y in range(Obt_dic[idf][0])] for x in range(len(Obt_dic[idf])-2)]
			for n in range(atom_num[typ]):
				idx = n + atom_sum
				line_tmp = [float(x) for x in dos[6+idx*ndos+i].split()[1:]]
				if ncl == 'T':
					dos_single = [[0] for x in range(len(Obt_dic[idf])-2)]
				else:
					dos_single = [[0 for y in range(Obt_dic[idf][0])] for x in range(len(Obt_dic[idf])-2)]
				line_idx = 0
				for orb in range(len(Obt_dic[idf])-2):
					if Obt_dic[idf][1] == 'sum':
						for spin_idx in range(Obt_dic[idf][0]):
							if ncl == 'T':
								if spin_idx == 0:
									dos_single[orb][spin_idx] = line_tmp[line_idx]
							else:
								dos_single[orb][spin_idx] = line_tmp[line_idx]
							line_idx = line_idx+1
					else:
						for orb_idx2 in range(Obt_num[Obt_dic[idf][orb+2]]):
							for spin_idx in range(Obt_dic[idf][0]):
								if ncl == 'T':
									if spin_idx == 0:
										dos_single[orb][spin_idx] = dos_single[orb][spin_idx]+line_tmp[line_idx]
								else:
									dos_single[orb][spin_idx] = dos_single[orb][spin_idx]+line_tmp[line_idx]
	
								line_idx = line_idx+1
				for orb in range(len(Obt_dic[idf])-2):
					dos_tmp[orb] = [x+y for x,y in zip(dos_tmp[orb],dos_single[orb])]
											
			par_dos[typ].append(dos_tmp) # [atom_type][energy][orbital][spin]
	return [Ene,Tot_dos,par_dos]

def write_tot_dos(Ene,Tot_dos,fermi,target):
	if not os.path.isdir(target+'/Pdos_dat'):
		os.mkdir(target+'/Pdos_dat',0o755)
	with open(target+'/Pdos_dat/Tot_dos.dat','w') as dos_out:
		dos_out.write('Energy')
		if len(Tot_dos[0]) == 1:
			dos_out.write('\ttot\n')
		else:
			dos_out.write('\ttot_up\ttot_down\n')
		for n in range(len(Ene)):
			dos_out.write(str(Ene[n]-fermi))
			for dos in Tot_dos[n]:
				dos_out.write('\t'+str(dos))
			dos_out.write('\n')

def write_par_dos(Ene,par_dos,atom_name,fermi,target):
	orb_name = ['s','p','d','f']
	if not os.path.isdir(target+'/Pdos_dat'):
		os.mkdir(target+'/Pdos_dat',0o755)
	for typ in range(len(atom_name)):
		with open(target+'/Pdos_dat/'+atom_name[typ]+'_dos.dat','w') as dos_out:
			dos_out.write('Energy')
			for orb in range(len(par_dos[typ][0])):
				dos_out.write('\t'+orb_name[orb])
				if len(par_dos[typ][0][0]) == 2:
					dos_out.write('_up\t'+orb_name[orb]+'_down')
			dos_out.write('\n')
			for n in range(len(Ene)):
				dos_out.write(str(Ene[n]-fermi))
				for dos_line in par_dos[typ][n]:
					for dos in dos_line:
						dos_out.write('\t'+str(dos))
				dos_out.write('\n')

# This function is for making input file for gnuplot
def make_dos_in(target,atom_name,spin,orb_len,plot_range):
	color_list = ['web-green','light-red','orange','web-blue','dark-violet','skyblue']
	with open(target+'/Pdos_dat/dos.in','w') as out:
		out.write("set terminal pdfcairo enhanced color font 'Arial, 14' size 3.6,5.4\n")
		out.write("set output 'dos_"+target.split('/')[-2]+".pdf'\n")
		out.write('set termoption dash\n')
		out.write("set ylabel 'Energy (eV)' font 'Arial, 14'\n")
		out.write("set xlabel 'Density of State' font 'Arial, 14'\n")
		out.write("set style fill transparent solid 0.25 noborder\n")
		out.write('set yr[-'+str(plot_range[0])+':'+str(plot_range[1])+']\n')
		out.write("plot 'Tot_dos.dat' u 2:1 w l lc rgb 'black' title 'Total'\n")
		for i in range(len(atom_name)):
			out.write("replot '"+atom_name[i]+"_dos.dat' u (")
			for k in range(orb_len-1):
				out.write("$"+str((k*int(spin)+2))+"+")
			out.write("$"+str(((orb_len-1)*int(spin)+2))+"):1 w filledcurves x1=0 lc rgb '"+color_list[i]+"' title '"+atom_name[i]+"'\n")

			out.write("replot '"+atom_name[i]+"_dos.dat' u (")
			for k in range(orb_len-1):
				out.write("$"+str((k*int(spin)+2))+"+")
			out.write("$"+str(((orb_len-1)*int(spin)+2))+"):1 w l lc rgb '"+color_list[i]+"' notitle\n")
		if spin == '2':
			out.write("replot 'Tot_dos.dat' u (0-$3):1 w l lc rgb 'black' notitle\n")
			for i in range(len(atom_name)):
				out.write("replot '"+atom_name[i]+"_dos.dat' u (0-")
				for k in range(orb_len-1):
					out.write("$"+str((k*2+3))+"-")
				out.write("$"+str(((orb_len-1)*2+3))+"):1 w filledcurves x1=0 lc rgb '"+color_list[i]+"' notitle\n")

				out.write("replot '"+atom_name[i]+"_dos.dat' u (0-")
				for k in range(orb_len-1):
					out.write("$"+str((k*2+3))+"-")
				out.write("$"+str(((orb_len-1)*2+3))+"):1 w l lc rgb '"+color_list[i]+"' notitle\n")

		out.write("\nset terminal pngcairo enhanced font 'Arial, 14' size 480,720\n")
		out.write("set output 'dos_"+target.split('/')[-2]+".png'\n")
		out.write("replot")

