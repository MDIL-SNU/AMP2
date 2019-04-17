###########################################
### Date: 2018-12-05			###
### yybbyb@snu.ac.kr			###
###########################################
import os,sys,subprocess,yaml,collections
def input_conf(conf):
	home = os.getcwd()
	with open(conf,'r') as inp:
		inp_yaml = yaml.load(inp)
	conf0 = inp_yaml['directory']['src_path']
	with open(conf0+'/config_def.yaml','r') as inp0:
		inp0_yaml = yaml.load(inp0)
	## default setting
#	for dir_key in inp0_yaml['directory'].keys():
#		inp0_yaml['directory'][dir_key] = home+inp0_yaml['directory'][dir_key]
	
	inp_override(inp0_yaml,inp_yaml)

	for dir_key in inp0_yaml['directory'].keys():
		inp0_yaml['directory'][dir_key] = os.path.expanduser(inp0_yaml['directory'][dir_key])	# home to absolute path
		inp0_yaml['directory'][dir_key] = os.path.abspath(inp0_yaml['directory'][dir_key])	# absolute path
		if dir_key == 'src_path' or dir_key == 'pot_path' or dir_key == 'submit':
			if not os.path.isdir(inp0_yaml['directory'][dir_key]):
				print('ERROR. check the config_directory')
				sys.exit()
		else:
			if not os.path.isdir(inp0_yaml['directory'][dir_key]):
				os.mkdir(inp0_yaml['directory'][dir_key])
	for dir_key in inp0_yaml['program'].keys():
		inp0_yaml['program'][dir_key] = os.path.expanduser(inp0_yaml['program'][dir_key])	# home to absolute path
		inp0_yaml['program'][dir_key] = os.path.abspath(inp0_yaml['program'][dir_key])	# absolute path

	if not os.path.isfile(inp0_yaml['program']['gnuplot']):	# gnuplot path check, if not, do not plot the figure
		inp0_yaml['calculation']['plot'] = 0


#	band_overlap_list = [['effective_mass','effective_mass']]
#	for calc_type in band_overlap_list:
#		if inp0_yaml['calculation'][calc_type[0]] == 1:
#			for pot_type in inp0_yaml[calc_type[1]]['potential_type']:
#				if not pot_type in inp0_yaml['band_calculation']['potential_type']:
#					inp0_yaml['band_calculation']['potential_type'].append(pot_type)

#	relax_overlap_list = [['band','band_calculation'],['density_of_states','density_of_states'],['dielectric','dielectric'],['effective_mass','effective_mass']]
#	for calc_type in relax_overlap_list:
#		if inp0_yaml['calculation'][calc_type[0]] == 1:
#			for pot_type in inp0_yaml[calc_type[1]]['potential_type']:
#				if not pot_type in inp0_yaml['relaxation']['potential_type']:
#					inp0_yaml['relaxation']['potential_type'].append(pot_type)


	conf_fin = home+'/config_fin.yaml'
	with open(conf_fin,'w') as inp_fin:
		yaml.dump(inp0_yaml,inp_fin,default_flow_style = False)
#	return inp0_yaml
	return conf_fin
		
def inp_override(source,override):
	for key in override.keys():
		if isinstance(source, collections.Mapping):
			if isinstance(override[key], collections.Mapping) and override[key]:
				returned = inp_override(source.get(key, {}),override[key])
				source[key] = returned
			else:
				source[key] = override[key]
		else:
			source = {key: override[key]}
	return source
