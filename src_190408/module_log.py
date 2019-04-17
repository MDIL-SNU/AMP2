###########################################
### Date: 2018-12-05			###
### yybbyb@snu.ac.kr			###
###########################################
import datetime, getpass, os
def make_amp2_log_default(path,src_path,comment,node,code_data):
	with open(path+'/amp2.log','a') as amp2_log:
		amp2_start_time = datetime.datetime.now().strftime("\tThe calculation is started at %Y-%m-%d : %H:%M \n")
		amp2_log.write('['+comment+']\n')
		user_name = getpass.getuser()
		amp2_log.write('\tThe calculation is performed by '+user_name+'\n')
		amp2_log.write(amp2_start_time)
		full_path = os.path.abspath(path)
		full_src_path = os.path.abspath(src_path)
		amp2_log.write('\tThe current running path is '+full_path+'\n')
		amp2_log.write('\tThe source path is '+full_src_path+'\n')
		amp2_log.write('\t'+code_data+'\n')
		amp2_log.write('\tRunning nodes are')
		for node_index in node:
			amp2_log.write(' '+node_index)
		amp2_log.write('\n')

def make_amp2_log(path,comment):
	lines = comment.splitlines()
	with open(path+'/amp2.log','a') as amp2_log:
		for line in lines:
			amp2_log.write('\t'+line+'\n')

def node_simple(node_file):
	node = []
	with open(node_file,'r') as nodefile:
		for ll in nodefile.readlines():
			if len(node) == 0:
				node.append(ll.split()[0])
			elif not ll.split()[0] in node:
				node.append(ll.split()[0])
	return node

def read_code_head(code,head_num):
	code_ver = subprocess.check_output(['head','-'+str(head_num),code])
	return code_ver

def write_log_in_outcar(outcar_file,log_file):
	if os.path.isfile(outcar_file) and os.path.isfile(log_file):
		with open(log_file,'r') as log_inp:
			log = log_inp.read()
		with open(outcar_file,'r') as out:
			outcar = out.read()
		with open(outcar_file,'w') as out:
			out.write(log+'\n')
			out.write(outcar)
