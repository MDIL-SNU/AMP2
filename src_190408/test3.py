from module_vasprun import *
from module_amp2_input import *
inp_file = sys.argv[1]
with open(inp_file,'r') as f:
	inp_yaml = yaml.load(f)

print inp_yaml['cif2vasp']['INCAR']

incar_from_yaml(sys.argv[2],inp_yaml['cif2vasp']['INCAR'])
