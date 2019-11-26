###########################################
### Date: 2018-12-05            ###
### yybbyb@snu.ac.kr            ###
###########################################
import shutil, os, sys, subprocess, yaml
from module_log import *
from module_vasprun import *
from module_band import *
from module_hse import *
from input_conf import set_on_off

code_data = 'Version 0.9.3. Modified at 2019-11-22'

dir = sys.argv[1]

inp_file = sys.argv[2]
with open(inp_file,'r') as f:
    inp_yaml = yaml.safe_load(f)
ERROR_path = inp_yaml['directory']['error']
src_path = inp_yaml['directory']['src_path']
vasp_std = inp_yaml['program']['vasp_std']
vasp_gam = inp_yaml['program']['vasp_gam']
vasp_ncl = inp_yaml['program']['vasp_ncl']
mpi = inp_yaml['program']['mpi_command']
gnuplot = inp_yaml['program']['gnuplot']
npar = inp_yaml['vasp_parallel']['npar']
kpar = inp_yaml['vasp_parallel']['kpar']
inp_band = inp_yaml['band_calculation']

node = node_simple(sys.argv[3])
nproc = sys.argv[4]

pot_type = sys.argv[5]
if pot_type == 'LDA':
    POT = 'LDA'
else:
    POT = 'GGA'

# Set directory for input structure and INCAR
dir_band = dir+'/band_'+pot_type

# Check existing data
if os.path.isdir(dir_band) and os.path.isfile(dir_band+'/Band_gap.log') :
    make_amp2_log_default(dir,src_path,'Band calculation with '+pot_type+' potential.',node,code_data)
    gap = pyhead(dir_band+'/Band_gap.log',1).split()[2]
#   gap = subprocess.check_output(['head','-n','1',dir_band+'/Band_gap.log']).split()[2]
#   print('Success! Band gap: '+gap)
    if gap == 'is':
        make_amp2_log(dir,'Band calculation is already done.\nIt is metallic.')
    else:
        make_amp2_log(dir,'Band calculation is already done.\nBand gap is '+gap+' eV.')
    print 1
    sys.exit()

if not os.path.isdir(dir_band):
    os.mkdir(dir_band,0755)

os.chdir(dir_band)

make_amp2_log_default(dir_band,src_path,'Band calculation with '+pot_type+' potential.',node,code_data)

# check relax calculation
no_rlx = 0
if not os.path.isdir(dir+'/relax_'+pot_type):
    make_amp2_log(dir_band,'Directory of relaxation does not exist.')
    no_rlx = 1
else:
    if not os.path.isfile(dir+'/relax_'+pot_type+'/CONTCAR'):
        make_amp2_log(dir_band,'CONTCAR file in relaxation does not exist.')
        no_rlx = 1
    elif count_line(dir+'/relax_'+pot_type+'/CONTCAR') < 9:
        make_amp2_log(dir_band,'CONTCAR file in relaxation is invalid.')
        no_rlx = 1

if no_rlx == 1 and set_on_off(inp_band['relax_check']) == 1:
    print 0
    sys.exit()
# Potential type
if pot_type == 'HSE':
    if no_rlx == 1:
        copy_input(dir+'/INPUT0',dir_band,POT)
    else:
        copy_input_cont(dir+'/relax_'+pot_type,dir_band)
    shutil.copy(dir+'/INPUT0/sym',dir_band+'/.')

    # make INCAR for CHGCAR
    incar_from_yaml(dir_band,inp_band['incar'])
    if no_rlx == 1:
        mag_on = 2
    else:
        mag_on = check_magnet(dir+'/relax_'+pot_type,inp_yaml['magnetic_ordering']['minimum_moment'])
    vasprun = make_incar_for_ncl(dir_band,mag_on,kpar,npar,vasp_std,vasp_gam,vasp_ncl)
    incar_for_hse(dir_band+'/INCAR')
    hse_algo = pygrep('ALGO',dir+'/relax_'+pot_type+'/INCAR',0,0).split()[2]
#   hse_algo = subprocess.check_output(['grep','ALGO',dir+'/relax_'+pot_type+'/INCAR']).split()[2]
    wincar(dir_band+'/INCAR',dir_band+'/INCAR',[['NSW','0'],['ALGO',hse_algo]],[])
    # Write k-points for band structure and band gap search
    make_sym_for_band(dir_band+'/POSCAR',dir_band+'/sym',dir_band)
    [symk,order,xticlabel,rec] = make_symk(dir_band+'/sym')
    nkpt = make_kp_for_band(symk,order,xticlabel,rec,inp_band['kspacing_for_band'],inp_band['type_of_kpt'],dir_band)
    make_kpts_for_hse(dir+'/relax_'+pot_type+'/IBZKPT',dir_band+'/KPOINTS_band',dir_band,'band')
    if no_rlx == 1:
        mag_on = 2
    else:
        mag_on = check_magnet(dir+'/relax_'+pot_type,inp_yaml['magnetic_ordering']['minimum_moment'])

    incar_from_yaml(dir_band,inp_band['incar'])
    vasprun = make_incar_for_ncl(dir_band,mag_on,kpar,npar,vasp_std,vasp_gam,vasp_ncl)

    # Band stucture calculation
    out = run_vasp(dir_band,nproc,vasprun,mpi)
    if out == 1:  # error in vasp calculation
        print 0
        sys.exit() 

# pot_type is LDA or GGA.
else:
    # check CHGCAR
    if os.path.isfile(dir_band+'/CHGCAR') and os.path.getsize(dir_band+'/CHGCAR') > 0 :
        make_amp2_log(dir_band,'Band calculation is performed by using existing CHGCAR file.')
    else:
        make_amp2_log(dir_band,'VASP calculation for CHGCAR file.')
        # Copy input data and write CHGCAR
        if no_rlx == 1:
            copy_input(dir+'/INPUT0',dir_band,POT)
        else:
            copy_input_cont(dir+'/relax_'+pot_type,dir_band)
        shutil.copy(dir+'/INPUT0/sym',dir_band+'/.')

        # make INCAR for CHGCAR
        incar_from_yaml(dir_band,inp_band['incar'])
        if no_rlx == 1:
            mag_on = 2
        else:
            mag_on = check_magnet(dir+'/relax_'+pot_type,inp_yaml['magnetic_ordering']['minimum_moment'])
        vasprun = make_incar_for_ncl(dir_band,mag_on,kpar,npar,vasp_std,vasp_gam,vasp_ncl)
        wincar(dir_band+'/INCAR',dir_band+'/INCAR',[['NSW','0'],['LCHARG','.T.']],[])

        # VASP calculation for CHGCAR
        out = run_vasp(dir_band,nproc,vasprun,mpi)
        if out == 1:  # error in vasp calculation
            print 0
            sys.exit() 


        make_amp2_log(dir_band,'CHGCAR file is generated successfully.')

    # Write k-points for band structure and band gap search
    make_sym_for_band(dir_band+'/POSCAR',dir_band+'/sym',dir_band)
    [symk,order,xticlabel,rec] = make_symk(dir_band+'/sym')
    nkpt = make_kp_for_band(symk,order,xticlabel,rec,inp_band['kspacing_for_band'],inp_band['type_of_kpt'],dir_band)
    shutil.copy(dir_band+'/KPOINTS_band',dir_band+'/KPOINTS')
    if no_rlx == 1:
        mag_on = 2
    else:
        mag_on = check_magnet(dir+'/relax_'+pot_type,inp_yaml['magnetic_ordering']['minimum_moment'])

    incar_from_yaml(dir_band,inp_band['incar'])
    vasprun = make_incar_for_ncl(dir_band,mag_on,kpar,npar,vasp_std,vasp_gam,vasp_ncl)
    wincar(dir_band+'/INCAR',dir_band+'/INCAR',[['ISTART','1'],['ICHARG','11'],['LCHARG','.F.'],['LWAVE','.T.']],[])

    # Band stucture calculation
    out = run_vasp(dir_band,nproc,vasprun,mpi)
    if out == 1:  # error in vasp calculation
        print 0
        sys.exit() 

    out = electronic_step_convergence_check(dir_band)

    while out == 1:
        make_amp2_log(dir_band,'Calculation options are changed. New calculation starts.')
        out = run_vasp(dir_band,nproc,vasprun,mpi)
        if out == 1:  # error in vasp calculation
            print 0
            sys.exit()
        out = electronic_step_convergence_check(dir_band)

    if out == 2:  # electronic step is not converged. (algo = normal)
        make_amp2_log(dir_band,'The calculation stops but electronic step is not converged.')
        print 0
        sys.exit()

# Extract data from VASP output files
spin = pygrep('ISPIN',dir_band+'/OUTCAR',0,0).split()[2]
#spin = subprocess.check_output(['grep','ISPIN',dir_band+'/OUTCAR']).split()[2]
ncl = pygrep('NONCOL',dir_band+'/OUTCAR',0,0).split()[2]
#ncl = subprocess.check_output(['grep','NONCOL',dir_band+'/OUTCAR']).split()[2]
[KPT,Band,nelect] = EIGEN_to_array(dir_band+'/EIGENVAL',spin)
fermi = get_fermi_level(Band,nelect,ncl)

# Error checking for deviated eigen value
warn = band_warning(Band,dir_band)

# Band gap calculation
gap = gap_estimation(dir_band,fermi,spin,ncl,KPT,Band,nelect) # gap is string

if gap == 'metal':
    make_amp2_log(dir,'Band calculation is already done.\nIt is metallic.')
    if os.path.isdir(dir+'/dos_'+POT):
        fermi = float(pyhead(dir+'/dos_'+POT+'/DOSCAR',6).splitlines()[-1].split()[3])
#       fermi = float(subprocess.check_output(['head',dir+'/dos_'+POT+'/DOSCAR','-n','6']).splitlines()[-1].split()[3])
    elif os.path.isdir(dir+'/relax_'+POT):
        fermi = float(pyhead(dir+'/relax_'+POT+'/DOSCAR',6).splitlines()[-1].split()[3])
#       fermi = float(subprocess.check_output(['head',dir+'/relax_'+POT+'/DOSCAR','-n','6']).splitlines()[-1].split()[3])
else:
    make_amp2_log(dir_band,'Band calculation is done.\nBand gap is '+gap+' eV.')
    os.remove(dir_band+'/WAVECAR')
with open(dir_band+'/KPOINTS_band','r') as kpt_inp:
    nkpt_for_band = int(kpt_inp.readlines()[1].split()[0])

# Drawing band structure (nkpt_for_band is used for hse band structure.)
plot_band_structure(spin,[x[-1*nkpt_for_band:] for x in Band],fermi,dir_band+'/xtic.dat',dir_band+'/xlabel.dat',[inp_band['y_min'],inp_band['y_max']],dir_band)

if inp_yaml['calculation']['plot'] == 1:
    os.chdir(dir_band)
    subprocess.call([gnuplot,dir_band+'/band.in'])

with open(dir_band+'/amp2.log','r') as amp2_log:
    with open(dir+'/amp2.log','a') as amp2_log_tot:
        amp2_log_tot.write(amp2_log.read())

print 1
#print('Success! Band gap: '+gap)
