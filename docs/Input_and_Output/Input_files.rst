Input files
===========

Configuration file
------------------
    Configuration file for AMP2. See :doc:`/Input_and_Output/Configuration/Configuration`

Structure file
--------------
    The valid formats for structure file are that for VASP and cif format. In the cif files, 
    symmetry operator (_space_group_symop_[] or _symmetry_equiv_[]), atomic label (_atom_site_label),
    occupancy (_atom_site_occupancy) and fractional positions (_atom_site_fract_[]) must be included.
    The name of structure files must be formatted as name.cif or POSCAR_name where tag is used for identification.
   
Etc.
----
    The user can optionally provide another input files for DFT calculation such as potential file (POTCAR),
    k-points file (KPOINTS) and initial input file (INCAR) when these files are placed at the same directory and 
    name and tag are the same (ex. POTCAR_name, KPOINTS_name, ...).