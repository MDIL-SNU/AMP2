======
Output
======

Output files
============

AMP\ :sup:`2`\  makes directory for each configuration file as its name (from name.cif or POSCAR_name). 
When the calculation is on progress, the directory is 
placed in output path in the configuration. If calculation is well finished, the calculation directory 
is moved to done path. If any error breaks out, it is move to error path.

Each directory includes several sub-directory as follow;

INPUT0
------
    Directory for input files for VASP calculation.

    - POSCAR_rlx_POT:
        Optimized structure file with POT functional.
    
    - KPOINTS:
        Converged k-points file

    - INCAR:
        VASP input file with converged cutoff energy and ground-state magnetic ordering 

kptest
------
    Directory for k-point convergence test.

    - kpoint.log: 
        Calculation log for k-points convergence test

encut
-----
    Directory for cutoff energy convergence test.

    - cutoff.log: 
        Calculation log for cutoff energy convergence test

relax_POT (POT = GGA or LDA)
---------------------------------
    Directory for structure relaxation.

magnetic_ordering
-----------------
    Directory for identifying magnetic spin ordering.

band_POT (POT = GGA or LDA)
--------------------------------
    Directory for band structure and band gap calculation.

dos_POT (POT = GGA or LDA)
-------------------------------
    Directory for density of states calculation.

dielectric_POT (POT = GGA or LDA)
--------------------------------------
    Directory for dielectric constant calculation.

hybrid_POT1_POT2 (POT = GGA or LDA)
----------------------------------------
    Directory for band gap calculation with hybrid oneshot scheme.

effm_POT (POT = GGA or LDA)
--------------------------------
    Directory for effective mass calculation.

Results
-------
    Directory for calculation results.

    - POSCAR_GGA:
        Optimized structure

    - Band_gap_GGA.log: 
        Information of band gap

    - band_GGA.png (band_GGA.pdf): 
        Band structure image

    - band_corrected.png (band_corrected.pdf): 
        Corrected band structure image

    - Band_gap_hybrid_GGA.log: 
        Information of band gap with HSE@PBE scheme
        
    - dos_GGA.png (dos_GGA.pdf): 
        Density of states image

    - dielectric_GGA.log: 
        Information of dielectric constant

    - effective_mass_hole_GGA.log: 
        Information of effective mass of hole

    - effective_mass_electron_GGA.log: 
        Information of effective mass of electron

    - Properties.json: 
        Summarized material properties

INPUT0_old
----------
    Directory for input files for VASP calculation with ferromagnetic ordering.
    If more stable magnetic spin ordering is obsevred, this directory is made.

relax_POT_old (POT = GGA or LDA)
-------------------------------------
    Directory for structure relaxation with ferromagnetic ordering.
    If more stable magnetic spin ordering is obsevred, this directory is made.

name_with_U
-----------
    Directory for AMP\ :sup:`2`\  calculation with DFT+U calculation.
    If the material is metallic and DFT+U calculation has been conducted,
    all of results move to this directory.

Additionally, AMP\ :sup:`2`\  provides log file as amp2.log for tracing the calculation.
