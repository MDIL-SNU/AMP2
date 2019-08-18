Output files
============

AMP2 makes directory for each configuration file. When the calculation is on progress, the directory is 
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

kptest
------
    Directory for k-point convergence test.

encut
-----
    Directory for cutoff energy convergence test.

relax_POT (POT = GGA, LDA or HSE)
---------------------------------
    Directory for structure relaxation.

magnetic_ordering
-----------------
    Directory for identifying magnetic spin ordering.

band_POT (POT = GGA, LDA or HSE)
--------------------------------
    Directory for band structure and band gap calculation.

dos_POT (POT = GGA, LDA or HSE)
-------------------------------
    Directory for density of states calculation.

dielectric_POT (POT = GGA, LDA or HSE)
--------------------------------------
    Directory for dielectric constant calculation.

hybrid_POT1_POT2 (POT = GGA, LDA or HSE)
----------------------------------------
    Directory for band gap calculation with hybrid oneshot scheme.

effm_POT (POT = GGA, LDA or HSE)
--------------------------------
    Directory for effective mass calculation.

Results
-------
    Directory for calculation results.

INPUT0_old
----------
    Directory for input files for VASP calculation with ferromagnetic ordering.
    If more stable magnetic spin ordering is obsevred, this directory is made.

relax_POT_old (POT = GGA, LDA or HSE)
-------------
    Directory for structure relaxation with ferromagnetic ordering.
    If more stable magnetic spin ordering is obsevred, this directory is made.

name_with_U
-----------
    Directory for AMP2 calculation with DFT+U calculation.
    If the material is metallic and DFT+U calculation has been conducted,
    all of results move to this directory.

Additionally, AMP2 provides log file as amp2.log for tracing the calculation.