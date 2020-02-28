Advanced configuration
======================

For advanced users, AMP\ :sup:`2`\  provides some additional configuration parameters written in
the default configuration file ('/src/cpnfig_def.yaml').

config_def.yaml:
::
    directory:
      submit: ./Submit                      # the path of structure file or the directory containg structure files
      output: ./Output                      # the path of the directory where calculation is conducted
      done: ./Done                          # the path of the directory where results are saved
      error: ./ERROR                        # the path of the directory where the materials with error are saved
      src_path: ./src                       # the path of the directory of AMP2 source codes
      pot_path_gga: ./pot/PB                # the path of directory for GGA pseudopotential
      pot_path_lda: ./pot/LDA               # the path of directory for LDA pseudopotential

    program:
      vasp_std: ./vasp_std                  # the path of standard version of VASP
      vasp_gam: ./vasp_gam                  # the path of gamma-only version of VASP
      vasp_ncl: ./vasp_ncl                  # the path of noncollinear version of VASP
      gnuplot: /usr/local/bin/gnuplot       # the path of executable file for gnuplot
      mpi_command: mpirun                   # mpi command (ex. mpirun, mpiexec, ...)

    vasp_parallel:
      npar: 2                               # the number of bands that are treated in parallel. It is same to NPAR tag in VASP.
      kpar: 2                               # the number of kpoints that are treated in parallel. It is same to NPAR tag in VASP.

    calculation:
      kp_test: T                            # On/Off for convergence test of k-points
      encut_test: T                         # On/Off for convergence test of cutoff energy
      relaxation: T                         # On/Off for structure optimization
      magnetic_ordering: T                  # On/Off for calculation to identify most stable magnetic spin ordering
      band: T                               # On/Off for the calculation for band structure and band gap
      density_of_states: T                  # On/Off for the calculation for density of states
      hse_oneshot: T                        # On/Off for the calculation for HSE@PBE
      dielectric: T                         # On/Off for the calculation for dielectric constant
      effective_mass: T                     # On/Off for the calculation for effective mass
      potential_type: GGA                   # calculation scheme (LDA or GGA)

    cif2vasp:
      pot_name:                             # the pseudopotential potential for element.
        GGA:                                # (Ex. GGA:\n    Ge:Ge_d\n    Cu:Cu_pv)
        LDA:
      soc_target:                           # the elements to carry out spin-orbit coupling calculation (Ex. soc_target:\n    - Bi\n    - Pb)
      u_value:                              # U values for PBE+U calculation (Ex. u_value:\n    La: 7.5\n    Ce: 8.5)
      max_nelm: 100                         # the maximum number of iteration for structure optimization.

    convergence_test:
      enconv: 0.01                          # convergence condition for energy (eV/atom). Negative value indicates that energy is not used as the condition.
      prconv: 10                            # convergence condition for pressure (bar). Negative value indicates that pressure is not used as the condition.
      foconv: -1                            # convergence condition for force (eV/angst). Negative value indicates that force is not used as the condition.
      initial_kpl: 1                        # Minimum value for the convergence test of k-points. It corresponds to the largest mesh grid in the three directions.
      max_kpl: 20                           # Maximum value for the convergence test of k-points. It corresponds to the largest mesh grid in the three directions.
      enstart: 200                          # Minimum value for the convergence test of cutoff energy
      enstep: 50                            # Interval for the convergence test of cutoff energy
      enmax: 1000                           # Maximum value for the convergence test of cutoff energy
      potential_type: GGA                   # Calculation scheme for convergence test. User have to choose one potential among the GGA, LDA and HSE.

    relaxation:
      potential_type:                        # Calculation scheme for structure optimization. User can choose one or more potential among the GGA, LDA and HSE.
        - GGA  
      max_iteration: 10                     # The maximum iteration number from previously optimized structure
      converged_ionic_step: -1              # The tolerance of steps for iteration. Until the relaxation finishes within the tolerance, we iterate the structure relaxation from previously optimized structure. In negative value, it is neglected.
      length_tolerance: 0.002               # The tolerance of length (ratio). Until the relaxation finishes within the tolerance, we iterate the structure relaxation from previously optimized structure. In negative value, it is neglected.
      angle_tolerance: 0.01                 # The tolerance of angle (degrees). Until the relaxation finishes within the tolerance, we iterate the structure relaxation from previously optimized structure. In negative value, it is neglected.  
      energy: -1                            # The energy tolerance (eV) to break the loop for structure optimization in VASP. In negative value, it is neglected.
      pressure: 10                          # The pressure tolerance (bar) to break the loop for structure optimization in VASP. In negative value, it is neglected.
      force: 0.02                           # The force tolerance (eV/angst) to break the loop for structure optimization in VASP. In negative value, it is neglected.

    band_calculation:
      kspacing_for_band: 0.01               # The distance between adjacent points in the band structure (2pi/ang).
      type_of_kpt: all                      # Set the lines to calculate the band gap. In the 'all', AMP2 calculates the eigenvalues along the lines connecting every combination of high symmetric k-points. In the 'band', AMP2 calculates the eigenvalue along the line to draw band structure.
      y_min: 3                              # The minimum energy range for band structure.
      y_max: 2                              # The maximum energy range from conduction band minimum for band structure.
      potential_type:                       # Calculation scheme for band structure. User can choose one or more potential among the GGA, LDA and HSE.
        - GGA

    density_of_states:
      kp_multiplier: all                    # Multiplier for k-points for smooth figure.
      y_min: 3                              # The minimum energy range for density of states.
      y_max: 2                              # The maximum energy range from conduction band minimum for density of states.
      potential_type:                       # Calculation scheme for density of states. User can choose one or more potential among the GGA, LDA and HSE.
        - GGA

    hybrid_oneshot:
      alpha: 0.25                           # Mixing parameter for hybrid calculation. If 'Auto' is set, the mixing parameter is set to be one of permittivity and PBE0 calualation is performed.
      fermi_width: 0.3                      # The energy range for DF
      vb_dos_min: 1                         # The energy range for DVB
      vb_dos_max: 3                         # The energy range for DVB
      cutoff_df_dvb: 0.3                    # DF/DVB used to classify semiconductor candidates. (See paper)  
      band_structure_correction: True       # On/Off for the band structure correction
      potential_type:                       # The potential used for lattice parameter optimization and for identifying the points at VBM and VBM. If one variable is inserted, AMP2 uses the lattice parameter and the points of VBM and CBM with that potential. If two variables are inserted, AMP2 uses the lattice parameter with above potential and the points of VBM and CBM with below potential. (Ex. potential_type:\n    - - HSE\n      - GGA)
        - GGA

    dielectric:
      kp_multiplier: all                    # Multiplier for k-points for dielectric constant.
      potential_type:                       # Calculation scheme for dielectric constant. User can choose one or more potential among the GGA and LDA
        - GGA

    effective_mass:
      carrier_type:                         # carrier type of effective mass to be estimated
        - hole
        - electron
      temperature_for_fermi: 300            # The temperature to estimate the Fermi distribution
      fermi_for_cutoff: 0.99                # Boundary condition for valid Fermi distribution (1-f)

To get more accurate band gap
-----------------------------

We suggest two approaches to get more accurate band gap.

- Band calculation with hybrid functional

  In the basic version, the band calculation is performed using PBE scheme.
  However, users can add the tags below to use hybrid functional for structure
  optimization and band calculation.
  ::
    relaxation:
      potential_type:
        - HSE
    band_calculation:
      potential_type:
        - HSE

- Using HSE@PBE scheme with hybrid structure

  Second approach is still using HSE@PBE method but the optimized structure is 
  calculated using hybrid functional. Since the band calculation with hybrid functional
  is too expensive, the k-points corresponding to the VBM and CBM are determined by using
  GGA method. For this calculation, users can use the commands below. Here, if potential_type
  in hybrid_oneshot is the main category, the method tags (HSE and GGA) are child subcategory 
  not parent subcategory. Please be careful.
  ::
    relaxation:
      potential_type:
        - GGA
        - HSE
    hybrid_oneshot:
      potential_type:
        - - HSE
          - GGA

Organic crystal
---------------
Organic crystals usually have lower Young's modulus than inorganic materials.
Thus, the error in the structural parameters can be substantial and they require
high precision for calculation. The tags below can control the precision of calculation.
::
  cif2vasp:
    INCAR:
      EDIFF: 1e-08

  convergence_test:
    enconv: 0.001
    prconv: 1

  relaxation:
    pressure: 1
    force: 0.002
