Input files
===========

Structure file
--------------
    The valid formats for structure file are that for VASP and cif format. In the cif files, 
    symmetry operator (_space_group_symop_[] or _symmetry_equiv_[]), atomic label (_atom_site_label),
    occupancy (_atom_site_occupancy) and fractional positions (_atom_site_fract_[]) must be included.
    The name of structure files must be formatted as name.cif or POSCAR_name where tag is used for identification.

    VASP structure file format:
    ::
        Primitive Cell
           1.000000000
              0.0    2.714895    2.714895
              2.714895    0.0    2.714895
              2.714895    2.714895    0.0
            Si
            2
        Selective dynamics
        Direct
            0.5    0.5    0.5  T  T  T ! Si1
            0.75    0.75    0.75  T  T  T ! Si1

Configuration
-------------

    All of parameters can be tuned in the configuration file as following. 
    The detail for each parameter is explained in :doc:`/Input/Configuration`.

    config.yaml:
    ::
        directory:
          submit: ./Submit                      # the path of structure file or the directory containg structure files
          output: ./Output                      # the path of the directory where calculation is conducted
          done: ./Done                          # the path of the directory where results are saved
          error: ./ERROR                        # the path of the directory where the materials with error are saved
          src_path: ./src                       # the path of the directory of AMP2 source codes
          pot_path_gga: ./pot/PBE               # the path of directory for GGA pseudopotential
          pot_path_lda: ./pot/LDA               # the path of directory for LDA pseudopotential

        program:
          vasp_std: ./vasp_std                  # the path of standard version of VASP
          vasp_gam: ./vasp_gam                  # the path of gamma-only version of VASP
          vasp_ncl: ./vasp_ncl                  # the path of noncollinear version of VASP
          gnuplot: /gnuplot                     # the path of executable file for gnuplot
          mpi_command: mpirun                   # mpi command (ex. mpirun, mpiexec, ...)

        calculation:
          magnetic_ordering: T                  # On/Off for the calculation to idetify most stable magnetic spin ordering
          band: T                               # On/Off for the calculation for band structure and band gap
          density_of_states: T                  # On/Off for the calculation for density of states
          hse_oneshot: T                        # On/Off for the calculation for HSE@PBE
          dielectric: T                         # On/Off for the calculation for dielectric constant
          effective_mass: T                     # On/Off for the calculation for effective mass
          potential_type: GGA                   # calculation scheme (LDA or GGA)

        vasp_parallel:
          npar: 2                               # the number of bands that are treated in parallel. It is same to NPAR tag in VASP.
          kpar: 2                               # the number of kpoints that are treated in parallel. It is same to NPAR tag in VASP.

        cif2vasp:
          pot_name:                             # the pseudopotential potential for element.
            GGA:                                # (Ex. GGA:\n    Ge:Ge_d\n    Cu:Cu_pv) 
            LDA: 
          soc_target:                           # the elements to carry out spin-orbit coupling calculation (Ex. soc_target:\n    - Bi\n    -Pb)
          u_value:                              # U values for PBE+U calculation (Ex. u_value:\n    La: 7.5\n    Ce: 8.5)

        hybrid_oneshot:
          alpha: 0.25                           # mixing parameter for hybrid calculation. If "Auto" is set, the mixing parameter is set to be one of permittivity and PBE0 calualation is performed.
          cutoff_df_dvb: 0.3                    # DF/DVB used to classify semiconductor candidates. (See paper)
          band_structure_correction: True       # On/Off for the band structure correction

        effective_mass:
          carrier_type:                         # carrier type of effective mass to be estimated
            - hole
            - electron   
