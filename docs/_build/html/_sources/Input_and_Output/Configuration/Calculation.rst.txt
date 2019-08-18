Calculation
-----------
The all tags in calculation determine whether the calculation is performed or not.

    - kp_test:
        kp_test tag determines whether the convergence test for k-points is performed or not.

        Usage:
        ::
            Calculation:
              kp_test: True | False
        Default:
        ::
            Calculation:
              kp_test: True

    - encut_test: 
        encut_test tag determines whether the Convergence test for cutoff energy is performed or not.

        Usage:
        ::
            Calculation:
              encut_test: True | False
        Default:
        ::
            Calculation:
              encut_test: True

    - relaxation:
        relaxation tag determines whether the structure optimization is performed or not.

        Usage:
        ::
            Calculation:
              relaxation: True | False
        Default:
        ::
            Calculation:
              relaxation: True
    
    - magnetic_ordering:
        magnetic_ordering tag determines whether to identify the most stable magnetic spin ordering or not.

        Usage:
        ::
            Calculation:
              magnetic_ordering: True | False
        Default:
        ::
            Calculation:
              magnetic_ordering: True

    - band:
        band tag determines whether to estimate the band gap and to draw band structure or not.

        Usage:
        ::
            Calculation:
              band: True | False
        Default:
        ::
            Calculation:
              band: True

    - density_of_states:
        density_of_states tag determines whether to estimate the density of states or not.

        Usage:
        ::
            Calculation:
              density_of_states: True | False
        Default:
        ::
            Calculation:
              density_of_states: True

    - hse_oneshot: 
        hse_oneshot tag determines whether to perform the hybrid calculation or not. This hybrid calculation
        is conducted without full band searching and structure optimization.

        Usage:
        ::
            Calculation:
              hse_oneshot: True | False
        Default:
        ::
            Calculation:
              hse_oneshot: True

    - dielectric: 
        dielectric tag determines whether to estimate the dielectric constant or not.

        Usage:
        ::
            Calculation:
              dielectric: True | False
        Default:
        ::
            Calculation:
              dielectric: True

    - effective_mass: 
        effective_mass tag determines whether to estimate the hole (and/or electron) effective mass or not.

        Usage:
        ::
            Calculation:
              effective_mass: True | False
        Default:
        ::
            Calculation:
              effective_mass: True


