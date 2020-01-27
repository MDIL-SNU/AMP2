Vasp_parallel
-------------

npar and kpar tags are used to enhance the efficiency of parallel computing calculation of VASP. 
    - npar:
        napr tag determines the number of bands that are treated in parallel. It is same to NPAR tag in VASP.

        Usage:
        ::
            vasp_parallel:
              npar: [integer]
        Default:
        ::
            vasp_parallel:
              npar: 2

    - kpar: 
        kpar tag determines the number of kpoints that are treated in parallel. It is same to NPAR tag in VASP.
        
        Usage:
        ::
            vasp_parallel:
              kpar: [integer]
        Default:
        ::
            vasp_parallel:
              kpar: 2


