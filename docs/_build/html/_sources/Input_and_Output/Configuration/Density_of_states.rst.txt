Density_of_states
-----------------

Density of states are estimated using twice of k-points used in optimization calculation.

- potential_type:

  potential_type tag determines the functional scheme (LDA, GGA or HSE) for density of states calculation. 
  Multiple functionals can be chosen.
  
  Usage:
  ::
    density_of_states:
      potential_type:
        - GGA | LDA | HSE
        - GGA | LDA | HSE

  Default:
  ::
    density_of_states:
      potential_type:
        - GGA

- relax_check:
  
  Relax_check tag determines whether the density of states calculation is not conducted without optimized structure or
  the density of states calculation is performed without optimized structure (using initial configuration). If relax_check is True
  and no optimization has been performed, AMP2 gives an error message.

  Usage:
  ::
    density_of_states:
      relax_check: True | False
  Default:
  ::
    density_of_states:
      relax_check: True


- y_min (y_max):

  y_min and y_max tags control the energy range of figure for density of states. The maximum energy range is set to be y_max + band gap.
  If band calculation has not been conducted, the maximum energy range is set to be y_max.

  Usage:
  ::
    density_of_states:
      y_min: [real]
      y_max: [real]
  Default:
  ::
    density_of_states:
      y_min: 3
      y_max: 2  

- incar:

  User can additionally modulate the INCAR for VASP calculation using this tag.
    
  Usage:
  ::
    density_of_states:
      incar:
        [INCAR tag in VASP] : [INCAR command in VASP]
  Default:
  ::
    density_of_states:
      incar:
