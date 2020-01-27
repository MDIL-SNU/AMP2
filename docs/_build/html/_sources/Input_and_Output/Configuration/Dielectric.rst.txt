Dielectric
----------

Dielectric tensor is calculated using density functional perturbation theory calculation 
with twice of k-points in structure optimization. Dielectic calculation for metallic system is unphysical.
Thus, AMP2 conducts the calculation for semi-conducting or insulating materials.
The hybrid scheme and spin-orbit coupling calculation can not be used for dielectric calculation.

- potential_type:

  potential_type tag determines the functional scheme (LDA or GGA) for dielectric calculation. 
  Multiple functionals can be chosen.
  
  Usage:
  ::
    dielectric:
      potential_type:
        - GGA | LDA
        - GGA | LDA

  Default:
  ::
    dielectric:
      potential_type:
        - GGA

- relax_check:
  
  Relax_check tag determines whether the dielectric calculation is not conducted without optimized structure or
  the dielectric calculation is performed without optimized structure (using initial configuration). If relax_check is True
  and no optimization has been performed, AMP2 gives an error message.

  Usage:
  ::
    dielectric:
      relax_check: True | False
  Default:
  ::
    dielectric:
      relax_check: True

- metal_check:
  
  metal_check tag determines whether the dielectric calculation is conducted without band gap or not. If metal_check is True
  and no band calculation has been performed, AMP2 gives an error message.

  Usage:
  ::
    dielectric:
      metal_check: True | False
  Default:
  ::
    dielectric:
      metal_check: True

- incar:

  User can additionally modulate the INCAR for VASP calculation using this tag.
    
  Usage:
  ::
    dielectric:
      incar:
        [INCAR tag in VASP] : [INCAR command in VASP]
  Default:
  ::
    dielectric:
      incar:
        EDIFF: 1e-08
