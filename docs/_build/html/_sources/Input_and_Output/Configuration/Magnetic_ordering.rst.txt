Magnetic_ordering
-----------------

AMP2 performs spin-polarized calculation for the materials satisfying any of conditions in the following.

  Materials include metal elements whose valence electrons partially occupy d 
  orbitals (Group 4 ~ Group 12) or f orbitals (lanthanides and actinides).

  The total number of electrons in the unit cell is odd.

For that materials, AMP\ :sup:`2`\ supports automatic protocol to identify magnetic ordering.
By default, AMP\ :sup:`2`\ performs the convergence test and relaxation calculations with ferromagnetic spin ordering.
Then, according to the option, this code finds ground-state collinear spin ordering using the Ising model. 
The detail algorithm is provided in the paper. 
By using the Ising model, AMP2 drmatically reduces computational cost to find the ground-state magnetic spin ordering. However, 
if the number of atoms are pretty large, it is still difficult to find most stable spin ordering trying all combinations. Thus,
AMP2 uses genetic algorithm to identify spin ordering for efficiency if the number of atoms with magnetic moment is larger than 9.
Please see this page for details. (:doc:`/Input_and_Output/Configuration/genetic_algorithm`)

- potential_type:
    potential_type tag determines the functional scheme (LDA, GGA or HSE) for convergence test. Only one of them should be chosen.

    Usage:
    ::
        magnetic_ordering:
          potential_type: GGA | LDA | HSE
    Default:
    ::
        magnetic_ordering:
          potential_type: GGA

- minimum_moment: 
    minimum_moment tag determines the cutoff value whether to impose spin moments for the ions or not. 
    If a magnetic moment is lower then this cutoff value in DFT calculation, AMP2 do not impose initial spin moment
    for that ions in further property calculations. If all magnetic moments are lower than the cutoff value, AMP2 carries out
    spin-unpolarized calculations.

    Usage:
    ::
        magnetic_ordering:
          minimum_moment: [real]
    Default:
    ::
        magnetic_ordering:
          minimum_moment: 0.5

- from_relax: 
    from_relax tag determines that magnetic ordering calculation is performed whether after the structure optimization or 
    before the structure optimization. 

    Usage:
    ::
        magnetic_ordering:
          from_relax: True | False
    Default:
    ::
        magnetic_ordering:
          from_relax: True

- incar:

    User can additionally modulate the INCAR for VASP calculation using this tag.
    
    Usage:
    ::
        magnetic_ordering:
          incar:
            [INCAR tag in VASP] : [INCAR command in VASP]
    Default:
    ::
        magnetic_ordering:
          incar:



