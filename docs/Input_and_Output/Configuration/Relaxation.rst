Relaxation
----------

To optimize the structure, AMP2 performs structure relaxation until satisfying the convergence condition or 100 iteration steps.
For the optimization, AMP2 uses conjugate-gradient method. If it is not converged within specific iteration step, optimization
calculation is restarted from converged geometry in previous step until max iteration steps. Sometimes, the lattice parameters
largely changes. In that case, we give an error message in amp2.log file.

- potential_type:

  potential_type tag determines the functional scheme (LDA, GGA or HSE) for band calculation. 
  Multiple functionals can be chosen.
  
  Usage:
  ::
    relaxation:
      potential_type:
        - GGA | LDA | HSE
        - GGA | LDA | HSE

  Default:
  ::
    relaxation:
      potential_type:
        - GGA

- energy:
  energy tag determines convergence condition of energy for structure optimization. Unit is 'eV'. 
  If user set pressure or force for convergence condition, it does not work.
  
  Usage:
  ::
    relaxation:
      energy: [real]

  Default:
  ::
    relaxation:
      energy: -1

- pressure:
  pressure tag determines convergence condition of pressure for structure optimization. Unit is 'kbar'. 
  If it is negative, it does not work.
  
  Usage:
  ::
    relaxation:
      pressure: [real]

  Default:
  ::
    relaxation:
      pressure: 1e-06

- force:
  force tag determines convergence condition of force for structure optimization. Unit is 'eV/Ang'. 
  If it is negative, it does not work.

  Usage:
  ::
    relaxation:
      force: [real]

  Default:
  ::
    relaxation:
      force: 1e-06

- incar:

  User can additionally modulate the INCAR for VASP calculation using this tag.
    
  Usage:
  ::
    relaxation:
      incar:
        [INCAR tag in VASP] : [INCAR command in VASP]
  Default:
  ::
    relaxation:
      incar:
