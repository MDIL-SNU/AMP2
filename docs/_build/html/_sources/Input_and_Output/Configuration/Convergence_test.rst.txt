Convergence_test
----------------

For the accuracy and efficiency, AMP\ :sup:`2`\ carries out a convergence test for kpoint and cutoff energy.
In AMP\ :sup:`2`\, converged kpoints and cutoff energy are searched by increasing the number of kpoints and cutoff energy.
Then, if the differences of energy, pressure and force between specific point and its next two points are within 
convergence criteria, the specific point is determined as converged kpoints or cutoff energy.
AMP\ :sup:`2`\ determine the k-points along each direction by an integer :math:`n_{\textrm{k}}` as follows;

.. math::
    \begin{align}
        N_{i} &= \textrm{ceil}[b_{i}/(b_{\textrm{max}}/n_{k})] && (1)\\
        k_{i} &= {{2m-N_{i}-1} \over {2N_{i}}} && (2)
    \end{align}

where :math:`i` is the axis index, :math:`b_{i}` is the length of ith reciprocal lattice, 
:math:`b_{\textrm{max}}` is the maximum among :math:`b_{i}`.

- potential_type:
    potential_type tag determines the functional scheme (LDA, GGA or HSE) for convergence test. Only one of them should be chosen.

    Usage:
    ::
        convergence_test:
          potential_type: GGA | LDA | HSE
    Default:
    ::
        convergence_test:
          potential_type: GGA

- enconv: 
    enconv tag is convergence criteria for energy. The unit is 'eV/atom'. 
    If a negative value is inserted, AMP2 does not check the energy convergence.

    Usage:
    ::
        convergence_test:
          enconv: real
    Default:
    ::
        convergence_test:
          enconv: 0.01

- prconv:
    enconv tag is convergence criteria for pressure. The unit is 'kbar'. 
    If a negative value is inserted, AMP2 does not check the pressure convergence.

    Usage:
    ::
        convergence_test:
          prconv: real
    Default:
    ::
        convergence_test:
          prconv: 10

- foconv:
    foconv tag is convergence criteria for maximum difference in atomic force. The unit is 'eV/:math:`\textrm{\AA}`'. 
    If a negative value is inserted, AMP2 does not check the pressure convergence.

    Usage:
    ::
        convergence_test:
          foconv: real
    Default:
    ::
        convergence_test:
          foconv: -1

- initial_kpl:
    initial_kpl tag determines the initial :math:`n_{k}` for convergence test in the equation (1).

    Usage:
    ::
        convergence_test:
          initial_kpl: [integer]
    Default:
    ::
        convergence_test:
          initial_kpl: 1

- max_kpl:
    max_kpl tag determines the maximum :math:`n_{k}` for convergence test in the equation (1).

    Usage:
    ::
        convergence_test:
          max_kpl: [integer]
    Default:
    ::
        convergence_test:
          max_kpl: 20

- enstart:
    enstart tag determines the initial cutoff energy for convergence test.

    Usage:
    ::
        convergence_test:
          enstart: [integer]
    Default:
    ::
        convergence_test:
          enstart: 200

- enstep:
    enstep tag determines the interval between cutoff energy for convergence test.

    Usage:
    ::
        convergence_test:
          enstep: [integer]
    Default:
    ::
        convergence_test:
          enstep: 50

- enmax:
    enmax tag determines the maximum cutoff energy for convergence test.

    Usage:
    ::
        convergence_test:
          enmax: [integer]
    Default:
    ::
        convergence_test:
          enmax: 1000

- incar:
    User can additionally modulate the INCAR for VASP calculation using this tag.

    Ex) If you want to change ALGO as normal,
    ::
      convergence_test:
        incar:
          ALGO: Normal

