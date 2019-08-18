cif2vasp
--------

In AMP\ :sup:`2`\, input files for VASP calculation are automatically generated from structure files.
These parameters can control the initial input files for VASP. The potential file (POTCAR) is built using 
the preset pseudopotential. (Preset pseudopotential: :doc:`/Input_and_Output/Configuration/potential`)

- soc_target: 
  soc_target tag determines the elements to carry out spin-orbit coupling calculation. In AMP2, spin-orbit coupling calculation
  is performed only for band structure and density of states.

  Usage:
  ::
    cif2vasp:
      soc_target:
        - [element name]
        - Bi
  Default:
  ::
    cif2vasp:
      soc_target:

- u_value:
  u_value tag controls :math:`U` values for PBE + Hubbard :math:`U` method. By default, AMP2 imposes :math:`U` parameters for 3d
  transition metal. If all tag is used instead of element name, every :math:`U` value is set to be the target value.

  Usage:
  ::
    cif2vasp:
      u_value:
        - [element name]: real     
  Default:
  ::
    cif2vasp:
      u_value:
        V: 3.1
        Cr: 3.5
        Mn: 4
        Fe: 4
        Co: 3.3
        Ni: 6.4
        Cu: 4
        Zn: 7.5


- max_nelm:
  max_nelm tag control the maximum electronic selfconsistency steps. It is same to NELM tag in VASP.
  For magnetic ordering calculation, AMP2 uses two times of given value.

  Usage:
  ::
    cif2vasp:
      max_nelm: [integer]
  Default:
  ::
    cif2vasp:
      max_nelm: 100

- incar:

  User can additionally modulate the INCAR for VASP calculation using this tag.
    
  Usage:
  ::
    cif2vasp:
      incar:
        [INCAR tag in VASP] : [INCAR command in VASP]
  Default:
  ::
    cif2vasp:
      incar:
 