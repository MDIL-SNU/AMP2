Hybrid_oneshot
--------------

Conventional density functional theory calculation like LDA and PBE underestimates band gap and somtimes it gives 
wrong results for small gap materials such as Ge and InAs. Thus, AMP2 performs hybrid calculation for accurate band gap.
In the previous study[], it is shown that accurate band gap can be obtained using extremum points (valence band maximum and
conduction band minimum) and optimized structure in PBE scheme. Since hybrid calculation demands high computational cost,
this approach is imposed in AMP2.

For the small gap materials with metallic band structure in PBE functionals, DOS (density of states) based correction scheme
is applied in AMP2. (See :doc:`/Input_and_Output/Configuration/small_gap_correction`)

Finally, AMP2 provides a method to select mixing parameter using permittivity since there is an inverse correlation between 
mixing parameter and permittivity. []

- potential_type: 
  potential_type tag determines the functional scheme (LDA, GGA or HSE) for optimized structure and for extremum points.
  If one functional scheme is inserted, same functional scheme is used for optimized structure and for extremum points.
  If two functional schemes are inserted, above one is used for optimized structure and below one is used for extremum points.

  Usage:
  ::
    hybrid_oneshot:
      potential_type:
        - GGA | LDA
        - - GGA | LDA | HSE
          - GGA | LDA | HSE
  Default:
  ::
    hybrid_oneshot:
      potential_type:
        - GGA

- alpha:
  alpha tag determines a mixing parameter for hybrid calculation. As we mentioned above,
  mixing parameter in PBE0 has a inverse correlation with permittivity. If alpha: auto is used,
  the mixing parameter is determined as one of permittivity.

  Usage:
  ::
    hybrid_oneshot:
      alpha: [real] | Auto
  Default:
  ::
    hybrid_oneshot:
      alpha: 0.25

- cutoff_df_dvb:
  cutoff_df_dvb tag controls :math:`D_{\textrm{VB}}/D_{\textrm{F}}` used to classify semiconductor candidates.

  Usage:
  ::
    hybrid_oneshot:
      cutoff_df_dvb: [real]
  Default:
  ::
    hybrid_oneshot:
      cutoff_df_dvb: 0.3

- band_structure_correction:
  band_structure_correction determines whether to conduct scissor-correction for band structure or not.

  Usage:
  ::
    hybrid_oneshot:
      band_structure_correction: True | False
  Default:
  ::
    hybrid_oneshot:
      band_structure_correction: True

- incar:

  User can additionally modulate the INCAR for VASP calculation using this tag.
    
  Usage:
  ::
    hybrid_oneshot:
      incar:
        [INCAR tag in VASP] : [INCAR command in VASP]
  Default:
  ::
    hybrid_oneshot:
      incar:
