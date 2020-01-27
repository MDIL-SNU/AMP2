small gap correction
--------------------

Conventional density functional theory calculation like LDA and PBE underestimates band gap and somtimes it gives 
wrong results for small gap materials such as Ge and InAs.
For these small gap materials with metallic band structure in PBE functionals, DOS (density of states) based correction scheme
is applied in AMP2. The details are explained in the paper.

As we mentioned, AMP2 uses the overlap integral between pseudo-wavefunctions to draw gap-corrected band structure.
However, to calculation the overlap integral with pseudo-wavefunction at nearest neighbor k-point is insufficient
if the material has degenerated band. Thus, AMP2 considers two adjacent points to track the band line with similar
wave character.