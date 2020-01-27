Effective_mass
--------------

In AMP2, effective mass tensor is estimated using semiclassical transport theory.
The details are explained in the paper.

- carrier_type:
    carrier_type tag determines the type of carrier (hole or electron) to be estimated.

    Usage:
    ::
        effective_mass:
          carrier_type:
            - hole | electron
    Default:
    ::
        effective_mass:
          carrier_type:
            - hole
            - electron
