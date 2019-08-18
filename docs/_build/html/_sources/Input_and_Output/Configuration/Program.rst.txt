Program
-------
The all tags in program determine the path of program except mpi_command.

    - vasp_std:
        vasp_std tag should be set to be the path for standard version of VASP.

        Usage:
        ::
            Program:
              vasp_std: [path]
        Default:
        ::
            Program:
              vasp_std: ./vasp_std

    - vasp_gam:
        vasp_gam tag should be set to be the path for gamma only version of VASP.

        Usage:
        ::
            Program:
              vasp_gam: [path]
        Default:
        ::
            Program:
              vasp_gam: ./vasp_gam

    - vasp_ncl:
        vasp_ncl tag should be set to be the path for non-collinear version of VASP.
        Though wrong path is set, most of calculations except spin-orbit coupling calculation can be conducted.

        Usage:
        ::
            Program:
              vasp_ncl: [path]
        Default:
        ::
            Program:
              vasp_ncl: ./vasp_ncl

    - gnuplot:
        gnuplot tag should be set to be the path for gnuplot.
        Though wrong path is set, most of calculations except drawing images can be conducted.

        Usage:
        ::
            Program:
              gnuplot: [path]
        Default:
        ::
            Program:
              gnuplot: /gnuplot

    - mpi_command:
        mpi_command tag should be set to be the operation command to conduct parallel computing calculation.

        Usage:
        ::
            Program:
              mpi_command: [command]
        Default:
        ::
            Program:
              mpi_command: mpirun
