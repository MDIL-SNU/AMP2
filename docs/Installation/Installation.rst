==========================
Installation and execution
==========================

Installation AMP\ :sup:`2`\
===========================

System requirements
-------------------
    AMP\ :sup:`2`\  supports Python 2.7 and 3. Currently, the package is not compatible
    with lower version than 2.7. AMP\ :sup:`2`\  is utilizes Python modules 
    in the following with link to each site.

    - numpy [https://www.numpy.org]
    - scipy [https://www.scipy.org]
    - spglib [https://atztogo.github.io/spglib]
    - PyYAML [https://pypi.org/project/PyYAML]

    These modules should be pre-installed. In addition, AMP\ :sup:`2`\  needs gnuplot to draw
    various figures.


Installation
------------

    To use AMP\ :sup:`2`\, please download the file from https://github.com/MDIL-SNU/AMP2 under the 
    working directory.

Essential setting
=================


Setting of configuration file
------------------------------

    AMP\ :sup:`2`\  uses YAML style configuration file. All setting parameters used in AMP\ :sup:`2`\  can
    be controlled in "config.yaml". Before using AMP\ :sup:`2`\, proper pathes and mpi program command 
    should be set to be suitable for your system. Following commands are the essential directories 
    and programs to be set.
    ::
        Directory:
          submit:
          src:
          pot_path_gga:
          pot_path_lda:
        Program:
          vasp_std:
          vasp_gam:
          vasp_ncl:
          gnuplot:
          mpi_command:

    Details for the commands are in :doc:`/Input/Configuration`.

Execution AMP\ :sup:`2`\
========================

You can execute AMP\ :sup:`2`\  using Python command as following.
::
    python [src_path]/main.py [path for configuration file] [path for nodefile] [the number of cores] 

- [src_path] is the path for directory of source codes for AMP\ :sup:`2`\.
- [path for configuration file] is the path for configuration file (config.yaml).
- [path for nodefile] is used to record the information of computing nodes such as PBS_nodefile in 
  Portable Batch System (PBS) and HOSTNAME in Sun Grid Engine (SGE). In the PBS system, we recommand to use
  the command, "echo $PBS_nodefile > nodefile". Also, users can save an arbitrary text
  by writing in the nodefile.
- [the number of cores] is the number of cores to be used in parallel computing.

For the convenience, we provide the shell script file (run.sh) as following.
::
    echo 'node information' > nodefile
    NPROC=16       # The number of cores for parallel computing

    ### set path of config.yaml ###
    conf=./config.yaml
    ###############################

    ### Do not change #############
    src_path=`grep 'src_path' $conf | tr -s ' ' | cut -d " " -f 3`
    ###############################

    python $src_path/main.py $conf nodefile $NPROC >& stdout.x

Before execution, you need to modify *'node information'*, *NPROC* and *conf*.
Then, you can execute AMP\ :sup:`2`\  using shell script as following.
::
    sh run.sh

The shell script file can be easily integrated with job scheduler program such 
as PBS.
