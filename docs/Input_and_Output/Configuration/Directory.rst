Directory
---------

All tags in directory define the path of directories used in AMP2.
If there is no directory in the path for Output, Done and ERROR, AMP2 makes new directories. 

- Submit:
    submit tag should be set to be a directory for target materials. In AMP2, user 
    can designate a specific material or a bunch of materials as target materials.
    To perform the AMP2 for a specific materials, submit path is set to be structure
    file or directory for continuous calculation. The valid formats for structure file are
    illustrated in :doc:`/Input_and_Output/Input_files`.
    The user can optionally provide another input files for DFT calculation such as potential file (POTCAR),
    k-points file (KPOINTS) and initial input file (INCAR) when these files are placed at the same directory and 
    name and tag are the same (ex. POTCAR_name, KPOINTS_name, ...).
    For calculating a bunch of materials, Submit path is set to be the directory where the valid strcture 
    format files and directories for continuous calculation are placed.

    Usage:
    ::
        directory:
          submit: [path of structure file] | [path of directory]
    Default:
    ::
        directory:
          submit: ./Submit

- Output:
    Output tag defines a path where the material on calculation is located.

    Usage:
    ::
        directory:
          output: [path of directory]
    Default:
    ::
        directory:
          output: ./Output

- Done:
    Done tag defines a path where calculated materials are saved.

    Usage:
    ::
        directory:
          done: [path of directory]
    Default:
    ::
        directory:
          done: ./Done

- Error:
    Output tag defines a path saving the materials in which calculation error broke out.

    Usage:
    ::
        directory:
          error: [path of directory]
    Default:
    ::
        directory:
          error: ./ERROR

- src_path:
    src_path tag should be set to be the directory for source codes.

    Usage:
    ::
        directory:
          src_path: [path of directory]
    Default:
    ::
        directory:
          src_path: ./src

- pot_path_GGA (pot_path_LDA):
    pot_path_GGA (pot_path_LDA) should be set to be the directory for pseudopotential provided by VASP. 

    Usage:
    ::
        directory:
          pot_path_GGA: [path of directory]
          pot_path_LDA: [path of directory]
    Default:
    ::
        directory:
          pot_path_GGA: ./PBE
          pot_path_LDA: ./LDA


