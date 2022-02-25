# PDB_generation_from_LAMMPS_Trajectory_file


Author: Amirhossein Saali


## Project Overview

This project is aimed to generate pdb file from Lammps trajectory file. 

Using a pdb file can facilitate researchers to visualize protein structure with more atomistic details. Moreover, pdb files contain information required for every type of structural studies.

This project is coded in Fortran as this programming language is one of the fastest of its kind. 

### Required Files

in order to run this program, two Fortran files (.f90) are provided. There is a [sample file](https://github.com/Amirsaali/PDB_generation_from_LAMMPS_Trajectory_file/blob/main/pdb_SAMPLE.pdb) which is necessary for generating pdb file. This sample file must be compatible with your lammps data file. [Trajectory_input.txt](https://github.com/Amirsaali/PDB_generation_from_LAMMPS_Trajectory_file/blob/main/Trajectory_input.txt) is Lammps trajectory file. if you prefer to preseve trajectory file format, you need to open [Trajectory_into_pdb.f90](https://github.com/Amirsaali/PDB_generation_from_LAMMPS_Trajectory_file/blob/main/Trajectory_into_pdb.f90) and type the name of your trajectory file in line 32.
