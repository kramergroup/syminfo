# syminfo
Utility for symmetry information of crystal structures

'syminfo' is a command-line utility to obtain symmetry information from crystal structures. Specifically, 'syminfo' can do the following:

- Find all space group operations
- Find symmetry equivalent atoms in the unit cell

*Warning*: This code is largely untested. Determination of symmetry operations is believed to be reliable. Other functionality should be thoroughly tested!

## Installation

### Compilation

Compiling the code should be as simple as

```bash
make
```

The code is completely self-contained and does not depend on libraries. The only requirement is a Fortran 90 compiler. The Intel compiler (ifort) has been used for testing.

### Installation

Compilation produces a self-contained executable in the main folder. Copy 'syminfo' to a location on your path (e.g., $HOME/bin). That should be it.

## Usage

syminfo reads crystal structure information from `stdin` and prints symmetry information to `stdout`. The simplest call sequence
```bash
syminfo < my_structure.str
```
will read the crystal structure defined in `my_structure.str` and output all symmetry operations.

### Input format

The expected format is similar to the POSCAR file format used by [VASP](http://www.vasp.at) (the Vienna Ab Initio Simulation Package). The example below describes Iron, a FCC crystal:

```
3.0000000E+00   0.0000000E+00   0.0000000E+00
0.0000000E+00   3.0000000E+00   0.0000000E+00
0.0000000E+00   0.0000000E+00   3.0000000E+00
4
0.0000000E+00   0.0000000E+00   0.0000000E+00 Fe 0.00
0.5000000E+00   0.5000000E+00   0.0000000E+00 Fe 0.00
0.0000000E+00   0.5000000E+00   0.5000000E+00 Fe 0.00
0.5000000E+00   0.0000000E+00   0.5000000E+00 Fe 0.00
```

The input file consists of three sections:
1. The first three lines give the lattice vectors a,b, and c in the cartesian reference frame
2. A single line contained the **total** number of atoms. Note: that is different from the POSCAR format, which expects the number of atoms per type.
3. Specification of atomic coordinates. Each row defines one atom with the first three numbers specifying the position in direct coordinates (i.e., as fractions of the lattice vectors), followed by the element and the nominal charge (or spin).

### Output formats

#### Space group operations

All space group operations are returned as default. The outout for an FCC structure
starts like:

```
48
  0.000000000E+00   0.000000000E+00  -0.100000000E+01
  0.000000000E+00   0.000000000E+00  -0.100000000E+01
  0.000000000E+00  -0.100000000E+01   0.000000000E+00
--------------------------------------------------------
  0.000000000E+00   0.000000000E+00   0.000000000E+00

  0.000000000E+00   0.000000000E+00  -0.100000000E+01
  0.000000000E+00   0.000000000E+00  -0.100000000E+01
  0.000000000E+00  -0.100000000E+01   0.000000000E+00
--------------------------------------------------------
  0.500000000E+01   0.500000000E+01   0.000000000E+00

...
```

The first line prints the number of space group operations. Followed by a block of
symmetry operations **W**. Each symetry operation

![equation](http://latex.codecogs.com/gif.latex?y={\bf R}\cdot x + {\bf t})

consists of a 3x3 rotation matrix **R** and a translation **t** below a line. The translation vector is printed in direct coordinates.
