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

A number of command line parameters are available for other behaviour.

```
syminfo - Determine pointgroup, spacegroup, and reduced cell information
           usage: syminfo [-p] [-r] [-v] [-h] [-2d] [-a spacegroup.dat] [-m spacegroup.dat] < structure
           Options:
              -h       Display help message
              -v       Display version
              -p       Output point-group
              -r       Output coordinates of reduced basis
              -a       Applies the spacegroup in file spacegroup.dat to the
                       supplied structure
              -m       Matrix of symmetry equivalent coordinates under
                       the given spacegroup
              -2d      2D mode; ignores symmetry in z-direction
                                apart from mirror-plane in x-y
              -cart    Use cartesian reference frame. Usually, space group
                       operations operate on direct coordinates according to
                       p' = P*p + t. With this switch, the space group is
                       referencing the  cartesian frame. This is achieved via
                       the transformations P->C*P*C^-1 and t->C*t.
              -cartin  Use cartesian frame for input
              -cartout Use cartesian frame for outout
```

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
192
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

The first line prints the number of space group operations. There are 48x4=192 operations, because each of the 48 point group operations can be combined with a translation between atoms in the cubic unit cell (see [definition of space group 227](http://img.chem.ucl.ac.uk/sgp/large/227az1.htm)). The rest of the output consists of blocks of
symmetry operations **W**. Each symmetry operation

![equation](http://latex.codecogs.com/gif.latex?x%60%3D%7B%5Cbf%20R%7D%5Ccdot%20x%20%2B%20%7B%5Cbf%20t%7D)

consists of a 3x3 rotation matrix **R** and a translation **t** below a line. The translation vector is printed in direct coordinates. An empty line separates operations.

#### Symmetry matrix (-m)

Using the `-m <spacegroup_operators.sym>` switch will find all symmetry equivalent atoms under the set of space group operations defined in the file `spacegroup_operators.sym` (see above for the format). The returned square matrix has the dimensions of the number of atom, which is printed in the first line, and contains a `T` if two atoms are equivalent. Otherwise, a `F` designates symmetry distinct pairs. Atom order is the same as in the input file.

For example,

```bash
syminfo -m example/fcc.sym < example/fcc.str
```

will produce a 4x4 matrix, indicating that all atoms are symmetry equivalent.

```
  4
 T T T T
 T T T T
 T T T T
 T T T T
```

#### Reduced basis (-r)

The `-r` command-line argument produces a structure file (using the input format) that contains only symmetry inequivalent coordinates followed by all symmetry operations (see output format above).

```bash
syminfo -r < example/fcc.str
```

will produce

```
0.3000000E+01     0.0000000E+00     0.0000000E+00
0.0000000E+00     0.3000000E+01     0.0000000E+00
0.0000000E+00     0.0000000E+00     0.3000000E+01
1
0.0000000E+00   0.0000000E+00   0.0000000E+00  Fe    0.00
192
0.000000000E+00   0.000000000E+00  -0.100000000E+01
0.000000000E+00  -0.100000000E+01   0.000000000E+00
-0.100000000E+01   0.000000000E+00   0.000000000E+00
--------------------------------------------------------
0.000000000E+00   0.000000000E+00   0.000000000E+00

0.000000000E+00   0.000000000E+00  -0.100000000E+01
0.000000000E+00  -0.100000000E+01   0.000000000E+00
-0.100000000E+01   0.000000000E+00   0.000000000E+00
--------------------------------------------------------
0.000000000E+00   0.500000000E+00   0.500000000E+00

...
```

#### Apply space group operations (-a)

This is the *inverse* of the `-r` switch and finds all atomic positions within the unit cell from a list of symmetry inequivalent positions.
