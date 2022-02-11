PROGRAM SYMINFO

  USE DATATYPES
  USE STRUCTURE
  USE SYMMETRY
  USE IODEF

  !VERSION
  character(LEN=*), parameter :: VERSION = "1.0.0"

  ! Mode
  integer, parameter :: MODE_SG = 1 , MODE_PG = 2 , MODE_RED = 3, MODE_APPLY = 4, MODE_SYMMTX = 5, MODE_PRIM = 6, MODE_NEIGH = 7
  integer :: mode = MODE_SG

  ! 2D Mode (ignores z-direction apart from mirror plane in x-y)
  logical :: planar = .false.

  ! Maximum neighbor distance for connectivity analysis
  double precision :: max_neighbor_distance = 2.3

  ! Output mode (Direct or cartesian)
  ! Internally, the transformation operations are stored as relative to the direct coordinate frame
  ! This means that p' = P*p + t needs p in direct coordinates. This is the default mode for output as well
  ! If the cartesian reference frame is needed, the following conversion is applied
  !
  !     Pc = C*P*C^-1 and tc = C*t    with C being the unit cell matrix
  ! 
  ! If cartesian output is requested, (Pc,tc) will be plotted instead of (P,t)
  logical :: cartframein = .false.
  logical :: cartframeout = .false.

  ! Data
  type(lattice) :: cell,prim
  type(coordinate), dimension(:), allocatable :: coords, redcoords
  type(pointgroup), dimension(:), allocatable :: pg
  type(spacegroup), dimension(:), allocatable :: sg
  logical, dimension(:,:), allocatable :: symmtx
  integer :: nion
  character(len=100) :: sgfile      ! Used for externally defined space groups (MODE_APPLY and MODE_SYMMTX)
  character(len=20) :: dummy        ! A character buffer used for reading commandline parameters

  ! Counter and dummies
  integer :: i,j,k

  ! IO
  integer, parameter :: stdin = 5
  integer, parameter :: stdout = 6

  CALL INIT_IO()

  ! Process input arguments
  CALL PROCESS_ARGUMENTS()
  
  ! Read the structure from file
  CALL READSTRUCTURE(stdin,cell,coords,nion,cartframein)
  
  ! Read space group information 
  IF ( mode == MODE_APPLY .OR. mode == MODE_SYMMTX ) THEN
     OPEN(10,FILE=sgfile)
     if ( cartframein ) THEN
        CALL READ_SPACEGROUP(10,sg,cell)
     ELSE
        CALL READ_SPACEGROUP(10,sg)
     END if
  END IF

  IF ( mode == MODE_SG ) THEN
     CALL FIND_SPACEGROUP(cell,coords,sg)
     CALL PRINT_SPACEGROUPS()
  ELSE IF ( mode == MODE_PG ) THEN
     CALL FIND_POINTGROUP(cell,coords,pg)
     CALL PRINT_POINTGROUP
  ELSE IF ( mode == MODE_RED ) THEN
     CALL FIND_SPACEGROUP(cell,coords,sg)
     CALL REDUCED_COORDINATES(cell,coords,sg,redcoords)
     CALL ROUND_COORDINATES(redcoords, ZERO_TOLERANCE)
     CALL WRITESTRUCTURE(IO_OUT,cell,redcoords,cartframeout)
     CALL PRINT_SPACEGROUPS
  ELSE IF ( mode == MODE_APPLY ) THEN
     CALL APPLY_SPACEGROUP(cell, coords, sg, redcoords)
     CALL ROUND_COORDINATES(redcoords, ZERO_TOLERANCE)
     CALL WRITESTRUCTURE(IO_OUT, cell, redcoords,cartframeout)
  ELSE IF ( mode == MODE_SYMMTX ) THEN
     CALL SYMMETRY_MATRIX(cell, coords, sg, symmtx)
     CALL WRITE_SYMMETRYMATRIX(IO_OUT,symmtx)
  ELSE IF ( mode == MODE_PRIM ) THEN
     prim = PRIMITIVE_CELL(cell,coords)
     prim = prim * ( .INVERSE. cell )
     CALL WRITELATTICE(IO_OUT,prim)
  END IF

CONTAINS
  
  SUBROUTINE PROCESS_ARGUMENTS()
      
    character(len=100) :: arg
    character(len=100) :: dummy
    integer :: narg

      narg = IARGC()


      ! Process args with higher priority that might influence others
      DO i=1,narg
         CALL GETARG(i,arg)        
         
         IF (arg == "-cart") THEN
            cartframein = .true.
            cartframeout = .true.
         ELSE IF (arg == "-cartin") THEN
            cartframein = .true.
         ELSE IF (arg == "-cartout") THEN
            cartframeout = .true.
         ELSE IF (arg == "-eps") THEN
            CALL GETARG(i+1,dummy)
            READ(dummy,*) ZERO_TOLERANCE
         ELSE IF (arg == "-debug") THEN
            CALL ENABLE_DEBUG()
         ELSE IF (arg == "-trace") THEN
            CALL ENABLE_DEBUG()
            CALL ENABLE_TRACE()
         END IF
      END DO

      DO i=1,narg
         CALL GETARG(i,arg)

         IF (arg == "-h") THEN
            CALL DISPLAY_HELPMESSAGE()
            STOP
         ELSE IF (arg == "-v") THEN
            CALL DISPLAY_VERSION()
            STOP
         ELSE IF (arg == "-2d") THEN
            planar = .true.
         ELSE IF (arg == "-p") THEN
            mode = MODE_PRIM
         ELSE IF (arg == "-pg") THEN
            mode = MODE_PG
         ELSE IF (arg == "-r") THEN
            mode = MODE_RED
         ELSE IF (arg == "-a") THEN
            mode = MODE_APPLY
            CALL GETARG(i+1,sgfile)
         ELSE IF (arg == "-m") THEN
            mode = MODE_SYMMTX
            CALL GETARG(i+1,sgfile)
         ELSE IF (arg == "-n") THEN
            mode = MODE_NEIGH
            CALL GETARG(i+1,dummy)
            READ(dummy , *, err=110) max_neighbor_distance
         END IF
      END DO

      RETURN

      ! Handle argument format errors
110   CONTINUE  
      WRITE(IO_ERR,*) "Distance argument to -n has to be a number"
      ERROR STOP 1

  END SUBROUTINE PROCESS_ARGUMENTS
  
  SUBROUTINE PRINT_POINTGROUP()
    
    integer :: i,j

    DO i=1,size(pg)
       WRITE(IO_OUT,"(3(F10.4))") (pg(i)%data(1,j),j=1,3)
       WRITE(IO_OUT,"(3(F10.4))") (pg(i)%data(2,j),j=1,3)
       WRITE(IO_OUT,"(3(F10.4))") (pg(i)%data(3,j),j=1,3)
       WRITE(IO_OUT,*)
    END DO

  END SUBROUTINE PRINT_POINTGROUP

  SUBROUTINE PRINT_SPACEGROUPS()

    integer :: i
    character(len=8) :: i_char

    WRITE(i_char, '(I8)') size(sg)
    WRITE(IO_OUT,*) adjustl(i_char) 

    DO i=1,size(sg)
       IF ( cartframeout ) THEN 
          CALL PRINT_SPACEGROUP(IO_OUT,sg(i),cell)
       ELSE
          CALL PRINT_SPACEGROUP(IO_OUT,sg(i))
       END IF

       WRITE(IO_OUT,*)
    END DO

  END SUBROUTINE PRINT_SPACEGROUPS

  SUBROUTINE DISPLAY_HELPMESSAGE()

    WRITE(IO_OUT,*) "syminfo - Determine pointgroup, spacegroup, and reduced cell information"
    WRITE(IO_OUT,*) "          usage: syminfo [-p] [-r] [-v] [-h] [-2d] [-n dist] [-a spacegrp.dat] [-m spacegrp.dat] < structure"
    WRITE(IO_OUT,*) "          Options: "
    WRITE(IO_OUT,*) "             -h       Display help message"
    WRITE(IO_OUT,*) "             -v       Display version"
    WRITE(IO_OUT,*) "             -pg      Output point-group"
    WRITE(IO_OUT,*) "             -p       Output transformation matrix to primitive cell"
    WRITE(IO_OUT,*) "             -r       Output coordinates of reduced basis"
    WRITE(IO_OUT,*) "             -a       Applies the spacegroup in file spacegroup.dat to the"
    WRITE(IO_OUT,*) "                      supplied structure"
    WRITE(IO_OUT,*) "             -m       Matrix of symmetry equivalent coordinates under"
    WRITE(IO_OUT,*) "                      the given spacegroup"
    WRITE(IO_OUT,*) "             -2d      2D mode; ignores symmetry in z-direction"
    WRITE(IO_OUT,*) "                               apart from mirror-plane in x-y"
    WRITE(IO_OUT,*) "             -n       Perform neighborhood analysis including sites max. distance away"
    WRITE(IO_OUT,*) "             -cart    Use cartesian reference frame. Usually, space group"
    WRITE(IO_OUT,*) "                      operations operate on direct coordinates according to"
    WRITE(IO_OUT,*) "                      p' = P*p + t. With this switch, the space group is"
    WRITE(IO_OUT,*) "                      referencing the  cartesian frame. This is achieved via"
    WRITE(IO_OUT,*) "                      the transformations P->C*P*C^-1 and t->C*t."
    WRITE(IO_OUT,*) "             -cartin  Use cartesian frame for input"
    WRITE(IO_OUT,*) "             -cartout Use cartesian frame for outout"

  END SUBROUTINE DISPLAY_HELPMESSAGE
  
  SUBROUTINE DISPLAY_VERSION()
    
    WRITE(IO_OUT,*) "syminfo - Version " // VERSION 

  END SUBROUTINE DISPLAY_VERSION

END PROGRAM
