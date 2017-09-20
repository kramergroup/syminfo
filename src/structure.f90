MODULE STRUCTURE

  USE DATATYPES
  USE IODEF

  IMPLICIT NONE

CONTAINS

  ! input of structures for ewald summation from file
  ! It assumes that the unit has already been opened to read from
  SUBROUTINE READSTRUCTURE(iin,cell,coords,nion,cartframe)

    ! Arguments
    integer, intent(in) :: iin                                             ! Index of fortran unit to read from
 
    type(lattice), intent(inout) :: cell                                    ! Lattice vectors
    type (coordinate), intent(out), allocatable, dimension(:) :: coords    ! Coordinates of atoms in direct coordinates
    integer, intent(out) :: nion                                           ! Number of ions
    logical :: cartframe                                                   ! Input file has cartesian coordinates

    ! Parameters
    integer, parameter :: nmaxval = 10                                     ! Number of maximum valences

    ! Temporary and dummy variables
    integer :: i,n,ii,jj
    character*255:: line                                                   ! the line buffer
    double precision,dimension(nmaxval) :: tval
    type(lattice) :: invcell                                               ! inverse lattice

    ! Read lattice vectors
    DO i=1,3
       READ(iin, *) (cell % e(i) % data(ii), ii=1,3)
    END DO
    
    invcell = .INVERSE. cell

    READ(iin, *) nion

    IF ( nion < 1 ) THEN
       WRITE(1, *) "No atomic coordinates given"
       STOP
    END IF

    ALLOCATE( coords(nion) )

    DO i=1,nion
       READ(iin, '(A)', end=999) line
       READ(line, *, end=888) (coords(i)%tau%data(ii), ii=1,3), coords(i)%symbol, (tval(ii), ii = 1, nmaxval)
777 CONTINUE
       ! This point is only reached if more valences are specified than we anticipate
       WRITE(1,*) "More than ", nmaxval, " valences not supported!"
       STOP
888 CONTINUE
       ALLOCATE( coords(i)%q(ii-1) )
       DO n=1,ii-1
          coords(i)%q(n) = tval(n)
       END DO

       ! Convert to direct coordiantes if nessesary
       IF ( cartframe ) coords(i)%tau = invcell * coords(i)%tau
       
    END DO

999 CONTINUE

    WRITE(IO_DEBUG,*) "READSTRUCTURE------------------------------------------------------"
    WRITE(IO_DEBUG,"(A,E16.7,E16.7,E16.7)") "  Lattice vector A: " , (cell%e(1)%data(ii),ii=1, 3)
    WRITE(IO_DEBUG,"(A,E16.7,E16.7,E16.7)") "  Lattice vector B: " , (cell%e(2)%data(ii),ii=1, 3)
    WRITE(IO_DEBUG,"(A,E16.7,E16.7,E16.7)") "  Lattice vector C: " , (cell%e(3)%data(ii),ii=1, 3)
    WRITE(IO_DEBUG,*)
    WRITE(IO_DEBUG,"(A,I4,A)") " Coordinates : (" , nion, " total)"
    WRITE(IO_DEBUG,*) "         X              Y               Z        Elem  Valences"
    DO i=1,nion
       WRITE(IO_DEBUG,"(' ',E16.7,E16.7,E16.7,(A5),*(F7.2))") (coords(i)%tau%data(ii), ii=1,3), "  " //coords(i)%symbol, coords(i)%q
    END DO
    WRITE(IO_DEBUG,*) "-------------------------------------------------------------------"

    RETURN

  END SUBROUTINE READSTRUCTURE
  
  SUBROUTINE WRITELATTICE(out,cell)

    ! Parameters
    integer, intent(in) :: out                                              ! Fortran unit for output    
    type(lattice), intent(in) :: cell                                       ! Lattice vectors

    ! Counter
    integer :: i,ii

    DO i=1,3
       WRITE(out,"(E18.7,E18.7,E18.7)") (cell%e(i)%data(ii),ii=1, 3)
    END DO

  END SUBROUTINE WRITELATTICE  

  SUBROUTINE WRITESTRUCTURE(out,cell,coords,cartframe)
    
    ! Parameters
    integer, intent(in) :: out                                              ! Fortran unit for output
    type(lattice), intent(in) :: cell                                       ! Lattice vectors
    type (coordinate), intent(in), allocatable, dimension(:) :: coords      ! Coordinates of atoms in direct coordinates
    logical, intent(in) :: cartframe                                        ! Print coordinates in cartesian frame if true

    ! Counter and local variables
    integer :: i,ii
    character(len=8) :: i_char
    type(vector) :: tau

    CALL WRITELATTICE(out,cell)
    WRITE(i_char, '(I8)') size(coords)
    WRITE(out,*) adjustl(i_char) 

    DO i=1,size(coords)
       IF ( cartframe ) THEN
          tau = cell * coords(i)%tau
       ELSE
          tau = coords(i)%tau
       END IF

       WRITE(out,"(' ',E16.7,E16.7,E16.7,(A5),*(F7.2))") &
            (tau%data(ii), ii=1,3), "  " //coords(i)%symbol, coords(i)%q
    END DO

  END SUBROUTINE WRITESTRUCTURE

  SUBROUTINE ROUND_COORDINATES(coords, prec)
    type(coordinate), dimension(:), intent(inout) :: coords
    double precision, intent(in) :: prec

    integer i,j

    DO i=1,size(coords)
       DO j=1,3
          coords(i)%tau%data(j) = NINT(coords(i)%tau%data(j)/prec)*prec
       END DO
    END DO
  END SUBROUTINE ROUND_COORDINATES

END MODULE STRUCTURE
