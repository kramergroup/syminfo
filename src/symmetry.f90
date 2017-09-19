MODULE SYMMETRY

  USE DATATYPES
  USE IODEF

  IMPLICIT NONE

  integer, parameter :: MAX_NUM_PG = 200
  double precision :: ZERO_TOLERANCE = 1E-4

CONTAINS 

  SUBROUTINE FIND_POINTGROUP(cell,coords,pg)
    type(lattice), intent(in) :: cell
    type(coordinate), dimension(:), intent(in) :: coords
    type(pointgroup), dimension(:), allocatable, intent(out) :: pg

    type(pointgroup), dimension(MAX_NUM_PG) :: storage
    double precision :: maxlen
    double precision, dimension(3) :: norm
    integer :: maxm,maxn,maxl
    integer :: nops

    type(vector), dimension(3) :: u                                  ! temporary vectors
    integer, dimension(3) :: m,n,l
    logical :: finishedM, finishedN, finishedL
    type(lattice) :: invcell

    ! Counter and dummies
    integer :: i,j,k,ii
    double precision :: len

    WRITE(IO_DEBUG,*) "FIND_POINTGROUP----------------------------------------------------"

    ! Calculate inverse of cell
    invcell = .INVERSE. cell

    ! Norm of lattice vectors
    DO i=1,3
       norm(i) = .NORM. cell%e(i)
    END DO

    ! Find maximum lattice extend to search through
    ! Maping the largest onto the smallest lattice vector
    maxlen = .NORM. cell % e(1)

    DO i=2,3
       len = .NORM. cell % e(i)
       if ( len > maxlen ) maxlen = len
    END DO
    
    maxm = CEILING( maxlen / norm(1) )
    maxn = CEILING( maxlen / norm(2) )
    maxl = CEILING( maxlen / norm(3) )
    
    WRITE(IO_DEBUG,*) " Search domain:", maxm, maxn, maxl

    nops = 0
    
    ! Iterate over all possible combinations of three lattice vectors and 
    ! find the combinations that have same lengths as the lattice vectors
    ! and span the same area between pairs
    m = -maxm
    
    finishedM = .false.
    DO WHILE (.NOT. finishedM)
       ! Calculate candidate lattice vectors and check length
       u(1) = cell%e(1)*m(1) + cell%e(2)*m(2) + cell%e(3)*m(3)
       IF ( ABS(.NORM. u(1) -norm(1))>ZERO_TOLERANCE ) GOTO 773
       n = -maxn
       finishedN = .false.

       DO WHILE (.NOT. finishedN)
          u(2) = cell%e(1)*n(1) + cell%e(2)*n(2) + cell%e(3)*n(3)
          IF ( ABS(.NORM. u(2) -norm(2))>ZERO_TOLERANCE ) GOTO 772

          l = -maxl
          finishedL = .false.

          DO WHILE (.NOT. finishedL)
             u(3) = cell%e(1)*l(1) + cell%e(2)*l(2) + cell%e(3)*l(3)
             IF ( ABS(.NORM. u(3) -norm(3))>ZERO_TOLERANCE ) GOTO 771
             
             ! Check area of faces (implied correct angles)
             IF ( abs(cell%e(2)*cell%e(3) - u(2)*u(3)) < ZERO_TOLERANCE .AND. &
                  abs(cell%e(1)*cell%e(2) - u(1)*u(2)) < ZERO_TOLERANCE .AND. &
                  abs(cell%e(1)*cell%e(3) - u(1)*u(3)) < ZERO_TOLERANCE ) THEN
                
                nops = nops+1
                DO ii=1,3
                   storage(nops)%data(ii,1) = m(ii)
                   storage(nops)%data(ii,2) = n(ii)
                   storage(nops)%data(ii,3) = l(ii)
                END DO
             END IF

771          CONTINUE
             CALL INCREMENT(l,-maxl,maxl,finishedL) 
          END DO
772       CONTINUE
          CALL INCREMENT(n,-maxn,maxn,finishedN)
       END DO
773    CONTINUE
       CALL INCREMENT(m,-maxm,maxm,finishedM)

       ! Transfer result 
       IF ( ALLOCATED(pg) ) THEN
          DEALLOCATE(pg)
       END IF
       ALLOCATE( pg(nops) )
       DO i=1,nops
          pg(i) = storage(i)
       END DO
    END DO
    WRITE(IO_DEBUG,*) " Number of point-group operations: ", nops
    WRITE(IO_DEBUG,*) "-------------------------------------------------------------------"

  END SUBROUTINE FIND_POINTGROUP
  
  SUBROUTINE FIND_SPACEGROUP(cell,coords,sg)
    type(lattice), intent(in) :: cell
    type(coordinate), intent(in), dimension(:) :: coords
    type(spacegroup), allocatable, intent(out), dimension(:) :: sg

    type(pointgroup), allocatable, dimension(:) :: pg
    type(coordinate), allocatable, dimension(:) :: minority
    type(vector), allocatable, dimension(:) :: translations
    
    integer :: nops,i,j,k
    integer :: ntrans
    logical :: occ
    type(coordinate) :: coord
    type(vector) :: point
    type(lattice) :: invcell
    type(spacegroup), allocatable, dimension(:) :: storage

    invcell = .INVERSE. cell
    
    WRITE(IO_DEBUG,*) " Inverse cell:"
    DO i=1,3
       WRITE(IO_DEBUG,*) invcell%e(i)
    END DO

    CALL FIND_POINTGROUP(cell,coords,pg)
    CALL FIND_MINORITY_COORDINATES(coords, minority)
        
    ! Build vector of possible translations
    ALLOCATE( translations(size(minority)*size(pg)) )
    ntrans = 0
    DO i=1,size(minority)
       DO k=1,size(pg)
          point = pg(k) * (minority(i) % tau) - (minority(1) % tau)
          CALL FLIPINCELL(point)
          DO j=1,ntrans
             IF ( (.NORM. (point-translations(j))) < ZERO_TOLERANCE) GOTO 700
          END DO
          ntrans = ntrans+1
       
          !translations(ntrans) = cell * point
          translations(ntrans) = point
700       CONTINUE
      END DO
    END DO

    WRITE(IO_DEBUG,*) " Translations:"
    DO i=1,ntrans
       WRITE(IO_DEBUG,"(3(E19.7))") translations(i) 
    END DO

    ! Check all permutations of point-groups and translations
    nops = size(pg) * ntrans
    ALLOCATE( storage(nops) )
    
    nops = 0
    DO i=1,size(pg)
       DO j=1,ntrans
          DO k=1,size(coords)
             point = coords(k) % tau
             WRITE(IO_TRACE,"(A2,A,3(I3),A,3(E19.7))") "P:", "(",i,j,k,")",point

!             point = pg(i) * point + invcell * translations(j)
             point = pg(i) * point + translations(j)

             WRITE(IO_TRACE,"(A2,A,3(I3),A,3(E19.7))") "S:", "(",i,j,k,")",point
             CALL FIND_OCCUPANT(cell, coords, point, occ, coord)
             IF ((.NOT. occ) .OR. (.NOT. MATCHING_OCCUPANTS(coords(k),coord))) GOTO 800
          END DO
          ! Only reach this point if all coords match
          ! Found a valid space-group operation
          nops = nops + 1
          storage(nops)%pointgroup = pg(i)
          storage(nops)%translation = translations(j)
          WRITE(IO_DEBUG,"(A,A,2(I3),A)") "SPACEGROUP OPERATION (PG,TRANS):", "(",i,j,")" 
800       CONTINUE ! next translation
       END DO
    END DO

    ALLOCATE( sg(nops) )
    DO i=1,nops
       sg(i) = storage(i)
    END DO

    WRITE(IO_DEBUG,*) "FIND_SPACEGROUP----------------------------------------------------"
    WRITE(IO_DEBUG,"(A40I6)") " Number of considered centres: ", ntrans
    WRITE(IO_DEBUG,"(A40I6)") " Number of operations: ", nops
    WRITE(IO_TRACE,"(A40I6)") " Number of coordinates: ", size(coords)
    WRITE(IO_DEBUG,*)
    DO i=1,nops
       CALL PRINT_SPACEGROUP(IO_DEBUG,sg(i))
       WRITE (IO_DEBUG,*) "  TYPE: ", .TYPE. sg(i), "            ORDER: ", .ORDER. sg(i)
       WRITE (IO_DEBUG,*)
    END DO
    WRITE(IO_DEBUG,*) "-------------------------------------------------------------------"

  END SUBROUTINE FIND_SPACEGROUP


  ! Find a primitive cell 
  type(lattice) FUNCTION PRIMITIVE_CELL(cell,coords,minority) result (red)
    type(lattice), intent(in) :: cell                      ! Lattice vectors 
    type(coordinate), intent(in), dimension(:) :: coords   ! Atomic positions in direct coordinates
    type(coordinate), allocatable, dimension(:),optional, intent(inout) :: minority

    type(vector), allocatable, dimension(:) :: translations

    logical :: assigned, occ
    type(vector) :: point
    type(lattice) :: invcell
    integer :: i,j,k,l
    type(coordinate) :: coord
    double precision m

    red = cell
    invcell = .INVERSE. cell
    
    IF ( .NOT. present(minority) ) THEN
       CALL FIND_MINORITY_COORDINATES(coords, minority)
    END IF

    ! Build vector of possible translations
    ALLOCATE( translations(size(minority)) )
    DO i=1,size(translations)
       translations(i) = minority(i) % tau - minority(1) % tau
       translations(i) = cell * translations(i)
    END DO

    DO i=1,size(translations)
       
       ! Apply each translation and check if this is still the same structure
       DO j=1,size(coords)
          point = cell * (coords(j) % tau) + translations(i)
          point = invcell * point
          CALL FIND_OCCUPANT(cell, coords, point, occ, coord)
          IF ((.NOT. occ) .OR. (.NOT. MATCHING_OCCUPANTS(coords(j),coord))) GOTO 810
       END DO
       ! Only reach this point if translation is
       ! lattice translation
       DO k=1,3
          IF (((.NORM. translations(i)+ZERO_TOLERANCE)-(.NORM. red%e(k))) < 0.0) THEN
             ! Smaller vector found, check for co-planarity with others
             point = red%e(k)           ! Save old vector in case of co-planarity
             red%e(k) = translations(i)
             m = (red%e(3)-red%e(1)) .DOT. ((red%e(2)-red%e(1)).CROSS.red%e(3))

             IF (ABS(m) < ZERO_TOLERANCE) THEN
                ! Co-planar - restore old vector and try next
                red%e(k) = point
             ELSE
                ! Found new - non-coplanar lattice vector, try next translation
                GOTO 810
             END IF
          END IF
       END DO
810    CONTINUE
    END DO
    
    ! Ensure right-handedness
    m  = (red%e(1).CROSS.red%e(2)).DOT.red%e(3)
    IF ( m < 0.0 ) THEN
       point = red%e(1)
       red%e(1) = red%e(2)
       red%e(2) = point
    END IF
    
  END FUNCTION PRIMITIVE_CELL

  ! Calculates a matrix of symmetry-equivalent atoms
  SUBROUTINE SYMMETRY_MATRIX(cell, coords, sg, matrix)
    type(lattice), intent(in) :: cell
    type(coordinate), dimension(:), intent(in) :: coords
    type(spacegroup), dimension(:), intent(in) :: sg
    logical, dimension(:,:), allocatable, intent(out) :: matrix
    
    ! Counter and dummy variables
    integer :: i,j,k
    type(coordinate) :: dummy
    type(vector) :: point
    type(lattice) :: invcell
    integer :: idx
    logical :: occ

    invcell = .INVERSE. cell
    
    ALLOCATE( matrix(size(coords),size(coords)) )
    
    matrix = .false.

    DO i=1,size(sg) 
       DO j=1,size(coords)
          point = coords(j)%tau
!          point = sg(i)%pointgroup * point
!          point = invcell * ((cell*point) + sg(i)%translation)
          point = sg(i)%pointgroup * point + sg(i)%translation

          CALL FIND_OCCUPANT(cell,coords,point,occ,dummy,idx)
          IF ( occ ) THEN
             matrix(j,idx) = .true.
          END IF

       END DO ! loop over all coordinates

    END DO ! loop over spacewgroup operations

  END SUBROUTINE SYMMETRY_MATRIX
    


  SUBROUTINE REDUCED_COORDINATES(cell, coords, sg, redcoords)
    type(lattice), intent(in) :: cell
    type(spacegroup), dimension(:), intent(in) :: sg
    type(coordinate), dimension(:), intent(in) :: coords
    type(coordinate), dimension(:), intent(out), allocatable :: redcoords
    
    integer, dimension(size(coords)) :: idx
    integer :: i,j,k,ii,nred
    logical, dimension(:,:), allocatable :: symmtx

    CALL SYMMETRY_MATRIX(cell,coords,sg,symmtx)
    
    nred = 1
    idx(1) = 1
    DO i=2,size(coords)
       IF ( ALL( .NOT. symmtx(i,1:i-1) ) ) THEN
          nred = nred + 1
          idx(nred) = i
       END IF
    END DO

!    type(vector) :: point
!    type(lattice) :: invcell
!    type(coordinate) :: coord
!    logical :: occ
!
!    invcell = .INVERSE. cell
!    idx = 1
!    nred = 0
!    DO i=1,size(coords)
!      DO j=1,nred
!          point = coords(i)%tau
!          DO k=1,size(sg)
!             point = sg(k)%pointgroup * point
!             point = invcell * ((cell*point) + sg(k)%translation)
!             CALL FIND_OCCUPANT(cell,coords(idx),point,occ,coord)
!             IF ( occ ) GOTO 820
!          END DO
!       END DO
!       nred = nred+1
!       idx(nred) = i
!820    CONTINUE
!    END DO

    ALLOCATE( redcoords(nred) )
    redcoords = coords(idx)
    
    WRITE(IO_DEBUG,*) "REDUCED_COORDINATES---------------------------------------------------"
    WRITE(IO_DEBUG,"(A50I4)") "Number of coordinates: ", nred
    DO i=1,nred
       WRITE(IO_DEBUG,"(' ',E16.7,E16.7,E16.7,(A5),*(F7.2))") (redcoords(i)%tau%data(ii), ii=1,3), &
                                                              "  " //redcoords(i)%symbol, redcoords(i)%q
    END DO
    WRITE(IO_DEBUG,*) "----------------------------------------------------------------------"

  END SUBROUTINE REDUCED_COORDINATES

  ! Apply a spacegroup to a set of (primitive) coordinates
  SUBROUTINE APPLY_SPACEGROUP(cell, coords, sg, full)
    type(lattice), intent(in) :: cell
    type(coordinate), intent(in), dimension(:) :: coords
    type(spacegroup), intent(in), dimension(:) :: sg
    type(coordinate), intent(out), allocatable, dimension(:) :: full

    integer :: nmax, n
    type(coordinate), allocatable, dimension(:) :: storage
    integer :: i,j,k
    type(vector) :: point
    type(lattice) :: invcell
    logical :: occ
    type(coordinate) coord

    invcell = .INVERSE. cell

    nmax = size(sg) * size(coords)
    ALLOCATE( storage(nmax) )
    
    storage(1) = coords(1)
    n = 1
    DO i=1,size(coords)
       DO j=1,size(sg)
          point = coords(i) % tau
          point = sg(j)%pointgroup*point + sg(j)%translation
          
          CALL FIND_OCCUPANT(cell, storage(1:n), point, occ, coord)
          IF ( .NOT. occ ) THEN
             n = n+1
             storage(n) = coords(i)
             storage(n)%tau = point
          END IF
       END DO
    END DO

    ALLOCATE( full(n) )

    DO i=1,n
       full(i) = storage(i)
    END DO

    DEALLOCATE( storage )
 
 END SUBROUTINE APPLY_SPACEGROUP

  ! Check if a particular point in space is occupied with an atom
  ! and return the corresponding coordinate
  SUBROUTINE FIND_OCCUPANT(cell, coords, point, occ, coord, idx)
    type(lattice), intent(in) :: cell                      ! Lattice vectors
    type(coordinate), intent(in), dimension(:) :: coords   ! Atomic positions in direct coordinates
    type(vector), intent(in) :: point                      ! Point in direct coordinates
    logical, intent(out) :: occ                            ! True if point matches a coordinate
    type(coordinate), intent(out), optional :: coord       ! Contains the matched coordinate
    integer, intent(out), optional :: idx                  ! Index of the coordinate
    
    integer :: i,j
    type(vector) :: p
    double precision :: delta
    type(coordinate) :: dummy

    occ = .false.
    if ( present(idx) ) idx = -1
    if ( present(coord) ) coord = dummy
    coord%symbol = "XX"

    p = point
    CALL FLIPINCELL(p)
    WRITE(IO_TRACE,*) p
    DO i=1,size(coords)
       dummy = coords(i)
       CALL FLIPINCELL(dummy%tau)
       delta = .NORM. (dummy%tau - p)
       WRITE(IO_TRACE,*) dummy%tau, "(", delta, ")"
       
       IF ( abs(abs(delta)/ZERO_TOLERANCE-1.0) <= 1E-1 ) &
            WRITE(IO_ERROR,*) "Warning : Difference close to precision. Operations might have been lost"

       IF ( delta < ZERO_TOLERANCE ) THEN
          IF ( present(coord) ) coord = coords(i)
          occ = .true.
          IF ( present(idx) ) idx = i
          RETURN
       END IF
    END DO
    
  END SUBROUTINE FIND_OCCUPANT
  

  ! Convert spacegroup from direct to cartesian representation
  type(spacegroup) FUNCTION TO_CARTESIAN_FRAME(sg,cell) result (cart)
    type(spacegroup), intent(in) :: sg
    type(lattice), intent(in) :: cell

    cart = sg
    
    ! Convert transformation matrix (Pc = C*P*C^-1)
    cart%pointgroup = cell * sg%pointgroup * (.INVERSE. cell)
    

    ! Convert translation vector (tc = C*t)
    cart%translation = cell * sg%translation

  END FUNCTION TO_CARTESIAN_FRAME



  ! Flip direct coordinates back into the cell
  SUBROUTINE FLIPINCELL(coord)
    type(vector), intent(inout) :: coord                  ! Position in direct coordinates
    
    double precision :: dummy
    integer :: i
    
    DO i=1,3
       dummy = mod(coord%data(i),1.0)
       IF( abs(1.0-dummy) < ZERO_TOLERANCE) dummy = 0.0
       IF ( dummy < -ZERO_TOLERANCE ) dummy = dummy + 1.0
       coord%data(i) = dummy
    END DO
  END SUBROUTINE FLIPINCELL


  ! Recursively find the minority species and corresponding coordinates
  RECURSIVE SUBROUTINE FIND_MINORITY_COORDINATES(coords, minority, rec)
    type(coordinate), intent(in), dimension(:) :: coords
    type(coordinate), allocatable, intent(inout), dimension(:) :: minority
    logical, intent(in), optional :: rec ! internal to indicate that this is a recursive call

    integer :: ncoords, i, j, k
    integer, allocatable, dimension(:) :: minidx, otheridx

    IF (.NOT. PRESENT(rec)) THEN
       WRITE(IO_DEBUG,*) "FIND_MINORITY_COORDINATES------------------------------------------"
    END IF
    
    ALLOCATE( minidx(size(coords)) ) 

    ! Count coordinates with same occupancy as first
    ncoords = 0
    minidx = 0
    DO i=1,size(coords)
       IF ( MATCHING_OCCUPANTS(coords(i),coords(1)) ) THEN
          ncoords = ncoords + 1
          minidx(ncoords) = i
       END IF
    END DO

    WRITE(IO_DEBUG,*) " Species: ", coords(1) % symbol, coords(1) % q
    WRITE(IO_DEBUG,*) " Number of coordinates: ", ncoords
    WRITE(IO_DEBUG,*) " Coordinates: "

    DO i=1,ncoords
       WRITE(IO_DEBUG,"(3(E18.7))") (coords(minidx(i)) % tau % data(j), j=1,3) 
    END DO
    WRITE(IO_DEBUG,*)

    IF ( (.NOT. ALLOCATED(minority)) .OR. (size(minority) > ncoords) ) THEN
       IF ( ALLOCATED(minority) ) DEALLOCATE(minority)
       ALLOCATE(minority(ncoords))

       minority = coords(minidx)
    END IF

    ! Nothing left, all coordinates assigned
    IF ( size(coords) == ncoords ) RETURN 

    ! Build index of remaining coordinates
    ALLOCATE( otheridx(size(coords)-ncoords) ) 
    j=0
    otheridx = 0
    DO i=1,size(coords)
       IF ( ANY( minidx == i ) ) GOTO 889
       j=j+1
       otheridx(j) = i
889    CONTINUE
    END DO
    WRITE(IO_TRACE,*) " Current indices  : ", minidx
    WRITE(IO_TRACE,*) " Remaining indices: ", otheridx

    ! Recursively check remaining coords
    CALL FIND_MINORITY_COORDINATES(coords(otheridx),minority,.true.)

    IF (.NOT. PRESENT(rec)) THEN
       WRITE(IO_DEBUG,*) " Minority species: ", minority(1) % symbol, minority(1) % q
       WRITE(IO_DEBUG,*) " Number of minority coordinates: ", size(minority)
       WRITE(IO_DEBUG,*) "------------------------------------------------------------------"
    END IF

  END SUBROUTINE FIND_MINORITY_COORDINATES


  ! Returns true if the two coordinates have the same occupant, meaning
  ! the same symbol, number of valence states, and valence states
  logical FUNCTION MATCHING_OCCUPANTS(a,b) result (same)
    type(coordinate), intent(in) :: a,b
    
    integer :: j,k
    
    same = .false.

    WRITE(IO_TRACE,*) " OCCUPANTS: ", a % symbol, b % symbol

    IF ( a % symbol /= b % symbol ) RETURN
    IF ( size(a % q) /= size(b % q)) RETURN
       
    DO j=1,size(a % q)
       DO k=1,size(b % q)
          IF ( (a % q(k) == b % q(j)) ) GOTO 900
       END DO
       WRITE(IO_DEBUG,*) " VALENCES DID NOT MATCH"
       RETURN    ! No matching valence
900    CONTINUE  ! Found a matching valence
    END DO
    
901 CONTINUE ! Only reach this point if occupants match
    same = .true.
    RETURN

  END FUNCTION MATCHING_OCCUPANTS

  ! Increments an array of integers in the range [min,max] to
  ! generate all permutations
  SUBROUTINE INCREMENT(idx, min, max,finished)
    integer, dimension(3), intent(inout) :: idx
    integer, intent(in) :: min, max
    logical, intent(inout) :: finished
    

    integer :: i
    ! Advance to first element that can be incremented
    i = 1
    DO WHILE(idx(i) == max) 
       i = i+1
       if ( i > size(idx) ) GOTO 999
    END DO
    
    idx(i) = idx(i) + 1
    DO WHILE (i > 1)
       i = i-1
       idx(i) = min
    END DO
    finished = .false.
    RETURN

999 CONTINUE
    finished = .true.
    RETURN

  END SUBROUTINE INCREMENT

  SUBROUTINE PRINT_SPACEGROUP(out, sg, cell)
    integer, intent(in) :: out
    type(spacegroup), intent(in) :: sg
    type(lattice), optional :: cell

    integer :: i,j
    type(spacegroup) :: tmpsg

    IF ( present(cell) ) THEN
       tmpsg = TO_CARTESIAN_FRAME(sg,cell)
    ELSE
       tmpsg = sg
    END IF

    DO i=1,3
       WRITE(out,"(3(E18.9))") (tmpsg%pointgroup%data(i,j), j=1,3)
    END DO
    WRITE(out,*) "--------------------------------------------------------"
    WRITE(out,"(3(E18.9))") (tmpsg%translation%data(i), i=1,3)
    
  END SUBROUTINE PRINT_SPACEGROUP

  SUBROUTINE READ_SPACEGROUP(in, sg, cell)
    integer, intent(in) :: in
    type(spacegroup), intent(inout), allocatable, dimension(:) :: sg
    type(lattice), intent(in), optional :: cell
   
    integer nops, i,j,k
    type(lattice) :: invcell
    
    IF ( present(cell) ) invcell = .INVERSE. cell

    READ(in, "(I)") nops
    
    IF ( nops < 1 ) RETURN

    IF ( allocated(sg) ) deallocate(sg)
    ALLOCATE( sg(nops) )
    
    DO i=1,nops
       DO j=1,3
          READ(in,"(3(E18.9))") (sg(i)%pointgroup%data(j,k), k=1,3)
       END DO
       READ(in,*) ! skip the divider line
       READ(in,"(3(E18.9))") (sg(i)%translation%data(j), j=1,3)
       READ(in,*) ! skip separator line
       
       ! If cell is present, input is expected to be in cartesian frame,
       ! convert to internal frame
       IF ( present(cell) ) THEN
          sg(i)%pointgroup = invcell * sg(i)%pointgroup * cell
          sg(i)%translation = invcell * sg(i)%translation
       END IF

    END DO
    
  END SUBROUTINE READ_SPACEGROUP

  SUBROUTINE WRITE_SYMMETRYMATRIX(out, matrix)
    integer, intent(in) :: out
    logical, dimension(:,:), intent(in) :: matrix

    integer i,j

    WRITE(out,"(I3)") SIZE(matrix,1)
    
    DO i=1,SIZE(matrix,2)
       WRITE(out,"(*(L2))") (matrix(i,j),j=1,SIZE(matrix,1))
    END DO
  END SUBROUTINE WRITE_SYMMETRYMATRIX

  SUBROUTINE READ_SYMMETRYMATRIX(in, matrix)
    integer, intent(in) :: in
    logical, dimension(:,:), allocatable, intent(out) :: matrix

    integer ncoord,i

    READ(in, *) ncoord

    ALLOCATE( matrix(ncoord,ncoord) )

    DO i=1,ncoord
       READ(in,*) matrix(i,:)
    END DO

  END SUBROUTINE READ_SYMMETRYMATRIX

END MODULE SYMMETRY
