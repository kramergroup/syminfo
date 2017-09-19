MODULE DATATYPES

  USE IODEF
  
  IMPLICIT NONE

  TYPE :: VECTOR
     double precision :: data(3)
  END TYPE VECTOR

  TYPE :: COORDINATE
     type(vector) :: tau                    ! atomic coordinates in direct coordinates
     double precision, dimension(:), allocatable :: q   ! vector of nval valences
     integer :: nval                        ! number of valences
     character (len=5) :: symbol            ! atomic symbol
  END TYPE COORDINATE

  TYPE :: LATTICE
     type(vector), dimension(3) :: e
  END TYPE LATTICE

  TYPE :: POINTGROUP
     double precision :: data(3,3)
  END type POINTGROUP

  TYPE :: SPACEGROUP
     type(pointgroup) :: pointgroup
     type(vector) :: translation
  END type SPACEGROUP

  INTERFACE OPERATOR(+)
     module procedure vadd
  END INTERFACE OPERATOR(+)

  INTERFACE OPERATOR(-)
     module procedure vsub
  END INTERFACE OPERATOR(-)
  
  INTERFACE OPERATOR(*)
     module procedure vdot, pdot, sdot, udot,vidot, lpdot, pldot, ppdot
  END INTERFACE OPERATOR(*)

  INTERFACE OPERATOR(.DOT.)
     module procedure vdot, pdot, sdot, udot,vidot
  END INTERFACE OPERATOR(.DOT.)

  INTERFACE OPERATOR(.NORM.)
     module procedure vnorm
  END INTERFACE OPERATOR(.NORM.)

  INTERFACE OPERATOR(.INVERSE.)
     module procedure linv
  END INTERFACE OPERATOR(.INVERSE.)

  INTERFACE OPERATOR(.CROSS.)
     module procedure vcross
  END INTERFACE OPERATOR(.CROSS.)

  INTERFACE OPERATOR(.ORDER.)
     module procedure sgorder
  END INTERFACE OPERATOR(.ORDER.)

  INTERFACE OPERATOR(.TYPE.)
     module procedure sgtype
  END INTERFACE OPERATOR(.TYPE.)

CONTAINS
  
  ! VECTOR OPERATORS
  FUNCTION vadd(v1,v2) result (v3)
    IMPLICIT NONE
    type (vector), intent(in)  :: v1,v2
    type (vector)              :: v3
    v3%data = v1%data + v2%data ! Note the use of array operations
  END FUNCTION vadd

  FUNCTION vsub(v1,v2) result (v3) 
    type (vector), intent(in)  :: v1,v2
    type (vector)              :: v3
    v3%data = v1%data - v2%data ! Note the use of array operations
  END FUNCTION vsub

  FUNCTION vdot(v1,v2) result (dot)
    type (vector), intent(in) :: v1,v2
    double precision          :: dot
    dot = DOT_PRODUCT(v1 % data, v2 % data)
  END FUNCTION vdot

  FUNCTION vcross(a,b) result (c)
    type(vector), intent(in)  :: a,b
    type(vector)              :: c
    
    c%data(1) = a%data(2) * b%data(3) - a%data(3) * b%data(2)
    c%data(2) = a%data(3) * b%data(1) - a%data(1) * b%data(3)
    c%data(3) = a%data(1) * b%data(2) - a%data(2) * b%data(1)
  END FUNCTION vcross

  FUNCTION vidot(v1,i) result (v2)
    type(vector), intent(in) :: v1
    integer, intent(in)      :: i
    type(vector)             :: v2

    v2%data = v1%data * i
  END FUNCTION vidot
  
  FUNCTION pdot(pg,v1) result (v2)
    type(pointgroup), intent(in)  :: pg
    type(vector), intent(in)      :: v1
    type(vector)                  :: v2
    integer                       :: i
    
    DO i=1,3
       v2%data(i) = DOT_PRODUCT(pg%data(i,:),v1%data)
    END DO
  END FUNCTION pdot
  
  
  FUNCTION sdot(sg,v1) result (v2)
    type(spacegroup), intent(in)  :: sg
    type(vector), intent(in)      :: v1
    type(vector)                  :: v2
    
    v2 = pdot(sg % pointgroup, v1)
    v2 % data = (v2 % data) + (sg % translation % data)

  END FUNCTION sdot

  FUNCTION udot(cell,v1) result (v2)
    type(lattice), intent(in) :: cell
    type(vector), intent(in)   :: v1
    type(vector)               :: v2
    
    integer :: i
    
    DO i=1,3
       v2%data(i) = cell%e(1)%data(i)*v1%data(1) + cell%e(2)%data(i)*v1%data(2) + cell%e(3)%data(i)*v1%data(3)
    END DO
  END FUNCTION udot

  FUNCTION lpdot(cell,pg) result (pgout)
    type(lattice), intent(in)    :: cell
    type(pointgroup), intent(in) :: pg
    type(pointgroup)             :: pgout

    integer i,j,k
    DO i=1,3
       DO j=1,3
          pgout%data(i,j) = 0.00
          DO k=1,3
             pgout%data(i,j) = pgout%data(i,j) + cell%e(k)%data(i)*pg%data(k,j)
          END DO
          
       END DO
    END DO

  END FUNCTION lpdot

  FUNCTION pldot(pg,cell) result (pgout)
    type(lattice), intent(in)    :: cell
    type(pointgroup), intent(in) :: pg
    type(pointgroup)             :: pgout

    integer i,j,k
    DO i=1,3
       DO j=1,3
          pgout%data(i,j) = 0.00
          DO k=1,3
             pgout%data(i,j) = pgout%data(i,j) + pg%data(i,k) * cell%e(j)%data(k)
          END DO
          
       END DO
    END DO

  END FUNCTION pldot

  FUNCTION ppdot(pg1,pg2) result (pgout)
    type(pointgroup), intent(in) :: pg1,pg2
    type(pointgroup)             :: pgout

    integer i,j,k
    DO i=1,3
       DO j=1,3
          pgout%data(i,j) = 0.00
          DO k=1,3
             pgout%data(i,j) = pgout%data(i,j) + pg1%data(i,k)*pg2%data(k,j)
          END DO
          
       END DO
    END DO

  END FUNCTION ppdot



  
  FUNCTION vnorm(v1) result (norm)
    type(vector), intent(in) :: v1
    double precision         :: norm
    
    norm = SQRT(DOT_PRODUCT(v1%data,v1%data))
  END FUNCTION vnorm

  FUNCTION linv(lat) result (inv)
    type(lattice), intent(in) :: lat
    type(lattice) :: inv

    double precision :: det

  det = lat%e(1)%data(1)*(lat%e(3)%data(3)*lat%e(2)%data(2)-lat%e(3)%data(2)*lat%e(2)%data(3)) &
      -lat%e(2)%data(1)*(lat%e(3)%data(3)*lat%e(1)%data(2)-lat%e(3)%data(2)*lat%e(1)%data(3)) &
      +lat%e(3)%data(1)*(lat%e(2)%data(3)*lat%e(1)%data(2)-lat%e(2)%data(2)*lat%e(1)%data(3));

  inv%e(1)%data(1) =  (lat%e(3)%data(3)*lat%e(2)%data(2)-lat%e(3)%data(2)*lat%e(2)%data(3))/det;
  inv%e(1)%data(2) = -(lat%e(3)%data(3)*lat%e(1)%data(2)-lat%e(3)%data(2)*lat%e(1)%data(3))/det;
  inv%e(1)%data(3) =  (lat%e(2)%data(3)*lat%e(1)%data(2)-lat%e(2)%data(2)*lat%e(1)%data(3))/det;
  inv%e(2)%data(1) = -(lat%e(3)%data(3)*lat%e(2)%data(1)-lat%e(3)%data(1)*lat%e(2)%data(3))/det;
  inv%e(2)%data(2) =  (lat%e(3)%data(3)*lat%e(1)%data(1)-lat%e(3)%data(1)*lat%e(1)%data(3))/det;
  inv%e(2)%data(3) = -(lat%e(2)%data(3)*lat%e(1)%data(1)-lat%e(2)%data(1)*lat%e(1)%data(3))/det;
  inv%e(3)%data(1) =  (lat%e(3)%data(2)*lat%e(2)%data(1)-lat%e(3)%data(1)*lat%e(2)%data(2))/det;
  inv%e(3)%data(2) = -(lat%e(3)%data(2)*lat%e(1)%data(1)-lat%e(3)%data(1)*lat%e(1)%data(2))/det;
  inv%e(3)%data(3) =  (lat%e(2)%data(2)*lat%e(1)%data(1)-lat%e(2)%data(1)*lat%e(1)%data(2))/det;

  END FUNCTION linv

  FUNCTION sgtype(sg) result (type)
    type(spacegroup), intent(in) :: sg
    integer :: type

    double precision :: det, tr

    det = sg%pointgroup%data(1,1)*(sg%pointgroup%data(3,3)*sg%pointgroup%data(2,2)-sg%pointgroup%data(3,2)*sg%pointgroup%data(2,3)) &
        - sg%pointgroup%data(2,1)*(sg%pointgroup%data(3,3)*sg%pointgroup%data(1,2)-sg%pointgroup%data(3,2)*sg%pointgroup%data(1,3)) &
        + sg%pointgroup%data(3,1)*(sg%pointgroup%data(2,3)*sg%pointgroup%data(1,2)-sg%pointgroup%data(2,2)*sg%pointgroup%data(1,3));
    
    tr = sg%pointgroup%data(1,1) + sg%pointgroup%data(2,2) + sg%pointgroup%data(3,3);

    IF ( .NOT. (abs(det) == 1.0) ) THEN
       WRITE(IO_ERROR,*) "Invalid pointgroup. Operation changes volume."
       RETURN
    END IF

    ! Indicate unknown type
    type = 99

    IF ( det ==  1.0 .AND. tr ==  3 ) type = 1
    IF ( det ==  1.0 .AND. tr ==  2 ) type = 6
    IF ( det ==  1.0 .AND. tr ==  1 ) type = 4
    IF ( det ==  1.0 .AND. tr ==  0 ) type = 3
    IF ( det ==  1.0 .AND. tr == -1 ) type = 2

    IF ( det == -1.0 .AND. tr == -3 ) type = -1
    IF ( det == -1.0 .AND. tr == -2 ) type = -6
    IF ( det == -1.0 .AND. tr == -1 ) type = -4
    IF ( det == -1.0 .AND. tr ==  0 ) type = -3
    IF ( det == -1.0 .AND. tr ==  1 ) type = -2

    IF ( type == 99 ) THEN
       WRITE(IO_ERROR,*) "Unknown spacegroup type. trace = ", tr, " det = ", det
       RETURN
    END IF

    RETURN

  END FUNCTION sgtype


  FUNCTION sgorder(sg) result (order)
    type(spacegroup), intent(in) :: sg
    integer :: order
    
    integer :: type

    ! Find spacegroup operation type
    type = .TYPE. sg

    ! Intialise order
    order = 99
    
    IF ( type > 0 .AND. type < 7 ) THEN
       order = type
       RETURN
    END IF
    
    IF ( type == -1 ) order = 2
    IF ( type == -6 ) order = 6
    IF ( type == -4 ) order = 4
    IF ( type == -3 ) order = 6
    IF ( type == -2 ) order = 2

    RETURN

  END FUNCTION sgorder

END MODULE DATATYPES
