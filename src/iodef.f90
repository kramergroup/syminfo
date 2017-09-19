MODULE IODEF

  IMPLICIT NONE
  
  INTEGER :: IO_IN = 5,IO_OUT = 6, IO_DEBUG = 99, IO_TRACE = 99, IO_ERROR = 0, IO_NULL = 99
  SAVE

CONTAINS
  
  SUBROUTINE INIT_IO(enable_debug, enable_trace)
    logical,optional :: enable_debug, enable_trace
    

    IF ( present(enable_debug) .AND. enable_debug ) THEN
       IO_DEBUG = IO_OUT
    END IF

    IF ( present(enable_trace) .AND. enable_trace ) THEN
       IO_TRACE = IO_OUT
    END IF
    
    OPEN(99,FILE="/dev/null")

  END SUBROUTINE INIT_IO

  SUBROUTINE ENABLE_TRACE()
    IO_TRACE = IO_OUT
  END SUBROUTINE ENABLE_TRACE
  
  SUBROUTINE ENABLE_DEBUG()
    IO_DEBUG = IO_OUT
  END SUBROUTINE ENABLE_DEBUG
END MODULE IODEF
