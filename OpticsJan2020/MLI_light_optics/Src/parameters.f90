MODULE parameters
  !--------- -------- --------- --------- --------- --------- --------- --------- -----
  ! Specify data types
  !--------- -------- --------- --------- --------- --------- --------- --------- -----
  IMPLICIT NONE
  INTEGER, PARAMETER :: rn = KIND(0.0d0)          ! Precision of real numbers
  INTEGER, PARAMETER :: is = SELECTED_INT_KIND(4) ! Data type of bytecode
END MODULE parameters


