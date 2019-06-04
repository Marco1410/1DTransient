Module sparseMod
  Use tools
  Use Print_All
  Use quickSortMod
  Implicit none
  Private
  Public :: sparseType, tripletType
  Type TripletType
     Real(dp), Dimension(:), Allocatable :: A
     Integer, Dimension(:), Allocatable :: row
     Integer, Dimension(:), Allocatable :: col
  End type TripletType
  Type sparseType
     Real(dp), Dimension(:), Allocatable :: A
     Integer, Dimension(:), Allocatable :: AI
     Integer, Dimension(:), Allocatable :: AJ
     Type(tripletType) :: triplet
   Contains
     Procedure, Public :: init
     Procedure, Public :: append
     Procedure, Public :: getSparse
  End type sparseType

  Real(dp), Dimension(:), Allocatable :: valueVector
  Integer, Dimension(:), Allocatable :: rowVector
  Integer, Dimension(:), Allocatable :: rowCounter
  Integer :: counter, repeats, i, j, k, n, m
Contains

  Subroutine init(this, nnz, rows)
    Implicit none
    Class(sparseType), Intent(InOut) :: this
    Integer, Intent(In) :: nnz, rows
    n = nnz
    m = rows
    Allocate(this%triplet%A(n))
    Allocate(this%triplet%row(n))
    Allocate(this%triplet%col(n))
    this%triplet%A = 0
    this%triplet%row = 0
    this%triplet%col = 0
    counter = 0
  End Subroutine init

  Subroutine append(this, value, row, col)
    Implicit none
    Class(sparseType), Intent(InOut) :: this
    Real(dp), Intent(In) :: value
    Integer, Intent(In) :: row, col
    counter = counter + 1
    this%triplet%A(counter) = value
    this%triplet%row(counter) = row
    this%triplet%col(counter) = col
  End Subroutine append

  Subroutine getSparse(this)
    Use st_to_cc_burkardt
    Implicit none
    Class(sparseType), Intent(InOut) :: this
    Integer :: ncc
    !John Burkardt's program
    Call st_to_cc_size(counter, this%triplet%row, this%triplet%col, ncc)
    Allocate(this%AJ(ncc),this%AI(m+1))
    Call st_to_cc_index(counter, this%triplet%row, this%triplet%col &
         , ncc, m, this%AJ, this%AI)
    Allocate(this%A(ncc))
    Call st_to_cc_values(counter, this%triplet%row, this%triplet%col &
         , this%triplet%A, ncc, m, this%AJ, this%AI, this%A)
    Deallocate(this%triplet%A, this%triplet%row, this%triplet%col)
  End Subroutine getSparse

  Subroutine handleDuplicates(this)
    Implicit none
    Class(sparseType), Intent(InOut) :: this
    Logical :: mask
    Allocate(rowVector(maxval(rowCounter)))
    Allocate(valueVector(maxval(rowCounter)))
    counter = 1
    repeats = 0
    Do i = 1, m
       rowVector = this%triplet%col(counter+repeats:counter+repeats+rowCounter(i)-1)
       valueVector = this%triplet%A(counter+repeats:counter+repeats+rowCounter(i)-1)
       j = 0
       Do while(j < rowCounter(i))
          j = j + 1
          this%A(counter) = valueVector(j)
          this%AJ(counter) = rowVector(j)
          k = j
          Do while(k < rowCounter(i))
             k = k + 1
             mask = rowVector(j).eq.rowVector(k)
             If(mask) then
                !add values from k to j
                this%A(counter) = this%A(counter) + valueVector(k)
                !move k to the back
                Call swap(rowVector(k), rowVector(rowCounter(i)))
                Call swap(valueVector(k), valueVector(rowCounter(i)))
                rowCounter(i) = rowCounter(i) - 1
                repeats = repeats + 1
             End If
          End Do
          counter = counter + 1
       End Do
    End Do
    DeAllocate(rowVector)
    DeAllocate(valueVector)
  End Subroutine HandleDuplicates
 
End Module sparseMod
