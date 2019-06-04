Module domainMod
  Use tools
  Use element1DMod
  Implicit none
  Private
  Public domainType, element1DPoissonType
  Type, extends(element1DType) :: element1DPoissonType
     Integer :: materialIndex
     Integer :: sourceIndex
  End type element1DPoissonType
  Type domainType
     Integer :: nElem, nNodes
     Real(dp), dimension(:), allocatable :: x
     Type(element1DPoissonType), Dimension(:), Allocatable :: element
  End type domainType
End Module domainMod
