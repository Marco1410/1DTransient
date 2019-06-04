Program main
  Use GidDataMod
  Use Poisson1DMod
  Implicit none
  Type(Poisson1DType) :: problem
  Call ReadGid(problem)
  Call problem%init
  Call problem%solve
  Call WriteGid(problem)
End Program main


