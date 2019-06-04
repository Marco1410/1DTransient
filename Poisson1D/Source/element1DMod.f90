Module element1DMod
  Use tools
  Implicit none
  Private
  Public :: element1DType, linearN, lineardN &
       , cuadraticN, cuadraticdN, cubicN, cubicdN
  Type element1DType
     Integer :: nNodes
     Integer, Dimension(:), Allocatable :: node
     Real(dp), Dimension(:), Allocatable :: jacobian
     Procedure(linearN), Pointer, nopass :: N => Null()
     Procedure(lineardN), Pointer, nopass :: dN => Null()
   Contains
     Procedure :: init
  End type element1DType

Contains

  Subroutine init(this, nodes)
    Implicit none
    Class(element1DType), Intent(InOut) :: this
    Integer, Intent(In) :: nodes
    this%nNodes = nodes
    Allocate(this%node(nodes))
    Select Case(nodes)
       Case(2)
          this%N => linearN
          this%dN => lineardN
       Case(3)
          this%N => cuadraticN
          this%dN => cuadraticdN
       Case(4)
          this%N => cubicN
          this%dN => cubicdN
       Case Default
          Print*, 'Orden de elemento incorrecto'
       End Select
     End Subroutine init
     
  Function linearN(u, n)
    Real(dp), Intent(In) :: u
    Integer, Intent(In) :: n
    Real(dp), Dimension(n) :: linearN
    linearN(1) = .5 - .5*u
    linearN(2) = .5 + .5*u
  End Function linearN
  Function lineardN(u, n)
    Real(dp), Intent(In) :: u
    Integer, Intent(In) :: n
    Real(dp), Dimension(n) :: lineardN
    lineardN(1) = -.5
    lineardN(2) = .5
  End Function lineardN
  Function cuadraticN(u, n)
    Implicit none
    Real(dp), Intent(In) :: u
    Integer, Intent(In) :: n
    Real(dp), Dimension(n) :: cuadraticN
    cuadraticN(1) = 0.5*u*(u-1)
    cuadraticN(2) = (1+u)*(1-u)
    cuadraticN(3) = 0.5*u*(1+u)
  End Function cuadraticN
  Function cuadraticdN(u, n)
    Implicit none
    Real(dp), Intent(In) :: u
    Integer, Intent(In) :: n
    Real(dp), Dimension(n) :: cuadraticdN
    cuadraticdN(1) = u-0.5
    cuadraticdN(2) = -2*u
    cuadraticdN(3) = u+0.5
  End Function cuadraticdN
  Function cubicN(u, n)
    Implicit none
    Real(dp), Intent(In) :: u
    Integer, Intent(In) :: n
    Real(dp), Dimension(n) :: cubicN
    cubicN(1) = -(9.d0/16.d0)*(u+1.d0/3.d0)*(u-1.d0/3.d0)*(u-1)
    cubicN(2) = (27.d0/16.d0)*(u+1)*(u-1.d0/3.d0)*(u-1)
    cubicN(3) = -(27.d0/16.d0)*(u+1)*(u+1.d0/3.d0)*(u-1)
    cubicN(4) = (9.d0/16.d0)*(u+1)*(u+1.d0/3.d0)*(u-1.d0/3.d0)
  End Function cubicN
  Function cubicdN(u, n)
    Implicit none
    Real(dp), Intent(In) :: u
    Integer, Intent(In) :: n
    Real(dp), Dimension(n) :: cubicdN
    cubicdN(1) = -(9.d0/16.d0)*(3*u**2-2*u-1.d0/9.d0)
    cubicdN(2) = (27.d0/16.d0)*(3*u**2-(2.d0/3.d0)*u-1)
    cubicdN(3) = -(27.d0/16.d0)*(3*u**2+(2.d0/3.d0)*u-1)
    cubicdN(4) = (9.d0/16.d0)*(3*u**2+2*u-1.d0/9.d0)
  End Function cubicdN
End Module element1DMod
