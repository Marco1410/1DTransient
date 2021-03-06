Module GaussLegendreFEM
  Use tools
  Implicit none
  Public

  Interface valueGauss
     Module Procedure valueGauss1
     Module Procedure valueGauss2
     Module Procedure valueGauss3
     Module Procedure valueGauss4
  End Interface valueGauss

  Interface gauss
     Module Procedure gauss1vect
     Module Procedure gauss2vect
     Module Procedure gauss3vect
     Module Procedure gauss4vect
     Module Procedure gauss5vect
  End Interface gauss

  Integer, private :: i, j, n
  Real(dp), Dimension(:,:), Allocatable, Private :: G
Contains
  
  Subroutine setGauss(nGauss)
    Implicit none
    Integer, Intent(In) :: nGauss
    n = nGauss
    Allocate(G(2,n))
    Call gaussQuad
  End Subroutine setGauss

  Subroutine valueGauss1(f, vect)
    Implicit none
    Real(dp), Dimension(:) :: vect
    Interface
       Real(8) Function f(x)
         Real(8), Intent(In) :: x
       End Function f
    End Interface
    Do i = 1, n
       vect(i) = f(G(1,i))
    End Do
  End Subroutine valueGauss1

  Subroutine valueGauss2(f, mat)
    Implicit none
    Real(dp), Dimension(:,:), Intent(Out) :: mat
    Integer :: size_of_f
    Interface
       Function f(x, n)
         Real(8), Intent(In) :: x
         Integer, Intent(In) :: n
         Real(8), Dimension(n) :: f
       End Function f
    End Interface
    size_of_f = size(mat,1)
    Do i = 1, n
       mat(1:size_of_f,i) = f(G(1,i),size_of_f)
    End Do
  End Subroutine valueGauss2

  Subroutine valueGauss3(f, index,  vect, lower, upper)
    Implicit none
    Real(dp), Dimension(:), Intent(Out) :: vect
    Real(dp), Intent(In) :: lower, upper
    Integer, Intent(In) :: index
    Interface
       Function f(i, x)
         Integer, Intent(In) :: i
         Real(8), Intent(In) :: x
         Real(8) :: f
       End Function f
    End Interface
    Do i = 1, n
       vect(i) = f(index, ((upper+lower)/2.d0)+((upper-lower)/2.d0)*G(1,i))
    End Do
  End Subroutine valueGauss3

  Subroutine valueGauss4(f, index,  vect1, vect2, lower, upper)
    Implicit none
    Real(dp), Dimension(:), Intent(Out) :: vect1, vect2
    Real(dp), Intent(In) :: lower, upper
    Integer, Intent(In) :: index
    Interface
       Function f(i, x)
         Integer, Intent(In) :: i
         Real(8), Intent(In) :: x
         Real(8) :: f
       End Function f
    End Interface
    Do i = 1, n
       vect1(i) = f(index, ((upper+lower)/2.d0)+((upper-lower)/2.d0)*G(1,i))
       vect2(i) = f(index, ((upper+lower)/2.d0)+((upper-lower)/2.d0)*G(1,i))
    End Do
  End Subroutine valueGauss4

  Real(dp) Function gauss1vect(vect1)
    Implicit none
    Real(dp), Dimension(n), Intent(In) :: vect1
    gauss1vect = 0
    Do i = 1, n
       gauss1vect = gauss1vect + G(2,i)*vect1(i)
    End Do
  End Function gauss1vect
  Real(dp) Function gauss2vect(vect1, vect2)
    Implicit none
    Real(dp), Dimension(n), Intent(In) :: vect1, vect2
    gauss2vect = 0
    Do i = 1, n
       gauss2vect = gauss2vect + G(2,i)*vect1(i)*vect2(i)
    End Do
  End Function gauss2vect
  Real(dp) Function gauss3vect(vect1, vect2, vect3)
    Implicit none
    Real(dp), Dimension(n), Intent(In) :: vect1, vect2, vect3
    gauss3vect = 0
    Do i = 1, n
       gauss3vect = gauss3vect + G(2,i)*vect1(i)*vect2(i)*vect3(i)
    End Do
  End Function gauss3vect
  Real(dp) Function gauss4vect(vect1, vect2, vect3, vect4)
    Implicit none
    Real(dp), Dimension(n), Intent(In) :: vect1, vect2, vect3, vect4
    gauss4vect = 0
    Do i = 1, n
       gauss4vect = gauss4vect + G(2,i)*vect1(i)*vect2(i)*vect3(i)*vect4(i)
    End Do
  End Function gauss4vect
  Real(dp) Function gauss5vect(vect1, vect2, vect3, vect4, vect5)
    Implicit none
    Real(dp), Dimension(n), Intent(In) :: vect1, vect2, vect3, vect4, vect5
    gauss5vect = 0
    Do i = 1, n
       gauss5vect = gauss5vect + G(2,i)*vect1(i)*vect2(i)*vect3(i)*vect4(i)*vect5(i)
    End Do
  End Function gauss5vect
  Subroutine gaussQuad
    Implicit none
    Real(dp), Parameter :: pi = dacos(-1.d0)
    Real(dp) :: f, df, dx, r
    Integer :: i, iter, k
    Real(dp), Dimension(:), Allocatable :: p0, p1, tmp
    p0 = [1.d0]
    p1 = [1.d0, 0.d0]
    Do k = 2, n
       tmp = ((2*k-1)*[p1,0.d0]-(k-1)*[0.d0, 0.d0,p0])/k
       p0 = p1; p1 = tmp
    End Do
    Do i = 1, n
       r = cos(pi*(i-0.25)/(n+0.5))
       Do iter = 1, 10
          f = p1(1); df = 0.
          Do k = 2, Size(p1)
             df = f + r*df
             f  = p1(k) + r * f
          End Do
          dx =  f / df
          r = r - dx
          If (Abs(dx)<10*Epsilon(dx)) Exit
       End Do
       G(1,i) = r
       G(2,i) = 2/((1-r**2)*df**2)
    End Do
  End Subroutine gaussQuad
End Module GaussLegendreFEM
