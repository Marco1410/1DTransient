!***************************************************
!       Instituto Universitario Aeronautico
!       Dpto. Mecanica Aeronautica
!***************************************************
! Filename      : linearSys.f90
! Version       : 1.2
! Date          : 23-08-2017
! Programmer(s) : G. Weht
!                 F. Airaudo(fairaudo574@alumnos.iua.edu.ar)
!***************************************************
! Description (brief):
!                     This module take matrix A and 
! vector b and return the solution in vector x. Uses
! Jacobi, Gauss Seidel or Gauss Jordan methods.
!                     Also allows the user to obtain
! determinant, inverse and LU decomposition of given
! matrix A.
!*************************************************** 
! Dependecies:
!             Use tools - Filename: Utilities.f90
!             Use Defines - Filename: defines.f90
!***************************************************
! Public procedures:
!                   Type(lSystem)
!                   Subroutine GaussJordan
!                   Subroutine GaussSeidel
!                   Subroutine Jacobi
!                   Subroutine LUdecom
!                   Subroutine getx
!                   Generic setData
!                   Subroutine cleanMem
!                   Subroutine inverseA
!                   Subroutine determinant
!***************************************************

Module LinearSystems
  !***********************************************
  !*                 EXTERNS                     *
  !***********************************************
  Use tools, Only: dp, sp
  Use Defines
  Implicit none
  Private
  Public :: lSystem
  ! Type lSystem give access to the subroutines 
  Type lSystem
   Contains
     Procedure, Public :: GaussJordan
     Procedure, Public :: GaussSeidel
     Procedure, Public :: Jacobi
     Procedure, Public :: LUdecom
     Procedure, Public :: getx
     Procedure, Public :: cleanMem
     Procedure, Public :: inverseA
     Procedure, Public :: determinant
     Generic, Public :: setData => setDataFromStructure, setDataFromOther
     Procedure, Private :: setDataFromStructure, setDatafromOther
     Generic, Public :: solve => solveFromStructure, solveFromOther
     Procedure, Private :: solveFromStructure, solveFromOther

  End type lSystem

  !***********************************************
  !*          LOCAL PRIVATE VARIABLES            *
  !*********************************************** 
  Integer ndim, totalSize, PivotCount
  Integer, parameter :: log=97
  Integer, Allocatable :: ipiv(:)
  Real(dp), Allocatable :: aInit(:,:), bInit(:), x(:), x_old(:)
  Real(dp), Allocatable :: a(:,:), b(:)
  Type(lData) :: data
  Logical :: checkMat, isDiag
  Logical :: DataIsLoaded = .false.
  Logical :: dataLogging = .true.
  Character(6) ::  errorType
  
  Procedure(l2Metric), Pointer :: pError=>null()

  !*************************************************
  !*          MODULE PROCEDURES DEFINITION         *
  !*************************************************
Contains

  !***************************************************
  ! LogInit:
  !     Opens and initiates linearSys.log for other
  !     methods to write over. Writes down current
  !     date.  
  !  
  ! Parameters:
  !     Input, none
  !     Output, LinearSys.log file
  !***************************************************
  Subroutine LogInit
    implicit none
    Integer date_time(8)
    Open(unit=log, file='linearSys.log',action='write',position="append")
    Call date_and_time(VALUES=date_time)
    Write(log,'(A)') ""
    Write(log,'(A)') ""
    Write(log,'(A)') "**********************************************"
    Write(log,1)"LinearSys executed on:",date_time(3),"/",date_time(2),"/", date_time(1)
    Write(log,'(A)') "**********************************************"
1   Format(A,I2,A,I2,A,I4)
  End Subroutine LogInit

  !***************************************************
  ! setDataFromOther:
  !     Takes dimensions, matrix A, vector b from
  !     main program and stores them in local
  !     variables for use. Option to not store
  !     log on this session.
  !
  ! Parameters:
  !     Input, dimension, matrix, vector,
  !            logging(.true. or .false., true by default)
  !***************************************************

  Subroutine setDataFromOther(sys, dim, matrix, vector, logging)
    Implicit none
    Class(lSystem) :: sys
    Integer dim
    Real(dp), dimension(dim), optional :: vector
    Real(dp), dimension(dim,dim) :: matrix
    Logical, optional :: logging
    ndim=dim
    totalSize = -1
    Call sys%cleanMem()    
    If(Present(logging)) dataLogging = logging
    Allocate(aInit(ndim, ndim), bInit(ndim), x(ndim) &
         ,x_old(ndim),a(ndim,ndim),b(ndim))
    aInit=matrix
    x=0.d0; x_old=x
    a=matrix
    If (Present(vector)) then
       bInit = vector
       b = vector
    End if
    CheckMat=.true.
    totalSize=sizeof(a)+sizeof(b)+sizeof(x)+sizeof(x_old)+sizeof(ndim)+sizeof(aInit)+sizeof(bInit)
    If(dataLogging) then
       Call LogInit
       Write(log,'(A)')"----------------------------------------------"
       Write(log,'(A)')"Memory allocated"
       Write(log,'(A,I0,A)')"Total size: ", totalSize, " bytes"
       Write(log,'(A)')"----------------------------------------------"
    End If
    
    DataisLoaded=.true.
  End subroutine setDataFromOther

  !***************************************************
  ! setDataFromStructure:
  !     Takes dimensions, matrix A, vector b from
  !     lData type and stores it in local variables
  !     for use. Option to not store log on this session,
  !
  ! Parameters:
  !     Input, intdata(ndim,A(:,:),b(:)),
  !            logging(.true. or .false., true by default)
  !***************************************************
  Subroutine setDataFromStructure(sys, intdata, logging)
    Implicit none
    Class(lSystem) :: sys
    Type(lData), Intent(in) :: intdata
    Logical, optional :: logging

    totalSize = -1
    Call sys%cleanMem()
    If(Present(logging)) dataLogging = logging
    data=intdata
    ndim=data%dim
    Allocate(aInit(ndim, ndim), bInit(ndim), x(ndim) &
         ,x_old(ndim),a(ndim,ndim),b(ndim))
    aInit=data%a; bInit=data%b
    x=0.d0; x_old=x
    a=data%a; b=data%b
    CheckMat=.true.
    totalSize=sizeof(a)+sizeof(b)+sizeof(x)+sizeof(x_old)+sizeof(ndim)+sizeof(aInit)+sizeof(bInit)
    Print*,"----------------------------------------------"
    Print*, "Memory allocated"
    Print*, "Total size: ", totalSize, " bytes"
    Print*,"----------------------------------------------"
    
    If(dataLogging) then
       Call LogInit
       Write(log,'(A)')"----------------------------------------------"
       Write(log,'(A)')"Memory allocated"
       Write(log,'(A,I0,A)')"Total size: ", totalSize, " bytes"
       Write(log,'(A)')"----------------------------------------------"
    End If
       
    DataisLoaded=.true.
  End subroutine setDataFromStructure


  !***************************************************
  !  CheckDataIfLoaded: If DataIsLoaded variable is
  !                     true, then returns. If it isn't
  !                     Displays an error message and
  !                     stops program.  
  !
  !  Parameters:
  !    Input, DataIsLoaded
  !
  !    Output, Error message
  !***************************************************
  Subroutine checkDataIfLoaded
    Implicit none
    If(DataIsLoaded) Then
       a=aInit
       b=bInit   
       Return
    Else
       Print*, "***************************************"
       Print*, "Data was not set up. Program terminated"
       Print*, "Use setData to give values to operate"
       Print*, "***************************************"
       If(dataLogging) then
          Write(log,'(A)') "***************************************"
          Write(log,'(A)') "Data was not set up. Program terminated"
          Write(log,'(A)') "Use setData to give values to operate"
          Write(log,'(A)') "***************************************"
       End if
       Stop
    End If

  End Subroutine checkDataIfLoaded

  Subroutine solveFromStructure(sys, dataInOut, method)
    Implicit none
    Class(lSystem) ::  sys
    Type(lData), Intent(InOut) :: dataInOut
    Character(*), intent(In), Optional :: method
    data = dataInOut

    Call sys%setDataFromOther(data%dim, data%a, data%b, logging=.false.)
    Call chooseMethod(method)
    Call sys%getx(dataInOut%x)

  End Subroutine solveFromStructure
  
  Subroutine solveFromOther(sys, dim, matrix, vector, vectorOut, method)
    Implicit none
    Class(lSystem) :: sys
    Integer, intent(In) :: dim
    Real(dp), dimension(dim,dim), intent(In) :: matrix
    Real(dp), dimension(dim), intent(In) :: vector
    Real(dp), dimension(dim), intent(Out) :: vectorOut
    Character(*), intent(In), Optional :: method

    Call sys%setDataFromOther(dim, matrix, vector, logging=.false.)
    Call chooseMethod(method)
    Call sys%getx(vectorOut)
       
  End Subroutine solveFromOther

  Subroutine chooseMethod(method)
    Implicit none
    Character(*), Intent(In) :: method
    Type(lSystem) :: sys
    Select Case(method)
    Case('GaussJordan')
       Call sys%GaussJordan('piv')
       Return
    Case('GJ')
       Call sys%GaussJordan('piv')
       Return
    Case('GaussSeidel')
       Call sys%GaussSeidel('l2',1d-7)
       Return
    Case('GS')
       Call sys%GaussSeidel('l2',1d-7)
       Return
    Case('GaussSeidelSQRT')
       Call sys%GaussSeidel('sqrt',1d-7)
       Return
    Case('Jacobi')
       Call sys%Jacobi('l2',1d-7)
       Return
    Case('J')
       Call sys%Jacobi('l2',1d-7)
       Return
    Case('JacobiSQRT')
       Call sys%Jacobi('sqrt',1d-7)
       Return
    Case default
       Call sys%GaussJordan('piv')
       Return
    End Select
  End Subroutine chooseMethod

  
  !***************************************************
  ! GaussJordan:
  !     Executes GaussJordan method with or
  !     without pivoting
  !
  ! Parameters:
  !     Input, char(whether or not do pivot)
  !
  !***************************************************
  Subroutine GaussJordan(sys, char)
    Implicit none
    Class(lSystem) :: sys
    Character(*), Optional:: char
    Call checkDataIfLoaded
    If(Present(char)) then
       Call GaussPivot
    Else
       Call GaussDirec
    End If

  End Subroutine GaussJordan

  Subroutine GaussPivot
    Implicit none
    Call TriangSupPivot
    Call RetroSustitution
  End Subroutine GaussPivot

  Subroutine GaussDirec
    Implicit none
    Call CheckMatrix
    If(CheckMat)Then
       Call TriangSup
       Call RetroSustitution
    End If
  End Subroutine GaussDirec

  !***************************************************
  ! TriangSup:
  !     Given matrix A and vector b return 
  !     upper triangular matrix.
  !
  ! Parameters:
  !     Input, A(:,:), b(:), ndim
  !     Output, A(:,:) , b(:) upper triangular matrix
  !***************************************************
  Subroutine TriangSup
    Implicit none
    Integer i, j, k
    Real(dp) m
    If(dataLogging) then
       write(log,'(A)') "TriangSup executing"
       Write(log,'(A)')"----------------------------------------------"
    End if
    Do i=1,ndim-1
       Do j=i+1,ndim
          m=a(j,i)/a(i,i)
          Do k=i+1,ndim
             a(j,k)=a(j,k)-m*a(i,k)
          End Do
          b(j)=b(j)-m*b(i)
       End Do
    End Do
  End Subroutine TriangSup

  !***************************************************
  ! TriangSupPivot:
  !     Given matrix A and vector b return 
  !     upper triangular matrix. Using Pivot
  !
  ! Parameters:
  !     Input, A(:,:), b(:), ndim
  !     Output, A(:,:) , b(:) upper triangular matrix
  !***************************************************
  Subroutine TriangSupPivot
    Implicit none
    Integer i, j, k
    Real(dp) m
    If(dataLogging) then
       write(log,'(A)') "TriangSup with Pivot executing"
       Write(log,'(A)')"----------------------------------------------"
    End if
    Do i=1,ndim-1
       Call Pivot(i)
       Do j=i+1,ndim
          m=a(j,i)/a(i,i)
          Do k=i+1,ndim
             a(j,k)=a(j,k)-m*a(i,k)
          End Do
          b(j)=b(j)-m*b(i)
       End Do
    End Do
  End Subroutine TriangSupPivot

  !***************************************************
  ! RetroSustitution:
  !     Given Upper triangular matrix A
  !     and vector b, return vector x
  !  
  ! Parameters:
  !     Input, A(:,:), b(:), ndim
  !     Output, x(:)
  !***************************************************  
  Subroutine RetroSustitution
    Implicit none
    Integer i, j
    Real(dp) Sum
    If(dataLogging) then
       write(log,'(A)') "Retrosustitution executing"
       Write(log,'(A)')"----------------------------------------------"
    End if
    Do i=ndim, 1,-1
       Sum=b(i)
       Do j=i+1,ndim
          Sum=Sum-a(i,j)*x(j)
       End Do
       x(i)=Sum/a(i,i)
    End Do
  End Subroutine RetroSustitution

  !***************************************************
  ! Pivot:
  !     Finds biggest coefficient in 'i' column.
  !     Pivots row with biggest coefficient with
  !     the 'i' row. 
  !
  ! Parameters:
  !     Input, A(:,:), b(:), ndim, i
  !     Output, A(:,:) , b(:), pivoted matrix
  !***************************************************  
  Subroutine Pivot(k)
    Implicit none
    Integer i, ind
    Integer, Intent(In) :: k
    Logical enter
    Real(dp), Allocatable :: a_aux(:)
    Real(dp) a_max, b_aux, a_ik
    If(.not. Allocated(a_aux)) Allocate(a_aux(ndim))
    a_max=abs(a(k,k))
    enter=.false.
    Do i=k+1,ndim
       a_ik=abs(A(i,k))
       If(a_ik > a_max)Then
          a_max= a_ik
          ind = i
          enter=.true.
       End If
    End Do
    If(enter .and. ind .ne. k) Then
       If(dataLogging) write(log,'(A,I0,A,I0)') "Pivoting row", k, " with", ind
       a_aux(:)= a(ind,:)
       a(ind,:)= a(k,:)
       a(k,:)  = a_aux(:)

       b_aux = b(ind)
       b(ind)= B(k)
       b(k)  = b_aux
       PivotCount=PivotCount+1
    End If
    Deallocate(a_aux)
  End Subroutine Pivot

  !***************************************************
  ! CheckMatrix:
  !     Checks whether matrix A has a zero
  !     in its diagonal 
  !
  ! Parameters:
  !     Input, A(:,:), ndim
  !     Output, CheckMat, error message
  !***************************************************
  Subroutine CheckMatrix
    Implicit none
    Integer i
    Real(dp) limInf
    limInf=1.d-13
    Do i=1, ndim
       If(abs(a(i,i)) < limInf)Then
          Print*, "Error..."
          Print"(a,i2,a,i2)", "Zero in diagonal row: ", i," col: ", i
          Print*, "Try with Pivot option"
          Print*, "----------------------------------------------"
          CheckMat=.false.
          If(dataLogging) then
             Write(log,'(A)') "Error..."
             Write(log,'(A,I0,A,I0)') "Zero in diagonal row: ", i," col: ", i
             write(log,'(A)')"Try with Pivot option"
             Write(log,'(A)')"----------------------------------------------"
          End if
          Return
       End If
    End Do
  End Subroutine CheckMatrix

  !***************************************************
  ! sqrtMatErr:
  !     With current x values obtains the
  !     difference between the left side and
  !     right side of each equation, adds the
  !     squares, then gets the square root of
  !     the sum. ||Ax-b||
  !
  ! Parameters:
  !     Input, A(:,:), b(:), x(:), ndim
  !     Output, sqrtMatErr
  !***************************************************
  Real(dp) Function sqrtMatErr()
    Implicit none
    Integer i, j
    Real(dp) Sum, dummy
    dummy=0.d0
    Do i=1, ndim
       Sum=0.d0
       Do j=1,ndim
          Sum=Sum+data%a(i,j)*x(j)
       End Do
       dummy=dummy+(Sum-data%b(i))**2
    End Do
    sqrtMatErr=sqrt(dummy)
  End Function sqrtMatErr

  !***************************************************
  ! l2Metric:
  !     Gets the difference between current
  !     value of x and old one. Adds the squares
  !     then gets the square root of the sum.
  !     |x-x_old| 
  !
  ! Parameters:
  !     Input, x(:), x_old(:), ndim
  !     Output, l2Metric
  !***************************************************
  Real(dp) Function  l2Metric()
    Implicit none
    Integer i, j
    Real(dp) Sum
    Sum=0.d0
    Do i=1, ndim
       Sum=Sum+(x(i)-x_old(i))**2
    End Do
    l2Metric=sqrt(Sum)
  End Function l2Metric

  !***************************************************
  ! errorChoose:
  !     Checks input char and decides 
  !     whether to get error from sqrtMatErr 
  !     or l2Metric. 
  !   
  ! Parameters:
  !     Input, char (type of error wanted)
  !     Output, pError (will point to the wanted subroutine) 
  !***************************************************
  Subroutine errorChoose
    Implicit none
    Select case(errorType)
    case("sqrt")
       pError=> sqrtMatErr
       If(dataLogging) then
          Write(log,'(A)') "Running iterative method with Square matrix error ||Ax-b||"
          Write(log,'(A)')"----------------------------------------------"
       End if
    Case("SQRT")
       pError=> sqrtMatErr
       If(dataLogging) then
          Write(log,'(A)') "Running iterative method with Square matrix error ||Ax-b||"
          Write(log,'(A)')"----------------------------------------------"
       End if
    Case("l2")
       pError=> l2Metric
       If(dataLogging) then
          Write(log,'(A)') "Running iterative method with L2 metric: ||x-x_old||"
          Write(log,'(A)')"----------------------------------------------"
       End if
    Case("L2") 
       pError=> l2Metric
       If(dataLogging) then
          Write(log,'(A)') "Running iterative method with L2 metric: ||x-x_old||"
          Write(log,'(A)')"----------------------------------------------"
       End if
    Case default
       pError=> l2Metric
       If(dataLogging) then
          Write(log,'(A)') "Running iterative method with L2 metric: ||x-x_old||"
          Write(log,'(A)')"----------------------------------------------"
       End if
    End Select
    Print*, "----------------------------------------------"
  End Subroutine errorChoose

  !***************************************************
  ! GaussSeidel:
  !     Given matrix A and vector b, return
  !     vector x by iteration using certain
  !     error and tolerance. Using
  !     GaussSeidel method. 
  !
  ! Parameters:
  !     Input, A(:,:), b(:), ndim, error, tol
  !     Output, x(:)
  !***************************************************
  Subroutine GaussSeidel(sys, error, tol)
    Implicit none
    Class(lSystem) :: sys
    Integer, Parameter :: cut=1000
    Integer iter, i, j
    Real(dp) err, tol, Sum
    Character(*) error
    errorType=trim(error)
    Call checkDataIfLoaded
    Call isDominant
    Call errorChoose !!choose the errorType by error var
    err=tol+1.d0
    iter=0 ; x=0.d0 ; x_old=x
    If(dataLogging) then
       Write(log,'(A)') "Initial X0 set to zero. Running Gauss-Seidel"
       Write(log,'(A)')"----------------------------------------------"
    End if
    Do While(err > tol .and. iter < cut)
       iter=iter+1
       Do i=1,ndim
          Sum=b(i)
          Do j=1, i-1
             Sum=Sum-a(i,j)*x(j)
          End Do
          Do j=i+1,ndim
             Sum=Sum-a(i,j)*x(j)
          End Do
          x(i)=Sum/a(i,i)
       End Do
       err=pError()
       x_old=x
    End Do
    If(iter==cut)then 
       Print*, "**Method cut by iterations**"
       Print*, "----------------------------------------------"
       If(dataLogging) then
          Write(log,'(A)') "**Method cut by iterations**"
          Write(log,'(A)')"----------------------------------------------"
       End if
    End If

    If(dataLogging) then
       Write(log,"(2x,A, e10.3,2x,A,i0)") "Error:",err, "Nº Iter:",iter
       Write(log,'(A)') "End of Gauss-Seidel"
       Write(log,'(A)')"----------------------------------------------"
    End if
  End Subroutine GaussSeidel

  !***************************************************
  ! Jacobi:
  !     Given matrix A and vector b, return
  !     vector x by iteration using certain
  !     error and tolerance. Using
  !     Jacobi method. 
  !
  ! Parameters:
  !     Input, A(:,:), b(:), ndim, error, tol
  !     Output, x(:)
  !***************************************************
  Subroutine Jacobi(sys, error, tol)
    Implicit none  
    Class(lSystem) :: sys
    Integer, Parameter :: cut=1000
    Integer iter, i, j
    Real(dp) err, tol, Sum
    Character(*) error
    errorType=trim(error)
    Call checkDataIfLoaded
    Call isDominant
    Call errorChoose !!choose the errorType by error var
    err=tol+1.d0
    iter=0 ; x=0.d0 ; x_old=x ! starting from x solution
    If(dataLogging) then
       Write(log,'(A)') "initial x0 set to zero. Running Jacobi"
       Write(log,'(A)') "----------------------------------------------"
    End if
    Do While(err > tol .and. iter < cut)
       iter=iter+1
       Do i=1,ndim
          Sum=b(i)
          Do j=1, i-1
             Sum=Sum-a(i,j)*x_old(j)
          End Do
          Do j=i+1,ndim
             Sum=Sum-a(i,j)*x_old(j)
          End Do
          x(i)=Sum/a(i,i)
       End Do
       err=pError()
       x_old=x
    End Do
    If(iter==cut) then
       Print*, "Method cut by iterations"
       Print*, "----------------------------------------------"
       If(dataLogging) then
          Write(log,'(A)') "Method cut by iterations"
          Write(log,'(A)')"----------------------------------------------"
       End if
    End If
    If(dataLogging) then
       Write(log,"(2x,A, e10.3,2x,A,i0)") "Error:",err, "Nº Iter:",iter
       Write(log,'(A)') "End of Jacobi"
       Write(log,'(A)')"----------------------------------------------"
    End if
  End Subroutine Jacobi

  !***************************************************
  ! Dominant:
  !     Checks whether matrix A is diagonally
  !     dominant. If it isn't, execute pivot
  !     and check again.  
  !
  ! Parameters:
  !     Input, A(:,:), b(:), ndim
  !     Output, isDiag(true/false)
  !***************************************************
  Subroutine Dominant(isDiag)
    Implicit none
    Integer i, j
    Real(dp) Sum
    logical isDiag
    If(dataLogging) then
       Write(log,'(A)') "Checking if matrix if diagonally dominant"
       Write(log,'(A)')"----------------------------------------------"
    End if
    isDiag=.true.
    Do i =1, ndim
       Sum=0.d0
       Do j=1, i-1
          Sum=Sum+abs(a(i,j))
       End Do
       Do j=i+1, ndim
          Sum=Sum+abs(a(i,j))
       End Do
       If(abs(a(i,i))<Sum) Then
          isDiag=.false.
          Call pivot(i)
       End If
    End Do


  End Subroutine Dominant

  !***************************************************
  ! isDominant:
  !     Uses dominant routine to check if a
  !     matrix if diagonally dominant, if
  !     it isn't, attempts to diagonalize
  !     it, then displays a message.
  !
  ! Parameters:
  !     Input, A(:,:), b(:), ndim
  !     Output, message
  !*************************************************** 
  Subroutine isDominant
    Implicit none
    logical isDiag
    Call checkDataIfLoaded
    call Dominant(isDiag)
    If (isDiag) then
       If(dataLogging) then
          Write(log,'(A)')"Matrix is diagonally dominant"
          Write(log,'(A)')"----------------------------------------------"
       End if
       Return
    Else
       If(dataLogging) then
          Write(log,'(A)')"**Warning, matrix non-diagonal**"
          Write(log,'(A)')"**attempting diagonalization**"
          Write(log,'(A)')"----------------------------------------------"
       End if
       call Dominant(isDiag)
    End If
    If (isDiag) then
       If(dataLogging) then
          Write(log,'(A)') "Matrix is diagonally dominant"
          Write(log,'(A)')"----------------------------------------------"
       End if
       Return
    Else
       If(dataLogging) then
          Write(log,'(A)') "**Could not diagonalize matrix**"
          Write(log,'(A)')"----------------------------------------------"
       End if
    End If
  End Subroutine isDominant

  !***************************************************
  !  CalcDet:
  !     Obtains the determinant value from a
  !     upper triangular matrix by adding the
  !     the values from its diagonal.
  !
  !  Parameters:
  !     Input, mat(:,:)
  !     Output, calcDet
  !***************************************************
  Real(dp) Function calcDet()
    Implicit none
    integer i
    calcDet=1
    Do i=1,Size(a,1)
       calcDet=CalcDet*a(i,i)
    End Do
    calcDet=CalcDet*(-1)**(pivotCount)
  End Function calcDet

  !***************************************************
  ! determinantSys:
  !     Reduces matrix A from a system
  !     to its upper triangular to
  !     obtain determinant from.
  !  
  ! Parameters:
  !     Input, A(:,:)
  !     Output, -
  !***************************************************
  Subroutine determinant(sys)
    Implicit none
    Class(lSystem) :: sys
    real(8) :: det
    integer i, j
    Call CheckDataIfLoaded
    PivotCount=0
    If(dataLogging) then
       Write(log,'(A)') "Obtaining determinant from matrix A"
       Write(log,'(A)')"----------------------------------------------"
    End if
    Call TriangSupPivot
    print*, "Determinant:", calcDet()
    Print*, "----------------------------------------------"
    If(dataLogging) then
       Write(log,'(A,E10.3)') "Determinant:", calcDet()
       Write(log,'(A)')"----------------------------------------------"
    End if
  End Subroutine determinant

  !***************************************************
  ! LUdecom:
  !     Runs LU decomposition on a given matrix
  !     A. Then solves it with retrosustitution.
  !  
  ! Parameters:
  !     Input, A(:,:), b(:), ndim
  !     Output, message
  !***************************************************
  Subroutine LUdecom(sys, LU)
    Implicit none
    Class(lSystem) :: sys
    Real(dp), Intent(Out) :: LU(ndim,ndim)
    If(dataLogging) then
       Write(*,'(A)') "Running LU decomposition"
       Write(*,'(A)') "----------------------------------------------"
       Write(log,'(A)') "Running LU decomposition"
       Write(log,'(A)') "----------------------------------------------"
    End if
    Call GaussLU
    Call GaussRES(b, x)
    !Print"(2x,A,e10.3,/)", "Error: ", sqrtMatErr()
    LU=a
  End Subroutine LUdecom

  !***************************************************
  ! Gauss LU:
  !     Given A matrix Obtains decomposition
  !     A=L.U
  !  
  ! Parameters:
  !     Input, A(:,:), ndim
  !     Output, A(:,:) = L.U
  !***************************************************
  Subroutine GaussLU
    Implicit none
    Integer i, k, j, j_max
    Real(dp) a_max, pivinv, c
    Call checkDataIfLoaded
    If(.not.allocated(ipiv)) allocate(ipiv(ndim))
    Do i=1,ndim
       ipiv(i)=i
    End Do
    Do k=1,ndim-1
       a_max=abs(a(ipiv(k),k))
       j_max=k        
       Do i=k+1,ndim
          If(abs(a(ipiv(i),k)).gt.a_max) Then
             a_max=abs(a(ipiv(i),k))
             j_max=i
          End If
       End Do
       i=ipiv(k)
       ipiv(k)=ipiv(j_max)
       ipiv(j_max)=i
       pivinv=1.D0/a(ipiv(k),k)
       Do i=k+1,ndim
          c=a(ipiv(i),k)*pivinv
          a(ipiv(i),k)=c
          Do j=k+1,ndim
             a(ipiv(i),j)=a(ipiv(i),j)-c*a(ipiv(k),j)
          End Do
       End Do
    End Do
  End Subroutine GaussLU

  !***************************************************
  ! GaussRES:
  !     From triangular matrix A=L.U solves
  !     system Ax=b  
  !
  ! Parameters:
  !     Input, A(:,:), b(:), ndim
  !     Output, x(:)
  !***************************************************
  Subroutine GaussRES(b, x)
    Implicit none
    Real(dp), Intent(In) :: b(:)
    Real(dp), Intent(Out) :: x(:)
    Real(dp) c
    Integer i, j
    Do i=1,ndim
       c=0.D0
       Do j=1,i-1
          c=c+a(ipiv(i),j)*x(j)
       End Do
       x(i)=b(ipiv(i))-c
    End Do
    DO i=ndim,1,-1
       c=0.D0
       Do j=i+1,ndim
          c=c+a(ipiv(i),j)*x(j)
       End Do
       x(i)=(x(i)-c)/a(ipiv(i),i)
    ENDDO
  End Subroutine GaussRES

  !***************************************************
  ! inverseA:
  !     Using routine GaussLU and GaussRES obtains
  !     the inverse of matrix A.  
  !
  ! Parameters:
  !     Input, A(:,:), ndim
  !     Output, Inverse(:,:)
  !***************************************************  
  Subroutine inverseA(sys, Inverse)
    Implicit none
    Class(lSystem) :: sys
    Real(dp), Dimension(ndim) :: b_aux, x_aux
    Real(dp), Intent(Out) :: Inverse(ndim,ndim)
    Integer i, j
    Print*, "Obtaining inverse of matrix A"
    Print*, "----------------------------------------------"
    Write(log,'(A)') "Obtaining inverse of matrix A"
    Write(log,'(A)') "----------------------------------------------"
    Call GaussLU
    b_aux=0.d0
    Do i=1,ndim
       b_aux(i)=1.d0
       Call GaussRES(b_aux, x_aux)
       b_aux(i)=0.d0
       Do j=1, ndim
          Inverse(j,i)=x_aux(j)
       End Do
    End Do

  End Subroutine inverseA

  !***************************************************
  ! getx:
  !     Asign value of vector x(Solution of linear system)
  !     to vector ret.
  !
  ! Parameters:
  !     Input, x(:), ndim
  !     Output, ret(:)
  !***************************************************
  Subroutine getx(sys, ret)
    Implicit none
    Class(lSystem) :: sys
    Real(dp), Intent(Out)::Ret(ndim)
    Call checkDataIfLoaded
    ret=x
  End Subroutine getx

  !***************************************************
  ! cleanMem:
  !     deallocates all vectors and matrices,
  !     shows ammount of memory cleaned. 
  !
  ! Parameters:
  !     Input, A(:,:), b(:), x(:), x_old(:)
  !     Output, message of cleared memory
  !***************************************************
  Subroutine cleanMem(sys)
    Implicit none
    Class(lSystem) :: sys

    dataLogging = .true.
    If(Allocated(a)) Deallocate(a)
    If(Allocated(b)) Deallocate(b)
    If(Allocated(x)) Deallocate(x)
    If(Allocated(aInit)) Deallocate(aInit)
    If(Allocated(bInit)) Deallocate(bInit)
    If(Allocated(x_old)) Deallocate(x_old)
    If(Allocated(ipiv)) Deallocate(ipiv)

    If (totalSize > 0) then
       Print*, "Memory cleaned"
       Print"(a,i4,a)", " Total size: ", totalSize, " bytes"
       Print*,"----------------------------------------------"
       If(dataLogging) then
          Write(log,'(A)')"Memory cleaned"
          Write(log,'(A,I0,A)') "Total size: ", totalSize, " bytes"
          Write(log,'(A)')"----------------------------------------------"
       End If
    End If
  End Subroutine cleanMem
        
  !***********************************************
  !*       END MODULE PROCEDURES DEFINITION      *
  !***********************************************
End Module LinearSystems


