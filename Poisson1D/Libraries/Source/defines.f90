!***************************************************
!       Instituto Universitario Aeronautico
!       Dpto. Mecanica Aeronautica
!***************************************************
! Filename      : defines.f90
! Version       : 1.1
! Date          : 28-08-2017
! Programmer(s) : G. Weht
!                 F. Airaudo(fairaudo574@alumnos.iua.edu.ar)
!***************************************************
! Description (brief):
!                     This module allows the user to
!                     store within a single structure
!                     a whole linear system allocated
!                     with its dimension.
!                     The elements of the system can
!                     be input from a datafile or from
!                     the main program. It also allows
!                     the user to not initialize a
!                     vector b.
!*************************************************** 
! Dependecies:
!             Tools - Filename: Utilities.f90
!***************************************************
! Public procedures:
!                   Type(ldata)
!                   Generic Init
!***************************************************
Module Defines
  !***********************************************
  !*                 EXTERNS                     *
  !***********************************************
  Use tools
  Implicit none
  
  Private
  Public :: lData
  Type lData
     Integer dim
     Real(dp), Allocatable:: a(:,:), b(:), x(:)
   Contains
     Procedure, Private, Pass :: InitAB
     Procedure, Private, Pass :: InitA
     Procedure, Private, Pass :: InitFile
     Generic, Public :: Init => InitAB, InitA, InitFile

  End type lData

!***********************************************
!*          LOCAL PRIVATE VARIABLES            *
!***********************************************
  Integer, parameter :: log=98

!*************************************************
!*          MODULE PROCEDURES DEFINITION         *
!*************************************************
Contains

!***************************************************
! InitAB:
!     When the input is the dimention and the
!     components A and b from a linear system Ax=b  
!     from the main program, this routine allocates
!     and stores their values in a structure 'data' 
!  
! Parameters:
!     Input, dim, mat, vec_b
!     Output, data%dim, data%a, data%x, data%b
!***************************************************
  Subroutine InitAB(data, dim, mat, vec_b)
    Implicit none
    Class(lData) :: data
    Integer :: dim
    Real(dp) :: mat(dim, dim), vec_b(dim), vec_x(dim)
    Integer :: totalSize

    data%dim=dim
    If(.not.allocated(data%a)) Call AllocateVars(data)
    ! assign data to structure
    data%a=mat
    data%b=vec_b
    data%x=0.d0
    Print*, "Data initialized..."

    Open(unit=log, file='Define.log',action='write',position="append")
    Write(log,*) "********************************************"
    Write(log,*) "Data succesfully initiallized from program"
    totalSize=sizeof(data%a)+sizeof(data%b)+sizeof(data%x)
    Write(log,'(A,I4,A)') "Memory Allocated:", totalSize, "bytes"
    write(log,*)"----------------------------------------------"
    write(log,'(A,I0,A,I0)')"System of dimentions:", data%dim, "x", data%dim
    Write(log,*) "********************************************"
    
  End Subroutine InitAB

!***************************************************
! InitA:
!     When the input is the dimention and the
!     matrix A from the main program, this routine
!     allocates and stores its values in a
!     structure 'data' 
!  
! Parameters:
!     Input, dim, mat
!     Output, data%dim, data%a, data%x
!***************************************************
  Subroutine InitA(data, dim, mat)
    Implicit none
    Class(lData) :: data
    Integer :: dim
    Real(dp) :: mat(dim, dim), vec_x(dim)
    Integer :: totalSize

    data%dim=dim
    If(.not.allocated(data%a)) Call AllocateVars(data)
    ! assign data to structure
    data%a=mat
    data%x=0.d0
    data%b=0.d0
    Print*, "Data initialized..."

    Open(unit=log, file='Define.log',action='write',position="append")
    Write(log,*) "********************************************"
    Write(log,*) "Data succesfully initiallized from program"
    totalSize=sizeof(data%a)+sizeof(data%x)
    Write(log,'(A,I4,A)') "Memory Allocated:", totalSize, "bytes"
    write(log,*)"----------------------------------------------"
    write(log,'(A,I0,A,I0)')"Matrix of dimentions", data%dim, "x", data%dim
    Write(log,*) "********************************************"
  End Subroutine InitA

!***************************************************
! InitFile:
!     When the input is the dimention and the
!     components A and b from a linear system Ax=b  
!     from a datafile, this routine allocates
!     and stores their values in a structure 'data' 
!  
! Parameters:
!     Input, datafile
!     Output, data%dim, data%a, data%x, data%b
!***************************************************
  Subroutine InitFile(data, fileData)
    Implicit none
    Class(lData) :: data
    Logical ex
    Character(*), Optional:: fileData
    Integer :: totalSize
    Inquire(file=fileData, EXIST=ex)
    If(ex) then
       Call ReadData(data, fileData)
       Open(unit=log, file='Define.log',action='write',position="append")
       Write(log,*) "********************************************"
       Write(log,'(3A)') "Data succesfully initiallized from file: '", fileData, "'"
       totalSize=sizeof(data%a)+sizeof(data%b)+sizeof(data%x)
       Write(log,'(A,I4,A)') "Memory Allocated:", totalSize, "bytes"
       write(log,*)"----------------------------------------------"
       write(log,'(A,I0,A,I0)')"System of dimentions: ", data%dim, "x", data%dim
       Write(log,*) "********************************************"
    Else
       Print*,"The file: ", "'",fileData,"'", " does not exist.."
       Print*, "Program STOP.."
       Print*,
       Write(log,*) "********************************************"
       Write(log,'(5A)') "The file: ", "'",fileData,"'", " does not exist"
       Write(log,'(A)')"Program STOP"
       Write(log,*) "********************************************"
       Call DemoFile(fileData)
       Stop
    End If
  End Subroutine InitFile

!***************************************************
! ReadData:
!     Reads linear system elements from datafile
!     and allocates them in 'data' structure  
!  
! Parameters:
!     Input, datafile
!     Output, data%dim, data%a, data%x, data%b
!***************************************************
  Subroutine ReadData(data, fileData)
    Implicit none
    Class(lData) :: data
    Integer i, j
    Character(*) fileData
    Print*, "Reading file: ",fileData
    Open(1,file=trim(fileData))
    Read(1,*)
    Read(1,*)data%dim
    If(.not.allocated(data%a)) Call AllocateVars(data)
    Read(1,*)
    Do i=1,data%dim
       Read(1,*)(data%a(i,j), j=1,data%dim)
    End Do
    Read(1,*)
    Do i=1,data%dim
       Read(1,*) data%b(i)
    End Do
    data%x = 0.d0
    Print"(2x,a,/)", "Data initialized..."

  End Subroutine ReadData

!***************************************************
! AllocateVars:
!     Allocates all elements of linear system
!     with given dimention dim.  
!***************************************************
  Subroutine AllocateVars(data)
    Implicit none
    Class(ldata) :: data
    Allocate(data%a(data%dim, data%dim), data%b(data%dim), data%x(data%dim))
  End Subroutine AllocateVars

!***************************************************
! DemoFile: Format example of linear system file.
!***************************************************
Subroutine DemoFile(file1)
  Implicit none
  Integer, Parameter:: u=99
  Character(*) file1
  Open(unit=u, file=file1//".example")
  Write(u,*)"dim"
  Write(u,*)"4"
  Write(u,*)"matrix A"
  Write(u,*)"4 8 4 0"
  Write(u,*)"1 5 4 -3"
  Write(u,*)"1 4 7 2"
  Write(u,*)"1 3 0 -2"
  Write(u,*)"vector B"
  Write(u,*)"8"
  Write(u,*)"-4"
  Write(u,*)"10"
  Write(u,*)"-4"
  Close(u)
  Print*, "Demo file created '", file1//".example'"
  Close (u)
End Subroutine DemoFile

End Module Defines

