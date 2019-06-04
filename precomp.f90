Program recomp
  Implicit none
  Integer, parameter :: functionsData=1 , Functions=2, projectData=3
  Integer :: i, iSource, imat, nSource, nmat, numberSource, numbermat
  Integer :: nElem, nNodes, iElem, iNode, date_time(8)
  Character (len = 250) :: SourceFun, kcoefFun, projectName, path, method
  Character (len = 250) :: rhoFun, cFun

  Open(projectData, file = 'projectData.dat')
  Read(projectData,'(250A)') projectName
  Read(projectData,'(250A)') path
  Close(projectData)

  Print*,'Starting Poisson 1D'
  Call date_and_time(VALUES=date_time)
  Print*,'Date:',date_time(3),"/",date_time(2),"/", date_time(1)
  Print*,'Hour:',date_time(5),":",date_time(6),":",date_time(7)
  Print*,
  Print*,'Writing functions...'
  Print*,
  
  Open(functionsData, file = trim(projectName)//'.dat')

  Do i = 1, 8
  Read(functionsData,*)
  End Do
  Read(functionsData,'(15x,250A)') method
  Select case(method)
  Case('Transient_State')
  Do i = 1, 5
  Read(functionsData,*)
  End Do
  End Select
  Do i = 1, 4
  Read(functionsData,*)
  End Do
  Read(functionsData,*) nElem
  Read(functionsData,*)
  Read(functionsData,*) nNodes
   Do i = 1, 7
  Read(functionsData,*)
  End Do
  Do iElem = 1, nElem
  Read(functionsData,*)
  End Do
  Do i = 1, 5
  Read(functionsData,*)
  End Do
  Do iNode = 1, nNodes
     Read(functionsData,*) 
  End Do
    
  Do i = 1, 4
  Read(functionsData,*)
  End Do
  Read(functionsData,*) nmat
  Read(functionsData,*)
  Read(functionsData,*)

  Open(Functions, file = trim(path)//'/Poisson1D/Source/functions.f90')
  write(Functions,*)    'Function kCoef(i,Xvalue)' 
  write(Functions,*)    'Implicit none'
  write(Functions,*)    'Real(8) :: vecfunkcoef(',nmat,')'
  write(Functions,*)    'Real(8), Intent(In) :: Xvalue'
  write(Functions,*)    'Real(8) :: kcoef,X'
  write(Functions,*)    'Integer, Intent(In) :: i'
  write(Functions,*)    'x=Xvalue'
  
  Do imat = 1, nmat
  Read(functionsData,'(i2,250A)') numbermat, kcoefFun
  write(Functions,*)    'vecfunkcoef(',numbermat,')='//trim(kCoefFun)
  Print*, 'vecfunkcoef(',numbermat,')='//trim(kCoefFun)
  End Do
  Read(functionsData,*)
  Read(functionsData,*)

  write(Functions,*)    'kCoef= vecfunkcoef(i)'  
  write(Functions,*)    'end function kCoef'
  write(Functions,*)
  
  write(Functions,*)
  write(Functions,*)    'Function rho(i,Xvalue)' 
  write(Functions,*)    'Implicit none'
  write(Functions,*)    'Real(8) :: vecfunrho(',nmat,')'
  write(Functions,*)    'Real(8), Intent(In) :: Xvalue'
  write(Functions,*)    'Real(8) :: rho,X'
  write(Functions,*)    'Integer, Intent(In) :: i'
  write(Functions,*)    'x=Xvalue'
  
  Do imat = 1, nmat
  Read(functionsData,'(i2,250A)') numbermat, rhoFun
  write(Functions,*)    'vecfunrho(',numbermat,')='//trim(rhoFun)
  Print*, 'vecfunrho(',numbermat,')='//trim(rhoFun)
  End Do
  Read(functionsData,*)
  Read(functionsData,*)

  write(Functions,*)    'rho= vecfunrho(i)'  
  write(Functions,*)    'end function rho'
  write(Functions,*)

  write(Functions,*)
  write(Functions,*)    'Function c(i,Xvalue)' 
  write(Functions,*)    'Implicit none'
  write(Functions,*)    'Real(8) :: vecfunc(',nmat,')'
  write(Functions,*)    'Real(8), Intent(In) :: Xvalue'
  write(Functions,*)    'Real(8) :: c,X'
  write(Functions,*)    'Integer, Intent(In) :: i'
  write(Functions,*)    'x=Xvalue'
  
  Do imat = 1, nmat
  Read(functionsData,'(i2,250A)') numbermat, cFun
  write(Functions,*)    'vecfunc(',numbermat,')='//trim(cFun)
  Print*, 'vecfunc(',numbermat,')='//trim(cFun)
End Do

  write(Functions,*)    'c= vecfunc(i)'  
  write(Functions,*)    'end function c'
  write(Functions,*)

  Do i = 1, 4
  Read(functionsData,*)
  End Do
  Read(functionsData,*) nSource
  Read(functionsData,*)
  Read(functionsData,*)
  
  write(Functions,*)
  write(Functions,*)   'Function Source(i,Xvalue)' 
  write(Functions,*)   'Implicit none'
  write(Functions,*)   'Real(8) :: vecfunSource(',nSource,')'
  write(Functions,*)   'Real(8), Intent(In) :: Xvalue'
  write(Functions,*)   'Real(8) :: Source,X'
  write(Functions,*)   'Integer, Intent(In) :: i'
  write(Functions,*)   'x=Xvalue'
  
     Do iSource = 1, nSource
  Read(functionsData,'(i2,250A)') numberSource, SourceFun
  write(Functions,*)   'vecfunSource(',numberSource,')='//trim(SourceFun)
  Print*,'vecfunSource(',numberSource,')='//trim(SourceFun)
End Do

  write(Functions,*)   'Source= vecfunSource(i)'
  write(Functions,*)   'End Function Source'
  close(functionsData)
  close(Functions)

  Print*, 'Compiling functions...'
  Call execute_command_line('cd;cd '//trim(path)//'/ ;make')
  Print*, 'Finish compiling functions'
End Program recomp
  
